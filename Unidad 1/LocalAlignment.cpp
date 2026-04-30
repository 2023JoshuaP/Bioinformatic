#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>

const int MATCH_SCORE = 2;
const int MISMATCH_PENALTY = -1;
const int GAP_PENALTY = -2;
const int DISPLAY = 100;

struct FastaRecordStructure {
    std::string identifier;
    std::string sequence;
};

struct AlignmentResult {
    std::string sequence_1;
    std::string sequence_2;
    int score;
    int start_pos_seq1, end_pos_seq1;
    int start_pos_seq2, end_pos_seq2;
    std::vector<std::vector<int>> score_matrix;
    std::vector<std::pair<int, int>> alignment_path;
};

void read_file_fasta(const std::string &file_fasta, std::vector<FastaRecordStructure> &records) {
    std::ifstream file(file_fasta);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_fasta << std::endl;
        return;
    }

    auto strip_cr = [](std::string &line) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
    };

    std::string line, sequence, header;
    while (std::getline(file, line)) {
        strip_cr(line);
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            if (!header.empty()) {
                records.push_back({header, sequence});
                sequence.clear();
            }
            std::string raw_header = line.substr(1);
            size_t first_space = raw_header.find(' ');
            header = (first_space != std::string::npos) ? raw_header.substr(0, first_space) : raw_header;
        }
        else {
            for (char &c : line) {
                c = std::toupper((unsigned char)c);
            }
            sequence += line;
        }
    }

    if (!header.empty()) {
        records.push_back({header, sequence});
    }

    file.close();
}

void export_sequences(const std::string &output_path, const std::vector<FastaRecordStructure> &records) {
    std::ofstream output(output_path);
    if (!output.is_open()) {
        std::cerr << "Error opening file for writing: " << output_path << std::endl;
        return;
    }

    for (const auto &record : records) {
        output << ">" << record.identifier << "\n";
        for (size_t i = 0; i < record.sequence.length(); i += DISPLAY) {
            output << record.sequence.substr(i, DISPLAY) << "\n";
        }
    }

    output.close();
    std::cout << "Sequences exported to: " << output_path << std::endl;
}

void export_aligned_sequences(const AlignmentResult &result, const std::string &output_path) {
    std::ofstream output(output_path);
    if (!output.is_open()) {
        std::cerr << "Error opening file for writing: " << output_path << std::endl;
        return;
    }

    output << "Secuencia 1:\n";
    for (size_t k = 0; k < result.sequence_1.length(); k += DISPLAY)
        output << result.sequence_1.substr(k, DISPLAY) << "\n";

    output << "\nSecuencia 2:\n";
    for (size_t k = 0; k < result.sequence_2.length(); k += DISPLAY)
        output << result.sequence_2.substr(k, DISPLAY) << "\n";

    output.close();
    std::cout << "Aligned sequences exported to: " << output_path << std::endl;
}

AlignmentResult smith_waterman(const std::string &sequence_1, const std::string &sequence_2, int match_score, int mismatch_penalty, int gap_penalty, const std::string &export_file = "") {
    int n = (int)sequence_1.length();
    int m = (int)sequence_2.length();
    int max_score = 0, end_i = 0, end_j = 0;

    std::vector<std::vector<int>> score_matrix(n + 1, std::vector<int>(m + 1, 0));

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int diagonal_score = score_matrix[i - 1][j - 1] + (sequence_1[i - 1] == sequence_2[j - 1] ? match_score : mismatch_penalty);
            int up_score = score_matrix[i - 1][j] + gap_penalty;
            int left_score = score_matrix[i][j - 1] + gap_penalty;
            score_matrix[i][j] = std::max({0, diagonal_score, up_score, left_score});

            if (score_matrix[i][j] > max_score) {
                max_score = score_matrix[i][j];
                end_i = i;
                end_j = j;
            }
        }
    }

    std::string aligned_seq1, aligned_seq2;
    std::vector<std::pair<int, int>> alignment_path;
    int i = end_i, j = end_j;

    while (i > 0 && j > 0 && score_matrix[i][j] > 0) {
        alignment_path.push_back({i, j});
        int diagonal = score_matrix[i - 1][j - 1] + (sequence_1[i - 1] == sequence_2[j - 1] ? match_score : mismatch_penalty);
        int up = score_matrix[i - 1][j] + gap_penalty;

        if (score_matrix[i][j] == diagonal) {
            aligned_seq1 += sequence_1[i - 1];
            aligned_seq2 += sequence_2[j - 1];
            i--;
            j--;
        }
        else if (score_matrix[i][j] == up) {
            aligned_seq1 += sequence_1[i - 1];
            aligned_seq2 += '-';
            i--;
        }
        else {
            aligned_seq1 += '-';
            aligned_seq2 += sequence_2[j - 1];
            j--;
        }
    }

    std::reverse(aligned_seq1.begin(), aligned_seq1.end());
    std::reverse(aligned_seq2.begin(), aligned_seq2.end());
    std::reverse(alignment_path.begin(), alignment_path.end());

    return {aligned_seq1, aligned_seq2, max_score, i, end_i - 1, j, end_j - 1, score_matrix, alignment_path};
}

void export_alignment(const AlignmentResult &result, const std::string &base_path) {
    int n = (int)result.score_matrix.size() - 1;
    int m = (int)result.score_matrix[0].size() - 1;

    {
        std::ofstream output(base_path + "_matrix.txt");
        if (output.is_open()) {
            for (int r = 0; r <= n; r++) {
                for (int c = 0; c <= m; c++) {
                    output << result.score_matrix[r][c];
                    if (c < m) {
                        output << " ";
                    }
                }
                output << "\n";
            }
            output.close();
            std::cout << "Alignment matrix exported to: " << base_path + "_matrix.txt" << std::endl;
        }
        else {
            std::cerr << "Error opening file for writing: " << base_path + "_matrix.txt" << std::endl;
        }
    }

    {
        std::ofstream output(base_path + "_path.txt");
        if (output.is_open()) {
            for (const auto &pos : result.alignment_path) {
                output << pos.first << " " << pos.second << "\n";
            }
            output.close();
            std::cout << "Alignment path exported to: " << base_path + "_path.txt" << std::endl;
        }
        else {
            std::cerr << "Error opening file for writing: " << base_path + "_path.txt" << std::endl;
        }
    }

    {
        std::ofstream output(base_path + "_info.txt");
        if (output.is_open()) {
            output << "Score: " << result.score << "\n";
            output << "Start/End Seq1: [" << result.start_pos_seq1 << " - " << result.end_pos_seq1 << "]\n";
            output << "Start/End Seq2: [" << result.start_pos_seq2 << " - " << result.end_pos_seq2 << "]\n";
            output.close();
            std::cout << "Alignment info exported to: " << base_path + "_info.txt" << std::endl;
        }
        else {
            std::cerr << "Error opening file for writing: " << base_path + "_info.txt" << std::endl;
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <file1.fasta> <file2.fasta>" << std::endl;
        return 1;
    }

    std::cout << "  Match score: " << MATCH_SCORE << "\n";
    std::cout << "  Mismatch: " << MISMATCH_PENALTY << "\n";
    std::cout << "  Gap penalty: " << GAP_PENALTY << "\n";

    std::vector<FastaRecordStructure> records1, records2;
    read_file_fasta(argv[1], records1);
    read_file_fasta(argv[2], records2);

    export_sequences("File txt/sequences1.txt", records1);
    export_sequences("File txt/sequences2.txt", records2);

    for (const auto &record1 : records1) {
        for (const auto &record2 : records2) {
            std::cout << "Aligning " << record1.identifier << " with " << record2.identifier << "..." << std::endl;

            auto start = std::chrono::high_resolution_clock::now();
            AlignmentResult result = smith_waterman(record1.sequence, record2.sequence, MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Time: " << std::fixed << std::setprecision(3) << elapsed.count() << " seconds\n\n";
            
            std::string base_path = "File txt/alignment_" + record1.identifier + "_" + record2.identifier;
            export_alignment(result, base_path);
            export_aligned_sequences(result, base_path + "_aligned.txt");
        }
    }

    return 0;
}