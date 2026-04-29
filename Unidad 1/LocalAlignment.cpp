#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

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
};

/*
Smith Waterman algorithm for local sequence alignment. It finds the best local alignment between two sequences by constructing a scoring matrix and performing a traceback to identify the optimal alignment. The algorithm uses dynamic programming to efficiently compute the scores and alignments, allowing for gaps and mismatches in the sequences.
*/

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

    if (!export_file.empty()) {
        std::ofstream output(export_file);
        if (output.is_open()) {
            output << n << " " << m << "\n";
            output << sequence_1 << "\n";
            output << sequence_2 << "\n";

            for (int r = 0; r <= n; r++) {
                for (int c = 0; c <= m; c++) {
                    output << score_matrix[r][c];
                    if (c < m) {
                        output << " ";
                    }
                }
                output << "\n";
            }

            output << "PATH\n";
            for (const auto &it : alignment_path) {
                output << it.first << " " << it.second << "\n";
            }
            
            output.close();
            std::cout << "Alignment details exported to: " << export_file << std::endl;
        }
        else {
            std::cerr << "Error opening file for writing: " << export_file << std::endl;
        }
    }

    return {aligned_seq1, aligned_seq2, max_score, i, end_i - 1, j, end_j - 1};
}

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

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <file1.fasta> <file2.fasta>" << std::endl;
        return 1;
    }

    std::vector<FastaRecordStructure> records1, records2;
    read_file_fasta(argv[1], records1);
    read_file_fasta(argv[2], records2);

    std::cout << "Read " << records1.size() << " records from " << argv[1] << std::endl;
    for (const auto &record : records1) {
        std::cout << "Identifier: " << record.identifier << ", Sequence Length: " << record.sequence.length() << std::endl;
    }
    
    std::cout << "Read " << records2.size() << " records from " << argv[2] << std::endl;
    for (const auto &record : records2) {
        std::cout << "Identifier: " << record.identifier << ", Sequence Length: " << record.sequence.length() << std::endl;
    }
    std::cout << std::endl;

    for (const auto &record1 : records1) {
        for (const auto &record2 : records2) {
            AlignmentResult result = smith_waterman(record1.sequence, record2.sequence, MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY, "Files/matrix_NC005089_vs_NC012920.txt");
            std::cout << "Alignment between " << record1.identifier << " and " << record2.identifier << ":\n";
            std::cout << "Score: " << result.score << "\n";
            std::cout << "Start/End Seq1: [" << result.start_pos_seq1 << " - " << result.end_pos_seq1 << "]\n";
            std::cout << "Start/End Seq2: [" << result.start_pos_seq2 << " - " << result.end_pos_seq2 << "]\n";
        }
    }

    return 0;
}