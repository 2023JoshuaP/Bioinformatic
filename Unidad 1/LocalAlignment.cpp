#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

const int MATCH_SCORE = 2;
const int MISMATCH_PENALTY = -1;
const int GAP_PENALTY = -2;

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

AlignmentResult smith_waterman(const std::string &sequence_1, const std::string &sequence_2, int match_score, int mismatch_penalty, int gap_penalty) {
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
    int i = end_i, j = end_j;

    while (i > 0 && j > 0 && score_matrix[i][j] > 0) {
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

    return {aligned_seq1, aligned_seq2, max_score, i, end_i, j, end_j};
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
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file.fasta>" << std::endl;
        return 1;
    }

    std::vector<FastaRecordStructure> records;
    std::string file_fasta = argv[1];
    read_file_fasta(file_fasta, records);

    std::cout << "Number of records read: " << records.size() << std::endl;
    for (const auto &record : records) {
        std::cout << "Identifier: " << record.identifier << "\n";
        // std::cout << "Sequence: " << record.sequence << "\n\n";
        // std::cout << "sequence_1: " << record.sequence.substr(0, 10) << "...\n";
    }

    // const int LIMIT = 100;

    // for (const auto &record : records) {
        // std::cout << record.identifier << "\n";
        // if (record.sequence.length() > LIMIT) {
            // std::cout << record.sequence.substr(0, LIMIT) << "...\n";
        // } else {
            // std::cout << record.sequence << "\n";
        // }
    // }

    for (size_t i = 0; i < records.size(); i++) {
        for (size_t j = i + 1; j < records.size(); j++) {
            AlignmentResult result = smith_waterman(records[i].sequence, records[j].sequence, MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY);
            std::cout << "Alignment between " << records[i].identifier << " and " << records[j].identifier << ":\n";
            std::cout << "Score: " << result.score << "\n";
            std::cout << "Aligned Sequence 1: " << result.sequence_1 << "\n";
            std::cout << "Aligned Sequence 2: " << result.sequence_2 << "\n";
        }
    }

    return 0;
}