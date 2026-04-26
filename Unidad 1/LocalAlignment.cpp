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
    std::vector<std::vector<int>> score_matrix_h(sequence_1.length() + 1, std::vector<int>(sequence_2.length() + 1, 0));

    for (int i = 1; i <= sequence_1.length(); i++) {
        for (int j = 1; j <= sequence_2.length(); j++) {
            int score_diagonal = score_matrix_h[i - 1][j - 1] + (sequence_1[i - 1] == sequence_2[j - 1] ? match_score : mismatch_penalty);
            int score_up = score_matrix_h[i - 1][j] + gap_penalty;
            int score_left = score_matrix_h[i][j - 1] + gap_penalty;
            score_matrix_h[i][j] = std::max({0, score_diagonal, score_up, score_left});
        }
    }

    int max_score = 0;

    /*Traceback max_score*/

    int end_seq1 = 0, end_seq2 = 0;
    for (int i = 1; i <= sequence_1.length(); i++) {
        for (int j = 1; j <= sequence_2.length(); j++) {
            if (score_matrix_h[i][j] > max_score) {
                max_score = score_matrix_h[i][j];
                end_seq1 = i;
                end_seq2 = j;
            }
        }
    }

    int i = end_seq1, j = end_seq2;
    std::string aligned_seq1, aligned_seq2;

    while (i > 0 && j > 0 && score_matrix_h[i][j] > 0) {
        if (score_matrix_h[i][j] == score_matrix_h[i - 1][j - 1] + (sequence_1[i - 1] == sequence_2[j - 1] ? match_score : mismatch_penalty)) {
            aligned_seq1 = sequence_1[i - 1] + aligned_seq1;
            aligned_seq2 = sequence_2[j - 1] + aligned_seq2;
            i--;
            j--;
        }
        else if (score_matrix_h[i][j] == score_matrix_h[i - 1][j] + gap_penalty) {
            aligned_seq1 = sequence_1[i - 1] + aligned_seq1;
            aligned_seq2 = '-' + aligned_seq2;
            i--;
        }
        else {
            aligned_seq1 = '-' + aligned_seq1;
            aligned_seq2 = sequence_2[j - 1] + aligned_seq2;
            j--;
        }
    }
    return {aligned_seq1, aligned_seq2, max_score, 0, 0, 0, 0};
}

void read_file_fasta(const std::string &file_fasta, std::vector<FastaRecordStructure> &records) {
    std::ifstream file(file_fasta);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_fasta << std::endl;
        return;
    }

    std::string line, sequence, header;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            if (!header.empty()) {
                records.push_back({header, sequence});
                sequence.clear();
            }
            header = line.substr(1);
        }
        else {
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