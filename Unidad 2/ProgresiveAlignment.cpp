#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <algorithm>

using namespace std;

struct FastaSequence {
    string header;
    string sequence;
};

struct AlignmentResult {
    string aligned_seq1;
    string aligned_seq2;
    int score;
};

vector<FastaSequence> readFileFasta(const string &filename) {
    ifstream file_fasta(filename);
    vector<FastaSequence> sequences;

    if (!file_fasta.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return sequences;
    }

    string line;
    FastaSequence current_sequence;

    while(getline(file_fasta, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            if (!current_sequence.header.empty()) {
                sequences.push_back(current_sequence);
                current_sequence = FastaSequence();
            }
            current_sequence.header = line.substr(1);
        }
        else {
            current_sequence.sequence += line;
        }
    }

    if (!current_sequence.header.empty()) {
        sequences.push_back(current_sequence);
    }

    return sequences;
}

AlignmentResult needle_wunsch(const string &seq1, const string &seq2, int match, int mismatch, int gap) {
    int n = seq1.size();
    int m = seq2.size();

    vector<vector<int>> score_matrix(n + 1, vector<int>(m + 1, 0));

    for (int i = 0; i <= n; i++) {
        score_matrix[i][0] = i * gap;
    }
    for (int j = 0; j <= m; j++) {
        score_matrix[0][j] = j * gap;
    }

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score_diagonal = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int score_up = score_matrix[i - 1][j] + gap;
            int score_left = score_matrix[i][j - 1] + gap;

            score_matrix[i][j] = max({score_diagonal, score_up, score_left});
        }
    }

    /* Backtracking process */
    string aligned_seq1, aligned_seq2;
    int i = n, j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int score_diagonal = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int score_up = score_matrix[i - 1][j] + gap;
            int score_left = score_matrix[i][j - 1] + gap;

            if (score_matrix[i][j] == score_diagonal) {
                aligned_seq1 += seq1[i - 1];
                aligned_seq2 += seq2[j - 1];
                i--;
                j--;
            }
            else if (score_matrix[i][j] == score_up) {
                aligned_seq1 += seq1[i - 1];
                aligned_seq2 += '-';
                i--;
            }
            else {
                aligned_seq1 += '-';
                aligned_seq2 += seq2[j - 1];
                j--;
            }
        }
        else if (i > 0) {
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += '-';
            i--;
        }
        else {
            aligned_seq1 += '-';
            aligned_seq2 += seq2[j - 1];
            j--;
        }
    }

    reverse(aligned_seq1.begin(), aligned_seq1.end());
    reverse(aligned_seq2.begin(), aligned_seq2.end());

    return {aligned_seq1, aligned_seq2, score_matrix[n][m]};
}

int main(int argc, char *argv[]) {
    return 0;
}