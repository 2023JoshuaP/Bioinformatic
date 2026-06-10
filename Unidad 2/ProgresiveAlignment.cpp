#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <limits>

using namespace std;

const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -2;

struct FastaSequence {
    string header;
    string sequence;
};

struct AlignmentResult {
    string aligned_seq1;
    string aligned_seq2;
    int score;
};

struct GuideTree {
    int sequence_index1;
    int sequence_index2;
    double distance;
    int size;
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

/* Distance Matrix */
vector<vector<double>> computeDistanceMatrix(const vector<FastaSequence> &sequences) {
    int n = sequences.size();
    vector<vector<double>> distance_matrix(n, vector<double>(n, 0.0));
    vector<int> self_scores(n, 0);

    for (int i = 0; i < n; i++) {
        self_scores[i] = needle_wunsch(sequences[i].sequence, sequences[i].sequence, MATCH, MISMATCH, GAP).score;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            AlignmentResult result = needle_wunsch(sequences[i].sequence, sequences[j].sequence, MATCH, MISMATCH, GAP);
            double distance = 1.0 - (static_cast<double>(result.score) / max(self_scores[i], self_scores[j]));

            distance_matrix[i][j] = distance;
            distance_matrix[j][i] = distance;
        }
    }

    return distance_matrix;
}

/* Calculate UPGMA */
vector<GuideTree> computeUPGMA(vector<vector<double>> &distance_matrix, int n) {
    vector<GuideTree> guide_tree;
    vector<int> cluster_sizes(n, 1);
    vector<bool> active(n, true);          // ✅ all start active

    for (int step = 0; step < n - 1; step++) {
        int min_i = -1, min_j = -1;
        double min_distance = numeric_limits<double>::max();

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (active[i] && active[j] && distance_matrix[i][j] < min_distance) {  // ✅
                    min_distance = distance_matrix[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        guide_tree.push_back({min_i, min_j});

        for (int k = 0; k < n; k++) {
            if (!active[k] || k == min_i || k == min_j) {
                continue;
            }
            distance_matrix[min_i][k] = distance_matrix[k][min_i] =
                (cluster_sizes[min_i] * distance_matrix[min_i][k] +
                 cluster_sizes[min_j] * distance_matrix[min_j][k]) /
                (cluster_sizes[min_i] + cluster_sizes[min_j]);
        }

        cluster_sizes[min_i] += cluster_sizes[min_j];
        active[min_j] = false;
    }

    return guide_tree;
}

/* Process gaps */
vector<string> applyGaps(const vector<string> &group, const string &aligned_rap) {
    vector<string> result;

    for (const string &sequence : group) {
        string new_sequence = "";
        int sequence_index = 0;

        for (char c : aligned_rap) {
            if (c == '-') {
                new_sequence += '-';
            }
            else {
                new_sequence += sequence[sequence_index];
                sequence_index++;
            }
        }
        result.push_back(new_sequence);
    }

    return result;
}

vector<string> mergeGroups(const vector<string> &group1, const vector<string> &group2) {
    string rep_a = group1[0];
    string rep_b = group2[0];

    AlignmentResult result = needle_wunsch(rep_a, rep_b, MATCH, MISMATCH, GAP);

    vector<string> aligned_group1 = applyGaps(group1, result.aligned_seq1);
    vector<string> aligned_group2 = applyGaps(group2, result.aligned_seq2);

    vector<string> merged;
    merged.insert(merged.end(), aligned_group1.begin(), aligned_group1.end());
    merged.insert(merged.end(), aligned_group2.begin(), aligned_group2.end());

    return merged;
}

/* Guide tree to drive merges */
vector<string> performMSA(const vector<FastaSequence> &sequences, const vector<GuideTree> &guide_tree) {
    vector<vector<string>> groups;
    for (const auto &seq : sequences) {
        groups.push_back({seq.sequence});
    }

    for (const auto &merge : guide_tree) {
        vector<string> merged_groups = mergeGroups(groups[merge.sequence_index1], groups[merge.sequence_index2]);
        groups[merge.sequence_index1] = merged_groups;
        groups[merge.sequence_index2] = {};
    }

    for (const auto &group : groups) {
        if (!group.empty()) {
            return group;
        }
    }

    return {};
}

int main(int argc, char *argv[]) {
    return 0;
}