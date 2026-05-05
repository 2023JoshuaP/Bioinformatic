#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>
using namespace std;

const int LINE_WIDTH = 60;
const int MATCH_SCORE = 2;
const int MISMATCH_PENALTY = -1;
const int GAP_PENALTY = -2;

struct FastaSequences {
    string header;
    string sequence;
};

struct AlignmentResult {
    string aligned_seq1;
    string aligned_seq2;
    int score;
    int start_pos_seq1, end_pos_seq1;
    int start_pos_seq2, end_pos_seq2;
    vector<vector<int>> score_matrix;
    vector<pair<int, int>> alignment_path;
};

void read_fasta(const string &filename, vector<FastaSequences> &sequences) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    auto strip_cr = [](string &line) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
    };

    string line, header, sequence;

    while (getline(file, line)) {
        strip_cr(line);
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            if (!header.empty()) {
                sequences.push_back({header, sequence});
                sequence.clear();
            }
            string raw_header = line.substr(1);
            size_t first_space = raw_header.find(' ');
            header = (first_space != string::npos) ? raw_header.substr(0, first_space) : raw_header;
        }
        else {
            for (char &c : line) {
                c = toupper((unsigned char)c);
            }
            sequence += line;
        }
    }

    if (!header.empty()) {
        sequences.push_back({header, sequence});
    }

    file.close();
}

AlignmentResult needleman_wunsch(const string &sequence1, const string &sequence2, int match_score, int mismatch_penalty, int gap_penalty) {
    int n = (int)sequence1.size();
    int m = (int)sequence2.size();

    vector<vector<int>> score_matrix(n + 1, vector<int>(m + 1, 0));

    for (int i = 0; i <= n; i++) {
        score_matrix[i][0] = i * gap_penalty;
    }
    for (int j = 0; j <= m; j++) {
        score_matrix[0][j] = j * gap_penalty;
    }

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score_diagonal = score_matrix[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? match_score : mismatch_penalty);
            int score_up = score_matrix[i - 1][j] + gap_penalty;
            int score_left = score_matrix[i][j - 1] + gap_penalty;

            if (score_diagonal >= score_up && score_diagonal >= score_left) {
                score_matrix[i][j] = score_diagonal;
            }
            else if (score_up >= score_left) {
                score_matrix[i][j] = score_up;
            }
            else {
                score_matrix[i][j] = score_left;
            }
        }
    }

    string aligned_sequence1, aligned_sequence2;
    int i = n, j = m;

    while (i > 0 && j > 0) {
        int score_diagonal = score_matrix[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? match_score : mismatch_penalty);
        int score_up = score_matrix[i - 1][j] + gap_penalty;
        int score_left = score_matrix[i][j - 1] + gap_penalty;

        if (score_matrix[i][j] == score_diagonal) {
            aligned_sequence1 += sequence1[i - 1];
            aligned_sequence2 += sequence2[j - 1];
            i--;
            j--;
        }
        else if (score_matrix[i][j] == score_up) {
            aligned_sequence1 += sequence1[i - 1];
            aligned_sequence2 += '-';
            i--;
        }
        else {
            aligned_sequence1 += '-';
            aligned_sequence2 += sequence2[j - 1];
            j--;
        }
    }

    while (i > 0) {
        aligned_sequence1 += sequence1[i - 1];
        aligned_sequence2 += '-';
        i--;
    }

    while (j > 0) {
        aligned_sequence1 += '-';
        aligned_sequence2 += sequence2[j - 1];
        j--;
    }

    reverse(aligned_sequence1.begin(), aligned_sequence1.end());
    reverse(aligned_sequence2.begin(), aligned_sequence2.end());

    return {aligned_sequence1, aligned_sequence2, score_matrix[n][m], 0, n - 1, 0, m - 1, score_matrix, {}};
}

void save_alignment(const AlignmentResult &result, const string &header1, const string &header2) {
    string base_path = "File txt/" + header1 + "_vs_" + header2 + "_alignment";
    
    {
        string csv_path = base_path + "_matrix.csv";
        ofstream csv_file(csv_path);
        
        if (!csv_file.is_open()) {
            cerr << "Error opening file for writing: " << csv_path << endl;
            return;
        }

        int n = (int)result.score_matrix.size();
        int m = (int)result.score_matrix[0].size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                csv_file << result.score_matrix[i][j];
                if (j < m - 1) {
                    csv_file << ",";
                }
            }
            csv_file << "\n";
        }
        csv_file.close();
        cout << "Matrix saved to: " << csv_path << endl;
    }

    {
        string txt_path = base_path + "_alignment.txt";
        ofstream txt_file(txt_path);

        if (!txt_file.is_open()) {
            cerr << "Error opening file for writing: " << txt_path << endl;
            return;
        }

        txt_file << "Score: " << result.score << "\n\n";
        
        int len = (int)result.aligned_seq1.size();

        for (int i = 0; i < len; i += LINE_WIDTH) {
            int end = min(i + LINE_WIDTH, len);

            string seq1 = result.aligned_seq1.substr(i, end - i);
            string seq2 = result.aligned_seq2.substr(i, end - i);

            string match_line;
            for (int k = 0; k < (int)seq1.size(); k++) {
                if (seq1[k] == seq2[k]) {
                    match_line += '|';
                } 
                else if (seq1[k] == '-' || seq2[k] == '-') {
                    match_line += ' ';
                } 
                else {
                    match_line += '.';
                }
            }

            txt_file << seq1 << "\n";
            txt_file << match_line << "\n";
            txt_file << seq2 << "\n\n";
        }

        txt_file.close();
        cout << "Alignment saved to: " << txt_path << endl;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <file1.fasta> <file2.fasta>\n";
        return 1;
    }
 
    vector<FastaSequences> sequences1, sequences2;
    read_fasta(argv[1], sequences1);
    read_fasta(argv[2], sequences2);
 
    if (sequences1.empty() || sequences2.empty()) {
        cerr << "Error opening files.\n";
        return 1;
    }
 
    auto start = chrono::high_resolution_clock::now();
    AlignmentResult result = needleman_wunsch(
        sequences1[0].sequence,
        sequences2[0].sequence,
        MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY
    );
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Time: " << fixed << setprecision(3) << elapsed.count() << " seconds\n";

    save_alignment(result, sequences1[0].header, sequences2[0].header);
 
    return 0;
}