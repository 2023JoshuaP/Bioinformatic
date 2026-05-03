#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

const int LINE_WIDTH = 60;
const int MATCH_SCORE = 1;
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

void print_alignment(const AlignmentResult &res,
                     const string &header1, const string &header2) {
    int len = (int)res.aligned_seq1.size();
 
    // Línea de identidad y estadísticas
    string identity_line(len, ' ');
    int matches = 0, mismatches = 0, gaps = 0;
    for (int k = 0; k < len; k++) {
        char c1 = res.aligned_seq1[k], c2 = res.aligned_seq2[k];
        if (c1 == '-' || c2 == '-')          { gaps++;       identity_line[k] = ' '; }
        else if (c1 == c2)                    { matches++;    identity_line[k] = '|'; }
        else                                  { mismatches++; identity_line[k] = '.'; }
    }
 
    double pct_id = 100.0 * matches / len;
 
    // ── Cabecera ──
    cout << "\n";
    cout << "------------------------------------------------------------\n";
    cout << "          NEEDLEMAN-WUNSCH  |  Alineamiento Global\n";
    cout << "------------------------------------------------------------\n";
    cout << "  Seq 1 : " << header1 << "\n";
    cout << "  Seq 2 : " << header2 << "\n";
    cout << "  Score : " << res.score << "\n";
    cout << "------------------------------------------------------------\n";
    cout << "  Identidad  : " << matches << "/" << len
         << " (" << fixed << setprecision(1) << pct_id << "%)\n";
    cout << "  Mismatches : " << mismatches << "\n";
    cout << "  Gaps       : " << gaps << "\n";
    cout << "------------------------------------------------------------\n\n";
 
    // ── Alineamiento por bloques ──
    int pos1 = 1, pos2 = 1;
    for (int start = 0; start < len; start += LINE_WIDTH) {
        int end = min(start + LINE_WIDTH, len);
        string block1 = res.aligned_seq1.substr(start, end - start);
        string block2 = res.aligned_seq2.substr(start, end - start);
        string blockI = identity_line.substr(start, end - start);
 
        // Contar posición real (sin gaps)
        int end_pos1 = pos1, end_pos2 = pos2;
        for (char c : block1) if (c != '-') end_pos1++;
        for (char c : block2) if (c != '-') end_pos2++;
 
        cout << "Seq1  " << setw(5) << pos1 << "  " << block1 << "  " << end_pos1 - 1 << "\n";
        cout << "           " << blockI << "\n";
        cout << "Seq2  " << setw(5) << pos2 << "  " << block2 << "  " << end_pos2 - 1 << "\n\n";
 
        pos1 = end_pos1;
        pos2 = end_pos2;
    }
    cout << "------------------------------------------------------------\n\n";
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cerr << "Uso: " << argv[0] << " <archivo1.fasta> <archivo2.fasta>\n";
        return 1;
    }
 
    vector<FastaSequences> sequences1, sequences2;
    read_fasta(argv[1], sequences1);
    read_fasta(argv[2], sequences2);
 
    if (sequences1.empty() || sequences2.empty()) {
        cerr << "Error: no se pudieron leer las secuencias.\n";
        return 1;
    }
 
    AlignmentResult result = needleman_wunsch(
        sequences1[0].sequence,
        sequences2[0].sequence,
        MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY
    );
 
    print_alignment(result, sequences1[0].header, sequences2[0].header);
 
    return 0;
}