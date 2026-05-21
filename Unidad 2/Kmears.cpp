#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

using namespace std;
using KmerMap = unordered_map<string, int>;
using KmerSet = unordered_set<string>;

struct FastaSequence {
    string header;
    string sequence;
};

vector<FastaSequence> readFastaFile(const string& filename) {
    ifstream file_fasta(filename);
    vector<FastaSequence> sequences;

    if (!file_fasta.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return sequences;
    }

    string line;
    FastaSequence current_sequence;

    while (getline(file_fasta, line)) {
        if (line.empty()) continue;

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

KmerMap countMers(string sequence, int k) {
    int n_threads = omp_get_max_threads();
    vector<KmerMap> local_kmer_maps(n_threads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& local_map = local_kmer_maps[tid];

        local_map.reserve(100000);

        #pragma omp for schedule(dynamic, 10000)
        for (size_t i = 0; i <= sequence.size() - k; ++i) {
            string kmer = sequence.substr(i, k);
            local_map[kmer]++;
        }
    }

    KmerMap global;
    for (auto& local : local_kmer_maps) {
        for (const auto& [kmer, freq] : local) {
            global[kmer] += freq;
        }
    }

    return global;
}

KmerSet buildKmerSet(const string& sequence, int k) {
    KmerSet kmers;
    kmers.reserve(sequence.size());

    for (size_t i = 0; i + k <= sequence.size(); ++i) {
        kmers.insert(sequence.substr(i, k));
    }

    return kmers;
}

void sharedKmers(const KmerSet& set1, const KmerSet& set2, int limit = 10) {
    int count = 0;

    for (const auto& kmer : set1) {
        if (set2.count(kmer)) {
            cout << "Shared kmer: " << kmer << endl;
            if (++count >= limit) {
                break;
            }
        }
    }

    cout << "Shared kmers: " << count << endl;
}

void uniqueKmers(const KmerSet& set1, const KmerSet& set2, const string& label, int limit = 10) {
    int count = 0;

    for (const auto& kmer : set1) {
        if (!set2.count(kmer)) {
            cout << "Unique kmer in " << label << ": " << kmer << endl;
            if (++count >= limit) {
                break;
            }
        }
    }
}

void compareKmears(const KmerSet& set1, const KmerSet& set2) {
    size_t intersection = 0;

    for (const auto& kmer : set1) {
        if (set2.count(kmer)) {
            intersection++;
        }
    }

    size_t union_size = set1.size() + set2.size() - intersection;
    double jaccard_index = (double)intersection / union_size;

    cout << "\n-------- COMPARISON --------\n";

    cout << "Set 1 size: " << set1.size() << endl;
    cout << "Set 2 size: " << set2.size() << endl;
    cout << "Shared k-mers: " << intersection << endl;
    cout << "Unique in Seq1: " << (set1.size() - intersection) << endl;
    cout << "Unique in Seq2: " << (set2.size() - intersection) << endl;
    cout << "Union: " << union_size << endl;

    cout << "\nJaccard Index: " << jaccard_index << endl;

    if (jaccard_index > 0.8) {
        cout << "Interpretation: Very high similarity (almost identical genomes)\n";
    }
    else if (jaccard_index > 0.5) {
        cout << "Interpretation: Moderate similarity (closely related organisms)\n";
    }
    else if (jaccard_index > 0.2) {
        cout << "Interpretation: Low similarity (related but different genomes)\n";
    }
    else {
        cout << "Interpretation: Very low similarity (distant organisms)\n";
    }

    sharedKmers(set1, set2, 10);
    uniqueKmers(set1, set2, "Sequence 1", 10);
    uniqueKmers(set2, set1, "Sequence 2", 10);
}

vector<pair<string, int>> sortKmers(const unordered_map<string,int>& kmers) {
    vector<pair<string,int>> vec(kmers.begin(), kmers.end());

    sort(vec.begin(), vec.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    return vec;
}

void statistics(const string& filename, const KmerMap& kmers, const vector<pair<string, int>> sorted, int k, size_t total_kmers) {
    ofstream outputfile(filename);

    if (!outputfile.is_open()) {
        cerr << "Error opening output file: " << filename << endl;
        return;
    }

    outputfile << "k value: " << k << "\n";
    outputfile << "Total k-mers: " << total_kmers << "\n";
    outputfile << "Unique k-mers: " << kmers.size() << "\n\n";

    outputfile << "Most frequent k-mer:\n";
    outputfile << sorted.front().first << " -> " << sorted.front().second << "\n\n";

    outputfile << "Least frequent k-mer:\n";
    outputfile << sorted.back().first << " -> " << sorted.back().second << "\n\n";

    outputfile << "===== TOP 20 K-MERS =====\n";

    for (int i = 0; i < 20 && i < sorted.size(); ++i) {
        outputfile << sorted[i].first << " : " << sorted[i].second << "\n";
    }

    outputfile.close();
}

void csvResults(const string& filename, const vector<pair<string, int>>& sorted) {
    ofstream outputfile(filename);

    if (!outputfile.is_open()) {
        cerr << "Error opening output file: " << filename << endl;
        return;
    }

    outputfile << "k-mer,Frequency\n";

    for (const auto& [kmer, freq] : sorted) {
        outputfile << kmer << "," << freq << "\n";
    }

    outputfile.close();
}

void worldCloud(const vector<pair<string, int>>& sorted, int top = 10) {
    int max_freq = sorted.front().second;

    for (int i = 0; i < top && i < sorted.size(); ++i) {
        int stars = (50 * sorted[i].second) / max_freq;
        cout << sorted[i].first << " ";
        for (int j = 0; j < stars; ++j) {
            cout << "*";
        }
        cout << " (" << sorted[i].second << ")";
        cout << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc == 3) {
        string file = "Files/" + string(argv[1]);
        int k = stoi(argv[2]);

        auto seqs = readFastaFile(file);

        if (seqs.empty()) {
            cerr << "Error reading FASTA file\n";
            return 1;
        }

        string seq = seqs[0].sequence;

        cout << "Sequence length: " << seq.size() << endl;
        cout << "Using k = " << k << endl;

        auto kmers = countMers(seq, k);
        auto sorted = sortKmers(kmers);

        cout << "\n----- TOP 10 K-MERS -----\n";

        for (int i = 0; i < 10 && i < sorted.size(); ++i) {
            cout << sorted[i].first << " -> " << sorted[i].second << endl;
        }

        worldCloud(sorted, 10);

        size_t total_kmers = seq.size() - k + 1;
        statistics("kmer_statistics.txt", kmers, sorted, k, total_kmers);
        csvResults("kmer_frequencies.csv", sorted);

        cout << "\nResults saved.\n";
    }
    else if (argc == 4) {
        string file1 = "Files/" + string(argv[1]);
        string file2 = "Files/" + string(argv[2]);
        int k = stoi(argv[3]);

        auto seqs1 = readFastaFile(file1);
        auto seqs2 = readFastaFile(file2);

        if (seqs1.empty() || seqs2.empty()) {
            cerr << "Error reading FASTA files\n";
            return 1;
        }

        string seq1 = seqs1[0].sequence;
        string seq2 = seqs2[0].sequence;

        cout << "Seq1 length: " << seq1.size() << endl;
        cout << "Seq2 length: " << seq2.size() << endl;
        cout << "Using k = " << k << endl;

        auto set1 = buildKmerSet(seq1, k);
        auto set2 = buildKmerSet(seq2, k);

        compareKmears(set1, set2);
    }
    else {
        cerr << "Usage:\n";
        cerr << "Single sequence:\n";
        cerr << "./Kmears <fasta> <k>\n\n";
        cerr << "Genome comparison:\n";
        cerr << "./Kmears <fasta1> <fasta2> <k>\n";
        return 1;
    }

    return 0;
}