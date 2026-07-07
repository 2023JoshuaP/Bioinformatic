#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <unordered_map>

using namespace std;
using KmerMap = unordered_map<string, int>;
using Graph = unordered_map<string, vector<string>>;

KmerMap countMers(const string &sequence, int k) {
    KmerMap frequence;

    if (sequence.size() < k) {
        cerr << "Error: k is larger than the sequence length\n";
        return frequence;
    }

    for (size_t i = 0; i < sequence.size() - k + 1; i++) {
        string kmer = sequence.substr(i, k);
        frequence[kmer]++;
    }

    return frequence;
}

KmerMap combineMers(const vector<string> &fragments, int k) {
    KmerMap global;

    for (const auto &fragment : fragments) {
        KmerMap local = countMers(fragment, k);
        for (const auto &[kmer, count] : local) {
            global[kmer] += count;
        }
    }

    return global;
}

Graph BuildDeBruijnGraph(const KmerMap &kmerCounts) {
    Graph graph;

    for (const auto &[kmer, count] : kmerCounts) {
        string prefix = kmer.substr(0, kmer.size() - 1);
        string suffix = kmer.substr(1);
        for (int i = 0; i < count; ++i) {
            graph[prefix].push_back(suffix);
        }
    }

    return graph;
}

void printGraph(const Graph &graph) {
    for (const auto &[node, edges] : graph) {
        cout << node << " -> ";
        for (const auto &edge : edges) {
            cout << edge << " ";
        }
        cout << endl;
    }
}

int main() {
    vector<string> fragments = { "ATGCGATGAC" };
    int k = 4;
    KmerMap kmers = combineMers(fragments, k);
    cout << "Generate Kmers:" << endl;
    for (const auto &[kmer, count] : kmers) {
        cout << kmer << ": " << count << endl;
    }

    cout << "\nBuild De Bruijn Graph:" << endl;
    Graph graph = BuildDeBruijnGraph(kmers);
    printGraph(graph);
    return 0;
}