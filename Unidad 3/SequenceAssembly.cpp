#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <unordered_map>

using namespace std;
using KmerMap = unordered_map<string, int>;

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

int main() {
    return 0;
}