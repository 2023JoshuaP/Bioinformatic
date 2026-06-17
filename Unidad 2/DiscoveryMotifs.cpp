#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <algorithm>
#include <unordered_map>

using namespace std;
using KmerMap = unordered_map<string, int>;

vector<string> sequences_test() {
    return {
        "ATCGTACGATGACCTGATCG",
        "GGTATACGATGACGTTACCA",
        "TTTCTACGATGACCATAGGT",
        "AACGTACGATGACGGGTTAA",
        "CGGATACGATGACTTCCGTA",
        "TACCTACGATGACAGGTACA",
        "GACTTACGATGACCGATAGC",
        "TCGATACGATGACTGGCAAT",
        "AGGCTACGATGACATTCGGA",
        "CCTATACGATGACGGAATTC"
    };
}

vector<string> getAnnexSequences() {
    return {
        "ATCGGCTAACGTAGCTAGCTTGACCGTACGATCGATCGGATCGTAGCTAGCATCGATCGTACGATCGATGCTAGCTAGCATCGATCGATACGATCGTAGCTAGCTACGTAGCTAGCTACGTAGCTTACGATGACGGTACCGATCGATCGTAGCTAACGTA",
 
        "GCTAGCTAGCATCGATCGTAGCTAGCTAGCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATACGTAGCTACGTACGATGACATCGTAGCTAGCTAACGTAGCTAGCTAGCGATCGTAGCTAGCTA",
 
        "CGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTACGTACGATGATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTA",
        
        "AACGTAGCTAGCTAGCATCGATCGTAGCTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTACGTATACGATGACGCTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGT",

        "TAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTAACGTATACGATGTCCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACG",

        "GATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTATACGATGACGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCT",

        "CTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTATACGATGCCGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTA",

        "CGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTATACGATGACCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAAC"
    };
}

struct KmerCandidate {
    string kmer;
    int total_frequency;
    int sequence_coverage;
};

struct ExtractRegion {
    int sequence_index;
    int start_position;
    string region;
};

struct ConservationResult {
    vector<double> percentPerPosition;
    vector<bool> fullyConserved;
    double overallConservation;
};

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

void builsKmerStats(const vector<string> &sequences, int k, KmerMap &frequency, KmerMap &coverage) {
    frequency.clear();
    coverage.clear();

    for (const auto &sequence : sequences) {
        KmerMap local = countMers(sequence, k);
        for (const auto &[kmer, count] : local) {
            frequency[kmer] += count;
            coverage[kmer]++;
        }
    }
}

vector<KmerCandidate> getTopCandidates(const KmerMap &frequency, const KmerMap &coverage, int top_n = 5) {
    vector<KmerCandidate> candidates;
    candidates.reserve(coverage.size());

    for (const auto &[kmer, coverage_count] : coverage) {
        int freq_count = frequency.at(kmer);
        candidates.push_back({kmer, freq_count, coverage_count});
    }

    sort(candidates.begin(), candidates.end(), [](const KmerCandidate &a, const KmerCandidate &b) {
        if (a.sequence_coverage != b.sequence_coverage) {
            return a.sequence_coverage > b.sequence_coverage;
        }
        return a.total_frequency > b.total_frequency;
    });

    if ((int)candidates.size() > top_n) {
        candidates.resize(top_n);
    }

    return candidates;
}

vector<vector<int>> findKmerPositions(const vector<string> &sequences, const string &candidate) {
    vector<vector<int>> positions(sequences.size());
    int k = (int)candidate.size();

    for (size_t i = 0; i < sequences.size(); i++) {
        const string &sequence = sequences[i];
        if (sequence.size() < k) {
            continue;
        }
        for (size_t j = 0; j <= sequence.size() - k; j++) {
            if (sequence.compare(j, k, candidate) == 0) {
                positions[i].push_back(j);
            }
        }
    }

    return positions;
}

vector<ExtractRegion> extractRegions(const vector<string> &sequences, const vector<vector<int>> &positions, int region_size) {
    vector<ExtractRegion> regions;
    for (size_t i = 0; i < sequences.size(); i++) {
        if (positions[i].empty()) {
            continue;
        }
        
        int start = positions[i][0];
        const string &sequence = sequences[i];

        if (start + region_size > (int)sequence.size()) {
            continue;
        }

        string region = sequence.substr(start, region_size);
        regions.push_back({(int)i, start, region});
    }

    return regions;
}

bool verifyAlignment(const vector<ExtractRegion> &regions) {
    if (regions.empty()) {
        return false;
    }

    size_t region_length = regions[0].region.size();
    for (const auto &region : regions) {
        if (region.region.size() != region_length) {
            return false;
        }
    }

    return true;
}

map<char, vector<int>> buildFrequencyMatrix(const vector<ExtractRegion> &regions) {
    map<char, vector<int>> frequency_matrix;
    frequency_matrix['A'] = vector<int>(regions[0].region.size(), 0);
    frequency_matrix['C'] = vector<int>(regions[0].region.size(), 0);
    frequency_matrix['G'] = vector<int>(regions[0].region.size(), 0);
    frequency_matrix['T'] = vector<int>(regions[0].region.size(), 0);
    
    for (const auto &region : regions) {
        for (size_t i = 0; i < region.region.size(); i++) {
            char nucleotide = region.region[i];
            if (frequency_matrix.count(nucleotide)) {
                frequency_matrix[nucleotide][i]++;
            }
        }
    }

    return frequency_matrix;
}

string buildConsensus(const map<char, vector<int>> &frequency_matrix, int region_size) {
    string consensus(region_size, 'N');

    for (int pos = 0; pos < region_size; pos++) {
        char best_base = 'N';
        int best_count = -1;
        for (char base : {'A', 'C', 'G', 'T'}) {
            int count = frequency_matrix.at(base)[pos];
            if (count > best_count) {
                best_count = count;
                best_base = base;
            }
        }
        consensus[pos] = best_base;
    }
    return consensus;
}

ConservationResult evaluateConservation(const map<char, vector<int>> &freqMatrix,const string &consensus,int numSequences,int regionLength) {
    ConservationResult result;
    result.percentPerPosition.resize(regionLength);
    result.fullyConserved.resize(regionLength);
 
    double sumPercent = 0.0;
 
    for (int pos = 0; pos < regionLength; ++pos) {
        char consensusBase = consensus[pos];
        int matchCount = freqMatrix.at(consensusBase)[pos];
        double percent = (numSequences > 0)
                          ? (100.0 * matchCount / numSequences)
                          : 0.0;
 
        result.percentPerPosition[pos] = percent;
        result.fullyConserved[pos] = (matchCount == numSequences);
        sumPercent += percent;
    }
 
    result.overallConservation = (regionLength > 0) ? (sumPercent / regionLength) : 0.0;
    return result;
}

void printFrequencyMatrix(const map<char, vector<int>> &freqMatrix, int regionLength) {
    cout << "\nFrequency Matrix:\n";
    cout << "Pos:  ";
    for (int p = 0; p < regionLength; ++p) cout << (p % 10) << " ";
    cout << "\n";
 
    for (char base : {'A', 'C', 'G', 'T'}) {
        cout << base << ":    ";
        for (int p = 0; p < regionLength; ++p) {
            cout << freqMatrix.at(base)[p] << " ";
        }
        cout << "\n";
    }
}
 
void reportMotif(const string &consensus, const vector<ExtractRegion> &regions, const ConservationResult &conservation, int totalSequenceCount, int k) {
    cout << "\nMotif Candidate: " << consensus << " with length " << consensus.size() <<  "\n";
    cout << "Total Sequences: " << totalSequenceCount << "\n";
    cout << "Sequences with Motif: " << regions.size() << "\n";
    cout << "Sequences containing the motif: " << regions.size() << " / " << totalSequenceCount << "\n\n";

    cout << "Conservation per Position:\n";
    for (const auto &region : regions) {
        cout << "Sequence " << region.sequence_index << " (Pos " << region.start_position << "): " << region.region << "\n";
    }

    cout << "\nConservation per Position:\n";
    for (size_t i = 0; i < conservation.percentPerPosition.size(); i++) {
        cout << "Position " << i << ": " << conservation.percentPerPosition[i] << "% conserved, " << (conservation.fullyConserved[i] ? "Fully Conserved" : "Not Fully Conserved") << "\n";
    }

    cout << "\nOverall Conservation: " << conservation.overallConservation << "%\n";
}

bool process(const vector<string> &sequences, int k,  int coveragePercent) {
    cout << "Processing with k = " << k << "\n";
    KmerMap frequency, coverage;
    builsKmerStats(sequences, k, frequency, coverage);

    vector<KmerCandidate> candidates = getTopCandidates(frequency, coverage, 5);

    if (candidates.empty()) {
        cout << "No candidates found for k = " << k << "\n";
        return false;
    }

    cout << "Top candidates:\n";

    for (const auto &candidate : candidates) {
        cout << candidate.kmer << " - Total Frequency: " << candidate.total_frequency << ", Sequence Coverage: " << candidate.sequence_coverage << "\n";
    }

    const KmerCandidate &bestCandidate = candidates[0];
    int requiredCoverage = (int)ceil(sequences.size() * coveragePercent / 100.0);

    if (bestCandidate.sequence_coverage < requiredCoverage) {
        cout << "Best candidate coverage (" << bestCandidate.sequence_coverage << ") is less than required (" << requiredCoverage << ") for k = " << k << "\n";
        return false;
    }

    cout << "Selected candidate: " << bestCandidate.kmer << " with coverage " << bestCandidate.sequence_coverage << "\n";

    vector<vector<int>> positions = findKmerPositions(sequences, bestCandidate.kmer);
    int region_size = k;
    vector<ExtractRegion> regions = extractRegions(sequences, positions, k);

    if (regions.empty()) {
        cout << "No regions extracted for candidate " << bestCandidate.kmer << "\n";
        return false;
    }

    if (!verifyAlignment(regions)) {
        cout << "Extracted regions are not of the same length for candidate " << bestCandidate.kmer << "\n";
        return false;
    }

    map<char, vector<int>> frequencyMatrix = buildFrequencyMatrix(regions);
    printFrequencyMatrix(frequencyMatrix, region_size);

    string consensus = buildConsensus(frequencyMatrix, region_size);

    ConservationResult conservation = evaluateConservation(frequencyMatrix, consensus, (int)sequences.size(), region_size);
    reportMotif(consensus, regions, conservation, (int)sequences.size(), k);

    return true;
}

void runDiscovery(const vector<string> &sequences) {
    int coveragePercent = 70;
    cout << "Starting motif discovery...\n";

    for (int k : {7, 8, 9}) {
        process(sequences, k, coveragePercent);
    }
}

int main(int argc, char *argv[]) {
    cout << "Motif Discovery Program\n";
    cout << "=======================\n";
    cout << "Select an option:" << endl;
    cout << "1. Use test sequences" << endl;
    cout << "2. Use annex sequences" << endl;
    cout << "Enter your choice: ";
    int choice;
    cin >> choice;

    vector<string> sequences;

    if (choice == 1) {
        sequences = sequences_test();
    }
    else if (choice == 2) {
        sequences = getAnnexSequences();
    }
    else {
        cout << "Invalid choice. Exiting." << endl;
        return 1;
    }

    runDiscovery(sequences);
    return 0;
}