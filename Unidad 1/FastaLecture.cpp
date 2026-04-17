#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#define ll long long

struct FastaRecordStructure {
    std::string identifier;
    std::string sequence;
};

struct FastaMetrics {
    ll maximum;
    ll minimum;
    double average;
};

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

FastaMetrics calculate_metrics(const std::vector<FastaRecordStructure>& records) {
    FastaMetrics metrics_fasta = {0, LLONG_MAX, 0.0};
    ll total_length = 0;

    for (const auto &record : records) {
        ll sequence_length = record.sequence.length();
        total_length += sequence_length;

        if (sequence_length > metrics_fasta.maximum) {
            metrics_fasta.maximum = sequence_length;
        }

        if (sequence_length < metrics_fasta.minimum) {
            metrics_fasta.minimum = sequence_length;
        }
    }

    metrics_fasta.average = records.empty() ? 0.0 : static_cast<double>(total_length) / records.size();

    return metrics_fasta;
}

double count_proportion_gc(const std::string &sequence) {
    ll gc_count = 0;

    for (char gc : sequence) {
        if (gc == 'G' || gc == 'C' || gc == 'g' || gc == 'c') {
            gc_count++;
        }
    }

    double proportion_gc = sequence.empty() ? 0.0 : (static_cast<double>(gc_count) / sequence.length()) * 100.0;

    return proportion_gc;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file.fasta>" << std::endl;
        return 1;
    }

    std::vector<FastaRecordStructure> records;
    std::string file_fasta = argv[1];
    read_file_fasta(file_fasta, records);

    FastaMetrics metrics = calculate_metrics(records);

    const int LIMIT = 100;

    for (const auto &record : records) {
        std::cout << record.identifier << "\n";
        if (record.sequence.length() > LIMIT) {
            std::cout << record.sequence.substr(0, LIMIT) << "...\n";
        } else {
            std::cout << record.sequence << "\n";
        }
    }

    std::cout << "------------------------------\n";
    std::cout << "Maximum sequence length: " << metrics.maximum << "\n";
    std::cout << "Minimum sequence length: " << metrics.minimum << "\n";
    std::cout << "Average sequence length: " << metrics.average << "\n";
    std::cout << "------------------------------\n";
    std::cout << "Proportion GC per sequence:\n";

    for (const auto &record : records) {
        std::cout << record.identifier << ": " << std::fixed << std::setprecision(2) << count_proportion_gc(record.sequence) << "%\n";
    }

    return 0;
}