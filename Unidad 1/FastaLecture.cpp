#include <iostream>
#include <fstream>
#include <string>
#include <vector>

struct FastaRecordStructure {
    std::string identifier;
    std::string sequence;
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

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file.fasta>" << std::endl;
        return 1;
    }

    std::vector<FastaRecordStructure> records;
    std::string file_fasta = argv[1];
    read_file_fasta(file_fasta, records);

    const int LIMIT = 100;

    for (const auto &record : records) {
        std::cout << record.identifier << "\n";
        if (record.sequence.length() > LIMIT) {
            std::cout << record.sequence.substr(0, LIMIT) << "...\n";
        } else {
            std::cout << record.sequence << "\n";
        }
    }

    return 0;
}