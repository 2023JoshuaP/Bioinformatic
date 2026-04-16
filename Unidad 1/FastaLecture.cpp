#include <iostream>
#include <fstream>
#include <string>
#include <vector>

void read_file_fasta(const std::string &file_fasta) {
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
                std::cout << "Header: " << header << std::endl;
                std::cout << "Sequence: " << sequence << std::endl;
                sequence.clear();
            }
            header = line.substr(1);
        }
        else {
            sequence += line;
        }
    }

    if (!header.empty()) {
        std::cout << "Header: " << header << std::endl;
        std::cout << "Sequence: " << sequence << std::endl;
    }

    file.close();
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file.fasta>" << std::endl;
        return 1;
    }

    std::string file_fasta = argv[1];
    read_file_fasta(file_fasta);

    return 0;
}