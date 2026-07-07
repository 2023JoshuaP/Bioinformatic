#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <deque>
#include <fstream>

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

/*
Assembly of DNA sequences using De Bruijn graph
*/
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

/*
Compute the in-degree of each node in the graph. The in-degree of a node is the number of edges directed into it.
*/
unordered_map<string, int> computeInDegree(const Graph &graph) {
    unordered_map<string, int> inDegree;

    for (const auto &[node, edges] : graph) {
        if (inDegree.find(node) == inDegree.end()) {
            inDegree[node] = 0;
        }
        for (const auto &edge : edges) {
            inDegree[edge]++;
        }
    }

    return inDegree;
}

/*
Compute the first node for the Eulerian path.
*/
string firstNode(const Graph &graph,  const unordered_map<string, int> &inDegree) {
    string start = "";

    for (const auto &[node, edges] : graph) {
        int outDeg = edges.size();
        int inDeg = inDegree.count(node) ? inDegree.at(node) : 0;

        if (outDeg - inDeg == 1) {
            return node;
        }

        if (start.empty() && outDeg > 0) {
            start = node;
        }
    }

    return start;
}

/*
Find an Eulerian path in the De Bruijn graph using Hierholzer's algorithm. The algorithm constructs the path by traversing the graph and backtracking when necessary.
*/
vector<string> hierholzer(Graph graph, const string &start) {
    deque<string> path;
    vector<string> stack;
    stack.push_back(start);

    while (!stack.empty()) {
        string &current = stack.back();

        if (graph.count(current) && !graph.at(current).empty()) {
            string next = graph[current].back();
            graph[current].pop_back();
            stack.push_back(next);
        }
        else {
            path.push_front(current);
            stack.pop_back();
        }
    }

    return vector<string>(path.begin(), path.end());
}

/*
Reconstruct the DNA sequence from the Eulerian path.
*/
string reconstructSequence(const vector<string> &path) {
    if (path.empty()) {
        return "";
    }

    string sequence = path[0];
    for (size_t i = 1; i < path.size(); ++i) {
        sequence += path[i].back();
    }

    return sequence;
}

/*
Functions to print the graph, degrees, and path for debugging purposes.
*/
void printGraph(const Graph &graph) {
    for (const auto &[node, edges] : graph) {
        cout << node << " -> ";
        for (const auto &edge : edges) {
            cout << edge << " ";
        }
        cout << endl;
    }
}

void printDegrees(const Graph &graph, const unordered_map<string, int> &inDegree) {
    cout << "Node\tout\tin\tout-in\n";
    unordered_map<string, int> outDegree;
    for (const auto &[node, edges] : graph) {
        outDegree[node] = edges.size();
    }

    unordered_map<string, bool> allNodes;
    for (const auto &[node, _] : outDegree) {
        allNodes[node] = true;
    }
    for (const auto &[node, _] : inDegree) {
        allNodes[node] = true;
    }

    for (const auto &[node, _] : allNodes) {
        int outDeg = outDegree.count(node) ? outDegree.at(node) : 0;
        int inDeg = inDegree.count(node) ? inDegree.at(node) : 0;
        cout << node << "\t" << outDeg << "\t" << inDeg << "\t" << (outDeg - inDeg) << "\n";
    }
}

void printPath(const vector<string> &path) {
    for (const auto &node : path) {
        cout << node << " ";
    }
    cout << endl;
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

    cout << "\nCompute In-Degree:" << endl;
    auto inDegree = computeInDegree(graph);
    printDegrees(graph, inDegree);

    cout << "\nFirst Node for Eulerian Path:" << endl;
    string startNode = firstNode(graph, inDegree);
    cout << "Start Node: " << startNode << endl;

    cout << "\nFind Eulerian Path using Hierholzer's Algorithm:" << endl;
    vector<string> eulerianPath = hierholzer(graph, startNode);
    printPath(eulerianPath);

    cout << "\nReconstruct Sequence from Eulerian Path:" << endl;
    string reconstructedSequence = reconstructSequence(eulerianPath);
    cout << reconstructedSequence << endl;

    return 0;
}