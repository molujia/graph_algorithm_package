#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>

using namespace std;

void dfs(const unordered_map<char, vector<char>>& graph, char node, unordered_set<char>& visited) {
    visited.insert(node);
    cout << node << " ";
    
    for (char neighbor : graph.at(node)) {
        if (visited.find(neighbor) == visited.end()) {
            dfs(graph, neighbor, visited);
        }
    }
}

int main() {
    unordered_map<char, vector<char>> graph = {
        {'A', {'B', 'C'}},
        {'B', {'A', 'D', 'E'}},
        {'C', {'A', 'F', 'G'}},
        {'D', {'B'}},
        {'E', {'B', 'H'}},
        {'F', {'C'}},
        {'G', {'C'}},
        {'H', {'E'}}
    };
    
    unordered_set<char> visited;
    char start_node = 'A';
    cout << "DFS traversal starting from node " << start_node << ": ";
    dfs(graph, start_node, visited);
    
    cout << endl;
    return 0;
}
