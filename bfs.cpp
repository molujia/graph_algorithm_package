#include <iostream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <unordered_map>

using namespace std;

void bfs(const unordered_map<char, vector<char>>& graph, char start) {
    unordered_set<char> visited;
    queue<char> q;
    
    q.push(start);
    visited.insert(start);
    
    while (!q.empty()) {
        char node = q.front();
        q.pop();
        
        cout << node << " ";
        
        for (char neighbor : graph.at(node)) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                q.push(neighbor);
            }
        }
    }
    cout << endl;
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
    
    char start_node = 'A';
    cout << "BFS traversal starting from node " << start_node << ": ";
    bfs(graph, start_node);
    
    return 0;
}
