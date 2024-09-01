#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>

using namespace std;

void tarjanDFS(int u, vector<vector<int>>& graph, vector<int>& ids, vector<int>& low, stack<int>& s, vector<bool>& onStack, int& id, vector<vector<int>>& sccs) {
    ids[u] = low[u] = id++;
    s.push(u);
    onStack[u] = true;
    
    for (int v : graph[u]) {
        if (ids[v] == -1) {
            tarjanDFS(v, graph, ids, low, s, onStack, id, sccs);
        }
        if (onStack[v]) {
            low[u] = min(low[u], low[v]);
        }
    }
    
    if (ids[u] == low[u]) {
        vector<int> scc;
        while (true) {
            int node = s.top();
            s.pop();
            onStack[node] = false;
            scc.push_back(node);
            if (node == u) break;
        }
        sccs.push_back(scc);
    }
}

bool tarjanAndRemoveCycles(vector<vector<int>>& graph) {
    int n = graph.size();
    vector<int> ids(n, -1), low(n, -1);
    vector<vector<int>> sccs;
    vector<bool> onStack(n, false);
    stack<int> s;
    int id = 0;
    
    // Tarjan's algorithm to find SCCs
    for (int i = 0; i < n; i++) {
        if (ids[i] == -1) {
            tarjanDFS(i, graph, ids, low, s, onStack, id, sccs);
        }
    }
    
    // Check for cycles in SCCs and remove them
    bool hasCycle = false;
    for (const auto& scc : sccs) {
        if (scc.size() > 1) {
            hasCycle = true;
            cout << "Removing edges in SCC: ";
            for (int u : scc) {
                cout << u << " ";
                graph[u].clear(); // 移除强连通分量内的所有边
            }
            cout << endl;
        }
    }
    
    return !hasCycle;
}

int main() {
    int n = 4;
    vector<vector<int>> graph(n);
    
    // 构建图
    graph[0].push_back(1);
    graph[1].push_back(2);
    graph[2].push_back(3);
    graph[3].push_back(1); // 环

    if (!tarjanAndRemoveCycles(graph)) {
        cout << "Cycle detected and removed using Tarjan's SCC." << endl;
    } else {
        cout << "No cycle detected." << endl;
    }
    
    return 0;
}
