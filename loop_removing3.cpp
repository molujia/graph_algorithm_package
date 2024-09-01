#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class Graph {
    int V;
    vector<pair<int, int>> edges;
    vector<vector<int>> adj;

    void DFS(int v, vector<bool> &visited);
    bool isConnected();

public:
    Graph(int V);
    void addEdge(int u, int v);
    void removeEdge(int u, int v);
    bool detectAndRemoveCycle();
};

Graph::Graph(int V) {
    this->V = V;
    adj.resize(V);
}

void Graph::addEdge(int u, int v) {
    edges.push_back({u, v});
    adj[u].push_back(v);
}

void Graph::removeEdge(int u, int v) {
    adj[u].erase(remove(adj[u].begin(), adj[u].end(), v), adj[u].end());
}

void Graph::DFS(int v, vector<bool> &visited) {
    visited[v] = true;
    for (int u : adj[v]) {
        if (!visited[u]) {
            DFS(u, visited);
        }
    }
}

bool Graph::isConnected() {
    for (auto &[u, v] : edges) {
        removeEdge(u, v);

        vector<bool> visited(V, false);
        DFS(u, visited);

        addEdge(u, v);

        if (!visited[v]) {
            return false;
        }
    }

    return true;
}

bool Graph::detectAndRemoveCycle() {
    sort(edges.rbegin(), edges.rend());

    for (auto &[u, v] : edges) {
        removeEdge(u, v);
        cout << "Trying to remove edge: " << u << " -> " << v << endl;

        if (!isConnected()) {
            cout << "Cycle detected! Removed edge: " << u << " -> " << v << endl;
            return true;
        } else {
            addEdge(u, v);
        }
    }

    cout << "No cycle detected." << endl;
    return false;
}

int main() {
    Graph g(6);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 0);
    g.addEdge(2, 3);
    g.addEdge(3, 4);
    g.addEdge(4, 5);

    g.detectAndRemoveCycle();

    return 0;
}
