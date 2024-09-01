#include <iostream>
#include <vector>
#include <algorithm>  // for std::remove
#include <queue>

using namespace std;

class Graph {
public:
    int V;
    vector<vector<int>> adj;

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
    adj[u].push_back(v);
}

void Graph::removeEdge(int u, int v) {
    auto it = std::remove(adj[u].begin(), adj[u].end(), v);
    adj[u].erase(it, adj[u].end());  //erase需要传入一个迭代器？
}

bool Graph::detectAndRemoveCycle() {
    vector<int> in_degree(V, 0);

    for (int u = 0; u < V; u++) {
        for (int v : adj[u]) {
            in_degree[v]++;
        }
    }

    // Topological sorting
    queue<int> q;
    for (int i = 0; i < V; i++) {
        if (in_degree[i] == 0) {
            q.push(i);
        }
    }

    int count = 0;
    vector<int> top_order;

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        top_order.push_back(u);

        for (int v : adj[u]) {
            if (--in_degree[v] == 0) {
                q.push(v);
            }
        }

        count++;
    }

    if (count != V) {
        cout << "Cycle detected! Removing a cycle edge..." << endl;
        for (int u = 0; u < V; u++) {
            for (int v : adj[u]) {
                if (in_degree[v] > 0) {
                    removeEdge(u, v);
                    cout << "Removed edge: " << u << " -> " << v << endl;
                    return true;
                }
            }
        }
    } else {
        cout << "No cycle detected." << endl;
        return false;
    }

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
