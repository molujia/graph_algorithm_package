#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>

using namespace std;

const int NIL = -1;

void BCCUtilMatrix(int u, int disc[], int low[], int parent[], stack<pair<int, int>> &st, vector<vector<int>> &adjMatrix, vector<vector<pair<int, int>>> &bcc) {
    static int time = 0;
    disc[u] = low[u] = ++time;
    int children = 0;

    for (int v = 0; v < adjMatrix.size(); ++v) {
        if (adjMatrix[u][v] != 0) { // 检查是否存在边
            if (disc[v] == -1) {
                children++;
                parent[v] = u;
                st.push({u, v});
                BCCUtilMatrix(v, disc, low, parent, st, adjMatrix, bcc);

                low[u] = min(low[u], low[v]);

                if ((parent[u] == NIL && children > 1) || (parent[u] != NIL && low[v] >= disc[u])) {
                    vector<pair<int, int>> component;
                    while (st.top() != make_pair(u, v)) {
                        component.push_back(st.top());
                        st.pop();
                    }
                    component.push_back(st.top());
                    st.pop();
                    bcc.push_back(component);
                }
            } else if (v != parent[u] && disc[v] < disc[u]) {
                low[u] = min(low[u], disc[v]);
                st.push({u, v});
            }
        }
    }
}

vector<vector<pair<int, int>>> findBCCMatrix(vector<vector<int>> &adjMatrix) {
    int V = adjMatrix.size();
    int *disc = new int[V];
    int *low = new int[V];
    int *parent = new int[V];
    stack<pair<int, int>> st;
    vector<vector<pair<int, int>>> bcc;

    fill(disc, disc + V, NIL);
    fill(low, low + V, NIL);
    fill(parent, parent + V, NIL);

    for (int i = 0; i < V; ++i) {
        if (disc[i] == NIL) {
            BCCUtilMatrix(i, disc, low, parent, st, adjMatrix, bcc);
        }

        while (!st.empty()) {
            vector<pair<int, int>> component;
            component.push_back(st.top());
            st.pop();
            bcc.push_back(component);
        }
    }

    delete[] disc;
    delete[] low;
    delete[] parent;

    return bcc;
}

void BCCUtilList(int u, int disc[], int low[], int parent[], stack<pair<int, int>> &st, vector<vector<int>> &adjList, vector<vector<pair<int, int>>> &bcc) {
    static int time = 0;
    disc[u] = low[u] = ++time;
    int children = 0;

    for (int v : adjList[u]) {
        if (disc[v] == -1) {
            children++;
            parent[v] = u;
            st.push({u, v});
            BCCUtilList(v, disc, low, parent, st, adjList, bcc);

            low[u] = min(low[u], low[v]);

            if ((parent[u] == NIL && children > 1) || (parent[u] != NIL && low[v] >= disc[u])) {
                vector<pair<int, int>> component;
                while (st.top() != make_pair(u, v)) {
                    component.push_back(st.top());
                    st.pop();
                }
                component.push_back(st.top());
                st.pop();
                bcc.push_back(component);
            }
        } else if (v != parent[u] && disc[v] < disc[u]) {
            low[u] = min(low[u], disc[v]);
            st.push({u, v});
        }
    }
}

vector<vector<pair<int, int>>> findBCCList(vector<vector<int>> &adjList) {
    int V = adjList.size();
    int *disc = new int[V];
    int *low = new int[V];
    int *parent = new int[V];
    stack<pair<int, int>> st;
    vector<vector<pair<int, int>>> bcc;

    fill(disc, disc + V, NIL);
    fill(low, low + V, NIL);
    fill(parent, parent + V, NIL);

    for (int i = 0; i < V; ++i) {
        if (disc[i] == NIL) {
            BCCUtilList(i, disc, low, parent, st, adjList, bcc);
        }

        while (!st.empty()) {
            vector<pair<int, int>> component;
            component.push_back(st.top());
            st.pop();
            bcc.push_back(component);
        }
    }

    delete[] disc;
    delete[] low;
    delete[] parent;

    return bcc;
}

int main() {
    vector<vector<int>> adjMatrix = {
        {0, 1, 0, 1},
        {1, 0, 1, 0},
        {0, 1, 0, 1},
        {1, 0, 1, 0}
    };

    vector<vector<int>> adjList = {
        {1, 3},
        {0, 2},
        {1, 3},
        {0, 2}
    };

    cout << "Biconnected components in the adjacency matrix graph:" << endl;
    vector<vector<pair<int, int>>> bccMatrix = findBCCMatrix(adjMatrix);
    for (const auto& component : bccMatrix) {
        for (const auto& edge : component) {
            cout << edge.first << "-" << edge.second << " ";
        }
        cout << endl;
    }

    cout << "Biconnected components in the adjacency list graph:" << endl;
    vector<vector<pair<int, int>>> bccList = findBCCList(adjList);
    for (const auto& component : bccList) {
        for (const auto& edge : component) {
            cout << edge.first << "-" << edge.second << " ";
        }
        cout << endl;
    }

    return 0;
}
