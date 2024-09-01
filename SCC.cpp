#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>
#include <queue>

class Graph {

public:
    Graph(int V);
    void addEdge(int v, int w);
    std::vector<std::vector<int>> findSCCs();
    int countSCCs();
    std::vector<std::pair<int, int>> findEdgesToMakeSCC();
    bool isConnected();

private:
    void SCCUtil(int u, std::vector<int>& disc, std::vector<int>& low, std::stack<int>& st, std::vector<bool>& inStack, std::vector<std::vector<int>>& sccs);

    int V; 
    std::vector<std::vector<int>> adj; 
    int time; 
};

Graph::Graph(int V) : V(V), adj(V), time(0) {}

void Graph::addEdge(int v, int w) {
    adj[v].push_back(w);
}

std::vector<std::vector<int>> Graph::findSCCs() {
    std::vector<int> disc(V, -1);
    std::vector<int> low(V, -1);
    std::vector<bool> inStack(V, false);
    std::stack<int> st;
    std::vector<std::vector<int>> sccs;

    for (int i = 0; i < V; ++i) {
        if (disc[i] == -1) {
            SCCUtil(i, disc, low, st, inStack, sccs);
        }
    }

    return sccs;
}

int Graph::countSCCs() {
    return findSCCs().size();
}

std::vector<std::pair<int, int>> Graph::findEdgesToMakeSCC() {
    std::vector<std::vector<int>> sccs = findSCCs();
    if (sccs.size() == 1) {
        return {}; // 图已经是强连通的
    }

    std::vector<int> sccInDegree(sccs.size(), 0);
    std::vector<int> sccOutDegree(sccs.size(), 0);
    std::vector<bool> sccVisited(sccs.size(), false);

    // 计算每个SCC的入度和出度
    for (size_t i = 0; i < sccs.size(); ++i) {
        for (int node : sccs[i]) {
            for (int adjNode : adj[node]) {
                int adjSccIndex = std::find(sccs.begin(), sccs.end(), std::vector<int>{adjNode}) - sccs.begin();
                if (i != adjSccIndex) {
                    sccOutDegree[i]++;
                    sccInDegree[adjSccIndex]++;
                }
            }
        }
    }

    std::vector<std::pair<int, int>> additionalEdges;
    std::queue<int> sources; // 存储没有出边的SCC的索引

    // 找出所有没有出边的SCC
    for (size_t i = 0; i < sccs.size(); ++i) {
        if (sccOutDegree[i] == 0) {
            sources.push(i);
        }
    }

    while (!sources.empty()) {
        int sourceIndex = sources.front();
        sources.pop();

        for (int targetIndex = 0; targetIndex < sccs.size(); ++targetIndex) {
            if (!sccVisited[targetIndex] && sccInDegree[targetIndex] == 0) {
                // 添加从sourceIndex到targetIndex的边
                additionalEdges.emplace_back(sccs[sourceIndex][0], sccs[targetIndex][0]);
                sccVisited[targetIndex] = true;
                sccInDegree[targetIndex]++; // 更新入度

                // 如果targetIndex的所有顶点都已访问，将其出度置为0
                if (std::all_of(sccs[targetIndex].begin(), sccs[targetIndex].end(), [&](int node) {
                    return sccVisited[std::find(sccs.begin(), sccs.end(), std::vector<int>{node}) - sccs.begin()];
                })) {
                    sccOutDegree[targetIndex] = 0;
                }

                break; // 由于每个SCC只有一个代表性顶点，添加一条边即可
            }
        }
    }

    return additionalEdges;
}



void Graph::SCCUtil(int u, std::vector<int>& disc, std::vector<int>& low, std::stack<int>& st, std::vector<bool>& inStack, std::vector<std::vector<int>>& sccs) {
    disc[u] = low[u] = ++time;
    st.push(u);
    inStack[u] = true;

    for (int v : adj[u]) {
        if (disc[v] == -1) {
            SCCUtil(v, disc, low, st, inStack, sccs);
            low[u] = std::min(low[u], low[v]);
        } else if (inStack[v]) {
            low[u] = std::min(low[u], disc[v]);
        }
    }

    if (low[u] == disc[u]) {
        std::vector<int> scc;
        int w;
        do {
            w = st.top();
            st.pop();
            inStack[w] = false;
            scc.push_back(w);
        } while (w != u);
        sccs.push_back(scc);
    }
}

bool Graph::isConnected() {
    
    std::vector<std::vector<int>> sccs = findSCCs();
    if (sccs.size() == 1) {
        return 1; // 图已经是强连通的
    }

    return 0;
}

// Functions for adjacency matrix
std::vector<std::vector<int>> findSCCsFromAdjMatrix(int V, const std::vector<std::vector<int>>& adjMatrix) {
    Graph g(V);
    
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j] != 0) {
                g.addEdge(i, j);
            }
        }
    }
    
    return g.findSCCs();
}

int countSCCsFromAdjMatrix(int V, const std::vector<std::vector<int>>& adjMatrix) {
    Graph g(V);

    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j] != 0) {
                g.addEdge(i, j);
            }
        }
    }

    return g.countSCCs();
}

std::vector<std::pair<int, int>> findEdgesToMakeSCCFromAdjMatrix(int V, const std::vector<std::vector<int>>& adjMatrix) {
    Graph g(V);

    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j] != 0) {
                g.addEdge(i, j);
            }
        }
    }

    return g.findEdgesToMakeSCC();
}

// Functions for adjacency list
std::vector<std::vector<int>> findSCCsFromAdjList(int V, const std::vector<std::vector<int>>& adjList) {
    Graph g(V);
    
    for (int u = 0; u < V; ++u) {
        for (int v : adjList[u]) {
            g.addEdge(u, v);
        }
    }
    
    return g.findSCCs();
}

int countSCCsFromAdjList(int V, const std::vector<std::vector<int>>& adjList) {
    Graph g(V);

    for (int u = 0; u < V; ++u) {
        for (int v : adjList[u]) {
            g.addEdge(u, v);
        }
    }

    return g.countSCCs();
}

std::vector<std::pair<int, int>> findEdgesToMakeSCCFromAdjList(int V, const std::vector<std::vector<int>>& adjList) {
    Graph g(V);

    for (int u = 0; u < V; ++u) {
        for (int v : adjList[u]) {
            g.addEdge(u, v);
        }
    }

    return g.findEdgesToMakeSCC();
}

bool isConnectedFromAdjMatrix(int V, const std::vector<std::vector<int>>& adjMatrix) {
    Graph g(V);

    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j] != 0) {
                g.addEdge(i, j);
            }
        }
    }

    return g.isConnected();
}

bool isConnectedFromAdjList(int V, const std::vector<std::vector<int>>& adjList) {
    Graph g(V);

    for (int u = 0; u < V; ++u) {
        for (int v : adjList[u]) {
            g.addEdge(u, v);
        }
    }

    return g.isConnected();
}

int main() {
    // 示例图的顶点数和边
    int V = 5;
    std::vector<std::vector<int>> adjList = {
        {1, 3}, // 0 连接到 1 和 3
        {2},     // 1 连接到 2
        {3},     // 2 连接到 3
        {4},     // 3 连接到 4
        {}        // 4 没有连接
    };

    // 创建图实例
    Graph g(V);

    // 添加边
    for (int i = 0; i < V; ++i) {
        for (int j : adjList[i]) {
            g.addEdge(i, j);
        }
    }

    std::cout << "Graph is " << (g.isConnected() ? "connected" : "not connected") << std::endl;

    // 测试 findSCCs 函数
    std::vector<std::vector<int>> sccs = g.findSCCs();
    std::cout << "Strongly Connected Components (SCCs):" << std::endl;
    for (const auto& scc : sccs) {
        std::cout << "SCC: ";
        for (int v : scc) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }

    // 测试 countSCCs 函数
    std::cout << "Number of SCCs: " << g.countSCCs() << std::endl;

    // 测试 findEdgesToMakeSCC 函数
    std::vector<std::pair<int, int>> edges = g.findEdgesToMakeSCC();
    std::cout << "Edges to make the graph strongly connected:" << std::endl;
    for (const auto& edge : edges) {
        std::cout << "(" << edge.first << ", " << edge.second << ")" << std::endl;
    }

    // 测试基于邻接矩阵和邻接表的函数
    std::vector<std::vector<int>> adjMatrix(V, std::vector<int>(V, 0));
    for (int i = 0; i < V; ++i) {
        for (int j : adjList[i]) {
            adjMatrix[i][j] = 1;
        }
    }

    auto sccsFromMatrix = findSCCsFromAdjMatrix(V, adjMatrix);
    auto sccsFromList = findSCCsFromAdjList(V, adjList);
    std::cout << "SCCs from adjacency matrix:" << std::endl;
    for (const auto& scc : sccsFromMatrix) {
        std::cout << "SCC: ";
        for (int v : scc) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "SCCs from adjacency list:" << std::endl;
    for (const auto& scc : sccsFromList) {
        std::cout << "SCC: ";
        for (int v : scc) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }

    // 测试 isConnected 函数
    std::cout << "Graph is " << (g.isConnected() ? "connected" : "not connected") << std::endl;

    return 0;
}

