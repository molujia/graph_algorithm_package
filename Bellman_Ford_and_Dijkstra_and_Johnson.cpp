#include <iostream>
#include <vector>
#include <queue>
#include <climits>
#include <tuple>
#include <functional>



/*

B算法 和 D算法 都是求单源最短路径的问题
B算法时间复杂度O（mn） 适合任意图
D算法时间复杂度为O（mlongm） 适合非负权图


J算法用于求每对节点之间的最短路
时间复杂度为O（nmlogm）


三种算法封装到Grapg类中
*/


using namespace std;

class Graph {
public:
    Graph(int vertices);

    // 添加边
    void addEdge(int u, int v, int weight);

    // Bellman-Ford 算法
    vector<int> bellmanFord(int src);

    // Dijkstra 算法
    vector<int> dijkstra(int src);

    // Johnson 算法
    vector<vector<int>> johnson();

private:
    int V;
    vector<tuple<int, int, int>> edges; // 用于 Bellman-Ford (u, v, weight)
    vector<vector<pair<int, int>>> adj; // 用于 Dijkstra (node, weight)
};

// 构造函数
Graph::Graph(int vertices) : V(vertices), adj(vertices) {}

// 添加边
void Graph::addEdge(int u, int v, int weight) {
    edges.emplace_back(u, v, weight);  // 用于 Bellman-Ford
    adj[u].emplace_back(v, weight);    // 用于 Dijkstra
}

// Bellman-Ford 算法
vector<int> Graph::bellmanFord(int src) {
    vector<int> dist(V, INT_MAX);
    dist[src] = 0;

    // 放松所有边 |V| - 1 次
    for (int i = 1; i < V; ++i) {
        for (const auto& edge : edges) {
            int u, v, weight;
            tie(u, v, weight) = edge;

            if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
            }
        }
    }

    // 检测负权回路
    for (const auto& edge : edges) {
        int u, v, weight;
        tie(u, v, weight) = edge;

        if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
            throw runtime_error("Graph contains negative weight cycle");
        }
    }

    return dist;
}

// Dijkstra 算法
vector<int> Graph::dijkstra(int src) {
    vector<int> dist(V, INT_MAX);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> minHeap;

    dist[src] = 0;
    minHeap.emplace(0, src); // (distance, vertex)

    while (!minHeap.empty()) {
        int u = minHeap.top().second;
        int u_dist = minHeap.top().first;
        minHeap.pop();

        // 如果当前距离大于已知最短距离，则跳过
        if (u_dist > dist[u]) continue;

        for (const auto& edge : adj[u]) {
            int v = edge.first;
            int weight = edge.second;

            if (dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                minHeap.emplace(dist[v], v);
            }
        }
    }

    return dist;
}

// Johnson 算法
vector<vector<int>> Graph::johnson() {
    // 添加一个虚拟源点
    Graph g(V + 1);

    // 复制原图的边
    for (const auto& edge : edges) {
        int u, v, weight;
        tie(u, v, weight) = edge;
        g.addEdge(u, v, weight);
    }

    // 从虚拟源点添加边
    for (int i = 0; i < V; ++i) {
        g.addEdge(V, i, 0);
    }

    // 使用 Bellman-Ford 计算从虚拟源点到所有顶点的最短路径
    vector<int> h = g.bellmanFord(V);

    // 重新调整图的权重
    vector<vector<pair<int, int>>> newAdj(V);
    for (int u = 0; u < V; ++u) {
        for (const auto& edge : adj[u]) {
            int v = edge.first;
            int weight = edge.second;
            int newWeight = weight + h[u] - h[v];
            newAdj[u].emplace_back(v, newWeight);
        }
    }

    vector<vector<int>> distances(V, vector<int>(V, INT_MAX));

    // 对每个顶点运行 Dijkstra 算法
    for (int i = 0; i < V; ++i) {
        Graph tempGraph(V);
        tempGraph.adj = newAdj;
        vector<int> dist = tempGraph.dijkstra(i);
        for (int j = 0; j < V; ++j) {
            if (dist[j] != INT_MAX) {
                distances[i][j] = dist[j] - h[i] + h[j];
            }
        }
    }

    return distances;
}

int main() {


    /*
    example
    
    */
    int V = 5; // 顶点数
    Graph g(V);

    // 添加边
    g.addEdge(0, 1, -1);
    g.addEdge(0, 2, 4);
    g.addEdge(1, 2, 3);
    g.addEdge(1, 3, 2);
    g.addEdge(1, 4, 2);
    g.addEdge(3, 2, 5);
    g.addEdge(3, 1, 1);
    g.addEdge(4, 3, -3);

    try {
        vector<vector<int>> distances = g.johnson();

        cout << "Johnson's Algorithm: Distance Matrix" << endl;
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                cout << (distances[i][j] == INT_MAX ? "INF" : to_string(distances[i][j])) << "\t";
            }
            cout << endl;
        }
    } catch (const runtime_error& e) {
        cout << e.what() << endl;
    }

    return 0;
}
