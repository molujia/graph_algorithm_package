#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>  // for iota
#include <tuple>
#include <queue>
#include <limits>

using namespace std;

class MSTAlgorithms {
private:
    int numVertices;
    vector<tuple<int, int, int>> edges;  // {weight, u, v}
    vector<vector<pair<int, int>>> adjList;  // {neighbor, weight}
    bool useKruskal;  // 自动选择算法

    // 并查集查找函数 (用于 Kruskal)
    int find(vector<int>& parent, int i) {
        if (parent[i] != i)
            parent[i] = find(parent, parent[i]);
        return parent[i];
    }

    // 并查集合并函数 (用于 Kruskal)
    void unionSets(vector<int>& parent, vector<int>& rank, int x, int y) {
        int rootX = find(parent, x);
        int rootY = find(parent, y);

        if (rank[rootX] < rank[rootY])
            parent[rootX] = rootY;
        else if (rank[rootX] > rank[rootY])
            parent[rootY] = rootX;
        else {
            parent[rootY] = rootX;
            rank[rootX]++;
        }
    }

    // Kruskal算法
    vector<tuple<int, int, int>> kruskalMST() {
        vector<int> parent(numVertices);
        vector<int> rank(numVertices, 0);
        iota(parent.begin(), parent.end(), 0);  // 初始化并查集

        vector<tuple<int, int, int>> mst;
        sort(edges.begin(), edges.end());  // 按权重排序

        for (const auto& [weight, u, v] : edges) {
            int setU = find(parent, u);
            int setV = find(parent, v);

            if (setU != setV) {  // 如果u和v属于不同集合
                mst.push_back({weight, u, v});
                unionSets(parent, rank, setU, setV);
            }
        }

        return mst;
    }

    // Prim算法
    vector<tuple<int, int, int>> primMST() {
        vector<bool> inMST(numVertices, false);
        vector<int> key(numVertices, numeric_limits<int>::max());
        vector<int> parent(numVertices, -1);
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;

        key[0] = 0;
        pq.push({0, 0});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            if (inMST[u]) continue;

            inMST[u] = true;

            for (const auto& [v, weight] : adjList[u]) {
                if (!inMST[v] && weight < key[v]) {
                    key[v] = weight;
                    pq.push({key[v], v});
                    parent[v] = u;
                }
            }
        }

        vector<tuple<int, int, int>> mst;
        for (int v = 1; v < numVertices; ++v) {
            if (parent[v] != -1) { // 确保有效边
                mst.push_back({key[v], parent[v], v});
            }
        }

        return mst;
    }

    // 计算MST总权重
    int calculateTotalWeight(const vector<tuple<int, int, int>>& mst) {
        int totalWeight = 0;
        for (const auto& [weight, u, v] : mst) {
            totalWeight += weight;
        }
        return totalWeight;
    }

public:
    // 构造函数：自动选择算法
    MSTAlgorithms(int vertices, const vector<tuple<int, int, int>>& edgeList, const vector<vector<pair<int, int>>>& adjacencyList) {
        numVertices = vertices;
        edges = edgeList;
        adjList = adjacencyList;
        useKruskal = !edges.empty() && adjList.empty();  // 如果有边列表且无邻接表，使用Kruskal；否则使用Prim
    }

    // 自动选择Kruskal或Prim算法计算最小生成树
    vector<tuple<int, int, int>> findMST() {
        if (useKruskal) {
            return kruskalMST();
        } else {
            return primMST();
        }
    }
};

// 测试函数
int main() {
    vector<tuple<int, int, int>> edgeList = {
        {1, 0, 1}, {2, 0, 2}, {3, 1, 2}, {4, 1, 3}, {5, 2, 3}, {6, 2, 4}, {7, 3, 4}
    };
    vector<vector<pair<int, int>>> adjList = {
        {{1, 1}, {2, 2}},
        {{0, 1}, {2, 3}, {3, 4}},
        {{0, 2}, {1, 3}, {4, 6}},
        {{1, 4}, {4, 7}, {2, 5}},
        {{2, 6}, {3, 7}}
    };

    MSTAlgorithms mstAlgorithms(5, edgeList, adjList);

    // 自动选择算法
    vector<tuple<int, int, int>> mst = mstAlgorithms.findMST();
    cout << "Edges in MST:\n";
    for (const auto& [weight, u, v] : mst) {
        cout << u << " - " << v << " (weight: " << weight << ")\n";
    }

    return 0;
}
