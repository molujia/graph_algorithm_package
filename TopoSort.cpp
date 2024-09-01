#include <iostream>
#include <vector>
#include <queue>
#include <stdexcept>

using namespace std;

class TopologicalSort {
public:
    // 邻接表接口
    vector<int> sortByAdjList(int numVertices, const vector<vector<int>>& adjList) {
        return topologicalSort(numVertices, adjList);
    }

    // 邻接矩阵接口
    vector<int> sortByAdjMatrix(int numVertices, const vector<vector<int>>& adjMatrix) {
        // 将邻接矩阵转换为邻接表
        vector<vector<int>> adjList(numVertices);
        for (int i = 0; i < numVertices; ++i) {
            for (int j = 0; j < numVertices; ++j) {
                if (adjMatrix[i][j] != 0) {
                    adjList[i].push_back(j);
                }
            }
        }
        return topologicalSort(numVertices, adjList);
    }

private:
    // 通用的拓扑排序实现（基于邻接表）
    vector<int> topologicalSort(int numVertices, const vector<vector<int>>& adjList) {
        vector<int> inDegree(numVertices, 0);  // 记录每个节点的入度
        vector<int> topoOrder;                 // 存储拓扑排序的结果
        queue<int> q;                          // 用于处理入度为0的节点

        // 计算每个节点的入度
        for (int i = 0; i < numVertices; ++i) {
            for (int neighbor : adjList[i]) {
                inDegree[neighbor]++;
            }
        }

        // 将所有入度为0的节点入队
        for (int i = 0; i < numVertices; ++i) {
            if (inDegree[i] == 0) {
                q.push(i);
            }
        }

        // 处理队列中的节点
        while (!q.empty()) {
            int node = q.front();
            q.pop();
            topoOrder.push_back(node);

            // 更新相邻节点的入度
            for (int neighbor : adjList[node]) {
                inDegree[neighbor]--;
                if (inDegree[neighbor] == 0) {
                    q.push(neighbor);
                }
            }
        }

        // 如果无法生成完整的拓扑序列，说明图中存在环
        if (topoOrder.size() != numVertices) {
            throw runtime_error("Graph is not a DAG! Cannot perform topological sort.");
        }

        return topoOrder;
    }
};

int main() {
    TopologicalSort topoSort;

    // 示例1：邻接表
    int numVertices1 = 6;
    vector<vector<int>> adjList = {
        {2, 3},  // 0 -> 2, 0 -> 3
        {3, 4},  // 1 -> 3, 1 -> 4
        {},      // 2 has no outgoing edges
        {5},     // 3 -> 5
        {5},     // 4 -> 5
        {}       // 5 has no outgoing edges
    };

    try {
        vector<int> topoOrder1 = topoSort.sortByAdjList(numVertices1, adjList);
        cout << "Topological Sort (Adjacency List): ";
        for (int node : topoOrder1) {
            cout << node << " ";
        }
        cout << endl;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
    }

    // 示例2：邻接矩阵
    int numVertices2 = 6;
    vector<vector<int>> adjMatrix = {
        {0, 0, 1, 1, 0, 0},  // 0 -> 2, 0 -> 3
        {0, 0, 0, 1, 1, 0},  // 1 -> 3, 1 -> 4
        {0, 0, 0, 0, 0, 0},  // 2 has no outgoing edges
        {0, 0, 0, 0, 0, 1},  // 3 -> 5
        {0, 0, 0, 0, 0, 1},  // 4 -> 5
        {0, 0, 0, 0, 0, 0}   // 5 has no outgoing edges
    };

    try {
        vector<int> topoOrder2 = topoSort.sortByAdjMatrix(numVertices2, adjMatrix);
        cout << "Topological Sort (Adjacency Matrix): ";
        for (int node : topoOrder2) {
            cout << node << " ";
        }
        cout << endl;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
    }

    return 0;
}
