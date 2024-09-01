#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <boost/math/distributions/chi_squared.hpp> // 引入 Boost 的卡方分布功能
#include <unordered_map>  //用于独立性检验
#include <unordered_set>  //优化查找
#include <utility>
//#include <boost/functional/hash.hpp>
#include <algorithm> 

//传入的数据本来就是0~n的，不需要映射，但是目前就这样了（主要是size在c++不好算），后面可以优化
std::vector<std::vector<int>> create_contingency_table(const std::vector<int>& data1, const std::vector<int>& data2) {
    std::unordered_map<int, int> row_mapping;
    std::unordered_map<int, int> col_mapping;

    int row_counter = 0;
    int col_counter = 0;

    for (int value : data1) {
        if (row_mapping.find(value) == row_mapping.end()) {
            row_mapping[value] = row_counter++;
        }
    }
    for (int value : data2) {
        if (col_mapping.find(value) == col_mapping.end()) {
            col_mapping[value] = col_counter++;
        }
    }

    std::vector<std::vector<int>> contingency_table(row_mapping.size(), std::vector<int>(col_mapping.size(), 0));

    for (size_t i = 0; i < data1.size(); ++i) {
        int row = row_mapping[data1[i]];
        int col = col_mapping[data2[i]];
        contingency_table[row][col]++;
    }

    return contingency_table;
}

double chi_square_statistic(const std::vector<std::vector<int>>& observed) {
    int rows = observed.size();
    int cols = observed[0].size();
    std::vector<double> row_sums(rows, 0.0);
    std::vector<double> col_sums(cols, 0.0);
    double total = 0.0;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            row_sums[i] += observed[i][j];
            col_sums[j] += observed[i][j];
            total += observed[i][j];
        }
    }

    double chi_square = 0.0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double expected = (row_sums[i] * col_sums[j]) / total;
            chi_square += std::pow(observed[i][j] - expected, 2) / expected;
        }
    }

    return chi_square;
}

int chi_test(const std::vector<int>& data1, const std::vector<int>& data2, const double alpha) {
    int n = data1.size();
    int m = data2.size();
    if (n != m) {
        std::cout << "Input scale doesn't match.";
        return -1;
    }
    std::vector<std::vector<int>> contingency_table = create_contingency_table(data1, data2);
    double chi_square = chi_square_statistic(contingency_table);

    // 计算自由度，通常是 (行数 - 1) * (列数 - 1)
    int df = (contingency_table.size() - 1) * (contingency_table[0].size() - 1);

    //自由度是0，直接判断为独立
    if (df == 0) {
        return 1;
    }

    // 使用 Boost 计算卡方分布的分位点（临界值）
    boost::math::chi_squared chi_squared_dist(df);
    double critical_value = boost::math::quantile(chi_squared_dist, 1 - alpha);

    if (chi_square > critical_value) {
        return 0; // 不独立
    } else {
        return 1; // 独立
    }
}

double g_squared_statistic(const std::vector<std::vector<int>>& observed) {
    int rows = observed.size();
    int cols = observed[0].size();
    std::vector<double> row_sums(rows, 0.0);
    std::vector<double> col_sums(cols, 0.0);
    double total = 0.0;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            row_sums[i] += observed[i][j];
            col_sums[j] += observed[i][j];
            total += observed[i][j];
        }
    }

    double g_squared = 0.0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double expected = (row_sums[i] * col_sums[j]) / total;
            if (observed[i][j] > 0) {
                g_squared += 2 * observed[i][j] * std::log(observed[i][j] / expected);
            }
        }
    }

    return g_squared;
}

int g_sq_test(const std::vector<int>& data1, const std::vector<int>& data2, const double alpha) {
    int n = data1.size();
    int m = data2.size();
    if (n != m) {
        std::cout << "Input scale doesn't match.";
        return -1;
    }
    std::vector<std::vector<int>> contingency_table = create_contingency_table(data1, data2);
    double g_squared = g_squared_statistic(contingency_table);

    // 计算自由度，通常是 (行数 - 1) * (列数 - 1)
    int df = (contingency_table.size() - 1) * (contingency_table[0].size() - 1);
    
    //自由度是0，直接判断为独立
    if (df == 0) {
        return 1;
    }

    // 使用 Boost 计算卡方分布的分位点（临界值）
    boost::math::chi_squared chi_squared_dist(df);
    double critical_value = boost::math::quantile(chi_squared_dist, 1 - alpha);

    if (g_squared > critical_value) {
        return 0; // 不独立
    } else {
        return 1; // 独立
    }
}

std::unordered_map<int, std::pair<std::vector<int>, std::vector<int>>> classifyRows
    (const std::vector<int> & data_i, const std::vector<int> & data_j, const std::vector<int> & data_k) {
    // 最终返回的哈希表
    std::unordered_map<int, std::pair<std::vector<int>, std::vector<int>>> classifications;
    // 遍历第 k 行数据
    for (size_t col = 0; col < data_k.size(); ++col) {
        int key = data_k[col];

        // 将第 i 行和第 j 行的数据根据第 k 行的值进行分类
        classifications[key].first.push_back(data_i[col]);
        classifications[key].second.push_back(data_j[col]);
    }
    return classifications;
}

bool k_in_sij(int i, int j, int k, const std::unordered_map<int, std::unordered_set<int>>& undirectedEdges) 
{
    std::unordered_map<int, std::unordered_set<int>> modifiedEdges = undirectedEdges;
    modifiedEdges.erase(k);
    for (auto& [node, neighbors] : modifiedEdges) 
    {
        neighbors.erase(k);
    }

    std::queue<int> toVisit;
    std::unordered_set<int> visited;

    toVisit.push(i);
    visited.insert(i);

    while (!toVisit.empty()) {
        int current = toVisit.front();
        toVisit.pop();

        if (current == j) {
            return false;  // k not in S(i, j)
        }

        for (int neighbor : modifiedEdges[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                toVisit.push(neighbor);
            }
        }
    }

    return true;  //k in S(i, j)
}

#define directed true
#define undirected false

// 示例 PC 算法函数, columns其实可以不用, const std::vector<std::string>& columns
// 默认是行主序的,所以在传参前需要把df转置一下，以行的方式传参data
// i.e. data_encoded.transpose().to_numpy().tolist()
std::vector<std::vector<int>> pc_algorithm(std::function<int(const std::vector<int>&, const std::vector<int>&, const double)> independence_test,const std::vector<std::vector<int>>& data) 
{
    //有向边
    std::unordered_map<int, std::unordered_set<int>> directedEdges; 
    //无向边
    std::unordered_map<int, std::unordered_set<int>> undirectedEdges;
    //独立性
    std::unordered_map<int, std::unordered_set<int>> independence;

    int numRows = data.size();
    float alpha = 0.05;

    // std::cout<<data[0][0]<<std::endl;
    // std::cout<<(independence_test(data[0],data[1],alpha)==1)<<std::endl;

    //先判断独立性，独立则无需继续判断（N^2），不独立还要继续深入（N^3）
    //启发式剪枝：我们的假设是“只传一跳”，也就是说条件集最大为1，如果要不断的扩大条件集，还需要修改算法
    for(int i=0;i<numRows;i++){
        for(int j=i+1;j<numRows;j++){
            //独立
            if(independence_test(data[i],data[j],alpha)==1)
            {
                independence[i].insert(j);  // i 和 j 独立
            }
            //不独立
            else
            {
                bool outer_flag = false; //是否条件独立
                //开始条件独立性判断
                for(int k=0;k<numRows;k++)
                {
                    if((k==i)||(k==j))continue;
                    std::unordered_map<int, std::pair<std::vector<int>, std::vector<int>>> classification = classifyRows(data[i],data[j],data[k]);
                    bool inner_flag = true; //在k下是否条件独立
                    for (const auto& [key, value] : classification) 
                    {   
                        //必须全独立，才算独立
                        if(independence_test(value.first,value.second,alpha)!=1)  //某一case不独立
                        {
                            inner_flag=false;
                            break;
                        }
                    }
                    if(inner_flag)  //在k下条件独立
                    {
                        outer_flag = inner_flag;
                        break;
                    }
                }
                if(outer_flag == false) //不条件独立
                {
                    undirectedEdges[i].insert(j);  //无向边默认只有ij，没有有ji，j一定比i大，如此实现优化
                    //undirectedEdges[j].insert(i);
                }
            }
        }
    }

    // std::cout<<"测试输出！！"<<std::endl;
    // std::unordered_map<int, std::unordered_set<int>> seenEdges;
    // for (const auto& [i, neighbors] : undirectedEdges) 
    // {
    //     for (const auto& j : neighbors) 
    //     {
    //         if (seenEdges[i].find(j) != seenEdges[i].end()) 
    //         {
    //             std::cout<<"有重复的无向边！！"<<std::endl;
    //             break;
    //         } 
    //         else 
    //         {
    //             seenEdges[i].insert(j);
    //             // std::cout<<i<<" "<<j<<" ";
    //         }
    //     }
    // }

    //通过j=i+1这样的声明来不重复的遍历三元组
    //在这个循环中，全部的边在遍历到时都是无向的
    /*
    使用 std::unordered_map 时，如果键 i 不存在于 undirectedEdges 中，
    调用 undirectedEdges[i] 会隐式地插入一个默认值（即一个空的 std::unordered_set<int>）并返回该默认值
    因此，undirectedEdges[i].find(k) 是安全的，不会导致运行时错误或异常。
    */
    for(int i = 0; i < numRows; i++) {
        for(int j = i + 1; j < numRows; j++) {

            //如果两个节点完全独立，按理说它们应该是不连通的，但是由于我们的算法对独立的判断不可能完美，所以它们还是有可能连同
            //if(independence[i].find(j) != independence[i].end()) continue;  //完全独立，理论上不可能有边，但是还是掩掉，避免多边

            for(int k = 0; k < numRows; k++) {
                if(k == i || k == j) continue;
                if(k_in_sij(i,j,k,undirectedEdges))continue; //k在S(i,j)内

                //符合条件，注意需要判断无向边的内容，比如(1,0)需要修正成(0,1)，但是对应的有向边还是(1,0)
                auto [edge1_head, edge1_tail] = std::minmax(i, k);
                auto [edge2_head, edge2_tail] = std::minmax(j, k);

                if((undirectedEdges[edge1_head].find(edge1_tail) != undirectedEdges[edge1_head].end()) && (undirectedEdges[edge2_head].find(edge2_tail) != undirectedEdges[edge2_head].end())
                && (directedEdges[edge1_head].find(edge1_tail) == directedEdges[edge1_head].end()) && (directedEdges[edge2_head].find(edge2_tail) == directedEdges[edge2_head].end()))
                {
                    undirectedEdges[edge1_head].erase(edge1_tail);  // 删除旧的无向边
                    undirectedEdges[edge2_head].erase(edge2_tail);
                    directedEdges[i].insert(k);   // 插入新的有向边
                    directedEdges[j].insert(k);
                }
            }
        }
    }

    //不能直接遍历这个哈希edges，如果要遍历，删改必须在它外部进行
    //遍历它的副本，然后在本体上删改，避免错误
    std::unordered_map<int, std::unordered_set<int>> undirectedEdges_copy = undirectedEdges;
    // std::unordered_map<int, std::unordered_set<int>> directedEdges_copy = directedEdges;

    //R1规则
    //遍历无向边
    //我不太确定R1规则是否避免了冲突，保险起见，我们对有向边的每个端点都执行一次
    // for (const auto& [i, neighbors] : undirectedEdges_copy) 
    // {
    //     for (const auto& j : neighbors) 
    //     {
    //         //情况一
    //         bool flag = false;
    //         for (const auto& directed_neighbor : directedEdges) 
    //         {
    //             if((directed_neighbor.second.find(i) != directed_neighbor.second.end()) && (directed_neighbor.second.find(j) == directed_neighbor.second.end())) 
    //             {
    //                 flag = true;
    //                 break;
    //             }
    //         }
    //         if(flag) 
    //         {
    //             undirectedEdges[i].erase(j);
    //             directedEdges[i].insert(j);
    //             continue;
    //         } 

    //         //情况二
    //         flag = false;
    //         for (const auto& directed_neighbor : directedEdges) 
    //         {
    //             if((directed_neighbor.second.find(j) != directed_neighbor.second.end()) && (directed_neighbor.second.find(i) == directed_neighbor.second.end())) 
    //             {
    //                 flag = true;
    //                 break;
    //             }
    //         }
    //         if(flag) 
    //         {
    //             undirectedEdges[j].erase(i);
    //             directedEdges[j].insert(i);
    //         } 
    //     }
    // }

    //处理 R2/R3 规则（这是暴力版本，下面有优化版本）
    // for(int i = 0; i < numRows; ++i) {
    //     for(int j = 0; j < numRows; ++j) {
    //         if(i == j) continue;
    //         if(directedEdges[i].find(j) != directedEdges[i].end()) continue; // 有向边已存在
    //         if(undirectedEdges[i].find(j) == undirectedEdges[i].end()) continue; // 无向边不存在

    //         for(int k = 0; k < numRows; ++k) {
    //             if(k == i || k == j) continue;
    //             if(directedEdges[i].find(k) != directedEdges[i].end() && directedEdges[k].find(j) != directedEdges[k].end()) {
    //                 undirectedEdges[i].erase(j);  // 删除无向边
    //                 directedEdges[i].insert(j);   // 插入有向边
    //             }
    //         }
    //     }
    // }

    //优化版本，如果我没理解错（应该没理解错），还是要分两种情况
    //R2/R3可以一起，因为R2是一种特殊的R3

    //注意！R1之后要重新赋值一遍！否则会出现两节点之间的环！
    // undirectedEdges_copy = undirectedEdges;

    // for (const auto& [i, neighbors] : undirectedEdges_copy) 
    // {
    //     for (const auto& j : neighbors) 
    //     {
    //         //情况一
    //         bool flag = false;
    //         for (int k=0;k<numRows;k++) 
    //         {
    //             if((directedEdges[i].find(k) != directedEdges[i].end()) && (directedEdges[k].find(j)!= directedEdges[k].end())) 
    //             {
    //                 flag = true;
    //                 break;
    //             }
    //         }
    //         if(flag) 
    //         {
    //             undirectedEdges[i].erase(j);
    //             directedEdges[i].insert(j);
    //             continue;
    //         } 

    //         //情况二
    //         flag = false;
    //         for (int k=0;k<numRows;k++) 
    //         {
    //             if((directedEdges[j].find(k) != directedEdges[j].end()) && (directedEdges[k].find(i) != directedEdges[k].end()))
    //             {
    //                 flag = true;
    //                 break;
    //             }
    //         }
    //         if(flag) 
    //         {
    //             undirectedEdges[j].erase(i);
    //             directedEdges[j].insert(i);
    //         } 
    //     }
    // }

    //最后一步如果还有无向边，那就稳妥起见直接删去

    //除此之外，为了稳妥，还要进行无环验证，如果有环需要修改代码


    //转换为 std::vector<std::vector<int>> 类型
    std::vector<std::vector<int>> result;
    for (const auto& [node, neighbors] : directedEdges) {
        for (const auto& neighbor : neighbors) {
            result.push_back({node, neighbor});
        }
    }
    return result;
}

// pybind11 模块定义
PYBIND11_MODULE(pc, m) {
    m.def("pc_algorithm", &pc_algorithm, "Perform PC Algorithm and return directed edges",
          pybind11::arg("independence_test"), pybind11::arg("data"));
    
    // 将 chi_test 封装到模块中
    m.def("chi_test", &chi_test, "Chi-square test for independence",
          pybind11::arg("data1"), pybind11::arg("data2"), pybind11::arg("alpha"));

    // 将 g_sq_test 封装到模块中
    m.def("g_sq_test", &g_sq_test, "G-squared test for independence",
          pybind11::arg("data1"), pybind11::arg("data2"), pybind11::arg("alpha"));
}

// int main()
// {
//     return 0;
// }