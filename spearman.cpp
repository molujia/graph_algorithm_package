#include <iostream>
#include <vector>
#include <algorithm>

// 计算数据排名
std::vector<int> calculateRanks(const std::vector<double>& data) {
    int n = data.size();
    std::vector<int> ranks(n);
    std::vector<std::pair<double, int>> dataWithIndex(n);

    for (int i = 0; i < n; ++i) {
        dataWithIndex[i] = { data[i], i };
    }

    std::sort(dataWithIndex.begin(), dataWithIndex.end());

    //排名通常从1开始而不是0
    for (int i = 0; i < n; ++i) {
        ranks[dataWithIndex[i].second] = i + 1;
    }

    return ranks;
}

// 计算斯皮尔曼相关系数
double spearmanCorrelation(const std::vector<double>& X, const std::vector<double>& Y) {
    int n = X.size();

    // 计算X和Y的排名
    std::vector<int> rankX = calculateRanks(X);
    std::vector<int> rankY = calculateRanks(Y);

    // 计算排名差的平方和
    double dSquaredSum = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = rankX[i] - rankY[i];
        dSquaredSum += d * d;
    }

    // 计算斯皮尔曼相关系数
    double ans = 1 - (6 * dSquaredSum) / (n * (n * n - 1));
    return ans;
}

int main() {
    std::vector<double> X = { 3,5,1,6,7,2,8,9,4 };
    std::vector<double> Y = { 5,3,2,6,8,1,7,9,4 };

    double ret = spearmanCorrelation(X, Y);
    std::cout << "spearman correlation: " << ret << std::endl;

    return 0;
}
