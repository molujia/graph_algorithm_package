#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <cassert>

// 函数：计算两个离散变量的频率
void calculateFrequencies(const std::vector<int>& X,
                          const std::vector<int>& Y,
                          std::map<int, int>& freqX,
                          std::map<int, int>& freqY,
                          std::map<std::pair<int, int>, int>& jointFreq) {
    assert(X.size() == Y.size());

    for (size_t i = 0; i < X.size(); ++i) {
        freqX[X[i]]++;
        freqY[Y[i]]++;
        jointFreq[{X[i], Y[i]}]++;
    }
}

// 函数：计算概率
std::map<int, double> calculateProbabilities(const std::map<int, int>& freq, int totalCount) {
    std::map<int, double> probabilities;
    for (const auto& entry : freq) {
        probabilities[entry.first] = static_cast<double>(entry.second) / totalCount;
    }
    return probabilities;
}

// 函数：计算互信息
double computeMutualInformation(const std::map<std::pair<int, int>, int>& jointFreq,
                                const std::map<int, int>& freqX,
                                const std::map<int, int>& freqY,
                                int totalCount) {
    double mutualInformation = 0.0;
    for (const auto& entry : jointFreq) {
        int x = entry.first.first;
        int y = entry.first.second;
        double pXY = static_cast<double>(entry.second) / totalCount;
        double pX = static_cast<double>(freqX.at(x)) / totalCount;
        double pY = static_cast<double>(freqY.at(y)) / totalCount;

        if (pXY > 0) { // 只计算非零概率
            mutualInformation += pXY * std::log2(pXY / (pX * pY));
        }
    }
    return mutualInformation;
}

int main() {
    // 示例数据
    std::vector<int> X = {0, 0, 1, 1, 2, 2, 0, 1, 2};
    std::vector<int> Y = {0, 1, 0, 1, 0, 1, 1, 0, 1};

    // 统计频率
    std::map<int, int> freqX;
    std::map<int, int> freqY;
    std::map<std::pair<int, int>, int> jointFreq;

    calculateFrequencies(X, Y, freqX, freqY, jointFreq);

    // 计算总数
    int totalCount = X.size();

    // 计算概率
    std::map<int, double> probX = calculateProbabilities(freqX, totalCount);
    std::map<int, double> probY = calculateProbabilities(freqY, totalCount);

    // 计算互信息
    double mi = computeMutualInformation(jointFreq, freqX, freqY, totalCount);

    // 输出结果
    std::cout << "Mutual Information: " << mi << std::endl;

    return 0;
}
