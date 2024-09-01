#include <iostream>
#include <vector>
#include <numeric>

// 计算均值
double mean(const std::vector<double>& data) {
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

// 计算皮尔逊相关系数
double pearsonCorrelation(const std::vector<double>& X, const std::vector<double>& Y) {
    int n = X.size();
    double meanX = mean(X);
    double meanY = mean(Y);

    double numerator = 0.0;
    double sumXSquare = 0.0;
    double sumYSquare = 0.0;

    for (int i = 0; i < n; ++i) {
        double deltaX = X[i] - meanX;
        double deltaY = Y[i] - meanY;
        numerator += deltaX * deltaY;
        sumXSquare += deltaX * deltaX;
        sumYSquare += deltaY * deltaY;
    }

    double denominator = std::sqrt(sumXSquare * sumYSquare);

    if (denominator == 0) {
        return 0;  // 避免除以0
    }

    return numerator / denominator;
}

int main() {
    std::vector<double> X = { 3,5,1,6,7,2,8,9,4 };
    std::vector<double> Y = { 5,3,2,6,8,1,7,9,4 };

    double r = pearsonCorrelation(X, Y);
    std::cout << "pearson correlation: " << r << std::endl;

    return 0;
}
