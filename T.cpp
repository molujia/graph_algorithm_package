#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <thread>
#include <mutex>
#include <future>

double mean(const std::vector<double>& data);
double variance(const std::vector<double>& data, double mean);
double tTest(const std::vector<double>& sample1, const std::vector<double>& sample2);
double tDistributionCDF(double t, double df);
double betaIncomplete(double a, double b, double x);
double betaContinuedFraction(double a, double b, double x);


double mean(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::invalid_argument("The data vector is empty");
    }
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}


double variance(const std::vector<double>& data, double mean) {
    if (data.size() < 2) {
        throw std::invalid_argument("Data vector size must be at least 2 to calculate variance");
    }
    double sum = 0.0;
    for (double value : data) {
        sum += (value - mean) * (value - mean);
    }
    return sum / (data.size() - 1);
}


double tTest(const std::vector<double>& sample1, const std::vector<double>& sample2) {
    //并行均值
    auto future_mean1 = std::async(std::launch::async, mean, std::ref(sample1));
    auto future_mean2 = std::async(std::launch::async, mean, std::ref(sample2));
    double mean1 = future_mean1.get();
    double mean2 = future_mean2.get();

    //并行方差
    auto future_var1 = std::async(std::launch::async, variance, std::ref(sample1), mean1);
    auto future_var2 = std::async(std::launch::async, variance, std::ref(sample2), mean2);
    double var1 = future_var1.get();
    double var2 = future_var2.get();

    double s1 = var1 / sample1.size();
    double s2 = var2 / sample2.size();

    double t = (mean1 - mean2) / std::sqrt(s1 + s2);
    double df = (s1 + s2) * (s1 + s2) / ((s1 * s1) / (sample1.size() - 1) + (s2 * s2) / (sample2.size() - 1));

    double pValue = 2.0 * (1.0 - tDistributionCDF(std::fabs(t), df));
    return pValue;
}


double tDistributionCDF(double t, double df) {
    double x = df / (t * t + df);
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    return betaIncomplete(0.5 * df, 0.5, 1.0 - x);
}


double betaIncomplete(double a, double b, double x) {
    const double EPS = std::numeric_limits<double>::epsilon();
    const int MAX_ITER = 100;

    double bt = (x == 0.0 || x == 1.0) ? 0.0 :
        std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + a * std::log(x) + b * std::log(1.0 - x));

    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * betaContinuedFraction(a, b, x) / a;
    } else {
        return 1.0 - bt * betaContinuedFraction(b, a, 1.0 - x) / b;
    }
}


double betaContinuedFraction(double a, double b, double x) {
    const double EPS = std::numeric_limits<double>::epsilon();
    const double FPMIN = std::numeric_limits<double>::min() / std::numeric_limits<double>::epsilon();
    const int MAX_ITER = 100;

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;

    for (int m = 1, m2 = 2; m <= MAX_ITER; ++m, m2 += 2) {
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < EPS) break;
    }
    return h;
}

int main() {
    std::vector<int> sample1_int, sample2_int;
    std::string line;
    
    //两组整型数据
    std::cout << "Enter the first set of integers (space-separated): ";
    std::getline(std::cin, line);
    std::istringstream iss1(line);
    int num;
    while (iss1 >> num) {
        sample1_int.push_back(num);
    }
    std::cout << "Enter the second set of integers (space-separated): ";
    std::getline(std::cin, line);
    std::istringstream iss2(line);
    while (iss2 >> num) {
        sample2_int.push_back(num);
    }
    
    std::vector<double> sample1(sample1_int.begin(), sample1_int.end());
    std::vector<double> sample2(sample2_int.begin(), sample2_int.end());
    
    try {
        double pValue = tTest(sample1, sample2);
        std::cout << "p-value: " << pValue << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
