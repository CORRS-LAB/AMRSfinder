#include <RcppEigen.h>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;

// Helper function to compute S
double computeS(int s, int t, const Eigen::Ref<const MatrixXd>& XS) {
    // R uses 1-based indexing, C++ uses 0-based
    // Adjust indices
    int s_idx = s - 1;
    int t_idx = t - 1;
    
    double Sst = XS(t_idx, 1);  // Column 2 in R is column 1 in C++ (0-based)
    if (s > 1) {
        Sst = Sst - XS(s_idx - 1, 1);  // s-1 in R is s_idx-1 in C++
    }
    return std::abs(Sst);
}

// [[Rcpp::export]]
IntegerVector findMaxZ(int s, int t, int mincpgs, const Eigen::MatrixXd& XS) {
    double max_val = 0.0;
    int a_max_id = 0;
    int b_max_id = 1;
    
    if (s <= t) {
        for (int a = s; a <= t; a++) {
            if (a + mincpgs - 1 <= t) {
                for (int b = a + mincpgs - 1; b <= t; b++) {
                    double Zab = 0.0;
                    
                    if ((b - a - mincpgs == 0) || ((a == s) && (b == t)) || (a == b)) {
                        Zab = 0.0;
                    } else {
                        int lab = b - a + 1;
                        int lst = t - s + 1;
                        double S_ab = computeS(a, b, XS);
                        double S_st = computeS(s, t, XS);
                        double u = std::pow(S_ab - (static_cast<double>(lab) * S_st / lst), 2);
                        Zab = u / (lab * (1.0 - static_cast<double>(lab) / lst));
                    }
                    
                    if (std::isfinite(Zab) && ((a != s) || (b != t)) && (max_val < Zab)) {
                        max_val = Zab;
                        a_max_id = a - s + 1;
                        b_max_id = b - s + 1;
                    }
                }
            }
        }
    }
    
    IntegerVector max_id = IntegerVector::create(a_max_id, b_max_id);
    return max_id;
}