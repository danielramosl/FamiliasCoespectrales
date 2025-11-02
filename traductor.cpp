#include <fstream>
#include <iostream>
#include <vector>

void imp_vec(const std::vector<int> &arr) {
    for(int i = 0; i < arr.size(); ++i) {
        std::cout << arr[i];
        if(i < arr.size() - 1) {
            std::cout << " ";
        }
    }
    std::cout << "\n";
}

std::vector<std::vector<int>> g6_mat(const std::string &s) {
    auto val = [](char c) { return int(c) - 63; };
    int n = val(s[0]);
    int p = 1;
    int m = n * (n - 1) / 2;
    std::vector<std::vector<int>> A(n, std::vector<int>(n));
    for (int k = 0; k < m; ++k) {
        int q = k / 6;
        int r = 5 - (k % 6);
        int bit = (val(s[p + q]) >> r) & 1;
        int kk = k, j = 1;
        while (kk >= j) { kk -= j; ++j; }
        int i = kk;
        if (bit) {
            A[i][j] = 1;
            A[j][i] = 1;
        }
    }
    return A;
}

std::vector<std::vector<int>> multi(const std::vector<std::vector<int>> &a, const std::vector<std::vector<int>> &b) {
    std::vector<std::vector<int>> vec(a.size(), std::vector<int>(b[0].size(), 0));
    for(int i = 0; i < vec.size(); ++i) {
        for(int j = 0; j < vec[0].size(); ++j) {
            for(int k = 0; k < a[0].size(); ++k) {
                vec[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return vec;
}

std::vector<std::vector<int>> toeplitz(const std::vector<std::vector<int>> &m, const int &id) {
    int tam = m.size() - id, a = m[id][id];
    std::vector<int> u(tam - 1), nu(tam - 1, 0);
    std::vector<std::vector<int>> t(tam + 1, std::vector<int>(tam));
    for(int i = 0; i < tam; ++i) {
        for(int j = i; j < tam; ++j) {
            if(i == j) {
                t[i][j] = 1;
                t[i + 1][j] = -a;
            } else {
                t[i][j] = 0;
            }
        }
    }
    for(int i = 0; i < tam - 1; ++i) {
        u[i] = m[id + i + 1][id];
    }
    for(int i = 0; i < tam - 1; ++i) {
        int alpha = 0;
        for (int j = 0; j < tam - 1; ++j) {
            alpha += m[id][id + 1 + j] * u[j];
        }
        for (int j = 0; j <= tam - (i + 2); ++j) {
            t[j + (i + 2)][j] = -alpha;
        }
        for(int j = 0; j < tam - 1; ++j) {
            for(int k = 0; k < tam - 1; ++k) {
                nu[j] += u[k] * m[id + 1 + j][id + 1 + k];
            }
        }
        u.swap(nu);
        std::fill(nu.begin(), nu.end(), 0);
    }
    return t;
}

std::vector<int> berkowitz(const std::vector<std::vector<int>> &a) {
    std::vector<std::vector<int>> vec;
    for (int i = 0; i < a.size(); ++i) {
        std::vector<std::vector<int>> t = toeplitz(a, i);
        vec = (i == 0) ? t : multi(vec, t);
    }
    std::vector<int> res(vec.size());
    for(int i = 0; i < res.size(); ++i) {
        res[i] = vec[i][0];
    }
    return res;
}
 
int main() {
    std::string s;
    std::cin >> s;
    auto mat = g6_mat(s);
    for(int i = 0; i < mat.size(); ++i) {
        imp_vec(mat[i]);
    }
    imp_vec(berkowitz(mat));
}