#include <bitset>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <vector>

void imp_vec1(const std::vector<double> &arr, std::ofstream &sal) {
    sal << std::fixed << std::setprecision(2);
    for(int i = 0; i < arr.size(); ++i) {
        sal << arr[i];
        if (i < arr.size() - 1) {
            sal << ",";
        }
    }
    sal << "\n";
}

void imp_vec2(const std::vector<double> &arr, std::ofstream &sal) {
    sal << std::fixed << std::setprecision(2);
    for(int i = 0; i < arr.size(); ++i) {
        sal << arr[i];
        sal << "\n";
    }
}

std::vector<double> jacobi(std::vector<std::vector<double>> &a) {
    int n = a.size();
    double tol = 1e-12;
    for (int iter = 0; iter < 2000; ++iter) {
        int p = 0, q = 1;
        double maxabs = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (std::fabs(a[i][j]) > maxabs) {
                    maxabs = std::fabs(a[i][j]);
                    p = i;
                    q = j;
                }
            }
        }
        if (maxabs < tol) {
            break;
        }
        double app = a[p][p], aqq = a[q][q], apq = a[p][q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0.0)
            ? 1.0 / (std::fabs(tau) + std::sqrt(1.0 + tau * tau))
            : -1.0 / (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
        double c = 1.0 / std::sqrt(1.0 + t * t);
        double s = t * c;
        for (int k = 0; k < n; ++k) {
            if (k != p && k != q) {
                double aip = a[k][p];
                double aiq = a[k][q];
                a[k][p] = c * aip - s * aiq;
                a[p][k] = a[k][p];
                a[k][q] = s * aip + c * aiq;
                a[q][k] = a[k][q];
            }
        }
        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
    }
    std::vector<double> vals(n);
    for (int i = 0; i < n; ++i) {
        vals[i] = a[i][i];
    }
    return vals;
}

std::vector<std::vector<double>> laplaciana(const std::vector<std::bitset<10>> &a) {
    std::vector<std::vector<double>> l(a.size(), std::vector<double>(a.size(), 0));
    for(int i = 0; i < a.size(); ++i) {
        int grado = 0;
        for(int j = 0; j < a.size(); ++j) {
            if(a[i][j] == 1) {
                grado += 1;
                l[i][j] = -1;
            }
        }
        l[i][i] = grado;
    }
    return l;
}

std::vector<std::bitset<10>> g6_mat(const std::string &s) {
    auto val = [](char c) { return int(c) - 63; };
    int n = val(s[0]);
    int p = 1;
    int m = n * (n - 1) / 2;
    std::vector<std::bitset<10>> A(n);  
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

int main() {
    
    for(int u = 1; u <= 10; ++u) {
        std::ifstream adj("n" + std::to_string(u) + "/lap" + std::to_string(u) + ".dat");
        std::ifstream cont("n" + std::to_string(u) + "/n" + std::to_string(u) + ".g6");
        std::ofstream sal ("n" + std::to_string(u) + "/lap" + std::to_string(u) + "_val.csv");
        std::ofstream sal2 ("n" + std::to_string(u) + "/lap" + std::to_string(u) + "_val_lista.csv");
        
        std::string matrices((std::istreambuf_iterator<char>(cont)),
            std::istreambuf_iterator<char>());
        int tam = 1 + ((u * (u - 1) / 2 + 5) / 6);
        std::string s;
        std::cout << "HOla\n";
        while(std::getline(adj, s)) {
            std::stringstream familia(s);
            int id_rep;
            familia >> id_rep;
            std::string s_rep(&matrices[(tam + 1) * id_rep], &matrices[(tam + 1) * id_rep + tam]);
            std::vector<std::vector<double>> mat_rep = laplaciana(g6_mat(s_rep));
            std::vector<double> valores = jacobi(mat_rep);
            imp_vec1(valores, sal);
            imp_vec2(valores, sal2);
        }
    }
}