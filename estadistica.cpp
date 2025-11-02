#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <float.h>
#include <cmath>

int main() {
    std::ofstream sal("salida_est_lap.txt");
    for(int u = 3; u <= 3; ++u) {
        std::ifstream ent("n" + std::to_string(u) + "/lap" + std::to_string(u) + "_val_lista.csv");
        double minv = DBL_MAX, maxv = -DBL_MAX;
        double M1 = 0.0, M2 = 0.0, M3 = 0.0, M4 = 0.0;
        double sum_abs = 0.0;
        long long n = 0;
        ent.tie(nullptr);
        double l;
        while (ent >> l) {
            if (l < minv) minv = l;
            if (l > maxv) maxv = l;
            sum_abs += std::fabs(l);
            n++;
            double delta = l - M1;
            double delta_n = delta / n;
            double delta_n2 = delta_n * delta_n;
            double term1 = delta * (n - 1) * delta_n;
            M4 += term1 * delta_n2 * (n * n - 3 * n + 3)
                + 6 * delta_n2 * M2
                + 4 * delta_n * M3;
            M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
            M2 += term1;
            M1 += delta_n;
        }
        double mean = M1;
        double var_pop = M2 / n;
        double std_pop = std::sqrt(var_pop);
        double m3 = M3 / n;
        double m4 = M4 / n;
        double skew = (var_pop > 0) ? m3 / std::pow(var_pop, 1.5) : 0.0;
        double kurt = (var_pop > 0) ? m4 / (var_pop * var_pop) : 0.0;
        double energy_norm = sum_abs / n;
        sal << std::fixed << std::setprecision(2);
        sal << mean << "\t";
        sal << std_pop << "\t";
        sal << skew << "\t";
        sal << kurt << "\t";
        sal << minv << "\t";
        sal << maxv << "\t";
        sal << energy_norm << "\t";
        sal << energy_norm * u << "\n";
    }
    return 0;
}