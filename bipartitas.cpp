#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

int main() {
    for(int u = 2; u <= 10; ++u) {
        std::ifstream adj("n" + std::to_string(u) + "/adj" + std::to_string(u) + ".dat");
        std::ifstream val("n" + std::to_string(u) + "/adj" + std::to_string(u) + "_val.csv");
        std::string valores, familia;
        int max = -1; std::string ids = "";
        while(std::getline(val, valores)) {
            std::getline(adj, familia);
            std::stringstream ss(valores);
            double suma = 0, actual;
            while(ss >> actual) {
                suma += actual;
            }
            if(std::abs(suma) <= 1e-4) {
                int cont = std::count(familia.begin(), familia.end(), ' ') + 1;
                if(cont > max) {
                    ids = familia + "\n";
                    max = cont;
                } else if(cont == max) {
                    ids += familia + "\n";
                }
            }
        }
        std::ofstream sal("n" + std::to_string(u) + "/adj" + std::to_string(u) + "_bipar.dat");
        sal << max << "\n" << ids;
    }
}