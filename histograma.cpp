#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>

int main() {
    for(int u = 1; u <= 10; ++u) {
        std::ifstream ent("n" + std::to_string(u) + "/lap" + std::to_string(u) + "_val_lista.csv");
        std::map<std::pair<double, double>, int> mapa;
        for(double i = -0.5; i < 20; i += 0.5) {
            mapa[{i, i + 0.5}] = 0;
        }
        double actual;
        while(ent >> actual) {
            for(auto p = mapa.begin(); p != mapa.end(); ++p) {
                if(actual > (*p).first.first && actual <= (*p).first.second) {
                    (*p).second += 1;
                }
            }
        }
        for(auto p = mapa.begin(); p != mapa.end(); ++p) {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << (*p).first.second << "\t" << (*p).second << "\n";
        }
        std::cout << "\n";
    }
}