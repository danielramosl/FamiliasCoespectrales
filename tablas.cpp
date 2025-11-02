#include <algorithm>
#include <fstream>
#include <string>

int main() {
    for(int u = 1; u <= 10; ++u) {
        std::ifstream ent("n" + std::to_string(u) + "/lap" + std::to_string(u) + ".dat");
        std::string fam;
        int arr[100000] = {};
        while(std::getline(ent, fam)) {
            int cont = std::count(fam.begin(), fam.end(), ' ') + 1;
            arr[cont] += 1;
        }
        std::ofstream sal("n" + std::to_string(u) + "/lap" + std::to_string(u) + "_dis.csv");
        sal << "familia,cantidad\n";
        for(int i = 1; i < 100000; ++i) {
            if(arr[i] != 0) {
                sal << i << "," << arr[i] << "\n";
            }
        }
    }
} 