#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <sstream>
#include <thread>
#include <vector>

struct par {
    std::vector<int> coef;
    int id;
};

void imp_vec(const std::vector<int> &arr, std::ofstream &sal) {
    for(int i = 0; i < arr.size(); ++i) {
        sal << arr[i] << " ";
    }
    sal << "\n";
}

void imp_mat(const std::vector<std::vector<int>> &mat, std::ofstream &sal) {
    for(int i = 0; i < mat.size(); ++i) {
        for(int j = 0; j < mat.size(); ++j) {
            sal << mat[i][j] << " ";
        }
        sal << "\n";
    }
    sal << "\n";
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

std::vector<std::vector<int>> laplaciana(const std::vector<std::bitset<10>> &a) {
    std::vector<std::vector<int>> l(a.size(), std::vector<int>(a.size(), 0));
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

void concurrente(int ini, int tam, const std::vector<std::vector<std::bitset<10>>> &mat, std::vector<par> &conten, std::mutex &mtx) {
    std::vector<par> aux;
    for(int i = 0; i < tam; ++i) {
        std::vector<int> poli = berkowitz(laplaciana(mat[ini + i]));
        aux.push_back({poli, ini + i});
    }
    mtx.lock();
    for(par &p : aux) {
        conten.push_back(p);
    }
    mtx.unlock();
}

bool comp(const par &a, const par &b) {
    return a.coef < b.coef;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    for(int u = 1; u <= 10; ++u) {
        std::ifstream ent("n" + std::to_string(u) + "/n" + std::to_string(u) + ".g6");
        std::vector<par> conten;
        std::string s;
        std::vector<std::vector<std::bitset<10>>> mat;
        while (std::getline(ent, s)) {
            mat.push_back(g6_mat(s));
        }

        clock_t t0 = clock();
        if(u < 9) {
            for(int i = 0; i < mat.size(); ++i) {
                std::vector<int> poli = berkowitz(laplaciana(mat[i]));
                conten.push_back({poli, i});
            }
        } else {
            int total = mat.size();
            int num_hilos = std::thread::hardware_concurrency();
            int tam = total / num_hilos;
            int resto = total % num_hilos;
            std::thread hilo[num_hilos];
            std::mutex mtx;
            for(int i = 0; i < num_hilos; ++i) {
                int bloque_tam = tam + (i == num_hilos - 1 ? resto : 0);
                hilo[i] = std::thread(&concurrente, tam * i, bloque_tam, std::ref(mat), std::ref(conten), std::ref(mtx));
            }
            for(int i = 0; i < num_hilos; ++i) {
                hilo[i].join();
            }
        }
        std::sort(conten.begin(), conten.end(), comp);
        clock_t t1 = clock();
        std::ofstream sal("n" + std::to_string(u) + "/lap" + std::to_string(u) + ".dat");
        sal << conten[0].id;
        for(int i = 1; i < conten.size(); ++i) {
            if(conten[i].coef == conten[i - 1].coef) {
                sal << " " << conten[i].id;
            } else {
                sal << "\n";
                sal << conten[i].id;
            }
        }
        std::cout << "n = " << u << "\nTiempo = " << double(t1 - t0) / CLOCKS_PER_SEC << "s\nMemoria: "
        << ( (mat.size() * sizeof(std::vector<std::bitset<10>>)
            + std::accumulate(mat.begin(), mat.end(), size_t(0),
                [](size_t acc, const auto& fila){
                    return acc + fila.capacity() * sizeof(std::bitset<10>);
                })
            ) / (1024.0 * 1024.0) )
        << "MB\n" << std::flush;
    }
    return 0;
}