#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include<bits/stdc++.h>
using namespace std;

template <typename T>
class LatticeBoltzmann {
private:
    int n;
    int m;
    T*** f;
    T*** feq;
    T** u;
    T** v;
    T** rho;
    T*** source;
    T* w;
    T* ex;
    T* ey;
    T nu;
    T omega;
    T dx;
    T dy;
    T dt;
    T rho0;
    T tau;
    T g;
    T t1;
    T t2;
    T t3;
    T t4;
    T s;
    T sx;
    T sy;
    T* uexact;
    T** shearStress;

public:
    LatticeBoltzmann(int nx, int ny)
        : n(nx), m(ny), dx(1.0), dy(1.0), dt(1.0), rho0(1.0), tau(0.80), g(9.8) {

        f = new T**[9];
        feq = new T**[9];
        source = new T**[9];
        for (int k = 0; k < 9; k++) {
            f[k] = new T*[n + 1];
            feq[k] = new T*[n + 1];
            source[k] = new T*[n + 1];
            for (int i = 0; i <= n; i++) {
                f[k][i] = new T[m + 1];
                feq[k][i] = new T[m + 1];
                source[k][i] = new T[m + 1];
            }
        }

        u = new T*[n + 1];
        v = new T*[n + 1];
        rho = new T*[n + 1];
        for (int i = 0; i <= n; i++) {
            u[i] = new T[m + 1];
            v[i] = new T[m + 1];
            rho[i] = new T[m + 1];
        }

        w = new T[9];
        ex = new T[9];
        ey = new T[9];
        uexact = new T[m + 1];

        w[0] = static_cast<T>(4.0 / 9.0);
        w[1] = static_cast<T>(1.0 / 9.0);
        w[2] = static_cast<T>(1.0 / 9.0);
        w[3] = static_cast<T>(1.0 / 9.0);
        w[4] = static_cast<T>(1.0 / 9.0);
        w[5] = static_cast<T>(1.0 / 36.0);
        w[6] = static_cast<T>(1.0 / 36.0);
        w[7] = static_cast<T>(1.0 / 36.0);
        w[8] = static_cast<T>(1.0 / 36.0);

        ex[0] = static_cast<T>(0.0);
        ex[1] = static_cast<T>(1.0);
        ex[2] = static_cast<T>(0.0);
        ex[3] = static_cast<T>(-1.0);
        ex[4] = static_cast<T>(0.0);
        ex[5] = static_cast<T>(1.0);
        ex[6] = static_cast<T>(-1.0);
        ex[7] = static_cast<T>(-1.0);
        ex[8] = static_cast<T>(1.0);

        ey[0] = static_cast<T>(0.0);
        ey[1] = static_cast<T>(0.0);
        ey[2] = static_cast<T>(1.0);
        ey[3] = static_cast<T>(0.0);
        ey[4] = static_cast<T>(-1.0);
        ey[5] = static_cast<T>(1.0);
        ey[6] = static_cast<T>(1.0);
        ey[7] = static_cast<T>(-1.0);
        ey[8] = static_cast<T>(-1.0);

        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                u[i][j] = static_cast<T>(0.0);
                v[i][j] = static_cast<T>(0.0);
                rho[i][j] = static_cast<T>(rho0);
                for (int k = 0; k < 9; k++) {
                    source[k][i][j] = static_cast<T>(0.0);
                }
            }
        }

        nu = (tau - static_cast<T>(0.50)) / static_cast<T>(3.0);
        omega = static_cast<T>(1.0) / ((static_cast<T>(3.0) * nu) + static_cast<T>(0.50));
    }

    void initializeF() {
        for (int j = 0; j <= m; j++) {
            for (int i = 0; i <= n; i++) {
                T t1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
                for (int k = 0; k <= 8; k++) {
                    T t2 = u[i][j] * ex[k] + v[i][j] * ey[k];
                    feq[k][i][j] = rho[i][j] * w[k] * (static_cast<T>(1.0) + (static_cast<T>(3.0) * t2) + (static_cast<T>(4.50) * t2 * t2) - (static_cast<T>(1.50) * t1));
                    f[k][i][j] = feq[k][i][j];
                }
            }
        }
    }

    void collision(int steps) {
        for (int t = 1; t <= steps; t++) {
            // Collision
            for (int j = 0; j <= m; j++) {
                for (int i = 0; i <= n; i++) {
                    T t1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
                    for (int k = 0; k <= 8; k++) {
                        T t2 = u[i][j] * ex[k] + v[i][j] * ey[k];
                        T t3 = static_cast<T>(3.0) * (ex[k] - u[i][j]);
                        T t4 = static_cast<T>(9.0) * ((ex[k] * u[i][j]) + (ey[k] * v[i][j])) * ex[k];
                        source[k][i][j] = (static_cast<T>(1.0) - static_cast<T>(0.50) / tau) * w[k] * g * (t3 + t4);
                        feq[k][i][j] = rho[i][j] * w[k] * (static_cast<T>(1.0) + (static_cast<T>(3.0) * t2) + (static_cast<T>(4.50) * t2 * t2) - (static_cast<T>(1.50) * t1));
                        f[k][i][j] = (omega * feq[k][i][j]) + ((static_cast<T>(1.0) - omega) * f[k][i][j]) + source[k][i][j];
                    }
                }
            }
            for (int j = 0; j <= m; j++) {
                for (int i = n; i >= 1; i--) {
                    f[1][i][j] = f[1][i - 1][j];
                }
                for (int i = 0; i <= n - 1; i++) {
                    f[3][i][j] = f[3][i + 1][j];
                }
            }
            // Streaming.
            for (int j = m; j >= 1; j--) {
                for (int i = 0; i <= n; i++) {
                    f[2][i][j] = f[2][i][j - 1];
                }
                for (int i = n; i >= 1; i--) {
                    f[5][i][j] = f[5][i - 1][j - 1];
                }
                for (int i = 0; i <= n - 1; i++) {
                    f[6][i][j] = f[6][i + 1][j - 1];
                }
            }
            for (int j = 0; j < m; j++) {
                for (int i = 0; i <= n; i++) {
                    f[4][i][j] = f[4][i][j + 1];
                }
                for (int i = 0; i < n; i++) {
                    f[7][i][j] = f[7][i + 1][j + 1];
                }
                for (int i = n; i >= 1; i--) {
                    f[8][i][j] = f[8][i - 1][j + 1];
                }
            }
            // Boundary Conditions
            // North and South boundary - wall
            for (int i = 0; i <= n; i++) {
                f[2][i][0] = f[4][i][0];
                f[5][i][0] = f[7][i][0];
                f[6][i][0] = f[8][i][0];
                f[4][i][m] = f[2][i][m];
                f[7][i][m] = f[5][i][m];
                f[8][i][m] = f[6][i][m];
            }
            for (int j = 0; j <= m; j++) {
                // West boundary - periodicity
                f[1][0][j] = f[1][n][j];
                f[5][0][j] = f[5][n][j];
                f[8][0][j] = f[8][n][j];
                // East boundary - periodicity
                f[3][n][j] = f[3][0][j];
                f[6][n][j] = f[6][0][j];
                f[7][n][j] = f[7][0][j];
            }
            for (int j = 0; j <= m; j++) {
                for (int i = 0; i <= n; i++) {
                    T s = static_cast<T>(0.0);
                    for (int k = 0; k <= 8; k++) {
                        s += f[k][i][j];
                    }
                    rho[i][j] = s;
                }
            }
            for (int j = 0; j <= m; j++) {
                for (int i = 0; i <= n; i++) {
                    T sx = static_cast<T>(0.0);
                    T sy = static_cast<T>(0.0);
                    for (int k = 0; k <= 8; k++) {
                        sx += f[k][i][j] * ex[k];
                        sy += f[k][i][j] * ey[k];
                    }
                    u[i][j] = (sx / rho[i][j]) + (g * static_cast<T>(0.50) / rho[i][j]);
                    v[i][j] = sy / rho[i][j];
                }
            }
        }
    }

    void compute_exact() {
        for (int j = 0; j <= m; j++) {
        
            T y = static_cast<T>(j) / static_cast<T>(m) * m;
            
            
            uexact[j] = (-g / (2 * nu)) * (m * m - y * y);
        }

        // for (int j = 0; j <= m; j++) {
        //     uexact[j] = static_cast<T>(-0.50) * dpdx * (static_cast<T>(j * j) - (static_cast<T>(m * j))) / nu;
        // }
    }

    void generate_data() {
        std::ofstream file1("field.txt");
        std::ofstream file2("comparison.txt");

        for (int j = 0; j <= m; j++) {
            for (int i = 0; i <= n; i++) {
                file1 << i << " " << j << " " << u[i][j] << " " << v[i][j] << " " << rho[i][j] << std::endl;
            }
            file1 << std::endl;
        }

        for (int j = 0; j <= m; j++) {
            file2 << static_cast<double>(j) / m << " " << u[60][j] << " " << v[60][j] << " " << uexact[j] << std::endl;
        }

        file1.close();
        file2.close();
    }
    void initializeShearStress() {
        shearStress = new T*[n + 1];
        for (int i = 0; i <= n; i++) {
            shearStress[i] = new T[m + 1];
            memset(shearStress[i], 0, sizeof(T) * (m + 1));
        }
    }
    void computeShearStress() {
        for (int j = 1; j < m; j++) {
            for (int i = 1; i < n; i++) {
                
                T dudx = static_cast<T>((u[i + 1][j] - u[i - 1][j]) / (2.0 * dx));


                shearStress[i][j] = nu * (dudx);
                
            }
        }
    }
    void saveShearStressData() {
    ofstream file("sheardata.txt");
    

    
    for (int j = 0; j <= m; j++) {
        for (int i = 0; i <= n; i++) {
            file << i << " " << j << " " << shearStress[i][j] << endl;
        }
    }

    file.close();
}

};

int main() {
    int nx = 100;
    int ny = 30;
    LatticeBoltzmann<long long> lbm(nx, ny);
    int time = 10000;

    lbm.initializeF();
    lbm.collision(time);
    lbm.compute_exact();
    lbm.generate_data();

    lbm.initializeShearStress(); 
    lbm.computeShearStress();
    lbm.saveShearStressData();

    return 0;
}
