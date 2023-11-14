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
    T* w;
    T* ex;
    T* ey;
    
    T omega;
    T dx;
    T dy;
    
    T rho0;
    
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
    double cs2 = 1.0/3.0;
    double cs = sqrt(cs2);
    double Ma = 0.01;
    
    double uMax = Ma*cs;
    double Re = 100;
    double Kn = Ma/Re;
    double refLen = 1;
    double nu = uBoundary*refLen/Re;
    double tau = kinVisc/cs2;
    double dt = refLen/n;
    double tauNdim = tau/dt;
    double beta = 1.0 / (2.0*tauNdim + 1.0);
    int simulationTime = 3 /nu;

public:
    LatticeBoltzmann(int nx, int ny)
        : N(nx), M(ny), dx(1.0), dy(1.0), dt(1.0) {
        int n = N+2;
        int m = M + 2;

        
        
        f = new T**[9];
        feq = new T**[9];
        source = new T**[9];
        for (int k = 0; k < 9; k++) {
            f[k] = new T*[n ];
            feq[k] = new T*[n ];
            
            for (int i = 0; i < n; i++) {
                f[k][i] = new T[m ];
                feq[k][i] = new T[m ];
                
            }
        }

        u = new T*[n ];
        v = new T*[n ];
        rho = new T*[n ];
        for (int i = 0; i < n; i++) {
            u[i] = new T[m ];
            v[i] = new T[m ];
            rho[i] = new T[m ];
        }

        w = new T[9];
        ex = new T[9];
        ey = new T[9];
        opp = new T[9];
        uexact = new T[m ];

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

        for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        u[i][j] = static_cast<T>(0.0);
                        v[i][j] = static_cast<T>(0.0);
                        rho[i][j] = static_cast<T>(rho0);
                        
                    }
                }
        
                nu = (tau - static_cast<T>(0.50)) / static_cast<T>(3.0);
                omega = static_cast<T>(1.0) / ((static_cast<T>(3.0) * nu) + static_cast<T>(0.50));
            }

    void initializeF() {
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
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
        for (int t = 1; t <= simulationTime; t++) {

            // fiLLING GHOST NODE WITH EQUILBRIUM VALUES.
            for(int i  =0 ; i < n;i++){
                for(int k = 0 ; k <=8;k++){
                    f[k][i][0] = (  w[k] )*(1.0 + 3.0*uBoundary*ex[k] + 4.5*(uBoundary*ex[k]*uBoundary*ex[k]) - 1.5*uBoundary*uBoundary);
                }
            }
            for(int i  =0 ; i < n;i++){
                for(int k = 0 ; k <=8;k++){
                    f[k][i][m] = (  w[k] )*(1.0 + 3.0*uBoundary*ex[k] + 4.5*(uBoundary*ex[k]*uBoundary*ex[k]) - 1.5*uBoundary*uBoundary);
                }
            }
            for(int j  =0 ; i < m;i++){
                for(int k = 0 ; k <=8;k++){
                    f[k][0][j] = (  w[k] )*(1.0 + 3.0*uBoundary*ex[k] + 4.5*(uBoundary*ex[k]*uBoundary*ex[k]) - 1.5*uBoundary*uBoundary);
                }
            }
            for(int j  =0 ; j < m;i++){
                for(int k = 0 ; k <=8;k++){
                    f[k][m][j] = (  w[k] )*(1.0 + 3.0*uBoundary*ex[k] + 4.5*(uBoundary*ex[k]*uBoundary*ex[k]) - 1.5*uBoundary*uBoundary);
                }
            }
            // Collision
            for (int j = 1; j < m-1; j++) {
                for (int i = 1; i < n-1; i++) {
                    T t1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
                    for (int k = 0; k <= 8; k++) {
                        T t2 = u[i][j] * ex[k] + v[i][j] * ey[k];
                        T t3 = static_cast<T>(3.0) * (ex[k] - u[i][j]);
                        T t4 = static_cast<T>(9.0) * ((ex[k] * u[i][j]) + (ey[k] * v[i][j])) * ex[k];
                        
                        feq[k][i][j] = rho[i][j] * w[k] * (static_cast<T>(1.0) + (static_cast<T>(3.0) * t2) + (static_cast<T>(4.50) * t2 * t2) - (static_cast<T>(1.50) * t1));
                        f[k][i][j] = (omega * feq[k][i][j]) + ((static_cast<T>(1.0) - omega) * f[k][i][j]) ;
                    }
                }
            }
            
            // Updating ghost nodes at 0th  column before advecting after every time step to maintain constant conditions.
            //( 1- EAST, 2- NORTH, 3- WEST, 4- SOUTH,5 - NE, 6-NW,7- SW,8- SE )
            //LEFT WALL GHOST UPDATE
            for(int i = 0; i <n;i++){
                f[3][i][0] = f[3][i][1];
                f[6][i][0] = f[6][i][1];
                f[7][i][0] = f[7][i][1];
            }
            //RIGHT WALL GHOST UPDATE

            for(int i = 0 ; i <n;i++){
                f[1][i][N+1] = f[1][i][N];
                f[5][i][N+1] = f[5][i][N];
                f[8][i][N+1] = f[8][i][N];

            }

            //Bottom WALL GHOST UPDATE

            for(int j = 0 ; j < m; j++){
                f[4][0][j] = f[4][1][j];
                f[7][0][j] = f[7][1][j];
                f[8][0][j] = f[8][1][j];

            }
            //TOP WALL GHJOST UPDATE

            for(int j = 0 ; j < m; j++){
                f[2][N+1][j] = f[2][N][j];
                f[5][N+1][j] = f[5][N][j];
                f[6][N+1][j] = f[6][N][j];

            }



            for (int j = 1; j < m-1; j++) {
                for (int i = n-1; i >= 1; i--) {
                    f[1][i][j] = f[1][i - 1][j];
                }
                for (int i = 1; i < n - 1; i++) {
                    f[3][i][j] = f[3][i + 1][j];
                }
            }
            // Streaming.
            for (int j =1; j <= m-1; j++) {
                for (int i = n-1; i >=1; i--) {
                    f[2][i][j] = f[2][i][j - 1];
                }
                for (int i = n-1; i >= 1; i--) {
                    f[5][i][j] = f[5][i - 1][j - 1];
                }
                for (int i = 1; i <= n - 1; i++) {
                    f[6][i][j] = f[6][i + 1][j - 1];
                }
            }
            for (int j = 1; j < m-1; j++) {
                for (int i = 1; i < n-1; i++) {
                    f[4][i][j] = f[4][i][j + 1];
                }
                for (int i = 1; i < n-1; i++) {
                    f[7][i][j] = f[7][i + 1][j + 1];
                }
                for (int i = n-1; i >= 1; i--) {
                    f[8][i][j] = f[8][i - 1][j + 1];
                }
            }
            // Boundary Conditions
            // North and South boundary - wall
            for (int i = 1; i < n-1; i++) {
                f[2][i][1] = f[4][i][1];
                f[5][i][1] = f[7][i][1];
                f[6][i][1] = f[8][i][1];
                f[4][i][m-1] = f[2][i][m-1];
                f[7][i][m-1] = f[5][i][m-1];
                f[8][i][m-1] = f[6][i][m-1];
            }
            for (int j = 1; j < m-1; j++) {
                // West boundary - periodicity
                f[1][1][j] = f[1][n-1][j];
                f[5][1][j] = f[5][n-1][j];
                f[8][1][j] = f[8][n-1][j];
                // East boundary - periodicity
                f[3][n-1][j] = f[3][1][j];
                f[6][n-1][j] = f[6][1][j];
                f[7][n-1][j] = f[7][1][j];
            }
            for (int j = 1; j < m-1; j++) {
                for (int i = 1; i < n-1; i++) {
                    T s = static_cast<T>(0.0);
                    for (int k = 0; k <= 8; k++) {
                        s += f[k][i][j];
                    }
                    rho[i][j] = s;
                }
            }
            for (int j = 1; j < m-1; j++) {
                for (int i = 1; i < n-1; i++) {
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
        for (int j = 1; j < m-1; j++) {
        
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

        for (int j = 1; j < m-1; j++) {
            for (int i = 1; i < n-1; i++) {
                int position = std::pow(i-circleX, 2)+std::pow(j-circleY, 2);
                if ( position > circleR2){
                    file1 << i << " " << j << " " << u[i][j] << " " << v[i][j] << " " << rho[i][j] << std::endl;

                }
                
            }
            file1 << std::endl;
        }

        for (int j = 0; j < m-1; j++) {
            file2 << static_cast<double>(j) / m << " " << u[60][j] << " " << v[60][j] << " " << uexact[j] << std::endl;
        }

        file1.close();
        file2.close();
    }
    void initializeShearStress() {
        shearStress = new T*[n ];
        for (int i = 0; i < n-1; i++) {
            shearStress[i] = new T[m ];
            memset(shearStress[i], 0, sizeof(T) * (m));
        }
    }
    void computeShearStress() {
        for (int j = 1; j < m-1; j++) {
            for (int i = 1; i < n-1; i++) {
                
                T dudx = static_cast<T>((u[i + 1][j] - u[i - 1][j]) / (2.0 * dx));


                shearStress[i][j] = nu * (dudx);
                
            }
        }
    }
    void saveShearStressData() {
    ofstream file("sheardata.txt");
    

    
    for (int j = 1; j < m-1; j++) {
        for (int i = 1; i < n-1; i++) {

            file << i << " " << j << " " << shearStress[i][j] << endl;
        }
    }

    file.close();
}

};

int main() {
    int nx = 100;
    int ny = 30;
    LatticeBoltzmann<int> lbm(nx, ny);
    

    lbm.initializeF();
    lbm.collision();
    lbm.compute_exact();
    lbm.generate_data();

    lbm.initializeShearStress(); 
    lbm.computeShearStress();
    lbm.saveShearStressData();


    return 0;
}
