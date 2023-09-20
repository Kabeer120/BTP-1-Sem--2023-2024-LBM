#include <iostream>
#include <cmath>
#include<fstream>
#include<bits/stdc++.h>
#include <cstdlib>
#include <ctime>
using namespace std;


int main() {
    int Nx = 30;                 // resolution x-dir
    int Ny = 30;                 // resolution y-dir
    int rho0 = 100;               // average density
    double tau = 0.6;  
    double PI = 3.14;
               // collision timescale
    int Nt = 100;                // number of timesteps
    bool plotRealTime = true; 
    double inv_tau = -1.0 / tau;    // switch on for plotting as the simulation goes along
    
    // Lattice speeds / weights
    int NL = 9;
    int idxs[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    int cxs[] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
    int cys[] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
    double weights[] = {4.0/9, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36};
    double F[Ny][Nx][NL];
    std::srand(std::time(nullptr)); // Seed the random number generator
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            for (int l = 0; l < NL; ++l) {
                F[y][x][l] = 1.0;  // Initialize with ones
                F[y][x][l] += 0.01 * (static_cast<double>(std::rand()) / RAND_MAX - 0.5); // Add random noise
            }
        }
    }
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            F[y][x][3] += 2.0 * (1.0 + 0.2 * std::cos(2.0 * PI * static_cast<double>(x) / Nx * 4.0));
        }
    }
    int X[Ny][Nx], Y[Ny][Nx];
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            X[y][x] = x;
            Y[y][x] = y;
        }
    }

    // Calculate rho for each grid point
    double rho[Ny][Nx];
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            double sum = 0.0;
            for (int i = 0; i < NL; ++i) {
                sum += F[y][x][i];
            }
            rho[y][x] = sum;
        }
    }

    // Normalize F based on rho
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            for (int i = 0; i < NL; ++i) {
                F[y][x][i] *= rho0 / rho[y][x];
            }
        }
    }
    

    bool cylinder[Ny][Nx];
    int centerX = Nx / 4;
    int centerY = Ny / 2;
    int radius = Ny / 4;
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            int dx = x - centerX;
            int dy = y - centerY;
            cylinder[y][x] = dx * dx + dy * dy < radius * radius;
        }
    }

    // fig = plt.figure(figsize=(4,2), dpi=80)

    for(int it  =0 ; it < Nt;it++){

        cout << it << endl;

        // # Drift
        // for i, cx, cy in zip(idxs, cxs, cys):
        //     F[:,:,i] = np.roll(F[:,:,i], cx, axis=1)
        //     F[:,:,i] = np.roll(F[:,:,i], cy, axis=0)
        
        std::vector<std::vector<double>> bndryF;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (cylinder[i][j]) {
                    std::vector<double> temp;
                    for (int k = 0; k < NL; k++) {
                        temp.push_back(F[i][j][k]);
                    }
                    bndryF.push_back(temp);
                }
            }
        }

        for (int i = 0; i < bndryF.size(); i++) {
            std::vector<double> reordered = {
                bndryF[i][0], bndryF[i][5], bndryF[i][6], bndryF[i][7],
                bndryF[i][8], bndryF[i][1], bndryF[i][2], bndryF[i][3], bndryF[i][4]
            };
            bndryF[i] = reordered;
        }

        

        std::vector<std::vector<double>> rho(Ny, std::vector<double>(Nx, 0.0));
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                double sum = 0.0;
                for (int k = 0; k < NL; k++) {
                    sum += F[k][i][j];
                }
                rho[i][j] = sum;
            }
        }

        // Calculate ux
        std::vector<std::vector<double>> ux(Ny, std::vector<double>(Nx, 0.0));
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                double sum = 0.0;
                for (int k = 0; k < NL; k++) {
                    sum += F[k][i][j] * cxs[k];
                }
                ux[i][j] = sum / rho[i][j];
            }
        }

        // Calculate uy
        std::vector<std::vector<double>> uy(Ny, std::vector<double>(Nx, 0.0));
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                double sum = 0.0;
                for (int k = 0; k < NL; k++) {
                    sum += F[k][i][j] * cys[k];
                }
                uy[i][j] = sum / rho[i][j];
            }
        }

        std::vector<std::vector<std::vector<double>>> Feq(Ny, std::vector<std::vector<double>>(Nx, std::vector<double>(NL, 0.0)));

        for (int i = 0; i < NL; i++) {
            int cx = cxs[i];
            int cy = cys[i];
            double w = weights[i];

            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nx; k++) {
                    Feq[j][k][i] = rho[j][k] * w * (1 + 3 * (cx * ux[j][k] + cy * uy[j][k]) +
                                                     9 * std::pow(cx * ux[j][k] + cy * uy[j][k], 2) / 2 -
                                                     3 * (ux[j][k] * ux[j][k] + uy[j][k] * uy[j][k]) / 2);
                }
            }
        }



        for (int i = 0; i < NL; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nx; k++) {
                    F[i][j][k] += inv_tau * (Feq[j][k][i] - F[i][j][k]);
                }
            }
        }

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (cylinder[i][j]) {
                    for (int k = 0; k < NL; k++) {
                        F[k][i][j] = bndryF[i][k];
                    }
                }
            }
        }


        // if (plotRealTime and (it % 10) == 0) or (it == Nt-1):
        //     plt.cla()
        //     ux[cylinder] = 0
        //     uy[cylinder] = 0
        //     vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) - (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
        //     vorticity[cylinder] = np.nan
        //     vorticity = np.ma.array(vorticity, mask=cylinder)
        //     plt.imshow(vorticity, cmap='bwr')
        //     plt.imshow(~cylinder, cmap='gray', alpha=0.3)
        //     plt.clim(-.1, .1)
        //     ax = plt.gca()
        //     ax.invert_yaxis()
        //     ax.get_xaxis().set_visible(False)
        //     ax.get_yaxis().set_visible(False)   
        //     ax.set_aspect('equal')  
        //     plt.pause(0.001)


    }




    
    
    
    return 0;
}






