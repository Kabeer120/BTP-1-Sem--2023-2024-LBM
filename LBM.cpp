#include <iostream>
#include <cmath>
#include<fstream>
#include<bits/stdc++.h>
using namespace std;

int main() {
    //  We initialized and set up the required variables.
    int n = 100;
    int m = 30;
    double f[9][n+1][m+1];
    double feq[9][n+1][m+1];
    double u[n+1][m+1];
    double v[n+1][m+1];
    double rho[n+1][m+1];
    double source[9][n+1][m+1];
    double w[9];
    double ex[9];
    double ey[9];
    double temp1, temp2, omega;
    int i, j, k, t, steps = 20000;
    double nu, dx = 1.0, dy = 1.0, dt = 1.0, rho0 = 1.0, tau = 0.80, dpdx = 1.0e-5, t1, t2, t3, t4, s, sx, sy;
    double uexact[m+1];

    // Setting up weights and directions for D2Q9 LaTtice.
    w[0] = 4.0/9.0;
    w[1] = 1.0/9.0;
    w[2] = 1.0/9.0;
    w[3] = 1.0/9.0;
    w[4] = 1.0/9.0;
    w[5] = 1.0/36.0;
    w[6] = 1.0/36.0;
    w[7] = 1.0/36.0;
    w[8] = 1.0/36.0;

    ex[0] = 0.0;
    ex[1] = 1.0;
    ex[2] = 0.0;
    ex[3] = -1.0;
    ex[4] = 0.0;
    ex[5] = 1.0;
    ex[6] = -1.0;
    ex[7] = -1.0;
    ex[8] = 1.0;

    ey[0] = 0.0;
    ey[1] = 0.0;
    ey[2] = 1.0;
    ey[3] = 0.0;
    ey[4] = -1.0;
    ey[5] = 1.0;
    ey[6] = 1.0;
    ey[7] = -1.0;
    ey[8] = -1.0;

    for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            rho[i][j] = rho0;
            source[0][i][j] = 0.0;
            source[1][i][j] = 0.0;
            source[2][i][j] = 0.0;
            source[3][i][j] = 0.0;
            source[4][i][j] = 0.0;
            source[5][i][j] = 0.0;
            source[6][i][j] = 0.0;
            source[7][i][j] = 0.0;
            source[8][i][j] = 0.0;
        }
    }

    nu = (tau - 0.50) / 3.0;
    omega = 1.0 / ((3.0 * nu) + 0.50);

    // Initialising distribution function
    for (int j = 0; j <= m; j++) {
        for (int i = 0; i <= n; i++) {
            double t1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
            for (int k = 0; k <= 8; k++) {
                double t2 = u[i][j] * ex[k] + v[i][j] * ey[k];
                feq[k][i][j] = rho[i][j] * w[k] * (1.0 + (3.0 * t2) + (4.50 * t2 * t2) - (1.50 * t1));
                f[k][i][j] = feq[k][i][j];
            }
        }
    }
    for (int t = 1; t <= steps; t++) {
    // Collision
        for (int j = 0; j <= m; j++) {
            for (int i = 0; i <= n; i++) {
                double t1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
                for (int k = 0; k <= 8; k++) {
                    double t2 = u[i][j] * ex[k] + v[i][j] * ey[k];
                    double t3 = 3.0 * (ex[k] - u[i][j]);
                    double t4 = 9.0 * ((ex[k] * u[i][j]) + (ey[k] * v[i][j])) * ex[k];
                    source[k][i][j] = (1.0 - 0.50 / tau) * w[k] * dpdx * (t3 + t4);
                    feq[k][i][j] = rho[i][j] * w[k] * (1.0 + (3.0 * t2) + (4.50 * t2 * t2) - (1.50 * t1));
                    f[k][i][j] = (omega * feq[k][i][j]) + ((1.0 - omega) * f[k][i][j]) + source[k][i][j];
                }
            }
        }
        for (int j = 0; j <= m; j++) {
            for (int i = n; i >= 1; i--) {
                f[1][i][j] = f[1][i-1][j];
            }
            for (int i = 0; i <= n-1; i++) {
                f[3][i][j] = f[3][i+1][j];
            }
        }                                       
                                                        // Streaming.
        for (int j = m; j >= 1; j--) {
            for (int i = 0; i <= n; i++) {
                f[2][i][j] = f[2][i][j-1];
            }
            for (int i = n; i >= 1; i--) {
                f[5][i][j] = f[5][i-1][j-1];
            }
            for (int i = 0; i <= n-1; i++) {
                f[6][i][j] = f[6][i+1][j-1];
            }
        }
        for (int j = 0; j < m; j++) {
            for (int i = 0; i <= n; i++) {
                f[4][i][j] = f[4][i][j+1];
            }
            for (int i = 0; i < n; i++) {
                f[7][i][j] = f[7][i+1][j+1];
            }
            for (int i = n; i >= 1; i--) {
                f[8][i][j] = f[8][i-1][j+1];
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
                double s = 0.0;
                for (int k = 0; k <= 8; k++) {
                    s += f[k][i][j];
                }
                rho[i][j] = s;
            }
        }
        for (int j = 0; j <= m; j++) {
            for (int i = 0; i <= n; i++) {
                double sx = 0.0;
                double sy = 0.0;
                for (int k = 0; k <= 8; k++) {
                    sx += f[k][i][j] * ex[k];
                    sy += f[k][i][j] * ey[k];
                }
                u[i][j] = (sx / rho[i][j]) + (dpdx * 0.50 / rho[i][j]);
                v[i][j] = sy / rho[i][j];
            }
        }
        

    }
    for (int j = 0; j <= m; j++) {
        uexact[j] = -0.50 * dpdx * ((j * j) - (m * j)) / nu;
    }

    ofstream file1("field.txt");
    ofstream file2("comparison.txt");

    for (int j = 0; j <= m; j++) {
        for (int i = 0; i <= n; i++) {
            file1 << i << " " << j << " " << u[i][j] << " " << v[i][j] << " " << rho[i][j] << endl;
        }
        file1 << endl;
    }

    for (int j = 0; j <= m; j++) {
        file2 << static_cast<double>(j) / m << " " << u[60][j] << " " << v[60][j] << " " << uexact[j] << endl;
    }


    file1.close();
    file2.close();

    // Task left = We need to use some library to plot the graph from 'comparison.txt' file.

    return 0;

}
