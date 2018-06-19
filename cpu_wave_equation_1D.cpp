#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double normal_pdf(double x, double m, double s) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

int main() {
    /* Problem parameters */
    double c = 1.0;  // Propagation speed of the wave.
    double L = 1.0;  // Length of the domain.
    int N = 10;     // Number of grid points.

    double dx = L/N;
    double dt = 0.01;

    double alpha = c*c * dt*dt / (2*dx*dx);

    /* Discretize the wave equation using the Crank-Nicolson method and
     * solve the resulting linear system A*u = b at every time step.
     */
    MatrixXf A(N,N);
    VectorXf b(N);
    vectorXf u_n(N);
    vectorXf u_nm1(N);

    // Initialize matrix of coefficients.
    A(0,0) = 1 + 2*alpha;
    A(1,1) = 1 + 2*alpha;
    A(0,1)   = -alpha;
    A(N,N-1) = -alpha

    for (i = 1; i < N-1; i++) {
        A(i,i) = 1 + 2*alpha;
        A(i,i+1) = -alpha;
        A(i,i-1) = -alpha;
    }

    // Set initial conditions.
    u_nm1(0) = 0;
    u_nm1(N) = 0;
    for (i = 1; i < N-1; i++)
        u_nm1(i) = normal_pdf(i*dx, 0.5, 0.1)

    // Take the first time step.
    for (i = 1; i < N-1; i++)
        u_n(i) = u_nm1(i) + (c*c/2) * (u_nm1(i+1) - 2*u_nm1(i) + u_nm1(i-1));

    cout << "u_n-1 = " << u_nm1 << endl;
    cout << "u_n = " << u_n << endl;

    // for (i = 1; i < N-1; i++) {
    //     b(i) = 2*(1-alpha)*u_n(i) - u_nm1(i) + alpha*(u_n(i+1) + u_n(i-1))
    //            - dt*dt*f_n(i);
    // }

    // A << 1,2,3,  4,5,6,  7,8,10;
    // b << 3, 3, 4;
    
    // Vector3f x = A.colPivHouseholderQr().solve(b);
}