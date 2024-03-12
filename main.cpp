#include <iostream>
#include <cmath>
#include "inmost.h"


using namespace INMOST;
using namespace std;

void
create_matrix_and_vector(int N, Sparse::Matrix& A, Sparse::Vector& b, Sparse::Vector& real_sol)
{

    int val = 0;

    for(int i = 0; i < (N - 1) * (N - 1); i++) {

        b[i] = 26 * std::sin((double)(i % (N - 1) + 1) / N) * std::sin((double)5 * (i / (N - 1) + 1) / N) / N / N;

        real_sol[i] = std::sin((double)(i % (N - 1) + 1) / N) * std::sin((double)5 * (i / (N - 1) + 1) / N);

        A[i][i] = 4;

        if(i % (N - 1) + 1 < (N - 1)) {
            A[i][i + 1] = -1;
        } else {
            b[i] += std::sin(1) * std::sin((double)5 * (i / (N - 1) + 1) / N);
        }

        if(i % (N - 1) - 1 >= 0) {
            A[i][i - 1] = -1;
        }

        if(i + N - 1 < (N - 1) * (N - 1)) {
            A[i][i + N - 1] = -1;
        } else {
            b[i] += std::sin(5) *  std::sin((double)(i % (N - 1) + 1) / N);
        }

        if(i - N + 1 >= 0) {
            A[i][i - N + 1] = -1;
        }

    }

}

void
test_function(int N, Sparse::Vector sol, Sparse::Vector real_sol)
{
    double val, l2_norm = 0, c_norm = 0;

    for(int i = 0; i < (N - 1) * (N - 1); i++) {
        val = abs(sol[i] - real_sol[i]);
        l2_norm += val * val;

        if (val > c_norm) {
            c_norm = val;
        }
    }

    std::cout << "l2_norm: " << sqrt(l2_norm) / N << std::endl;
    std::cout << "c_norm: " << c_norm<< std::endl;
}

int
main()
{
    int N;
    std::cin >> N;

    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;
    Sparse::Vector real_sol;

    A.SetInterval(0, (N - 1) * (N - 1));
    b.SetInterval(0, (N - 1) * (N - 1));
    sol.SetInterval(0, (N - 1) * (N - 1));
    real_sol.SetInterval(0, (N - 1) * (N - 1));

    create_matrix_and_vector(N, A, b, real_sol);

    Solver S(Solver::INNER_MPTILUC);
    S.SetMatrix(A);
    
    bool solved = S.Solve(b, sol);

    cout << "num.iters: " << S.Iterations() << endl;
    cout << "prec.time: " << S.PreconditionerTime() << endl;
    cout << "iter.time: " << S.IterationsTime() << endl;
    if(!solved){
        cout << "Linear solver failure!" << endl;
        cout << "Reason: " << S.ReturnReason() << endl;
    }

    test_function(N - 1, sol, real_sol);
    return 0;
}

