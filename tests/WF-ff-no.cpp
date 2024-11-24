#include <vector>
#include <thread>
#include <cmath>
#include <iostream>

#include <ff/ff.hpp>

#if 0
#if !defined(USE_A2A)
#define USE_A2A
#endif
#endif

#if !defined(PRODUCT)
#define PRODUCT fake_product
#endif

using namespace ff;

struct task_t {
    int start;
    int stop;
    int iter;
};

// Compute the product Row * Columns of the given cell as Row * Transpose(Column)
double product(double **matrix, const int row, const int col) {
    double prod = 0.0;
    int len = col - row + 1;
    printf("(%d, %d) = (%d,%d..%d)*(%d,%d..%d))\n", row, col, row, col-1, col-len+1, col, row+1, row+len-1);
    for(int i = 1; i < len; i++) {
        // prod += matrix[row][col - i] * matrix[col][row + i];
        prod += matrix[row][col - i] * matrix[row + i][col];
    }
    return prod;
}

double fake_product(double **matrix, const int row, const int col) {
    return matrix[row][col-1] * matrix[col][row+1] + 3.14;
}

void initMatrix(double **matrix, int n) {
    for (int i = 0; i < n; i++)
        matrix[i][i] = ((double) (i+1)) / n;
}

void compute_iter(double **matrix, int start,int stop, int iter) {
    for(int m=start;m<stop;++m) {
        const int row=m;
        const int col=row+iter;
        double prod = PRODUCT(matrix, row, col);

        // Save the result in the corresponding cell
        matrix[row][col] = std::cbrt(prod);
        // Save the result also in the transposed cell
        // to speedup next computations using cached row accesses
        matrix[col][row] = matrix[row][col];
    }
}

struct Master:ff_monode_t<task_t> {
    Master(double **matrix, long n, int pardegree):
            matrix(matrix), n(n),pardegree(pardegree),taskv(pardegree) {
    }

    void setNextIteration() {
        ++iter;
    }

    // sends tasks to workers
    void send_tasks(int &start, int &stop) {
        const size_t chunksize = (n-iter) / pardegree;
        ssize_t more           = (n-iter) - chunksize*pardegree;

        for(int i=0; i<(pardegree-1); ++i) {
            start = stop;
            stop  = start + chunksize + (more>0 ? 1:0);
            --more;

            taskv[i] = {start, stop, iter};

            ff_send_out_to(&taskv[i], i);
        }
        if ((n-iter)==pardegree)  {
            ff_send_out_to(EOS, pardegree-1);
            --pardegree;
        }
        start = stop;
        stop  = start + chunksize + (more>0 ? 1:0);
    }

    task_t* svc(task_t*) {
        int start=0;
        int stop=0;

        send_tasks(start,stop);
        compute_iter(matrix, start,stop,iter);

        return EOS;
    }

    double **matrix; // the matrix
    int n;           // number of columns/rows
    int pardegree;   // initial number of workers
    int iter=0;      // current iterations (from 1 to n)

    std::vector<task_t> taskv;	 // memory needed for the tasks
};
struct Worker:ff_minode_t<task_t> {
    Worker(double **matrix, long n):matrix(matrix),n(n) {
    }
    task_t* svc(task_t* in) {
        compute_iter(matrix, in->start,in->stop,in->iter);
        return in; // notify, we've completed the task
    }
    double **matrix;  // the matrix
    long n;           // n. of columns/rows
};

void wavefront(double **matrix, int n, int pardegree) {
    if (pardegree==1) {  // sequential computation
        for (int k = 1; k < n; ++k) {
            for (int m = 0; m < n - k; ++m) {
                int row = m;
                int col = row + k;

                double prod = PRODUCT(matrix, row, col);

                // Save the result in the corresponding cell
                matrix[row][col] = std::cbrt(prod);
                // Save the result also in the transposed cell to speedup next computations using cached row accesses
                matrix[col][row] = matrix[row][col];
            }
        }
        return;
    }

#if defined(USE_A2A)
    ff_a2a sk;
    Master master(matrix, n, pardegree);

    std::vector<ff_node*> LW{&master};
    std::vector<ff_node*> RW;
    for(int i=0;i<pardegree;++i)
        RW.push_back(new Worker(matrix,n));

    sk.add_firstset(LW);
    sk.add_secondset(RW,true);

#else
    ff_farm sk;
	Master master(matrix, n, pardegree);
	sk.add_emitter(&master);
    std::vector<ff_node*> W;
	for(int i=0;i<pardegree;++i)
		W.push_back(new Worker(matrix,n));
	sk.add_workers(W);
	sk.cleanup_workers();
#endif
    // we consider the case with NO_DEFAULT_MAPPING
    // sk.no_mapping();
    //sk.blocking_mode();

    for (int k = 1; k < n; ++k) {
        master.setNextIteration();
        sk.run_and_wait_end();
    }

}

void printMatrix(double **matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%.2lf ",  matrix[i][j]);
        }
        std::cout << std::endl;
    }
}


int main(int argc, char **argv) {
    // Take matrix size as input
    /*if (argc < 3) {
        std::cout << "use: %s size pardegree" << std::endl;
        exit(-1);
    }*/
    long n = 2;
    long pardegree = 48;
    try {
        n = std::stol(argv[1]);
        pardegree = std::stol(argv[2]);

        if (pardegree <=0) {
            std::cout << "pardegree should be greater than zero\n";
            return -1;
        }

        // Allocate the matrix
        auto **matrix = new double *[n];
        for (int i = 0; i < n; ++i) {
            matrix[i] = new double[n]();
        }

        // Intialize the diagonal with random values
        initMatrix(matrix, n);
        // Compute the wavefront
        auto start = std::chrono::high_resolution_clock::now();
        wavefront(matrix, n, std::min(pardegree,n));
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Matrix size: " << n << " Time taken: " << duration.count() / 1000.0 << "ms Final cell: " << matrix[0][n-1] << "\n";
        std::cout << "Time per iteration: " << duration.count() / 1000.0 / n << std::endl;
        //printMatrix(matrix, n);
    } catch (std::exception &e) {
        std::cout << "Arguments must be positive integer numbers!" << std::endl;
    }
    return 0;
}
