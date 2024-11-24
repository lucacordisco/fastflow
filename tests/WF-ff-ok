#include <vector>
#include <thread>
#include <cmath>
#include <iostream>

#include <ff/ff.hpp>

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
    for(int i = 1; i < len; i++) {
        prod += matrix[row][col - i] * matrix[col][row + i];
    }
    return prod;
}

void initMatrix(double **matrix, int n) {
	for (int i = 0; i < n; i++)
		matrix[i][i] = ((double) (i+1)) / n;
}

void compute_iter(double **matrix, int start,int stop, int iter) {
	for(int m=start;m<stop;++m) {
		const int row=m;
		const int col=row+iter;
		double prod = product(matrix, row, col);
		
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

	// sends tasks to workers 
	void send_tasks(int &start, int &stop) {
		const size_t chunksize = (n-iter) / pardegree;
		ssize_t more           = (n-iter) - chunksize*pardegree;
		
		for(int i=0; i<(pardegree-1); ++i) {
			start = stop;
			stop  = start + chunksize + (more>0 ? 1:0);
			--more;
			
			taskv[i] = {start, stop, iter};

			//std::printf("Emitter sending to Worker%d task {%d,%d,%d}\n", i, taskv[i].start,taskv[i].stop,taskv[i].iter);
			ff_send_out_to(&taskv[i], i);
		}
		nsent=(pardegree-1);
		if ((n-iter)==pardegree)  {
			ff_send_out_to(EOS, pardegree-1);
			--pardegree;
		}		
		start = stop;
		stop  = start + chunksize + (more>0 ? 1:0);
	}
		
	task_t* svc(task_t* in) {
		int start=0;
		int stop=0;
		if (!in) {
			send_tasks(start,stop);
			//std::printf("Emitter(START) computing task {%d,%d,%d}\n", start,stop,iter);
			compute_iter(matrix, start,stop,iter);			
			return pardegree?GO_ON:EOS;
		}

		if (++nreceived == nsent) {
			++iter;
			nreceived=0;
			nsent=0;

			if (pardegree==0) {
				assert(iter == (n-1));
				compute_iter(matrix, 0, n-iter, iter);				
				return EOS;
			}
				
			send_tasks(start,stop);
			//std::printf("Emitter computing task {%d,%d,%d} nw=%d\n",start,stop,iter,pardegree);
			compute_iter(matrix, start,stop,iter);
		}
		return GO_ON;
	}
	
	double **matrix; // the matrix
	int n;           // number of columns/rows
	int pardegree;   // initial number of workers
	int iter=1;      // current iterations (from 1 to n)

	int nsent=0;     // n. of tasks sent to workers
	int nreceived=0; // n. of tasks got back
	std::vector<task_t> taskv;	 // memory needed for the tasks
};
struct Worker:ff_minode_t<task_t> {
	Worker(double **matrix, long n):matrix(matrix),n(n) {	
	}
	task_t* svc(task_t* in) {
		//std::printf("Worker%ld received {%d,%d,%d}\n", get_my_id(), in->start,in->stop,in->iter);
		compute_iter(matrix, in->start,in->stop,in->iter);
		return in; // notify, we've completed the task
	}
	double **matrix;  // the matrix
	long n;           // n. of columns/rows
};

void wavefront(double **matrix, int n, int pardegree, bool blocking, bool nomapping) {
	if (pardegree==1) {  // sequential computation
		for (int k = 1; k < n; ++k) {
			for (int m = 0; m < n - k; ++m) {
				int row = m;
				int col = row + k;
				
				double prod = product(matrix, row, col);
				
				// Save the result in the corresponding cell
				matrix[row][col] = std::cbrt(prod);
				// Save the result also in the transposed cell to speedup next computations using cached row accesses
				matrix[col][row] = matrix[row][col];
			}
		}
		return;
	}
	
	ff_farm farm;
	Master m(matrix, n, pardegree);
	farm.add_emitter(&m);
    std::vector<ff_node*> W;
	for(int i=0;i<pardegree;++i)
		W.push_back(new Worker(matrix,n));
	farm.add_workers(W);
	farm.wrap_around();
	farm.cleanup_workers();
	farm.blocking_mode(blocking);
	if (nomapping) farm.no_mapping();
	farm.run_and_wait_end();
}

void printMatrix(double **matrix, int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}


int main(int argc, char **argv) {
    // Take matrix size as input
    if (argc < 3) {
		std::printf("use: %s size pardegree [blocking] [nomapping]\n", argv[0]);
		std::printf("example: %s 1024 10 1 1    # means BLOCKING_MODE and NO_DEFAULT_MAPPING\n", argv[0]); 
		return -1;
    }
	bool blocking  = false;
	bool nomapping = false;
	long n         = std::stol(argv[1]);
	long pardegree = std::stol(argv[2]);

	if (argc == 4) {
		blocking = std::stoi(argv[3])?true:false;
	}
	if (argc == 5) {
		nomapping = std::stoi(argv[4])?true:false;
	}
	
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
	wavefront(matrix, n, std::min(pardegree,n), blocking, nomapping);
	auto end = std::chrono::high_resolution_clock::now();
	
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "Matrix size: " << n << " Time taken: " << duration.count() / 1000.0 << "ms Final cell: " << matrix[0][n-1] << "\n";
	
    return 0;
}
