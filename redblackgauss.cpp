#include "mpi.h"
#include "stdio.h"
#include <random>
#include <time.h>
#include <fstream>

#define TAG_PUSHBACK 2910

struct Matrix
{
	int w;
	int h;
	double* val;
};
int bigH = 0;
void print(double* matrix, int H, int W)
{
	printf("Ma tran [%dx%d]\n", H, W);
	for (auto i = 0; i < H; i++)
	{
		for (auto j = 0; j < W; j++)
		{
			printf("%4.2f ", matrix[i* W + j]);
		}
		printf("\n");
	}
}

double* random(int w, int h)
{
	//srand(time(NULL));
	auto matrix = new double[w*h * 10];

	for (auto i = 0; i < h; i++)
	{
		for (auto j = 0; j < w; j++)
		{
			matrix[i*h + j] = ((double)(rand() % 50)) / (rand() % 50 + 1) * (rand() % 50);
		}
	}
	return matrix;
}

Matrix read()
{
	std::ifstream in;
	int height, width;
	double val = -1;

	in.open("input.txt");
	in >> height >> width;
	Matrix mt;
	double* matrix = new double[height*width];
	for (auto i = 0; i < height; i++)
		for (auto j = 0; j < width; j++)
		{
			in >> val;
			matrix[i* width + j] = val;
		}
	in.close();
	mt.val = matrix;
	mt.w = width;
	mt.h = height;
	return mt;
}

//split data from process 0 to another process
Matrix init(int rank, int size)
{
	Matrix matrix;

	if (rank == 0)
	{
		matrix = read();// random(W, H);// read();
		auto m = (matrix.h - 2) % size;
		auto d = (matrix.h - 2) / size;
		bigH = matrix.h;
		//make sure d > 2
		auto r0 = m + d + 1;
		matrix.h = r0 + 1;
		Matrix matrixsize;
		matrixsize.h = d + 2;
		matrixsize.w = matrix.w;

		for (int i = 1; i < size; i++) {
			MPI_Send(&matrixsize, sizeof(Matrix), MPI_BYTE, i, 1, MPI_COMM_WORLD);
			MPI_Send(matrix.val + (r0 - 1 + (i - 1)*d) * matrix.w, matrixsize.h * matrix.w, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Status status;
		MPI_Recv(&matrix, sizeof(Matrix), MPI_BYTE, 0,1, MPI_COMM_WORLD, &status);
		matrix.val = new double[matrix.w * matrix.h];
		MPI_Recv(matrix.val, matrix.w* matrix.h, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
	}
	return matrix;
}

void gather(int rank, int size, double* mt, int H, int W, int bigH )
{
	if (rank == 0)
	{
		// W = H, because big matrix is a square matrix
		auto d = (bigH - 2) / size;
		auto r0 = H - 1;
		MPI_Status status;
		for (int i = 1; i < size; i++)
			MPI_Recv(mt + (r0 - 1 + (i - 1)*d) * W, (d + 2)*W, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &status);
	}
	if (rank > 0)
	{
		MPI_Send(mt, H * W, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
	}
}

void exchange_boundary(int rank, int size, double* matrix, int W, int H)
{
	MPI_Status status;
	if (rank < size - 1)
		MPI_Sendrecv(matrix + (H - 2) * W, W, MPI_DOUBLE, rank + 1, 5, matrix + (H - 1)*W, W, MPI_DOUBLE, rank + 1, 4, MPI_COMM_WORLD, &status);

	if (rank > 0)
		//SEND TO I-1
		//RECV FORM I-1
		MPI_Sendrecv(matrix + W, W, MPI_DOUBLE, rank - 1, 4, matrix, W, MPI_DOUBLE, rank - 1, 5, MPI_COMM_WORLD, &status);

}

int calc(int r0, int rank, int isblack, double* matrix, int W, int H)
{
	int maxc = 0;
	int colorfactor = rank == 0 ? isblack : (isblack + r0 - 1 + (rank - 1)*H) % 2;
	for (auto j = 1; j < H - 1; j++)  {
		for (auto i = 1; i < W - 1; i += 2){
			auto colorfactori = (colorfactor + j + i) % 2;
			auto columni = j*W + i + colorfactori;
			auto newvalue = 0.25 * (matrix[columni - W] + matrix[columni + W] + matrix[columni - 1] + matrix[columni + 1]);
			if (fabs(matrix[columni] - newvalue) > maxc) maxc = matrix[columni] - newvalue;
			matrix[columni] = newvalue;
		}
	}
	return maxc;
}

int gauss(int rank, int size, int W, int H, double expeps, double* matrix, int MAXSTEPS, int r0)
{
	double  eps;
	int steps;
	//if (rank == 0)
	//	printf("____1-%d:", rank),print(matrix, W,H);
	auto maxc = 0;
	for (steps = 0; steps < MAXSTEPS; steps++)
	{
		//tinh red
		//MPI_Barrier(MPI_COMM_WORLD);
		//if (rank == 0)printf("____2-%d:", rank), print(matrix, W);
		//assign(rank, size, W, H, matrix);
		//if (rank == 0)printf("____3-%d:", rank), print(matrix, W);
		//MPI_Barrier(MPI_COMM_WORLD);
		//auto redcureps = update(0, rank, W, H, matrix);

		maxc = calc(r0, rank, 0, matrix, W, H);

		exchange_boundary(rank, size, matrix, W, H);

		auto maxblackc = calc(r0, rank, 1, matrix, W, H);
		if (maxblackc > maxc) maxc = maxblackc;

		//get max epsilon
		MPI_Allreduce(&maxc, &eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (eps <= expeps) break;
	}
	//if (rank == 0)
	//	printf("____1-%d:", rank), print(matrix, W, H);
	return steps;
}



//mpiexec /np {numberofprecoess} ltss
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//initializing
	int MAXSTEPS = 60000;
	double EPSILON = 0.00001;

	//split data from process 0 to another process
	Matrix matrix = init(rank, size);

	//start running
	double starttime = MPI_Wtime();

	Matrix bigmatrix = matrix;
	MPI_Bcast(&bigmatrix, sizeof(Matrix), MPI_BYTE, 0, MPI_COMM_WORLD);

	auto r0 = bigmatrix.h - 1;

	int steps = gauss(rank, size, matrix.w, matrix.h, EPSILON, matrix.val, MAXSTEPS, r0);

	gather(rank, size, matrix.val, matrix.h, matrix.w, bigH);
	MPI_Barrier(MPI_COMM_WORLD);
	double deltaT = MPI_Wtime() - starttime;

	if (rank == 0)
	{
		print(matrix.val, bigH, matrix.w);
		printf("time  : %f\r\n", deltaT);
	}
	//freeing resource
	delete[] matrix.val;
	MPI_Finalize();
	return 0;
}
