#include <stdio.h>


void distribute(int rank, int size, double* mt, int H, int W)
{
  auto m = (H-2)%size;
  auto d = (H-2)/size;

  //make sure d > 2
  
  auto r0 = m + d + 1; // number of row skipped
  Matrix matrix;
  if(0==rank) {
    matrix.val = mt;
    Matrix matrixsize;
    matrixsize.h = d + 2;
    matrixsize.w = W;
    
    for( int i = 1 ; i < size; i++) {
      MPI_Send(&matrixsize, sizeof(Matrix), MPI_BYTE, i, 1, MPI_COMM_WORLD);
      MPI_Send(mt + (r0-1 +(i-1)*d) * W, matrixsize.h * W, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
    }
  }
  else {
  
    Matrix matrix;
    MPI_Status  status;
    MPI_Recv(&matrix, sizeof(Matrix), MPI_BYTE, 0, 1, MPI_COMM_WORLD, &status);
    matrix.val = new double[matrix.w * matrix.h];
    MPI_Recv(matrix.val, matrix.w * matrix.h, MPI_DOUBLE, i, 2, MPI_COMM_WORLD); 
  }
  return matrix;
}
void gather(int rank, int size, double* mt, int H, int W)
{
  
}

  void exchange_boundary(int rank, int size, double* matrix, int W, int H)
  {
    if(rank < size - 1){
      MPI_Sendrecv(matrix + (H - 2) * W, W, MPI_DOUBLE, rank + 1, 5, matrix + (H - 1)*W, W, MPI_DOUBLE, rank + 1, 4, MPI_COMM_WORLD, &status);
    }

    if(rank > 0) {

      //SEND TO I-1
      //RECV FORM I-1
      MPI_Sendrecv(matrix + W, W, MPI_DOUBLE, rank - 1, 4, matrix, W, MPI_DOUBLE, rank - 1, 5, MPI_COMM_WORLD, &status);
    }
  }

  void calc(int rank, int isblack, double* matrix, int W, int H)
  {
    int colorfactor = rank == 0 ? isblack : (isblack + r0 -1 + (rank -1)*H) % 2;
    for(auto j = 1; i < H - 1; j++)  {
      for(auto i = 1; i< W - 1; i+=2){
	auto colorfactori = (colorfactor + j +i ) % 2;
	auto columni = j*W + i + colorfactori ; 
	matrix[columni] = 0.25 * (
				  matrix[ columni - W ] + matrix[columni + W ] + matrix[ columni -1] + matrix[columnti +1]);
      }
    }
  }
