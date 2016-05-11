#include <stdio.h>


void distribute(int rank, int size, double* mt, int H, int W)
{
  auto m = (H-2)%size;
  auto d = (H-2)/size;
  
  auto r0 = m + d + 1; // number of row skipped
  Matrix matrix;
  if(0==rank) {
    matrix.val = mt;
    Matrix matrixsize;
    matrixsize.h = d + 2;
    matrixsize.w = W;
    
    for( int i = 1 ; i < size; i++) {
      MPI_Send(&matrixsize, sizeof(Matrix), MPI_BYTE, i, 1, MPI_COMM_WORLD);
      MPI_Send(mt + ((r0-1 +(i-1)*d) * W, matrixsize.h * W, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
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
