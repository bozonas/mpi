#include <stdio.h>
#include "mpi.h"

#define MAX_ORDER 100

typedef float MATRIX_T[MAX_ORDER][MAX_ORDER];

main(int argc, char *argv[])
{
  int my_rank; /* My process rank           */
  int p;       /* The number of processes   */
  MATRIX_T A;
  MATRIX_T B;
  int m;
  float local_x[MAX_ORDER];
  float global_x[MAX_ORDER];
  float local_y[MAX_ORDER][MAX_ORDER];

  /* Let the system do what it needs to start up MPI */
  MPI_Init(&argc, &argv);

  /* Get my process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  void CreateMatrix(MATRIX_T A, int m, int local_m);
  void Print_matrix(float A[MAX_ORDER][MAX_ORDER], int n);
  void Get_data(int *m_ptr, int my_rank);
  void Parallel_matrix_matrix_prod(MATRIX_T local_A,
                                   MATRIX_T local_B, int m, int local_m,
                                   float local_x[], float global_x[],
                                   float local_y[][MAX_ORDER],
                                   int my_rank, int p);

  Get_data(&m, my_rank);

  int local_m = m / p;
  float tempA[MAX_ORDER][MAX_ORDER];
  float tempB[MAX_ORDER][MAX_ORDER];

  // Print Input
  CreateMatrix(A, m, local_m);
  CreateMatrix(B, m, local_m);
  MPI_Gather(A, local_m * MAX_ORDER, MPI_FLOAT, tempA, local_m * MAX_ORDER, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(B, local_m * MAX_ORDER, MPI_FLOAT, tempB, local_m * MAX_ORDER, MPI_FLOAT, 0, MPI_COMM_WORLD);
  if (my_rank == 0)
  {
    Print_matrix(tempA, m);
    Print_matrix(tempB, m);
  }

  Parallel_matrix_matrix_prod(A, B, m, local_m, local_x, global_x, local_y, my_rank, p);

  float tempC[MAX_ORDER][MAX_ORDER];
  MPI_Gather(local_y, local_m * MAX_ORDER, MPI_FLOAT, tempC, local_m * MAX_ORDER, MPI_FLOAT, 0, MPI_COMM_WORLD);
  if (my_rank == 0)
  {
    printf("Result:\n");
    Print_matrix(tempC, m);
  }

  /* Shut down MPI */
  MPI_Finalize();
}

void Parallel_matrix_matrix_prod(MATRIX_T local_A,
                                 MATRIX_T local_B, int m, int local_m,
                                 float local_x[], float global_x[],
                                 float local_y[][MAX_ORDER],
                                 int my_rank, int p)
{
  int i, j, l;

  for (l = 0; l < m; l++)
  {
    for (i = 0; i < local_m; i++)
      local_x[i] = local_B[i][l];
    MPI_Allgather(local_x, local_m, MPI_FLOAT, global_x, local_m, MPI_FLOAT, MPI_COMM_WORLD);
    for (i = 0; i < local_m; i++)
    {
      local_y[i][l] = 0.0;
      for (j = 0; j < m; j++)
        local_y[i][l] = local_y[i][l] + local_A[i][j] * global_x[j];
    }
  }
}

void Get_data(
    int *m_ptr /* out */,
    int my_rank /* in  */)
{
  if (my_rank == 0)
  {
    printf("Enter the order of the matrix (m x m) and vector\n");
    scanf("%d", m_ptr);
  }
  MPI_Bcast(m_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void CreateMatrix(MATRIX_T A, int m, int local_m)
{
  int i;
  int j;
  MATRIX_T temp = {
      {0},
  };
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      temp[i][j] = 1;

  MPI_Scatter(temp, local_m * MAX_ORDER, MPI_FLOAT, A,
              local_m * MAX_ORDER, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void Print_matrix(
    float A[MAX_ORDER][MAX_ORDER],
    int n)
{
  int i, j;
  printf("Matrix is \n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
      printf("%4.1f ", A[i][j]);
    printf("\n");
  }

  printf("\n");
}
