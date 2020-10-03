#include <stdio.h>
#include "mpi.h"

#define MAX_ORDER 100

typedef float MATRIX_T[MAX_ORDER][MAX_ORDER];

main(int argc, char *argv[])
{
  int my_rank; /* My process rank           */
  int p;       /* The number of processes   */
  MATRIX_T A;
  float x[MAX_ORDER];
  float y[MAX_ORDER];
  int m;
  float val;

  /* Let the system do what it needs to start up MPI */
  MPI_Init(&argc, &argv);

  /* Get my process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  void CreateMatrix(MATRIX_T A, int m);
  void CreateVector(float v[], int m);
  void Print_matrix(MATRIX_T A, int n);
  void Print_vector(float y[], int n);
  void Get_data(int *m_ptr, int my_rank);
  void Global_sum(float *val, int my_rank, int p);

  val = 1;
  printf("my_rank: %d val: %f\n", my_rank, val);

  Global_sum(&val, my_rank, p);

  /* Shut down MPI */
  MPI_Finalize();
}

void Global_sum(float *val, int my_rank, int p)
{
  float tmp_val;

  MPI_Reduce(val, &tmp_val, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_rank == 0)
  {
    *val = tmp_val;
    printf("%f\n", *val);
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

void CreateMatrix(MATRIX_T A, int m)
{
  int i;
  int j;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      A[i][j] = 1;
}

void CreateVector(float v[], int m)
{
  int i;
  for (i = 0; i < m; i++)
    v[i] = 1;
}

void Print_matrix(
    MATRIX_T A,
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

void Print_vector(
    float y[],
    int n)
{
  int i;

  printf("Vector is \n");
  for (i = 0; i < n; i++)
    printf("%4.1f ", y[i]);
  printf("\n");
}
