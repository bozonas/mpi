#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <stdlib.h>

typedef struct
{
  int p;             /* Total number of processes    */
  MPI_Comm comm;     /* Communicator for entire grid */
  MPI_Comm row_comm; /* Communicator for my row      */
  MPI_Comm col_comm; /* Communicator for my col      */
  int q;             /* Order of grid                */
  int my_row;        /* My row number                */
  int my_col;        /* My column number             */
  int my_rank;       /* My rank in the grid comm     */
} GRID_INFO_T;

#define MAX 65536
typedef struct
{
  int n_bar;
#define Order(A) ((A)->n_bar)
  float entries[MAX];
#define Entry(A, i, j) (*(((A)->entries) + ((A)->n_bar) * (i) + (j)))
} LOCAL_MATRIX_T;

/*********************************************************/
void Print_matrix(
    char *title /* in  */,
    LOCAL_MATRIX_T *local_A /* out */,
    GRID_INFO_T *grid /* in  */,
    int n /* in  */)
{
  int mat_row, mat_col;
  int grid_row, grid_col;
  int source;
  int coords[2];
  float *temp;
  MPI_Status status;

  if (grid->my_rank == 0)
  {
    temp = (float *)malloc(Order(local_A) * sizeof(float));
    printf("%s\n", title);
    for (mat_row = 0; mat_row < n; mat_row++)
    {
      grid_row = mat_row / Order(local_A);
      coords[0] = grid_row;
      for (grid_col = 0; grid_col < grid->q; grid_col++)
      {
        coords[1] = grid_col;
        MPI_Cart_rank(grid->comm, coords, &source);
        if (source == 0)
        {
          for (mat_col = 0; mat_col < Order(local_A); mat_col++)
            printf("%4.1f ", Entry(local_A, mat_row, mat_col));
        }
        else
        {
          MPI_Recv(temp, Order(local_A), MPI_FLOAT, source, 0,
                   grid->comm, &status);
          for (mat_col = 0; mat_col < Order(local_A); mat_col++)
            printf("%4.1f ", temp[mat_col]);
        }
      }
      printf("\n");
    }
    free(temp);
  }
  else
  {
    for (mat_row = 0; mat_row < Order(local_A); mat_row++)
      MPI_Send(&Entry(local_A, mat_row, 0), Order(local_A),
               MPI_FLOAT, 0, 0, grid->comm);
  }

} /* Print_matrix */

/*********************************************************/
LOCAL_MATRIX_T *Local_matrix_allocate(int local_order)
{
  LOCAL_MATRIX_T *temp;

  temp = (LOCAL_MATRIX_T *)malloc(sizeof(LOCAL_MATRIX_T));
  return temp;
} /* Local_matrix_allocate */

/*********************************************************/
void Free_local_matrix(
    LOCAL_MATRIX_T **local_A_ptr /* in/out */)
{
  free(*local_A_ptr);
} /* Free_local_matrix */

void Fill_matrix(
    LOCAL_MATRIX_T *local_A /* out */,
    GRID_INFO_T *grid /* in  */,
    int n /* in  */)
{
  int mat_row, mat_col;
  int grid_row, grid_col;
  int dest;
  int coords[2];
  float *temp;
  MPI_Status status;

  if (grid->my_rank == 0)
  {
    temp = (float *)malloc(Order(local_A) * sizeof(float));
    fflush(stdout);
    for (mat_row = 0; mat_row < n; mat_row++)
    {
      grid_row = mat_row / Order(local_A);
      coords[0] = grid_row;
      for (grid_col = 0; grid_col < grid->q; grid_col++)
      {
        coords[1] = grid_col;
        MPI_Cart_rank(grid->comm, coords, &dest);
        if (dest == 0)
        {
          for (mat_col = 0; mat_col < Order(local_A); mat_col++)
          {
            float *pntr = (local_A->entries) + mat_row * Order(local_A) + mat_col;
            *pntr = 2;
          }
        }
        else
        {
          for (mat_col = 0; mat_col < Order(local_A); mat_col++)
          {
            float *pntr = temp + mat_col;
            *pntr = 2;
          }
          MPI_Send(temp, Order(local_A), MPI_FLOAT, dest, 0,
                   grid->comm);
        }
      }
    }
    free(temp);
  }
  else
  {
    for (mat_row = 0; mat_row < Order(local_A); mat_row++)
      MPI_Recv(&Entry(local_A, mat_row, 0), Order(local_A),
               MPI_FLOAT, 0, 0, grid->comm, &status);
  }
}

/*********************************************************/
void Setup_grid(
    GRID_INFO_T *grid /* out */)
{
  int old_rank;
  int dimensions[2];
  int wrap_around[2];
  int coordinates[2];
  int free_coords[2];

  /* Set up Global Grid Information */
  MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

  /* We assume p is a perfect square */
  grid->q = (int)sqrt((double)grid->p);
  dimensions[0] = dimensions[1] = grid->q;

  /* We want a circular shift in second dimension. */
  /* Don't care about first                        */
  wrap_around[0] = wrap_around[1] = 1;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions,
                  wrap_around, 1, &(grid->comm));
  MPI_Comm_rank(grid->comm, &(grid->my_rank));
  MPI_Cart_coords(grid->comm, grid->my_rank, 2,
                  coordinates);
  grid->my_row = coordinates[0];
  grid->my_col = coordinates[1];

  /* Set up row communicators */
  free_coords[0] = 0;
  free_coords[1] = 1;
  MPI_Cart_sub(grid->comm, free_coords,
               &(grid->row_comm));

  /* Set up column communicators */
  free_coords[0] = 1;
  free_coords[1] = 0;
  MPI_Cart_sub(grid->comm, free_coords,
               &(grid->col_comm));
} /* Setup_grid */

void Fill_vector(
    float local_x[] /* out */,
    int n /* in  */,
    GRID_INFO_T *grid /* in  */)
{
  int i;
  float temp[MAX];
  int n_bar;

  n_bar = n / grid->p;

  if (grid->my_rank == 0)
  {
    for (i = 0; i < n; i++)
    {
      temp[i] = 2;
    }
  }
  MPI_Scatter(temp, n_bar, MPI_FLOAT, local_x, n_bar, MPI_FLOAT,
              0, grid->comm);
}

/**********************************************************************/
void Print_vector(
    char *title /* in */,
    float local_x[] /* in */,
    int n /* in */,
    GRID_INFO_T *grid)
{

  int i;
  float temp[MAX];
  int n_bar;

  n_bar = n / grid->p;

  MPI_Gather(local_x, n_bar, MPI_FLOAT, temp, n_bar, MPI_FLOAT,
             0, grid->comm);

  if (grid->my_rank == 0)
  {
    printf("%s\n", title);
    for (i = 0; i < n; i++)
      printf("%4.1f ", temp[i]);
    printf("\n");
  }
} /* Print_vector */

void Multiply(int n /* in  */,
              GRID_INFO_T *grid /* in  */,
              LOCAL_MATRIX_T *local_A /* in  */,
              float local_x[],
              float local_out[])
{
  float global_x[MAX];
  int i;
  int j;
  int n_bar;

  n_bar = n / grid->p;

  MPI_Allgather(local_x, n_bar, MPI_FLOAT, global_x,
                n_bar, MPI_FLOAT, grid->comm);

  for (i = 0; i < n; i++)
  {
    local_out[i] = 0.0;
    for (j = 0; j < n; j++)
    {
      local_out[i] += Entry(local_A, i, j) * global_x[j];
    }
  }
}

main(int argc, char *argv[])
{
  int p;
  int my_rank;
  GRID_INFO_T grid;
  LOCAL_MATRIX_T *local_A;
  float local_out[MAX];
  float result[MAX];
  int n;
  int n_bar;
  int i;
  int j;


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  Setup_grid(&grid);
  if (my_rank == 0)
  {
    printf("What's the order of the matrices?\n");
    scanf("%d", &n);
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  n_bar = n / grid.q;

  local_A = Local_matrix_allocate(n_bar);
  Order(local_A) = n_bar;
  Fill_matrix(local_A, &grid, n);
  Print_matrix("Generated matrix:", local_A, &grid, n);

  float local_x[MAX];
  Fill_vector(local_x, n, &grid);
  Print_vector("Generated vector:", local_x, n, &grid);

  Multiply(n, &grid, local_A, local_x, local_out);

  MPI_Gather(local_out, n_bar/2, MPI_FLOAT, result, n_bar/2, MPI_FLOAT,
             0, grid.comm);

  if (my_rank == 0)
  {
    printf("Result:\n");
    for (i = 0; i < n; i++)
      printf("%4.1f ", result[i]);
    printf("\n");
  }

  Free_local_matrix(&local_A);

  MPI_Finalize();
} /* main */
