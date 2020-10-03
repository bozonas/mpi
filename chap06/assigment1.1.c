/* send_col_to_row.c -- send column 1 of a matrix on process 0 to row 1
 *     on process 1.
 *
 * Input: none
 * Output: The row received by process 1.
 *
 * Note:  This program should only be run with 2 processes
 *
 * See Chap 6., pp. 98 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"

main(int argc, char *argv[])
{
    int p;
    int my_rank;
    float A[10][10];
    MPI_Status status;
    MPI_Datatype column_mpi_t;
    int i, j;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int n = 10;
    int batch = n*n/p;

    MPI_Type_vector(batch, 1, 1, MPI_FLOAT, &column_mpi_t);
    MPI_Type_commit(&column_mpi_t);

    if (my_rank == 0)
    {
        int source;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                A[i][j] = 0.0;

        for (source = 1; source < p; source++)
        {
            MPI_Send(&(A[source][0]), 1, column_mpi_t, source, 0,
                     MPI_COMM_WORLD);
        }

        for (source = 1; source < p; source++)
        {
            MPI_Recv(&(A[source][0]), 10, MPI_FLOAT, source, 0,
                     MPI_COMM_WORLD, &status);
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                printf("%3.1f ", A[i][j]);
            printf("\n");
        }
        printf("\n");
    }
    else
    {
        MPI_Recv(&(A[my_rank][0]), 10, MPI_FLOAT, 0, 0,
                     MPI_COMM_WORLD, &status);

        for (j = 0; j < n; j++)
            A[my_rank][j] = (float)j;

        // for (j = 0; j < 10; j++)
        //     printf("%3.1f ", A[my_rank][j]);
        // printf("Sending %3.1f - %d\n", A[my_rank][1], my_rank);
        MPI_Send(&(A[my_rank][0]), 1, column_mpi_t, 0, 0,
                 MPI_COMM_WORLD);
        // }
    }

    MPI_Finalize();
} /* main */
