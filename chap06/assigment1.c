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

    MPI_Type_vector(10, 1, 10, MPI_FLOAT, &column_mpi_t);
    MPI_Type_commit(&column_mpi_t);

    if (my_rank == 0)
    {
        int source;
        for (i = 0; i < 10; i++)
            for (j = 0; j < 10; j++)
                A[i][j] = 0.0;

        for (source = 1; source < 10; source++)
        {
            MPI_Recv(&(A[0][source]), 10, MPI_FLOAT, source, 0,
                     MPI_COMM_WORLD, &status);
        }

        for (i = 0; i < 10; i++)
        {
            for (j = 0; j < 10; j++)
                printf("%3.1f ", A[i][j]);
            printf("\n");
        }
        printf("\n");
    }
    else
    {
        for (j = 0; j < 10; j++)
            A[j][my_rank] = (float)5;

        printf("Sending %3.1f - %d\n", A[0][my_rank], my_rank);
        MPI_Send(&(A[0][my_rank]), 1, column_mpi_t, 0, 0,
                 MPI_COMM_WORLD);
    }

    MPI_Finalize();
} /* main */
