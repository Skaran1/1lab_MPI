#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <chrono>
#include <vector>
#include <cmath>

//globals
const int p = 2;
const int Tag = 0;


int main(int argc, char* argv[])
{

    int rc;
    MPI_Status status;

    if (rc = MPI_Init(&argc, &argv))
    {
        std::cout << "Ошибка запуска, выполнение остановлено " << std::endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    int rank;
    int numprocs;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); \

    int size = 3;

    int* array = new int[size];

    if (rank == 1)
        for (int i = 0; i < size; i++)
            array[i] = i + 3;


    auto start = std::chrono::high_resolution_clock::now();

    if (rank == 1)
        MPI_Send(array, size, MPI_INT, 0, Tag, MPI_COMM_WORLD);

    if (rank == 0)
        MPI_Recv(array, size, MPI_INT, 1, Tag, MPI_COMM_WORLD, &status);


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> delay = end - start;

    if (rank == 0)
    {
        std::cout << "Delay = " << delay.count() << "sec" << std::endl;
        for (int i = 0; i < size; i++)
            std::cout << "Array[" << i << "] = " << array[i] << std::endl;
    }

    MPI_Finalize();
    return 0;
}
