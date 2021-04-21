#include <stdio.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <mpi.h>
#include <cmath>
#include <iomanip>

//globals
const int Tag = 0;
const int itter_x = 6;
const double X = 1;
const int itter_t = 6;
const double T = 2;
const double a = 1;

const int show_proc = 0;


double f_function(double t, double x);            //x+t
double time_initial_function(double t);           //t
double coord_initial_function(double x);          //x

void show_table(double** u);


int main(int argc, char* argv[])
{

    int rc;
    MPI_Status status;
    MPI_Request request;

    if (rc = MPI_Init(&argc, &argv))
    {
        std::cout << "Ошибка запуска, выполнение остановлено " << std::endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    int rank;
    int procs;

    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double h = X / (itter_x - 1);
    double tau = T / (itter_t - 1);

    double** u = new double* [itter_t];
    for (int counter = 0; counter < itter_t; counter++)
        u[counter] = new double[itter_x]{};


    for (int counter = 0; counter < itter_x; counter++)
        u[0][counter] = time_initial_function(counter * h);

    for (int counter = 0; counter < itter_t; counter++)
        u[counter][0] = coord_initial_function(counter * tau);


    double* buf = new double[itter_x]{0};
    int time = rank;
    while (time < itter_t - 1)
    {
        for (int coord = 0; coord < itter_x - 1; coord++)
         {
            if (time != 0)
            {
                MPI_Recv(buf, 1, MPI_DOUBLE, (time - 1) % procs, Tag, MPI_COMM_WORLD, &status);
                u[time][coord + 1] = buf[0];
            }
            u[time + 1][coord + 1] = 1.0 / (1.0 / 2 / h + 1.0 / 2 / tau) * (f_function((time + 1) * tau
                - 1.0 / 2 * tau, (coord + 1) * h - 1.0 / 2 * h) - 1.0 / 2 / tau * (u[time + 1][coord] - u[time][coord] - u[time][coord + 1]) - 1.0 / 2 / h * (u[time][coord + 1] - u[time + 1][coord] - u[time][coord]));

            buf[0] = u[time + 1][coord + 1];
            //MPI_Send(buf, 1, MPI_DOUBLE, (rank + 1) % procs, Tag, MPI_COMM_WORLD);
            MPI_Isend(buf, 1, MPI_DOUBLE, (rank + 1) % procs, Tag, MPI_COMM_WORLD, &request);
         }
        time += procs;
    }

    for (time = 2; time < itter_t - 1; time++)
    {
        if (time % procs != 0)
        {
            if (rank == time % procs)
                //MPI_Send(u[time], itter_x, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD);
                MPI_Isend(u[time], itter_x, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, &request);

            if (rank == 0)
               MPI_Recv(u[time], itter_x, MPI_DOUBLE, time % procs, Tag, MPI_COMM_WORLD, &status);
        }
    }
    if ((itter_t - 2) % procs != 0)
    {
        if (rank == (itter_t - 2) % procs)
            //MPI_Send(u[itter_t - 1], itter_x, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD);
            MPI_Isend(u[itter_t - 1], itter_x, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, &request);
        if (rank == 0)
            MPI_Recv(u[itter_t - 1], itter_x, MPI_DOUBLE, (itter_t - 2) % procs, Tag, MPI_COMM_WORLD, &status);
    }


    if (rank == 0)
        show_table(u);

    MPI_Finalize();
    return 0;
}


double f_function(double t, double x)
{
    return x + t;
}


double time_initial_function(double t)
{
    return t;
}


double coord_initial_function(double x)
{
    return x;
}


void show_table(double** u)
{
    int time;
    std::cout << std::endl << "\t\tTable u[time][coord]" << std::endl << std::endl << "time" << std::endl << std::endl;
    for (time = 0; time < itter_t; time++)
    {
        for (int coord = 0; coord < itter_x; coord++)
            printf("%5.2lf ", u[itter_t - 1 - time][coord]);
        if (time == itter_t - 1)
            std::cout << "\t coord";
        std::cout << std::endl;
    }
}
