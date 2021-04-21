#include <stdio.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <iomanip>

//globals
const int Tag = 0;
const int itter_x = 6;
const double X = 1;
const int itter_t = 6;
const double T = 2;
const double a = 1;


double f_function(double t, double x);            //x+t
double time_initial_function(double t);           //t
double coord_initial_function(double x);          //x

void show_table(double** u);


int main(int argc, char* argv[])
{
    double h = X / (itter_x - 1);
    double tau = T / (itter_t - 1);

    double** u = new double*[itter_t];
    for (int counter = 0; counter < itter_t; counter++)
        u[counter] = new double [itter_x];


    for (int counter = 0; counter < itter_x; counter++)
        u[0][counter] = time_initial_function(counter * h);

    for (int counter = 0; counter < itter_t; counter++)
        u[counter][0] = coord_initial_function(counter * tau);


    for (int time = 1; time < itter_t; time++)
    {
        for (int coord = 1; coord < itter_x; coord++)
        {
            u[time][coord] = 1.0 / (1.0 / 2 / h + 1.0 / 2 / tau) * (f_function(time * tau - 1.0 / 2 * tau, coord * h - 1.0 / 2 * h) - 1.0 / 2 / tau * (u[time][coord - 1]-u[time - 1][coord - 1] - u[time - 1][coord]) - 1.0 / 2 / h * (u[time - 1][coord - 1]-u[time][coord - 1] - u[time - 1][coord - 1]));
        }
    }

    show_table(u);


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
    return x ;
}


void show_table(double** u)
{
    int time;
    std::cout << std::endl << "\t\tTable u[time][coord]" << std::endl << std::endl << "time" << std::endl << std::endl;
    for (time = 0; time < itter_t; time++)
    {
        for (int coord = 0; coord < itter_x; coord++)
            printf("%5.2lf ", u[itter_t - 1 -time][coord]);
        if (time == itter_t - 1)
            std::cout << "\t coord";
        std::cout << std::endl;
    }
}
