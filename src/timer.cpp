#include "timer.h"
#include <fstream>
#include<iostream>

void parallel_timer::start(PetscInt node_rank)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if(node_rank==0)
    {
        time = MPI_Wtime();
    }
}
void parallel_timer::stop(PetscInt node_rank)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if(node_rank==0)
    {
        time = MPI_Wtime() - time ;
    }
}
void parallel_timer::display_result(PetscInt node_rank)
{
    if(node_rank==0)
    {
        std::cout<<timing_description<<' '<<time<<'\n';
    }
}


PetscScalar parallel_timer::get_time(PetscInt node_rank)
{
	if(node_rank == 0)
	{
		return this->time;
	}
	else
	{
		return 0.0; //dummy int..
	}

}

parallel_timer::parallel_timer(std::string timing_description_in,PetscInt node_rank)
{
    if(node_rank ==0)
    {
        timing_description = std::move(timing_description_in);
    }
}
