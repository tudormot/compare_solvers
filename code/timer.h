#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED

#include<string>
#include"petscksp.h"


//a wrapper of the MPI timing fuction, object to be used from a global context
//each process will create one parallel_timer object when called from a global context, however the timing_description string will only be stored inthe object of process with rank=0
//slight waste of memory in the sense that an empty string and an integer are created on the other processes for nothing, however this shouldn't be a problem
//MPI barriers used internally to ensure a good timing approximation
class parallel_timer{
private:

    PetscScalar time; //after start holds time of start, after stop holds tstart-tstop
    std::string timing_description; //this is populated via the constructor, will be used as a description of the timing


public:

    void start(PetscInt node_rank);
    void stop(PetscInt node_rank);
    void display_result(PetscInt node_rank);
    PetscScalar get_time(PetscInt node_rank);

    parallel_timer(std::string timing_description_in, PetscInt node_rank);
};

#endif // TIMER_H_INCLUDED
