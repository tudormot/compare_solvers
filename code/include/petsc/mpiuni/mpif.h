!
!     Trying to provide as little support for fortran code in petsc as needed
!
#include "mpiunifdef.h"
!
!     External objects outside of MPI calls
       integer MPI_COMM_WORLD
       parameter (MPI_COMM_WORLD = 2)
       integer MPI_COMM_SELF
       parameter (MPI_COMM_SELF = 1)
       integer  MPI_COMM_NULL
       parameter (MPI_COMM_NULL = 0)
       integer MPI_IDENT
       parameter (MPI_IDENT = 0)
       integer MPI_UNEQUAL
       parameter (MPI_UNEQUAL = 3)
       integer MPI_KEYVAL_INVALID
       parameter (MPI_KEYVAL_INVALID = 0)
       integer MPI_SUCCESS
       parameter (MPI_SUCCESS = 0)
       integer MPI_ERR_OTHER
       parameter (MPI_ERR_OTHER = 17)
       integer MPI_ERR_UNKNOWN
       parameter (MPI_ERR_UNKNOWN = 18)
       integer MPI_ERR_INTERN
       parameter (MPI_ERR_INTERN = 21)

       integer MPI_PACKED
       parameter (MPI_PACKED=0)
       integer MPI_ANY_SOURCE
       parameter (MPI_ANY_SOURCE=2)
       integer MPI_ANY_TAG
       parameter (MPI_ANY_TAG=-1)
       integer MPI_UNDEFINED
       parameter (MPI_UNDEFINED=-32766)
       INTEGER MPI_INFO_NULL
       PARAMETER (MPI_INFO_NULL=0)


       integer MPI_REQUEST_NULL
       parameter (MPI_REQUEST_NULL=0)

       integer MPI_STATUS_SIZE
       parameter (MPI_STATUS_SIZE=3)
       INTEGER MPI_SOURCE,MPI_TAG,MPI_ERROR
       PARAMETER(MPI_SOURCE=1,MPI_TAG=2,MPI_ERROR=3)


!     Data Types. Same Values used in mpi.c
       integer MPI_INTEGER,MPI_LOGICAL
       integer MPI_REAL,MPI_DOUBLE_PRECISION
       integer MPI_COMPLEX, MPI_CHARACTER
       integer MPI_2INTEGER
       integer MPI_DOUBLE_COMPLEX
       integer MPI_INTEGER4
       integer MPI_INTEGER8
       integer MPI_2DOUBLE_PRECISION
       integer MPI_REAL4,MPI_REAL8

!
!  These should match the values in mpi.h many below are wrong
!
       parameter (MPI_INTEGER=x'400104')
       parameter (MPI_LOGICAL=x'400104')
       parameter (MPI_REAL=x'100104')
       parameter (MPI_REAL4=x'100104')
       parameter (MPI_DOUBLE_PRECISION=x'100108')
       parameter (MPI_REAL8=x'100108')
       parameter (MPI_COMPLEX=x'200108')
       parameter (MPI_CHARACTER=x'300101')
       parameter (MPI_2INTEGER=x'e00108')
       parameter (MPI_DOUBLE_COMPLEX=x'200110')
       parameter (MPI_INTEGER4=x'400104')
       parameter (MPI_INTEGER8=x'400108')
       parameter (MPI_2DOUBLE_PRECISION=x'100208')

       integer MPI_SUM
       parameter (MPI_SUM=1)
       integer MPI_MAX
       parameter (MPI_MAX=2)
       integer MPI_MIN
       parameter (MPI_MIN=3)
       integer MPI_MAXLOC
       parameter (MPI_MAXLOC=12)
       integer MPI_MINLOC
       parameter (MPI_MINLOC=13)

       integer MPI_MAX_PROCESSOR_NAME
       parameter (MPI_MAX_PROCESSOR_NAME=128-1)

!
!  some parameters require common blocks?
!
       integer MPI_IN_PLACE
       common /MPIUNIPRIV/ MPI_IN_PLACE
       save /MPIUNIPRIV/

