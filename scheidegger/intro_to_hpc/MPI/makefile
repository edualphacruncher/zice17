
### Fortran compiler wrapper
FC=mpif90

### CPP compiler wrapper
CXX=mpicxx

### Optimization flags
FFLAGS=-O3

#######################################################################	

all: 1.hello_world_mpi.exec 2.ping_pong.exec 3.MPI_reduce.exec  \
     1a.hello_world_mpi_cpp.exec \
     2a.integrate.exec  \


#######################################################################	

## Example 1
1.hello_world_mpi.exec: 1.hello_world_mpi.f90
	$(FC) $(FFLAGS) 1.hello_world_mpi.f90 -o 1.hello_world_mpi.exec

## Example 1a
1a.hello_world_mpi_cpp.exec: 1a.hello_world_mpi.cpp
	$(CXX) $(FFLAGS) 1a.hello_world_mpi.cpp -o 1a.hello_world_mpi_cpp.exec
	
## Example 2
2.ping_pong.exec: 2.ping_pong.f90
	$(FC) $(FFLAGS) 2.ping_pong.f90 -o 2.ping_pong.exec

## Example 2a
2a.integrate.exec: 2a.integrate.f90
	$(FC) $(FFLAGS) 2a.integrate.f90 -o 2a.integrate.exec
	
## Example 3
3.MPI_reduce.exec: 3.MPI_reduce.f90
	$(FC) $(FFLAGS) 3.MPI_reduce.f90 -o 3.MPI_reduce.exec


#######################################################################	

clean :
	rm -f *.exec

#######################################################################	
