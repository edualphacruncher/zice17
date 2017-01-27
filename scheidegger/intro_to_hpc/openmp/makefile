#######################################################################	

all: 1.hello_world.exec 1a.hello_world.exec 2.racing_cond.exec 3.racing_cond_fix.exec \
	4.integration_pi.exec 4a.integration_pi.exec \
	4b.integration_pi_reduction.exec
	
#######################################################################	

## Example 1
1.hello_world.exec : 1.hello_world.f90
	gfortran 1.hello_world.f90 -fopenmp -o 1.hello_world.exec

## Example 1a
1a.hello_world.exec : 1a.hello_world.cpp
	g++ 1a.hello_world.cpp -fopenmp -o 1a.hello_world.exec	
	
## Example 2
2.racing_cond.exec : 2.racing_cond.f90
	gfortran 2.racing_cond.f90 -fopenmp -o 2.racing_cond.exec	

## Example 3
3.racing_cond_fix.exec : 3.racing_cond_fix.f90
	gfortran 3.racing_cond_fix.f90 -fopenmp -o 3.racing_cond_fix.exec	
	
## Example 4
4.integration_pi.exec: 4.integration_pi.f90
	gfortran 4.integration_pi.f90 -fopenmp -o 4.integration_pi.exec

## Example 4a
4a.integration_pi.exec: 4a.integration_pi.cpp
	g++ 4a.integration_pi.cpp -fopenmp -o 4a.integration_pi.exec	
	
## Example 4b
4b.integration_pi_reduction.exec: 4b.integration_pi_reduction.f90
	gfortran 4b.integration_pi_reduction.f90 -fopenmp -o 4b.integration_pi_reduction.exec	

#######################################################################	
	
clean :
	rm -f *.exec
	
#######################################################################	
	