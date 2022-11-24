#include <stdio.h>
#include <mpi.h>
#include <math.h>

double Trap(double a,double h,int k); 
double f(double x);

int main(void){
	int my_rank, comm_sz,local_n;
	double a= 0.0,b=1.0,local_a,local_b;
	double local_int,total_int;
	int source;
	double h;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);

	h = (b-a)/comm_sz;

	local_a = a + my_rank*comm_sz*h;
	local_b = local_a + comm_sz*h;
	total_int = (f(a) + f(b))/2.0;
	local_int = Trap(local_a,h,my_rank);

	if (my_rank != 0){
		MPI_Send(&local_int,1,MPI_DOUBLE,0,0,
			MPI_COMM_WORLD);
	} else {
		total_int = local_int;
		for (source = 1; source < comm_sz; source++){
			MPI_Recv(&local_int,1,MPI_DOUBLE,source,0,
				MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			total_int += local_int;
		}
	}
	
	if (my_rank == 0){
		printf("With n= %d trapezoids, our estimate\n", comm_sz);
		printf("of the integral from %f to %f = %.15e\n",
			a,b,total_int);
	}
	MPI_Finalize();
	return 0;
}


double Trap(double a,double h,int k) {
 double estimate, x;
 
 estimate = 2*(f(a+k*h));
 estimate *= h;
 
 return estimate;
}

double f(double x){
 double res;
 
 res = (x*x*x) + 2*(x*x) - x +1;
 return res;
 }
