#include <stdio.h>
#include <algorithm>
#include "mpi.h"
#include <time.h> 
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


//Task 1

int one()
{
	int rank, n, i, message; 
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n);
	printf("Hello MPI from process = %d, total number of process = %d\n", rank, n);

	MPI_Finalize();

	return 0;
}

// Task 2

int two()
{
	int rank;
	float a, b;
	double t01, t00, t10, t11, dt0 = 0, dt1 = 0;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	a = 0.0;
	b = 0.0;

	if (rank == 0)
	{
		b = 1.0;
		t00 = MPI_Wtime();
		MPI_Send(&b, 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&a, 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &status);
		t01 = MPI_Wtime();
		dt0 = t01 - t00;
		printf("\n Time : %f", dt0);
	}
	if (rank == 1)
	{
		a = 2.0;
		t10 = MPI_Wtime();
		MPI_Recv(&b, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Send(&a, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
		t11 = MPI_Wtime();
		dt1 = t11 - t10;
		printf("\n Time: %f", rank, dt1);
	}
	MPI_Finalize();
	return (0);

}

// Task 3

#define processes 4
void three(int schema) {
	int size, rank, recvRank;
	double start_time, end_time;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	switch (schema) {
		// кольцо - "эстафетная палочка"
	case 1:
	{
		if (rank == 0) {
			MPI_Send(&rank, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("Process %d sent message '%d'.\n", rank, rank);
			MPI_Recv(&recvRank, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Process %d received message '%d'.\n", rank, recvRank);
		}
		else {
			MPI_Recv(&recvRank, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Process %d received message '%d'.\n", rank, recvRank);
			rank = recvRank + 1;
			if (rank == size - 1) {
				MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				printf("Process %d sent message '%d'.\n", rank, rank);
			}
			else {
				MPI_Send(&rank, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
				printf("Process %d sent message '%d'.\n", rank, rank);
			}
		}
		break;
	}

	// колько - "сдвиг"
	case 2:
	{
		MPI_Request request;

		int prev = rank - 1;
		int next = rank + 1;
		if (rank == 0)
			prev = size - 1;
		if (rank == size - 1)
			next = 0;

		MPI_Isend(&rank, 1, MPI_INT, next, 0, MPI_COMM_WORLD, &request);
		MPI_Irecv(&recvRank, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
		for (int i = 0; i < size; i++) {
			if (rank == i) {
				printf("Process %d sent message '%d' to %d.\n", rank, rank, next);
				printf("Process %d received message '%d' from %d.\n", rank, recvRank, prev);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		break;
	}

	// master-slave
	case 3:
	{
		if (rank == 0)
		{
			for (int i = 1; i < size; i++)
			{
				MPI_Recv(&recvRank, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				printf("Process %d received a message '%d' from %d.\n", rank, recvRank, recvRank);
			}
		}
		else
			MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		break;
	}


	// каждый -> каждому
	case 4:
	{
		int send[processes];
		int recv[processes];
		int i;
		for (i = 0; i < processes; i++) {
			send[i] = rank;
		}
		for (i = 0; i < processes; i++) {
			if (i != rank) {
				MPI_Send(&send[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Recv(&recv[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				printf("Process %d received message %d from process %d.\n", rank, recv[i], i);
			}
		}
		break;
	}
	}
	MPI_Finalize();
}

// Lab
	
double f(double x, double y) {
	double arg;
	arg = -(x*x+y*y);
	return exp(arg);
}

void integral(double a1, double b1, double a2, double b2) {
	int n, rank;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n);
	double a, b;
	a = a1 + rank * ((b1 - a1) / n);
	b = a1 + ((rank + 1) * ((b1 - a1) / n));
	double step = pow(10, -4);
	double count_x = (b - a) / step;
	double count_y = (b2 - a2) / step;
	double sum = 0;
	for (int i=0; i < count_x; i++) {
		for (int j = 0; j < count_y; j++) {
			sum += f(a + i*step + step / 2.0, a2+ j*step + step / 2.0);
		}
	}
	sum = sum * step * step;
	if (rank!=0) {
		MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	else {
		double tmp = 0;
		for (int i = 1; i < n; i++) {
			MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			sum += tmp;
		}
		printf("\nIntegral = %.12f", sum);
	}
	MPI_Finalize();
}



int main(int *argc, char **argv) {
	MPI_Init(argc, &argv);
	//one();
	//two();
	//three(4);
	//integral(0, 3 / 2.0, 0, 5 / 3.0);
	//integral(0, 3 / 2.0, 0, 1);
	//integral(-3/2.0, 5/2.0 , -1, 1);
	integral(-4, 4, -4, 4);
}