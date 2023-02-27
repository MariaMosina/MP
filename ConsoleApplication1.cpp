#include "pch.h"
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <time.h>       

// Task 1
int main() {
#pragma omp parallel
	{
		printf("Num threads: %d, Thread %d is working now: Hello, World!\n", omp_get_num_threads(), omp_get_thread_num());
	}
	return 0;
}


// Task 2
#define n 16000
int main()
{
	double time_spent = 0.0;
	clock_t begin = clock();
	
	int i;
	double a[n]; 
	double b[n - 2];
	for (i = 0; i < n; i++) {
		a[i] = i;
	}

#pragma omp parallel shared(a, b) 
	{
#pragma omp for schedule(static)
//#pragma omp for schedule(static, 1000)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 4)
//#pragma omp for schedule(guided, 4)
		for (i = 1; i < n-1; i++)
		{
			b[i-1] = (a[i - 1] + a[i] + a[i + 1]) / 3.0;
		}
	}

	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("The time: %f  seconds\n", time_spent);
	return 0;
}

// Task 4
#define m 100
#define n 100
#define k 1000
int main()
{
	double time_spent = 0.0;
	clock_t begin = clock();
	int i, j;
	double** A, ** B, ** C;
	A = new double*[m];
	B = new double*[k];
	for (int i = 0; i < m; i++) {
		A[i] = new double[k];
	}
	for (int i = 0; i < k; i++) {
		B[i] = new double[n];
	}
	for (i = 0; i < m; i++) {
		for (j = 0; j < k; j++) {
			A[i][j] = i + j;
		}
	}
	for (i = 0; i < k; i++) {
		for (j = 0; j < n; j++) {
			B[i][j] = i * (j+1);
		}
	}
	C = new double*[m];
	int number_of_threads = 4;
#pragma omp parallel for num_threads(number_of_threads)
	for (int i = 0; i < m; i++)
	{
		C[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			C[i][j] = 0;
			for (int g = 0; g < k; g++)
				C[i][j] += A[i][g] * B[g][j];
		}
	}

	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Number of threads: %d \nThe time: %f  seconds\n", number_of_threads, time_spent);
	return 0;
}