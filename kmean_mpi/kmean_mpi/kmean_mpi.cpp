// kmean_mpi.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.
//

#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <winsock.h>
#include <fstream>

extern "C" {
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
}

typedef struct
{
	double _r;
	double _g;
	double _b;
	double _m;
	double _n;

} Point;

//Reads the image dimensions a x b pixels stream
void readImageSize(FILE* ifp, int* K, int* a, int* b)
{
	fscanf_s(ifp, "%d\n", K);
	printf("%d\n", *K);

	fscanf_s(ifp, "%d\n", a);
	printf("%d\n", *a);

	fscanf_s(ifp, "%d\n", b);
	printf("%d\n", *b);
}

//reads the ifp file and stores in structure
void readPoints(FILE* ifp, Point* points, int num_points)
{
	int i;
	for (i = 0; i < num_points; i++)
	{
		fscanf_s(ifp, "%lf,%lf,%lf,%lf,%lf", &points[i]._r, &points[i]._g, &points[i]._b, &points[i]._m, &points[i]._n);
		//printf("%lf,%lf,%lf,%lf,%lf", points[i]._r, points[i]._g, points[i]._b, points[i]._m, points[i]._n);
	}
}

//Initialize random points as assumed means
void initialize(Point* mean, int K, int num_points, Point* points)
{
	int i, a, p = 2;
	srand(time(NULL));
	for (i = 0; i < K; i++)
	{
		a = num_points / p;
		//printf("\n num_points: %d\n", num_points/p);
		mean[i]._r = points[a]._r;
		mean[i]._g = points[a]._g;
		mean[i]._b = points[a]._b;
		mean[i]._m = points[a]._m;
		mean[i]._n = points[a]._n;
		/*mean[i]._r=((double)(rand()%1000))/1000;
		mean[i]._g=((double)(2*rand()%1000))/1000;
		mean[i]._b=((double)(3*rand()%1000))/1000;
		mean[i]._m=((double)(4*rand()%1000))/1000;
		mean[i]._n=((double)(5*rand()%1000))/1000;*/
		//printf("%lf,%lf,%lf,%lf,%lf\n",mean[i]._r,mean[i]._g,mean[i]._b,mean[i]._m,mean[i]._n);
		p++;
		/*mean[i]._r=((double)(rand()%1000))/1000;
		mean[i]._g=((double)(2*rand()%1000))/1000;
		mean[i]._b=((double)(3*rand()%1000))/1000;
		mean[i]._m=((double)(4*rand()%1000))/1000;
		mean[i]._n=((double)(5*rand()%1000))/1000;*/
	}
}
//All points having no clusters
void IntClusterMem(int* cluster, int num_points)
{
	int i;
	for (i = 0; i < num_points; i++)
	{
		cluster[i] = -1;
	}
}

//Distance
double calculateDistance(Point point1, Point point2)
{
	return sqrt((pow((point1._r - point2._r), 2) + pow((point1._g - point2._g), 2) + pow((point1._b - point2._b), 2) + pow((point1._m - point2._m), 2) + pow((point1._n - point2._n), 2)));
}

//to calculate which cluster is the point belonging to.
int pointsCluster(Point point, Point* mean, int K)
{
	int parent = 0;
	double dist = 0;
	double minDist = calculateDistance(point, mean[0]);
	int i;
	for (i = 1; i < K; i++)
	{
		dist = calculateDistance(point, mean[i]);
		if (minDist >= dist)
		{
			parent = i;
			minDist = dist;
		}
		//printf("\n%lf,%lf,%lf,%lf,%lf\n",mean[i]._r,mean[i]._g,mean[i]._b,mean[i]._m,mean[i]._n);
	}
	//printf("\nMin dist = %lf",minDist);

	//printf("\n%lf,%lf,%lf,%lf,%lf\n",point._r,point._g,point._b,point._m,point._n);
	return parent;
}

//calculate new mean
void calcNewMean(Point* points, int* cluster, Point* mean, int K, int num_points)
{
	Point* newMean = new Point[K];
	int* members = new int[K];
	int i;
	for (i = 0; i < K; i++)
	{
		members[i] = 0;
		newMean[i]._r = 0;
		newMean[i]._g = 0;
		newMean[i]._b = 0;
		newMean[i]._m = 0;
		newMean[i]._n = 0;
	}
	for (i = 0; i < num_points; i++)
	{
		members[cluster[i]]++;
		newMean[cluster[i]]._r += points[i]._r;
		newMean[cluster[i]]._g += points[i]._g;
		newMean[cluster[i]]._b += points[i]._b;
		newMean[cluster[i]]._m += points[i]._m;
		newMean[cluster[i]]._n += points[i]._n;
	}
	for (i = 0; i < K; i++)
	{
		if (members[i] != 0.0)
		{
			newMean[i]._r /= members[i];
			newMean[i]._g /= members[i];
			newMean[i]._b /= members[i];
			newMean[i]._m /= members[i];
			newMean[i]._n /= members[i];
		}
		else
		{
			newMean[i]._r = 0;
			newMean[i]._g = 0;
			newMean[i]._b = 0;
			newMean[i]._m = 0;
			newMean[i]._n = 0;
		}
	}
	for (i = 0; i < K; i++)
	{
		mean[i]._r = newMean[i]._r;
		mean[i]._g = newMean[i]._g;
		mean[i]._b = newMean[i]._b;
		mean[i]._m = newMean[i]._m;
		mean[i]._n = newMean[i]._n;
	}
}

//check for convergence
// it checks that is each points cluster remaining the same
int chkConvrg(int* before_clusters, int* after_cluster, int num_points, float tol)
{
	int i;
	tol = num_points * tol;
	for (i = 0; i < num_points; i++)
		if ((before_clusters[i] - after_cluster[i]) > tol)
			return -1;
	return 0;
}

void createImage(Point* points, int width, int height)
{
	std::ofstream img("picture.ppm");
	img << "P3" << std::endl;
	img << width << " " << height << std::endl;
	img << "255" << std::endl;

	for (int y = 0; y < width; y++) {
		for (int x = 0; x < height; x++) {

			int r = points[y * height + x]._r;
			int g = points[y * height + x]._g;
			int b = points[y * height + x]._b;

			img << r << " " << g << " " << b << std::endl;
		}
	}
}

void createImage2(int* after_cluster, int K, int num_points, Point* points, int width, int height)
{
	Point* colors = new Point[K];
	srand(time(NULL));
	for (int i = 0; i < K; i++)
	{
		colors[i]._m = 0;
		colors[i]._n = 0;
		colors[i]._r = rand() % 255;
		colors[i]._g = rand() % 255;
		colors[i]._b = rand() % 255;
	}

	std::ofstream img("picture1.ppm");

	img << "P3" << std::endl;
	img << width << " " << height << std::endl;
	img << "255" << std::endl;

	for (int y = 0; y < width; y++) {
		for (int x = 0; x < height; x++) {

			//std::cout << points[x * width + y].centroide << std::endl;
			int r = colors[after_cluster[y * height + x]]._r;
			int g = colors[after_cluster[y * height + x]]._g;
			int b = colors[after_cluster[y * height + x]]._b;

			img << r << " " << g << " " << b << std::endl;
		}
	}
}

int main(int argc, char* argv[])
{
	int rank;
	int size;
	struct timespec start_t, stop_t;
	//clock_t start_t, stop_t;
	double span_t;
	//start_t = clock();
	char hostname[256];
	int K = 0;
	int num_points;
	int i, j, l = 1;
	int job_size;
	int job_done = 0;
	int x, y;
	float tol=0;
	//int f = 0;
	float check;
	int tam;
	clock_t tic, toc;
	double tspan = 0.0, tspantemp = 0.0;

	Point* mean;
	Point* points;
	Point* get_points;
	int* formed_clusters;
	int* before_clusters;
	int* after_cluster;
	tic = clock();
	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	tam = size;

	//printf("Hello world!  I am process number: %d on host %s started at time %f\n", rank, hostname,0);

   //stop_t = clock();
   //span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;
   //printf("(clock)Hello world!  I am process number: %d on host %s started at time %f\n", rank, hostname,span_t);

   //creating of derived MPI structure
	MPI_Datatype MPI_POINT;
	MPI_Datatype type = MPI_DOUBLE;
	int blocklen = 2;
	MPI_Aint disp = 0; //C type that holds any valid address
	MPI_Type_create_struct(1, &blocklen, &disp, &type, &MPI_POINT);
	MPI_Type_commit(&MPI_POINT);

	//Master

	if (rank != 0)
	{
		//Receiving the cluster
		//printf("Receiving data...\n");
		tic = clock();
		MPI_Barrier(MPI_COMM_WORLD);
		
		//printf("Start time for recieving Init data from master %f:  I am process number: %d on host %s\n",span_t, rank, hostname);
		MPI_Recv(&job_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&K, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

		mean = new Point[K];
		MPI_Recv(mean, K, MPI_POINT, 0, 0, MPI_COMM_WORLD, &status);
		//printf("part_size =%d\n",job_size);
		points = new Point[job_size];
		after_cluster = new int[job_size];

		for (i = 0; i < job_size; i++)
		{
			MPI_Recv(&points[i]._r, 1, MPI_POINT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&points[i]._g, 1, MPI_POINT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&points[i]._b, 1, MPI_POINT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&points[i]._m, 1, MPI_POINT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&points[i]._n, 1, MPI_POINT, 0, 0, MPI_COMM_WORLD, &status);
		}
		//printf("Received successfully [%d]\n",rank);
		

		MPI_Barrier(MPI_COMM_WORLD);
		//printf("Stop time for recieving Init data from master %f:  I am process number: %d on host %s\n",span_t, rank, hostname);
		toc = clock();
		tspantemp = (double)(toc - tic) / (double)CLOCKS_PER_SEC;
		tspan += tspantemp;

		while (1)
		{
			//printf("Calculating new clusters [%d]\n",rank);
			//printf("part_size =%d\n",job_size);
			tic = clock();
			for (i = 0; i < job_size; i++)
			{
				after_cluster[i] = pointsCluster(points[i], mean, K);
				//printf("\n Mean = %lf",mean);
				//printf("\n%lf,%lf,%lf,%lf,%lf\n",points[i]._r,points[i]._g,points[i]._b,points[i]._m,points[i]._n);
			}
			
			//stop_t = clock();
			//span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;
			MPI_Barrier(MPI_COMM_WORLD);
			
			

			//printf("Start time for sending index to master %f:  I am process number: %d on host %s\n",span_t, rank, hostname);

			//printf("sending to master [%d]\n",rank);
			MPI_Send(after_cluster, job_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);


			//stop_t = clock();
			//span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;

			//printf("Stop time for sending index to master %f:  I am process number: %d on host %s\n",tocPC5, rank, hostname);

			//if(rank == 1)
	

			MPI_Bcast(&job_done, 1, MPI_INT, 0, MPI_COMM_WORLD);

			if (job_done == 1) //No more work to be done
				break;
			MPI_Barrier(MPI_COMM_WORLD);


			//stop_t = clock();
			//span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;
			//printf("Start time for recieving mean from master %f:  I am process number: %d on host %s\n",span_t, rank, hostname);

			//Receiving recently created mean from master
			//MPI_Recv(mean,K,MPI_POINT,0,0, MPI_COMM_WORLD,&status);
			MPI_Bcast(mean, K, MPI_POINT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			toc = clock();
			tspantemp = (double)(toc - tic) / (double)CLOCKS_PER_SEC;
			tspan += tspantemp;

			//stop_t = clock();
			//span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;

			//printf("Stop time for recieving mean from master %f:  I am process number: %d on host %s\n",tocPC6, rank, hostname);

		}

	}
	else
	{

		//f = 1;
		//Readinf file
		tic = clock();
		FILE* ifp;
		fopen_s(&ifp, "input1.txt", "r");

		readImageSize(ifp, &K, &x, &y);
		num_points = x * y;
		points = (Point*)malloc(sizeof(Point) * num_points);
		readPoints(ifp, points, num_points);
		fclose(ifp);

		//Allocates memory
		before_clusters = (int*)malloc(sizeof(int) * num_points);
		after_cluster = (int*)malloc(sizeof(int) * num_points);
		mean = new Point[K];
		//cuts job depending on number of threads
		check = num_points % (size - 1);
		if (check == 0.00)
		{
			job_size = num_points / (size - 1);
		}
		else
		{
			printf("\n Enter no. of Processes as n+1 (where n divides  %d  in equal parts)\n\n", num_points);
			exit(1);
		}
		//initializing to default values
		initialize(mean, K, num_points, points);
		IntClusterMem(before_clusters, num_points);
		IntClusterMem(after_cluster, num_points);
		

		printf("Tolerance = %f\n", tol);
		//printf("Enter tolerence: ");
		//scanf("%f",&tol);
		//Ttic = MPI_Wtime();
		//MPIticC1 = MPI_Wtime();
		MPI_Barrier(MPI_COMM_WORLD);
		toc = clock();
		tspantemp = (double)(toc - tic) / (double)CLOCKS_PER_SEC;
		tspan += tspantemp;
		//printf("(M)Start time for sending Init data to slaves %f:  I am process number: %d on host %s\n",span_t, rank, hostname);

		//Sending the essential cluster to other processors
		tic = clock();
		for (i = 1; i < size; i++)
		{
			//printf("Sending data to [%d]\n",i);
			MPI_Send(&job_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&K, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(mean, K, MPI_POINT, i, 0, MPI_COMM_WORLD);
			for (j = 0; j < job_size; j++)
			{
				MPI_Send(&points[j]._r + (i - 1) * job_size, 1, MPI_POINT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&points[j]._g + (i - 1) * job_size, 1, MPI_POINT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&points[j]._b + (i - 1) * job_size, 1, MPI_POINT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&points[j]._m + (i - 1) * job_size, 1, MPI_POINT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&points[j]._n + (i - 1) * job_size, 1, MPI_POINT, i, 0, MPI_COMM_WORLD);
			}
		}
		//printf("Sent Successfully!\n");
		

		MPI_Barrier(MPI_COMM_WORLD);
		// printf("(M)Stop time for sending Init data to slaves %f:  I am process number: %d on host %s\n",span_t, rank, hostname);

		toc = clock();
		tspantemp = (double)(toc - tic) / (double)CLOCKS_PER_SEC;
		tspan += tspantemp;
		//printf("\n2 : %f\n",MPI_Wtime());
		//master processor job
		while (1)
		{
			//ticC2 = gettimeofday(&stop_t, NULL);
			tic = clock();
			MPI_Barrier(MPI_COMM_WORLD);
			
			//stop_t = clock();
			//span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;
			//printf("(M) Start time for recieving index to master %f:  I am process number: %d on host %s\n",ticC2, rank, hostname);

			//printf("\n3 : %f\n",MPI_Wtime());
			//printf("Master Receiving data...\n");
			for (i = 1; i < size; i++)
				MPI_Recv(after_cluster + (job_size * (i - 1)), job_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			//printf("Master Received data successfully\n");
			MPI_Barrier(MPI_COMM_WORLD);
			//tocC2 = gettimeofday(&stop_t, NULL);
			
			//stop_t = clock();
			//span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;
	//		printf("(clock M)Stop time for recieving index to master %f:  I am process number: %d on host %s\n",span_t, rank, hostname);

			
			//printf("(M) Stop time for recieving index to master %f:  I am process number: %d on host %s\n",tocC2, rank, hostname);
			
			//printf("\n4 : %f\n",toc);
			calcNewMean(points, after_cluster, mean, K, num_points);
			toc = clock();
			tspantemp = (double)(toc - tic) / (double)CLOCKS_PER_SEC;
			tspan += tspantemp;

			//printf("New Centroids are calculated!\n");

			if (chkConvrg(after_cluster, before_clusters, num_points, tol) == 0)
			{
				printf("K-mean algorithm Converged at iteration %d\n", l);
				job_done = 1;

			}
			else
			{
				l++;
				//printf("Iteration = %d\n",l);
				//printf("Not converged!\n");
				for (i = 0; i < num_points; i++)
					before_clusters[i] = after_cluster[i];
			}









			//Informing slaves that no more job to be done
			//clock_gettime(CLOCK_MONOTONIC,&stop_t);
					//span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));

			//ticJobdone = span_t;
			//for(i=1;i<size;i++)
			tic = clock();
			MPI_Bcast(&job_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
			//clock_gettime(CLOCK_MONOTONIC,&stop_t);
					//span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			//tocJobdone = span_t;
			//tspanJobdone += (double)(tocJobdone - ticJobdone);
			if (job_done == 1)
				break;
			MPI_Barrier(MPI_COMM_WORLD);
			//ticC3 = span_t;
		   //stop_t = clock();
		   //span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;
   //		printf("(clock M)Start time for sending mean from master %f:  I am process number: %d on host %s\n",span_t, rank, hostname);


			//printf("(M) Start time for sending mean from master %f:  I am process number: %d on host %s\n",tic, rank, hostname);

			//Sending the recently created mean
			//for(i=1;i<size;i++)
			MPI_Bcast(mean, K, MPI_POINT, 0, MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			toc = clock();
			tspantemp = (double)(toc - tic) / (double)CLOCKS_PER_SEC;
			tspan += tspantemp;

			//stop_it = clock();
			//span_t = (double)(stop_t-start_t)/(double)CLOCKS_PER_SEC;
			//sizeof
			//printf("(clock M)Stop time for sending mean from master %f:  I am process number: %d on host %s\n",span_t, rank, hostname);

			//printf("(M) Stop time for sending mean from master %f:  I am process number: %d on host %s\n",tocC3, rank, hostname);

		}

		//printf("Total time elapsed in forming clusters : %f msec\n", (tspan)*1000);
		//Outputting to the ofp file
		tic = clock();
		createImage2(after_cluster, K, num_points, points, x, y);
		createImage(points, x, y);
		toc = clock();
		tspantemp = (double)(toc - tic) / (double)CLOCKS_PER_SEC;
		tspan += tspantemp;
		printf("Total time elapsed in forming clusters : %f sec\n", tspan);
	}

	if (rank == 0)
	{
		//printf("Input time : %f sec\n", tocIn - ticIn);
		//printf("send initial data time : %f sec\n", tocC1 - ticC1);

		//printf("Recieve results from slave time : %f sec\n", tspanC2);
		//printf("NewMean calculation time : %f sec\n", tspanNewMean);
		//printf("Convergence check time : %f sec\n", tspanChCon);
		//printf("send newmeans to slaves time : %f sec\n", tspanC3);
		//printf("Total time for iterative communication data (mean, after cluster) : %f\n", tspanC2 + tspanC3);
		//printf("Output time : %f sec\n", tocOut - ticOut);
		//printf("\tMaster time (Communication) : %f sec\n", ((tocC1 - ticC1) + tspanC2 + tspanC3));
		//printf("\tMaster time (Calculation) : %f sec\n", (tspanNewMean + tspanChCon));
		//printf("\tMaster time (File operation) : %f sec\n", ((tocIn - ticIn) + (tocOut - ticOut)));
		//printf("\tTotal Master time : %f sec\n", ((tocIn - ticIn) + (tocOut - ticOut)) + (tspanNewMean + tspanChCon) + ((tocC1 - ticC1) + tspanC2 + tspanC3));
		//printf("\nTotal Time : %f sec\n", TOC - TIC);
	}
	else if (rank == 1)
	{
		//	printf("Recieve initial data time : %f sec\n", tspanPC4);
		//printf("Calculate After_Cluster time : %f sec\n", tspanPCal);
		//	printf("send results to master time : %f sec\n", tspanPC5);
		//	printf("Recieve 'job_done' status and act accordingly time : %f sec\n", tspanPC6);
		//	printf("\tSlave time (Communication) : %f sec\n", (tspanPC4 +  tspanPC5 +tspanPC6));
		//	printf("\tSlave time (Calculation) : %f sec\n", tspanPCal);
		//	printf("\tTotal Slave time : %f sec\n", (tspanPCal)+tspanPC4 +  tspanPC5+tspanPC6);
		//	printf("\tTotal time : %f sec\n", (tspanPCal)+tspanPC4 +  tspanPC5+tspanPC6 + ((tocIn - ticIn)+  (tocOut - ticOut)) + (tspanNewMean+  tspanChCon) + ((tocC1 - ticC1) +  tspanC2 + tspanC3));
	
		
	}


	//End of all

	MPI_Finalize();

	
	
	/*
	free (get_points);
	free (formed_clusters);
	free (mean);
	free (before_clusters);
	free (points);
	free (after_cluster);*/

	

	return 0;
}