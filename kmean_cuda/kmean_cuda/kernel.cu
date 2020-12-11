
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda.h>

#include <iostream>
#include <fstream>

typedef struct
{
	double _r;
	double _g;
	double _b;
	double _m;
	double _n;

} Point;

#define CUDA_CALL(x) {if((x) != cudaSuccess){ \
  printf("CUDA error at %s:%d\n",__FILE__,__LINE__); \
  printf("  %s\n", cudaGetErrorString(cudaGetLastError())); \
  exit(EXIT_FAILURE);}} 

void readImageSize(FILE* ifp, int* K, int* a, int* b)
{
	fscanf(ifp, "%d\n", K);
	printf("%d\n", *K);

	fscanf(ifp, "%d\n", a);
	printf("%d\n", *a);

	fscanf(ifp, "%d\n", b);
	printf("%d\n", *b);
}

//reads the ifp file and stores in structure
void readPoints(FILE* ifp, Point* points, int num_points)
{
	int i;
	for (i = 0; i < num_points; i++)
	{
		fscanf(ifp, "%lf,%lf,%lf,%lf,%lf", &points[i]._r, &points[i]._g, &points[i]._b, &points[i]._m, &points[i]._n);
		//printf("%lf,%lf,%lf,%lf,%lf\n", points[i]._r, points[i]._g, points[i]._b, points[i]._m, points[i]._n);
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
		mean[i]._r = points[a]._r;
		mean[i]._g = points[a]._g;
		mean[i]._b = points[a]._b;
		mean[i]._m = points[a]._m;
		mean[i]._n = points[a]._n;
		p++;
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

//to calculate which cluster is the point belonging to.
__global__ void pointsCluster(int* after_cluster_d, Point* point_d, Point* Dmean, int K, int x, int y)
{
	//__shared__ Point Dmean[105];
	//printf("\n%d\t%d\t%d\n",K,x,y);
	int j, k, i;
	j = blockIdx.x * blockDim.x + threadIdx.x;
	k = blockIdx.y * blockDim.y + threadIdx.y;

	if (j < x && k < y)
	{
		int parent = 0;
		double dist = 0;
		int t = (k * (x)+j);
		//if(t>204790)
		//printf("t = %d\n",t);
		double minDist = sqrt((pow((point_d[t]._r - Dmean[0]._r), 2) + pow((point_d[t]._g - Dmean[0]._g), 2) + pow((point_d[t]._b - Dmean[0]._b), 2) + pow((point_d[t]._m - Dmean[0]._m), 2) + pow((point_d[t]._n - Dmean[0]._n), 2)));
		for (i = 1; i < K; i++)
		{
			dist = sqrt((pow((point_d[t]._r - Dmean[i]._r), 2) + pow((point_d[t]._g - Dmean[i]._g), 2) + pow((point_d[t]._b - Dmean[i]._b), 2) + pow((point_d[t]._m - Dmean[i]._m), 2) + pow((point_d[t]._n - Dmean[i]._n), 2)));
			if (minDist >= dist)
			{
				parent = i;
				minDist = dist;

			}
		}
		after_cluster_d[t] = parent;
	}
}


//calculate new mean
void calcNewMean(Point* points, int* cluster, Point* mean, int K, int num_points)
{
	Point* newMean = (Point*)malloc(sizeof(Point) * K);
	int* members = (int*)malloc(sizeof(int) * (K));
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
int chkConvrg(int* before_clusters, int* after_cluster, int num_points, float tol)
{
	int i;
	tol = num_points * tol;
	for (i = 0; i < num_points; i++)
	{
		if (abs(before_clusters[i] - after_cluster[i]) > tol)
		{
			return -1;
		}
	}
	return 0;
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
	//cpu variables
	int K;
	int num_points;
	int* before_clusters;
	int i;
	int job_done = 0;
	int x, y, iter = 0, iterations;

	Point* mean;
	Point* points;

	int* after_cluster;
	float tol;

	//gpu variables
	Point* points_d;
	Point* mean_d;
	int* after_cluster_d;
	int* before_cluster_d;

	cudaEvent_t startinit, endinit, startmean, endmean, startcal, endcal, startindex, endindex;
	cudaEvent_t start1, end1;
	float timeinit, timemean, timecal, timeindex;
	float time1;
	//float totTime = 0;
	tol = 0.0;

	cudaEventCreate(&start1);
	cudaEventCreate(&end1);
	cudaEventRecord(start1, 0);

	//Readinf file
	FILE* ifp;
	fopen_s(&ifp, "input2.txt", "r");
	readImageSize(ifp, &K, &x, &y);
	num_points = x * y;
	int blockX = 32;
	int blockY = 32;
	//allocate CPU memory
	points = (Point*)malloc(sizeof(Point) * num_points);
	readPoints(ifp, points, num_points);
	fclose(ifp);
	//printf("Input Read Successfully \n");
	before_clusters = new int[num_points];
	after_cluster = new int[num_points];
	mean = new Point[K];

	//initializing to default values
	initialize(mean, K, num_points, points);
	IntClusterMem(before_clusters, num_points);
	IntClusterMem(after_cluster, num_points);

	//printf("points = %lf",points[0]._r);
	//allocate gpu memory
	//printf("No problem till Here1\n");
	CUDA_CALL(cudaMalloc((void**)&after_cluster_d, sizeof(int) * num_points));
	CUDA_CALL(cudaMalloc((void**)&before_cluster_d, sizeof(int) * num_points));
	CUDA_CALL(cudaMalloc((void**)&points_d, sizeof(Point) * num_points));
	CUDA_CALL(cudaMalloc((void**)&mean_d, sizeof(Point) * K));
	//printf("No problem till Here2\n");

	cudaEventCreate(&startinit);
	cudaEventCreate(&endinit);
	cudaEventRecord(startinit, 0);

	//copy data points to device
	CUDA_CALL(cudaMemcpy(points_d, points, sizeof(Point) * num_points, cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(after_cluster_d, after_cluster, sizeof(int) * num_points, cudaMemcpyHostToDevice));

	cudaEventRecord(endinit, 0);
	cudaEventSynchronize(endinit);
	cudaEventElapsedTime(&timeinit, startinit, endinit);

	while (1)
	{
		//printf("No problem till Here4\n");
		iter++;
		cudaEventCreate(&startmean);
		cudaEventCreate(&endmean);
		cudaEventRecord(startmean, 0);
		//copy initial centroids to device
		CUDA_CALL(cudaMemcpy(mean_d, mean, sizeof(Point) * K, cudaMemcpyHostToDevice));
		cudaEventRecord(endmean, 0);
		cudaEventSynchronize(endmean);
		cudaEventElapsedTime(&timemean, startmean, endmean);
		//cuda memory copy
		//CUDA_CALL(cudaMemcpy(after_cluster_d, after_cluster, sizeof(int)*num_points, cudaMemcpyHostToDevice));
		//CUDA_CALL(cudaMemcpy(before_cluster_d, before_clusters, sizeof(int)*num_points, cudaMemcpyHostToDevice));
		//CUDA_CALL(cudaMemcpy(x_d, &x, sizeof(int), cudaMemcpyHostToDevice));
		//CUDA_CALL(cudaMemcpy(y_d, &y, sizeof(int), cudaMemcpyHostToDevice));
		//CUDA_CALL(cudaMemcpy(K_d, &K, sizeof(int), cudaMemcpyHostToDevice));
		cudaEventCreate(&startcal);
		cudaEventCreate(&endcal);
		cudaEventRecord(startcal, 0);

		dim3 block(blockX, blockY);
		dim3 grid((x + blockX - 1) / blockX, (y + blockY - 1) / blockY);

		pointsCluster <<< grid, block >>> (after_cluster_d, points_d, mean_d, K, x, y);

		//printf("Time taken by parallel portion: %f\n",time);
		//totTime +=time;
		//printf("No problem till Here5\n");
		cudaDeviceSynchronize();
		cudaEventRecord(endcal, 0);
		cudaEventSynchronize(endcal);
		cudaEventElapsedTime(&timecal, startcal, endcal);

		cudaEventCreate(&startindex);
		cudaEventCreate(&endindex);
		cudaEventRecord(startindex, 0);

		CUDA_CALL(cudaMemcpy(after_cluster, after_cluster_d, sizeof(int) * num_points, cudaMemcpyDeviceToHost));
		
		cudaEventRecord(endindex, 0);
		cudaEventSynchronize(endindex);
		cudaEventElapsedTime(&timeindex, startindex, endindex);
		calcNewMean(points, after_cluster, mean, K, num_points);
		//printf("New Centroids are calculated!\n");

		if (chkConvrg(after_cluster, before_clusters, num_points, tol) == 0)
		{
			printf("K-mean algorithm Converged!\n");
			job_done = 1;
		}
		else
		{
			//printf("Not converged!\n");
			for (i = 0; i < num_points; i++)
			{
				//printf("1 after_cluster[%d]=%d, before_clusters[%d]=%d\n",i,after_cluster[i],i,before_clusters[i]);

				before_clusters[i] = after_cluster[i];

				//printf("after_cluster[%d]=%d, before_clusters[%d]=%d\n",i,after_cluster[i],i,before_clusters[i]);
			}
		}

		if (job_done == 1)
			break;

	}
	
	//Outputting to the ofp file
	//createImage2(after_cluster, K, num_points, points, x, y);

	cudaEventRecord(end1, 0);
	cudaEventSynchronize(end1);
	cudaEventElapsedTime(&time1, start1, end1);
	printf("Time for sending initial data from host to device : %f\t sec\n", timeinit / 1000);
	printf("Time for sending new means from host to device : %f\t sec\n", timemean / 1000);
	printf("Time for calculation : %f\t sec\n", timecal / 1000);
	printf("Time for sending new index from device to host : %f\t sec\n", timeindex / 1000);
	printf("Total Time : %f\t sec\n", time1 / 1000);
	printf("Total Iterations = %d\n", iter);
	CUDA_CALL(cudaFree(after_cluster_d));
	CUDA_CALL(cudaFree(mean_d));
	CUDA_CALL(cudaFree(points_d));
	free(before_clusters);
	free(mean);
	free(points);
	free(after_cluster);

	return 0;
}

