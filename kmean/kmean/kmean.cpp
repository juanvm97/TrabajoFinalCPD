#include <stdio.h>
#include <math.h>
#define _CRT_SECURE_NO_DEPRECATE
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <iostream>
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

//Reads the image dimensions a x b pixels
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
		//printf("\n num_points: %d\n", num_points/p);
		mean[i]._r = points[a]._r;
		mean[i]._g = points[a]._g;
		mean[i]._b = points[a]._b;
		mean[i]._m = points[a]._m;
		mean[i]._n = points[a]._n;
		//mean[i]._r=((double)(rand()%1000))/1000;
		//mean[i]._g=((double)(2*rand()%1000))/1000;
		//mean[i]._b=((double)(3*rand()%1000))/1000;
		//mean[i]._m=((double)(4*rand()%1000))/1000;
		//mean[i]._n=((double)(5*rand()%1000))/1000;
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
	}
	return parent;
}

//calculate new mean
void calcNewMean(Point* points, int* cluster, Point* mean, int K, int num_points)
{
	Point* newMean = new Point[K];
	int* members = new int[K];
	//Point* newMean = malloc(sizeof(Point) * K);
	//int* members = malloc(sizeof(int) * K);
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
		if (abs(before_clusters[i] - after_cluster[i]) > tol)
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
	int K;
	int num_points;
	int i;
	int job_done = 0;
	int x, y;
	clock_t tic, toc;
	double tspan = 0.0, tspantemp = 0.0;
	int iter = 0;
	Point* mean;
	Point* points;
	Point* get_points;
	int* formed_clusters;
	int* before_clusters;
	int* after_cluster;
	float tol = 0.000001;

	//tol = atof(argv[3]);

	//printf("Enter Tolerance:  ");
	//scanf("%f",&tol);
	printf("Tolerance = %.10f\n", tol);
	//Readinf file
	FILE* ifp;
	fopen_s(&ifp, "input4.txt", "r");
	readImageSize(ifp, &K, &x, &y);
	num_points = x * y;
	points = new Point[num_points];
	readPoints(ifp, points, num_points);
	fclose(ifp);

	before_clusters = (int*)malloc(sizeof(int) * num_points);
	after_cluster = new int[num_points];

	mean = new Point[K];
	//mean = malloc(sizeof(Point) * K);

	//initializing to default values
	initialize(mean, K, num_points, points);
	IntClusterMem(before_clusters, num_points);
	IntClusterMem(after_cluster, num_points);

	while (1)
	{

		tic = clock();
		iter++;
		for (i = 0; i < num_points; i++)
		{
			after_cluster[i] = pointsCluster(points[i], mean, K);
		}
		calcNewMean(points, after_cluster, mean, K, num_points);
		//printf("New Centroids are calculated!\n");
		toc = clock();
		tspantemp = (double)(toc - tic) / (double)CLOCKS_PER_SEC;
		tspan += tspantemp;
		if (chkConvrg(after_cluster, before_clusters, num_points, tol) == 0)
		{
			printf("K-mean algorithm Converged!\n");
			job_done = 1;
		}
		else
		{
			//printf("Not converged!\n");
			for (i = 0; i < num_points; i++)
				before_clusters[i] = after_cluster[i];
		}

		if (job_done == 1)
			break;

	}
	printf("Total Iterations = %d\n", iter);
	printf("Total time elapsed in forming clusters : %f msec\n", tspan * 1000);
	//Outputting to the ofp file

	createImage2(after_cluster, K, num_points, points, x, y);
	createImage(points, x, y);

	//End of all
	return 0;
}