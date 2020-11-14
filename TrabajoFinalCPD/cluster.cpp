
#include <iostream>

#include "cluster.h"

Cluster::Cluster(Point point)
{
	this->mean = Point(point.x, point.y, point.r, point.g, point.b);
}

void maxAndMin(int& max, int& min, int value)
{
	if (max < value)
	{
		//std::cout << max << " " << value << std::endl;
		max = value;
	}
	if (min > value)
	{
		//std::cout << min << " " << value << std::endl;
		min = value;
	}
}

void Cluster::updateMean()
{
	int maxX = 0;
	int minX = INT_MAX;
	int maxY = 0;
	int minY = INT_MAX;
	int maxR = 0;
	int minR = INT_MAX;
	int maxG = 0;
	int minG = INT_MAX;
	int maxB = 0;
	int minB = INT_MAX;

	for (int i = 0; i < this->elements.size(); i++)
	{
		maxAndMin(maxX, minX, this->elements[i]->x);
		maxAndMin(maxY, minY, this->elements[i]->y);
		maxAndMin(maxR, minR, this->elements[i]->r);
		maxAndMin(maxG, minG, this->elements[i]->g);
		maxAndMin(maxB, minB, this->elements[i]->b);
	}
	//std::cout << "->> " <<  maxX << " " << minX << std::endl;
	this->mean.x = (maxX + minX) / 2;
	this->mean.y = (maxY + minY) / 2;
	this->mean.r = (maxR + minR) / 2;
	this->mean.g = (maxG + minG) / 2;
	this->mean.b = (maxB + minB) / 2;
}

bool Cluster::FindAndRemove(Point* point)
{
	std::vector<Point*>::iterator toErease;

	toErease = std::find(this->elements.begin(), this->elements.end(), point);

	// And then erase if found
	if (toErease != this->elements.end())
	{
		this->elements.erase(toErease);
		return true;
	}
	return false;
}
