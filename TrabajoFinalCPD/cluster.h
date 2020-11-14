#pragma once

#include <vector>

#include "point.h"

class Cluster
{
public:
	Point mean;
	std::vector<Point*> elements;

	Cluster(Point point);
	void updateMean();
	bool FindAndRemove(Point* point);
};

