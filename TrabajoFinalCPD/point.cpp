
#include <iostream>

#include "point.h"

Point::Point(int x, int y, int r, int g, int b)
{
	this->x = x;
	this->y = y;
	this->r = r;
	this->g = g;
	this->b = b;
	this->centroide = -1;
}

Point::Point()
{
	int x = 0;
	int y = 0;
	int r = 0;
	int g = 0;
	int b = 0;
}

void Point::print()
{
	std::cout << "---------------------" << std::endl;
	std::cout << this->x << " " << this->y << std::endl;
	std::cout << this->r << " " << this->g << " " << this->b << std::endl;
	std::cout << "---------------------" << std::endl;
}
