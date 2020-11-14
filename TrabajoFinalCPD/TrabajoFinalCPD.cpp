// TrabajoFinalCPD.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.
//

#include <time.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <fstream>

#include "cluster.h"

extern "C" {
    #define STB_IMAGE_IMPLEMENTATION
    #include "stb_image.h"
}

bool load_image(std::vector<unsigned char>& image, const std::string& filename, int& x, int& y)
{
    int n;
    unsigned char* data = stbi_load(filename.c_str(), &x, &y, &n, 4);
    if (data != nullptr)
    {
        image = std::vector<unsigned char>(data, data + x * y * 4);
    }
    stbi_image_free(data);
    return (data != nullptr);
}

void getPoints(std::vector<Point>& points, std::vector<unsigned char>& image, int width, int height)
{
    const size_t RGBA = 4;
    size_t index;
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            index = RGBA * (j * width + i);
            points.push_back(Point(i, j, static_cast<int>(image[index + 0]), static_cast<int>(image[index + 1]), static_cast<int>(image[index + 2])));
        }
    }
}

void createImage(std::vector<Point>& points, int width, int height)
{
    std::ofstream img("picture.ppm");
    img << "P3" << std::endl;
    img << width << " " << height << std::endl;
    img << "255" << std::endl;

    for (int y = 0; y < width; y++) {
        for (int x = 0; x < height; x++) {

            int r = points[x * width + y].r;
            int g = points[x * width + y].g;
            int b = points[x * width + y].b;

            img << r << " " << g << " " << b << std::endl;
        }
    }
}

void createImage2(std::vector<Point>& points, int width, int height, int k)
{
    std::vector<Point> colors;
    srand(time(NULL));
    for (int i = 0; i < k; i++)
    {
        colors.push_back(Point(0, 0, rand() % 255, rand() % 255, rand() % 255));
    }


    std::ofstream img("picture.ppm");
    img << "P3" << std::endl;
    img << width << " " << height << std::endl;
    img << "255" << std::endl;

    for (int y = 0; y < width; y++) {
        for (int x = 0; x < height; x++) {

            //std::cout << points[x * width + y].centroide << std::endl;
            int r = colors[points[x * width + y].centroide].r;
            int g = colors[points[x * width + y].centroide].g;
            int b = colors[points[x * width + y].centroide].b;

            img << r << " " << g << " " << b << std::endl;
        }
    }
}

void getMeans(int k, std::vector<Cluster>& clusters, std::vector<Point>& points)
{
    srand(time(NULL));
    for (int i = 0; i < k; i++)
    {
        int temp = rand() % points.size();
        clusters.push_back(Cluster(points[temp]));
    }
}

bool comparePoint(std::vector<Cluster>& clusters, Point& point)
{
    double distance = INFINITY;
    double tempDistance = 0;
    bool res = false;

    for (int i = 0; i < clusters.size(); i++)
    {
         /*
        tempDistance = sqrt(pow(clusters[i].mean.x - point.x, 2) +
            pow(clusters[i].mean.y - point.y, 2) +
            pow(clusters[i].mean.r - point.r, 2) +
            pow(clusters[i].mean.g - point.g, 2) +
            pow(clusters[i].mean.b - point.b, 2));
        */
        tempDistance = sqrt(pow(clusters[i].mean.x - point.x, 2) +
            pow(clusters[i].mean.y - point.y, 2));


        if (tempDistance < distance) {
            
            if (point.centroide != -1)
            {
                clusters[point.centroide].FindAndRemove(&point);
            }
            
            point.centroide = i;
            clusters[i].elements.push_back(&point);
            res = true;
            distance = tempDistance;
        }
    }
    //std::cout << point.centroide << std::endl;
    return res;
}

void updateMeans(std::vector<Cluster>& clusters)
{
    for (int i = 0; i < clusters.size(); i++)
    {
        clusters[i].updateMean();
    }
}

void kmeansSecuencial(int k, double tolerance, std::vector<Point>& points)
{
    std::vector<Cluster> clusters;

    getMeans(k, clusters, points);
    bool temp;
    double converge;
    int count1;
    int count0 = points.size();
    //contador temporal para las iteraciones 
    int countt = 0;

    do {
        count1 = 0;
        for (int i = 0; i < points.size(); i++)
        {
            temp = comparePoint(clusters, points[i]);
            //std::cout << points[i].centroide << std::endl;
            if (temp)
            {
                count1++;
            }
        }
        updateMeans(clusters);
        clusters[0].mean.print();

        double adsCount = (double)(count1 + count0) / 2.0;
        double sizePoints = (double)(points.size() * 1.0);
        converge = (double)(count1 + count0) / 2.0 / (double)(points.size() * 1.0);
        count0 = count1;
        std::cout << adsCount << " / " << sizePoints << " = " << converge << std::endl;

        //aumento del contador temporal
        countt++;
    } while (countt < 30);//(converge > tolerance);
    std::cout << "end kmeans" << std::endl;
}

int main()
{
    std::string filename = "image.jpg";

    int width, height;
    std::vector<unsigned char> image;
    bool success = load_image(image, filename, width, height);
    if (!success)
    {
        std::cout << "Error loading image\n";
        return 1;
    }

    std::cout << "Image width = " << width << '\n';
    std::cout << "Image height = " << height << '\n';

    const size_t RGBA = 4;

    int x = 3;
    int y = 4;
    size_t index = RGBA * (y * width + x);

    std::cout << "RGBA pixel @ (x=3, y=4): "
        << static_cast<int>(image[index + 0]) << " "
        << static_cast<int>(image[index + 1]) << " "
        << static_cast<int>(image[index + 2]) << " "
        << static_cast<int>(image[index + 3]) << '\n';

    std::vector<Point> points;
    getPoints(points, image, width, height);

    kmeansSecuencial(20, 0.5, points);
    std::cout << "-----" << std::endl;
    createImage2(points, width, height, 20);

    return 0;
}

// Ejecutar programa: Ctrl + F5 o menú Depurar > Iniciar sin depurar
// Depurar programa: F5 o menú Depurar > Iniciar depuración

// Sugerencias para primeros pasos: 1. Use la ventana del Explorador de soluciones para agregar y administrar archivos
//   2. Use la ventana de Team Explorer para conectar con el control de código fuente
//   3. Use la ventana de salida para ver la salida de compilación y otros mensajes
//   4. Use la ventana Lista de errores para ver los errores
//   5. Vaya a Proyecto > Agregar nuevo elemento para crear nuevos archivos de código, o a Proyecto > Agregar elemento existente para agregar archivos de código existentes al proyecto
//   6. En el futuro, para volver a abrir este proyecto, vaya a Archivo > Abrir > Proyecto y seleccione el archivo .sln
