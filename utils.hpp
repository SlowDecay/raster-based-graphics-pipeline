#ifndef UTIL // #include guards
#define UTIL

#include<bits/stdc++.h>
using namespace std;
#include "point.hpp"
#include "matrix.hpp"

// a is the point to be rotated
// axis is the rotation axis. It should be normalized
// angle is the rotation angle (in degrees)
Point rodrigues(Point a, Point axis, double angle)
{
    angle *= PI / 180; // taking angle to radian

    Point res = cos(angle) * a + (1 - cos(angle)) * axis.dot(a) * axis + sin(angle) * axis.cross(a);
    return res;
}

void init()
{
    for (int i = 0; i < K; i++)
        idt[i][i] = 1;
}

Matrix getTransMatrix(Point t)
{
    Matrix res = idt;

    for (int i = 0; i < K; i++)
        res[i][K - 1] = t[i];
    res[K - 1][K - 1] = 1;

    return res;
}

Matrix getScaleMatrix(Point s)
{
    Matrix res = idt;

    for (int i = 0; i < K; i++)
        res[i][i] = s[i];
    res[3][3] = 1;

    return res;
}

Matrix getRotMatrix(double degree, Point axis)
{
    Matrix res = idt;

    axis.normDist();

    Point i(1, 0, 0);
    Point j(0, 1, 0);
    Point k(0, 0, 1);

    Point c1 = rodrigues(i, axis, degree);
    Point c2 = rodrigues(j, axis, degree);
    Point c3 = rodrigues(k, axis, degree);

    for (int i = 0; i < K - 1; i++)
        res[i][0] = c1[i], res[i][1] = c2[i], res[i][2] = c3[i];

    return res;
}
#endif