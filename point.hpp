#ifndef POINT // #include guards
#define POINT

#include<bits/stdc++.h>
#include "constants.hpp"
using namespace std;
typedef vector<double> vd;

struct Point
{
    vd cords;

    Point() : cords(K) { cords[K - 1] = 1; }
    Point(double x, double y, double z) : cords(K)
    {
        cords[0] = x;
        cords[1] = y;
        cords[2] = z;
        cords[3] = 1;
    }

    Point &normDist()
    {
        if (fabs(cords[0]) < EPS && fabs(cords[1]) < EPS && fabs(cords[2]) < EPS)
            return *this;

        double d = 0;
        for (int i = 0; i < K - 1; i++)
            d += cords[i] * cords[i];
        d = sqrt(d);

        for (int i = 0; i < K - 1; i++)
            cords[i] /= d;

        return *this;
    }

    Point &normWeight()
    {
        for (int i = 0; i < K; i++)
            cords[i] /= cords[K - 1];
        return *this;
    }

    double &operator[](int idx) { return cords[idx]; }

    double dot(Point &rhs)
    {
        double ans = 0;
        for (int i = 0; i < K - 1; i++)
            ans += cords[i] * rhs[i];

        return ans;
    }

    Point cross(Point &rhs)
    {
        Point p;

        p[0] = cords[1] * rhs[2] - cords[2] * rhs[1];
        p[1] = cords[2] * rhs[0] - cords[0] * rhs[2];
        p[2] = cords[0] * rhs[1] - cords[1] * rhs[0];

        return p;
    }

    Point operator-() { return Point(-cords[0], -cords[1], -cords[2]); }
};

Point operator*(Point p, double m) { return Point(p[0] * m, p[1] * m, p[2] * m); }
Point operator*(double m, Point p) { return Point(p[0] * m, p[1] * m, p[2] * m); }
Point operator+(Point p1, Point p2)
{
    return Point(p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]);
}
Point operator-(Point p1, Point p2)
{
    return Point(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
}

ostream &operator<<(ostream &dout, Point p)
{
    dout << setprecision(7) << fixed << p[0] << " " << p[1] << " " << p[2];
    return dout;
}

istream &operator>>(istream &din, Point &p)
{
    din >> p[0] >> p[1] >> p[2];
    return din;
}



#endif