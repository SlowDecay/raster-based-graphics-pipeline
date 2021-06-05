#ifndef TRIANGLE // #include guards
#define TRIANGLE

#include <bits/stdc++.h>
#include "point.hpp"

using namespace std;
typedef vector<double> vd;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

struct Segment
{
    Point p1, p2;

    Segment(Point p1, Point p2) : p1(p1), p2(p2) {}

    bool isHorz() { return fabs(p1[1] - p2[1]) < EPS; }

    double yxSlope() { return (p1[0] - p2[0]) / (p1[1] - p2[1]); }
    double yzSlope() { return (p1[2] - p2[2]) / (p1[1] - p2[1]); }

    bool contains(double ys)
    {
        if(p1[1] < ys && p2[1] < ys) return false;
        if(p1[1] > ys && p2[1] > ys) return false;

        return true;
    }

    Point interpolate(double ys)
    {
        double x = p1[0] + yxSlope() * (ys - p1[1]);
        double y = ys;
        double z = p1[2] + yzSlope() * (ys - p1[1]);

        return Point(x, y, z);
    }
};

struct Triangle
{
    Point points[3];
    int color[3];

    void setColor()
    {
        for (int i = 0; i < 3; i++)
            color[i] = uniform_int_distribution<int>(0, 255)(rng);
    }

    Point &operator[](int idx) { return points[idx]; }

    double max_y() { return max(points[0][1], max(points[1][1], points[2][1])); }
    double min_y() { return min(points[0][1], min(points[1][1], points[2][1])); }

    Point ched(double ys, bool daan)
    {
        Segment sides[] = {Segment(points[0], points[1]), 
                            Segment(points[1], points[2]),
                            Segment(points[0], points[2])};
        
        Point res = Point((daan? -1e3: 1e3), 0, 0);

        for(int i = 0; i < 3; i++)
        {
            if(sides[i].isHorz()) continue;
            if(!sides[i].contains(ys)) continue;

            Point p = sides[i].interpolate(ys);

            if(daan && p[0] > res[0]) res = p;
            else if(!daan && p[0] < res[0]) res = p;
        }

        return res;
    }
};

ostream &operator<<(ostream &dout, Triangle t)
{
    dout << t[0] << "\n"
         << t[1] << "\n"
         << t[2] << endl;
    return dout;
}

istream &operator>>(istream &din, Triangle &t)
{
    din >> t[0] >> t[1] >> t[2];
    return din;
}

#endif