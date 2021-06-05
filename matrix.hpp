#ifndef MATRIX // #include guards
#define MATRIX

#include<bits/stdc++.h>
#include "point.hpp"
#include "triangle.hpp"
using namespace std;

struct Matrix
{
    vector<vd> vals;

    Matrix() : vals(K)
    {
        for (int i = 0; i < K; i++)
            vals[i].resize(K);
    }

    vd &operator[](int idx) { return vals[idx]; }

    Matrix operator*(Matrix &rhs)
    {
        Matrix res;

        for (int i = 0; i < K; i++)
        {
            for (int j = 0; j < K; j++)
            {
                for (int k = 0; k < K; k++)
                    res[i][j] += vals[i][k] * rhs[k][j];
            }
        }

        return res;
    }

    Point operator*(Point &rhs)
    {
        Point res;
        for (int i = 0; i < K; i++)
        {
            res[i] = 0;
            for (int j = 0; j < K; j++)
                res[i] += vals[i][j] * rhs[j];
        }

        return res.normWeight();
    }

    Triangle operator*(Triangle &rhs)
    {
        Triangle res;
        for(int i = 0; i < 3; i++) res[i] = (*this)*rhs[i];

        return res;
    }
} idt;

ostream &operator<<(ostream &dout, Matrix m)
{
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < K; j++)
            dout << setprecision(7) << fixed << m[i][j] << (j == K - 1 ? "\n" : " ");
    }

    return dout;
}

#endif