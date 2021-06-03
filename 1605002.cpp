#include<bits/stdc++.h>
using namespace std;
typedef vector<double> vd;

const int K = 4;

struct Point
{
    vd cords;

    Point(): cords(K) {}
    Point(vd cords): cords(cords) {}

    double& operator[](int idx) { return cords[idx]; }
};

ostream& operator<<(ostream& dout, Point& p)
{
    dout << p[0] << " " << p[1] << " " << p[2];
    return dout;
}

struct Matrix
{
    vector<vd> vals;

    Matrix(): vals(K)
    {
        for(int i = 0; i < K; i++) vals[i].resize(K);
    }

    vd& operator[](int idx) { return vals[idx]; }

    Matrix operator*(Matrix& rhs)
    {
        Matrix res;

        for(int i = 0; i < K; i++)
        {
            for(int j = 0; j < K; j++)
            {
                for(int k = 0; k < K; k++) res[i][j] += vals[i][k]*rhs[k][j];
                return res;
            }
        }
    }
};

ostream& operator<<(ostream& dout, Matrix& m)
{
    for(int i = 0; i < K; i++)
    {
        for(int j = 0; j < K; j++) dout << m[i][j] << (j == K-1? "\n": " ");
    }

    return dout;
}

int main()
{
    Point p;
    p[2] = 4, p[1] = 7;

    cout << p << endl;

    Matrix m;
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1;

    cout << m << endl;

    return 0;
}