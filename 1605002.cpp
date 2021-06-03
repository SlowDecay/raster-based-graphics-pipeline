#include<bits/stdc++.h>
using namespace std;
typedef vector<double> vd;

const int K = 4;
const double EPS = 1e-6;

struct Point
{
    vd cords;

    Point(): cords(K) { cords[K-1] = 1; }
    Point(vd cords): cords(cords) {}

    Point& normDist()
    {
        if(fabs(cords[0]) < EPS && fabs(cords[1]) < EPS && fabs(cords[2]) < EPS) return *this;

        double d = 0;
        for(int i = 0; i < K-1; i++) d += cords[i]*cords[i];
        d = sqrt(d);

        for(int i = 0; i < K-1; i++) cords[i] /= d;

        return *this;
    }

    Point& normWeight()
    {
        for(int i = 0; i < K; i++) cords[i] /= cords[K-1];
        return *this;
    }

    double& operator[](int idx) { return cords[idx]; }
};

ostream& operator<<(ostream& dout, Point& p)
{
    dout << p[0] << " " << p[1] << " " << p[2];
    return dout;
}

istream& operator>>(istream& din, Point& p)
{
    din >> p[0] >> p[1] >> p[2];
    return din;
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
            }
        }

        return res;
    }

    Point operator*(Point& rhs)
    {
        Point res;
        for(int i = 0; i < K; i++)
        {
            for(int j = 0; j < K; j++) res[i] += vals[i][j]*rhs[j];
        }

        return res;
    }
} idt;

ostream& operator<<(ostream& dout, Matrix& m)
{
    for(int i = 0; i < K; i++)
    {
        for(int j = 0; j < K; j++) dout << m[i][j] << (j == K-1? "\n": " ");
    }

    return dout;
}

void init()
{
    for(int i = 0; i < K; i++) idt[i][i] = 1;
}

Matrix getTransMatrix(Point t)
{
    Matrix res = idt;

    for(int i = 0; i < K; i++) res[i][3] = t[i];
    res[3][3] = 1;

    return res;
}

Matrix getScaleMatrix(Point s)
{
    Matrix res = idt;

    for(int i = 0; i < K; i++) res[i][i] = s[i];
    res[3][3] = 1;

    return res;
}

Matrix getRotMatrix(double degree, Point axis)
{
    Matrix res = idt;

    axis.normDist();

    Point i; i[0] = 1;
    Point j; j[1] = 1;
    Point k; k[2] = 1;



    return res;
}

void transform_model()
{
    Point eye, look, up;

    ifstream fin("data/scene.txt");
    ofstream fout("data/stage1.txt");

    fin >> eye >> look >> up;

    double fovY, aspectRatio, near, far;
    fin >> fovY >> aspectRatio >> near >> far;

    stack<Matrix> bases;
    stack<int> szs;

    bases.push(idt);

    while (true)
    {
        string cmd;
        fin >> cmd;

        if(cmd == "triangle")
        {
            Point p1, p2, p3;
            fin >> p1 >> p2 >> p3;

            Matrix gunHobe = bases.top();

            p1 = gunHobe*p1;
            p2 = gunHobe*p2;
            p3 = gunHobe*p3;

            fout << p1 << endl;
            fout << p2 << endl;
            fout << p3 << endl << endl;
        }
        else if(cmd == "translate")
        {
            Point t;
            fin >> t;

            Matrix m = getTransMatrix(t);
            m = m*bases.top();
            bases.push(m);
        }
        else if(cmd == "scale")
        {
            Point s;
            fin >> s;

            Matrix m = getScaleMatrix(s);
            m = m*bases.top();
            bases.push(m);
        }
        else if(cmd == "rotate")
        {
            double angle;
            Point axis;
            fin >> angle >> axis;

            Matrix m = getRotMatrix(angle, axis);
            m = m*bases.top();
            bases.push(m);
        }
        else if(cmd == "push") szs.push(bases.size());
        else if(cmd == "pop")
        {
            if(szs.empty()) continue;

            while(bases.size() > szs.top()) bases.pop();
            szs.pop();
        }
        else if(cmd == "end") break;
        else assert(false);
    }

    fin.close();
    fout.close();
}

int main()
{
    init();
    transform_model();

    return 0;
}