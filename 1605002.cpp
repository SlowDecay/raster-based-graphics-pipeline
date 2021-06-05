#include <bits/stdc++.h>
#include "bitmap_image.hpp"
using namespace std;
typedef vector<double> vd;

const int K = 4;
const double EPS = 1e-6;
const double PI = acos(-1);

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

// a is the point to be rotated
// axis is the rotation axis. It should be normalized
// angle is the rotation angle (in degrees)
Point rodrigues(Point a, Point axis, double angle)
{
    angle *= PI / 180; // taking angle to radian

    Point res = cos(angle) * a + (1 - cos(angle)) * axis.dot(a) * axis + sin(angle) * axis.cross(a);
    return res;
}

struct Triangle
{
    Point points[3];
    int color[3];

    Triangle()
    {
        for (int i = 0; i < 3; i++)
            color[i] = rand() % 256;
    }

    Point &operator[](int idx) { return points[idx]; }
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

        if (cmd == "triangle")
        {
            Point p1, p2, p3;
            fin >> p1 >> p2 >> p3;

            Matrix gunHobe = bases.top();

            p1 = gunHobe * p1;
            p2 = gunHobe * p2;
            p3 = gunHobe * p3;

            fout << p1 << endl;
            fout << p2 << endl;
            fout << p3 << endl
                 << endl;
        }
        else if (cmd == "translate")
        {
            Point t;
            fin >> t;

            Matrix m = getTransMatrix(t);
            m = bases.top() * m;
            bases.push(m);
        }
        else if (cmd == "scale")
        {
            Point s;
            fin >> s;

            Matrix m = getScaleMatrix(s);
            m = bases.top() * m;
            bases.push(m);
        }
        else if (cmd == "rotate")
        {
            double angle;
            Point axis;
            fin >> angle >> axis;

            Matrix m = getRotMatrix(angle, axis);
            m = bases.top() * m;
            bases.push(m);
        }
        else if (cmd == "push")
            szs.push(bases.size());
        else if (cmd == "pop")
        {
            if (szs.empty())
                continue;

            while (bases.size() > szs.top())
                bases.pop();
            szs.pop();
        }
        else if (cmd == "end")
            break;
        else
            assert(false);
    }

    fin.close();
    fout.close();
}

void transform_view()
{
    Point eye, look, up;

    ifstream fin("data/scene.txt");
    ofstream fout("data/stage2.txt");

    fin >> eye >> look >> up;

    fin.close();

    fin.open("data/stage1.txt");

    Point l = look - eye;
    l.normDist();
    Point r = l.cross(up);
    r.normDist();
    Point u = r.cross(l);

    Matrix tr = getTransMatrix(-eye);

    Matrix rot = idt;
    for (int j = 0; j < K - 1; j++)
        rot[0][j] = r[j], rot[1][j] = u[j], rot[2][j] = -l[j];

    Matrix v = rot * tr;

    Point p1, p2, p3;
    while (fin >> p1 >> p2 >> p3)
    {
        p1 = v * p1;
        p2 = v * p2;
        p3 = v * p3;

        fout << p1 << endl;
        fout << p2 << endl;
        fout << p3 << endl
             << endl;
    }

    fin.close();
    fout.close();
}

void transform_projection()
{
    Point eye, look, up;

    ifstream fin("data/scene.txt");
    ofstream fout("data/stage3.txt");

    fin >> eye >> look >> up;

    double fovY, aspectRatio, near, far;
    fin >> fovY >> aspectRatio >> near >> far;

    fovY *= PI / 180;

    fin.close();

    fin.open("data/stage2.txt");

    double fovX = fovY * aspectRatio;
    double t = near * tan(fovY / 2);
    double r = near * tan(fovX / 2);

    Matrix pr;

    pr[0][0] = near / r, pr[1][1] = near / t;
    pr[2][2] = -(far + near) / (far - near), pr[2][3] = -(2 * far * near) / (far - near);
    pr[3][2] = -1;

    Point p1, p2, p3;
    while (fin >> p1 >> p2 >> p3)
    {
        p1 = pr * p1;
        p2 = pr * p2;
        p3 = pr * p3;

        fout << p1 << endl;
        fout << p2 << endl;
        fout << p3 << endl
             << endl;
    }

    fin.close();
    fout.close();
}

void remove_hidden_surface()
{
    ifstream fin("data/config.txt");
    ofstream fout("data/z_buffer.txt");

    int screen_width, screen_height; // might need to change to double
    double box_left, box_right, box_up, box_down;
    double z_max, z_min;

    fin >> screen_width >> screen_height >> box_left >> box_down;
    fin >> z_min >> z_max;

    box_right = -box_left, box_up = -box_down;

    fin.close();

    fin.open("data/stage3.txt");

    vector<vd> z_buffer(screen_width);
    for(int i = 0; i < screen_width; i++) z_buffer[i].resize(screen_height, z_max);

    bitmap_image image(screen_width, screen_height);
    for(int j = 0; j < screen_height; j++)
    {
        for(int i = 0; i < screen_width; i++) image.set_pixel(i, j, 0, 0, 0);
    }

    // do things

    image.save_image("data/out.bmp");

    fin.close();
    fout.close();
}

int main()
{
    init();

    transform_model();
    transform_view();
    transform_projection();

    remove_hidden_surface();

    return 0;
}