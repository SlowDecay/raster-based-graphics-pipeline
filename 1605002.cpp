#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#include "triangle.hpp"
#include "matrix.hpp"
#include "utils.hpp"
using namespace std;

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
            Triangle t;
            fin >> t;

            Matrix gunHobe = bases.top();
            t = gunHobe * t;

            fout << t << endl;
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

    Triangle t;
    while (fin >> t)
    {
        t = v * t;
        fout << t << endl;
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

    Triangle triangle;
    while (fin >> triangle)
    {
        triangle = pr * triangle;
        fout << triangle << endl;
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
    for (int i = 0; i < screen_height; i++)
        z_buffer[i].resize(screen_width, z_max);

    bitmap_image image(screen_width, screen_height);
    for (int j = 0; j < screen_height; j++)
    {
        for (int i = 0; i < screen_width; i++)
            image.set_pixel(i, j, 0, 0, 0);
    }

    double dx = (box_right - box_left) / screen_width;
    double dy = (box_up - box_down) / screen_height;

    double top_y = box_up - dy / 2;
    double bottom_y = box_down + dy / 2;
    double left_x = box_left + dx / 2;
    double right_x = box_right - dx / 2;

    Triangle t;
    while (fin >> t)
    {
        t.setColor();

        double dcord = max(t.min_y(), bottom_y), ucord = min(t.max_y(), top_y);
        int upor = round((top_y - ucord) / dy), nich = round((top_y - dcord) / dy);

        double ys = top_y - upor * dy;
        for (int row = upor; row <= nich; row++)
        {
            Point lp = t.ched(ys, false);
            Point rp = t.ched(ys, true);

            ys -= dy;

            if (lp[0] > rp[0])
                continue;

            double dhaal = (rp[2] - lp[2]) / (rp[0] - lp[0]);
            double zp = lp[2];

            if (lp[0] < left_x)
                zp += (left_x - lp[0]) * dhaal, lp[0] = left_x;
            if (rp[0] > right_x)
                rp[0] = right_x;

            int baam = round((lp[0] - left_x) / dx), daan = round((rp[0] - left_x) / dx);

            for (int col = baam; col <= daan; col++)
            {
                if (zp >= z_min && z_buffer[row][col] > zp)
                {
                    z_buffer[row][col] = zp;
                    image.set_pixel(col, row, t.color[0], t.color[1], t.color[2]);
                }

                zp += dx * dhaal;
            }
        }
    }

    for (int i = 0; i < screen_height; i++)
    {
        for (int j = 0; j < screen_width; j++)
        {
            if (z_buffer[i][j] < z_max)
                fout << setprecision(6) << fixed << z_buffer[i][j] << "    ";
        }
        fout << endl;
    }
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