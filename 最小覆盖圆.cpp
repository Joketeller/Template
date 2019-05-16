/*最小圆覆盖*/
/*给定n个点, 让求半径最小的圆将n个点全部包围,可以在圆上*/
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#define EPS 1e-8
using namespace std;
const int maxn = 550;
struct point{
    double x, y;
};
int sgn(double x)
{
    if (fabs(x) < EPS)
        return 0;
    return x < 0 ? -1 : 1;
}
double get_distance(const point a, const point b)//两点之间的距离
{
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}
point get_circle_center(const point a, const point b, const point c)//得到三角形外接圆的圆心
{
    point center;
    double a1 = b.x - a.x;
    double b1 = b.y - a.y;
    double c1 = (a1 * a1 + b1 * b1) / 2.0;
    double a2 = c.x - a.x;
    double b2 = c.y - a.y;
    double c2 = (a2 * a2 + b2 * b2) / 2.0;
    double d = a1 * b2 - a2 * b1;
    center.x = a.x + (c1 * b2 - c2 * b1) / d;
    center.y = a.y + (a1 * c2 - a2 * c1) / d;
    return center;
}
//p表示定点, n表示顶点的个数, c代表最小覆盖圆圆心, r是半径
void min_cover_circle(point *p, int n, point &c, double &r)//找最小覆盖圆(这里没有用全局变量p[], 因为是为了封装一个函数便于调用)
{
    random_shuffle(p, p + n);//随机函数,使用了之后使程序更快点,也可以不用
    c = p[0];
    r = 0;
    for (int i = 1; i < n; i++)
    {
        if (sgn(get_distance(p[i], c) - r) > 0)//如果p[i]在当前圆的外面, 那么以当前点为圆心开始找
        {
            c = p[i];//圆心为当前点
            r = 0;//这时候这个圆只包括他自己.所以半径为0
            for (int j = 0; j < i; j++)//找它之前的所有点
            {
                if (sgn(get_distance(p[j], c) - r) > 0)//如果之前的点有不满足的, 那么就是以这两点为直径的圆
                {
                    c.x = (p[i].x + p[j].x) / 2.0;
                    c.y = (p[i].y + p[j].y) / 2.0;
                    r = get_distance(p[j], c);
                    for (int k = 0; k < j; k++)
                    {
                        if (sgn(get_distance(p[k], c) - r) > 0)//找新作出来的圆之前的点是否还有不满足的, 如果不满足一定就是三个点都在圆上了
                        {
                            c = get_circle_center(p[i], p[j], p[k]);
                            r = get_distance(p[i], c);
                        }
                    }
                }
            }
        }
    }
}
int main()
{
    int n;
    point p[maxn];
    point c; double r;
    while (~scanf("%d", &n) && n)
    {
        for (int i = 0; i < n; i++)
            scanf("%lf %lf", &p[i].x, &p[i].y);
        min_cover_circle(p, n, c, r);
        printf("%.2lf %.2lf %.2lf\n", c.x, c.y, r);
    }
    return 0;
}
