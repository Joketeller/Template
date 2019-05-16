// point Plus
int sgn(const double &x){ return x < -eps? -1 : (x > eps);}
inline double sqr(const double &x){ return x * x;}
struct Point
{
    double x, y;
    Point(const double &x = 0, const double &y = 0):x(x), y(y){}
    Point operator -(const Point &a)const{ return Point(x - a.x, y - a.y);}
    Point operator +(const Point &a)const{ return Point(x + a.x, y + a.y);}
    Point operator *(const double &a)const{ return Point(x * a, y * a);}
    Point operator /(const double &a)const{ return Point(x / a, y / a);}
    bool operator < (const Point &a)const{ return sgn(x - a.x) < 0 || (sgn(x - a.x) == 0 && sgn(y - a.y) < 0);}
    bool operator == (const Point &a)const{ return sgn(sgn(x - a.x) == 0 && sgn(y - a.y) == 0);}
    friend double det(const Point &a, const Point &b){ return a.x * b.y - a.y * b.x;}
    friend double dot(const Point &a, const Point &b){ return a.x * b.x + a.y * b.y;}
    friend double dist(const Point &a, const Point &b){ return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));}
    void in(){ scanf("%lf %lf", &x, &y); }
    void out()const{ printf("%lf %lf\n", x, y); }
};