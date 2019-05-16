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

struct Poly //多边形类
{
    vector<Point>p; //顺时针凸包
    vector<Point>tb;// 逆时针凸包
    void in(const int &r)
    {
        p.resize(r);  //不早凸包的时候可以把p改为a
        for(int i = 0; i < r; i++) p[i].in();
    }
    //判断点集是否为凸包(返回m-1==n),或者用凸包点算出凸包顶点tb(本题即是)
    void isCanHull()
    {
        sort(p.begin(), p.end());
        p.erase(unique(p.begin(), p.end()), p.end());
        int n = p.size();
        tb.resize(n * 2 + 5);
        int m = 0;
        for(int i = 0; i < n; i++)
        {
            while(m > 1 && sgn(det(tb[m - 1] - tb[m - 2], p[i] - tb[m - 2])) <= 0)m--;
            tb[m++] = p[i];
        }
        int k = m;
        for(int i = n - 2; i >= 0; i--)
        {
            while(m > k && sgn(det(tb[m - 1] - tb[m -2], p[i] - tb[m - 2])) <= 0)m--;
            tb[m++] = p[i];
        }
        tb.resize(m);
        if(m > 1)tb.resize(m - 1);
        //for(int i = 0; i < m - 1; i++) tb[i].out();
    }  
    //旋转卡壳算法：输入凸包，输入最远两个点的坐标及最远距离
    int maxdist()//int &first,int &second,若要输出坐标可放入
    {
        int n=tb.size(),first,second;
        int Max=-INF;
        if(n==1) return Max;//first=second=0;
        for(int i=0,j=1;i<n;i++)
        {
            while(sgn(det(tb[(i+1)%n]-tb[i],tb[j]-tb[i])-det(tb[(i+1)%n]-tb[i],tb[(j+1)%n]-tb[i]))<0)
                j=(j+1)%n;
            int d=dist(tb[i],tb[j]);
            if(d>Max) Max=d;//first=i,second=j;
            d=dist(tb[(i+1)%n],tb[(j+1)%n]);
            if(d>Max) Max=d;//first=i,second=j;
        }
        return Max;
    }
}poly;