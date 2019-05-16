#include<bits/stdc++.h>
using namespace std;
inline void read(int &x)
{
    x=0;
    static int p;p=1;
    static char c;c=getchar();
    while(!isdigit(c)){if(c=='-')p=-1;c=getchar();}
    while(isdigit(c)) {x=(x<<1)+(x<<3)+(c-48);c=getchar();}
    x*=p;
}
double a,b,c,d,e,f;
double ansx,ansy,ansz;
const double eps=1e-10;
double dis(double x,double y,double z)
{
    return sqrt(x*x+y*y+z*z);
}
double calc(double x,double y)
{
    double A=c;
    double B=d*y+e*x;
    double C=a*x*x+b*y*y+f*x*y-1.0;
    double delta=B*B-4.0*A*C;
    if(delta<0)return 210000000.0;
    double x1=(-B+sqrt(delta))/(2.0*A);
    double x2=(-B-sqrt(delta))/(2.0*A);
    if(dis(x,y,x1)>dis(x,y,x2))return x2;
    return x1;
}
void MNTH()
{
    for(int times=1;times<=1;times++)
    {
        double T=10000;
        while(T>eps)
        {
            double nowx=ansx+(rand()*2-RAND_MAX)*T;
            double nowy=ansy+(rand()*2-RAND_MAX)*T;
            double nowz=calc(nowx,nowy);
            if(nowz==210000000.0){T*=0.99;continue;}
            double delta=dis(nowx,nowy,nowz)-dis(ansx,ansy,ansz);
            if(delta<0)ansx=nowx,ansy=nowy,ansz=nowz;
            else if(exp(delta/T)*RAND_MAX<rand())ansx=nowx,ansy=nowy,ansz=nowz;
            T*=0.99;
        }
    }
    printf("%.6lf\n",dis(ansx,ansy,ansz));
}
int main()
{
    srand(19890604);
    while(scanf("%lf%lf%lf%lf%lf%lf",&a,&b,&c,&d,&e,&f)!=EOF)
    {
        ansx=0;
        ansy=0;
        ansz=sqrt(1.0/c);
        MNTH();
    }
    return 0;
}
