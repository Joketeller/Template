//a 为 余数， m为模数
#include <cstdio>
#include <algorithm>
#define maxn 100000
#define ll long long
using namespace std;
ll a[maxn], m[maxn], N;
void exgcd(ll a, ll b, ll &x, ll &y)
{
	if(!b){x=1,y=0;return;}
	ll xx, yy;
	exgcd(b,a%b,xx,yy);
	x=yy, y=xx-a/b*yy;
}
ll gcd(ll a, ll b){return !b?a:gcd(b,a%b);}
ll lcm(ll a, ll b){return a/gcd(a,b)*b;}
ll CRT(ll *a, ll *m, ll n)
{
	ll a1, m1, a2, m2, d, x1, x2;
	a1=a[1], m1=m[1];
	for(ll i=2;i<=n;i++)
	{
		a2=a[i], m2=m[i];
		d=gcd(m1,m2);
		if((a2-a1)%d!=0)return -1;
		exgcd(m1,m2,x1,x2);
		x1=(a2-a1)/d*x1%m2;
		a1=(a1+x1*m1)%lcm(m1,m2);
		m1=lcm(m1,m2);
	}
	return (a1+m1)%m1;
}
int main()
{
	ll i, N;
	while(~scanf("%lld",&N))
	{
		for(i=1;i<=N;i++)scanf("%lld%lld",m+i,a+i);
		printf("%lld\n",CRT(a,m,N));
	}
	return 0;
}

