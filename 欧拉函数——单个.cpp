// 欧拉函数——单个,O(n^(1/2))
#include <bits/stdc++.h>
#define me0(x) memset(x,0,sizeof(x))
#define ll long long
#define ull unsigned long long
#define scan(x) scanf("%d",&x)
#define scanll(x) scanf("%lld",&x)
#define dscan(x,y) scanf("%d%d",&x,&y)
#define rep(x,be,en) for (x=be;x<=en;x++)
#define fr1(n) for (int i=1;i<=n;i++)
#define fr(n) for (int i=0;i<n;i++)
using namespace std;
ll euler(ll n)
{
	ll ans=n;
	for (int i=2;i*i<=n;i++)
	{
		if (n%i==0)
		{
			ans=ans/i*(i-1);
			while (n%i==0)
			 n/=i;
		}
	}
	if (n>1) ans=ans/n*(n-1);
	return ans;
 } 
const ll mod=1e9+7;
ll num[100];
int main()
{
	int i,j;
	int top=1;
	num[top]=mod;
	for (i=2;i<=  ;  i++)
	 num[i]=euler(num[i-1]);
    for (i=1;i<=100;i++) 
    cout<<num[i]<<" ";
}

