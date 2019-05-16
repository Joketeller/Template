#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<queue>
#include<vector>
#include<iomanip>

#define ll long long
#define scan(x) scanf("%d",&x)
#define dscan(x,y) scanf("%d%d",&x,&y)
using namespace std;
const int  MOD=1e9+7;
const int MAXN=1e5+100;
struct node
{
	int m,n,id;
	bool operator < (const node& b)
	{
	   	return n<b.n;
	}
} query[MAXN];
vector<node> left_area[MAXN];
ll lad_and[MAXN];
ll inv_and[MAXN];
int ans_to_query[MAXN];
int in_chunk[MAXN];
const int MAXX=1e5;
ll qpow(ll base, ll po)
{
	ll ans=1;
	for (;po;base*=base,base%=MOD,po>>=1)
	{
		if (po&1) 
		{
			ans*=base;
			ans%=MOD;
		}
	}
	return ans;
}

ll C(int n,int m)
{
	return (lad_and[n] * inv_and[n-m] %MOD *inv_and[m]) %MOD;
}


int main(int argc, char const *argv[])
{
	lad_and[0]=1;
	for (int i=1;i<=MAXX;i++)
	{
		lad_and[i]=lad_and[i-1]*(ll)i % MOD;
	}
	inv_and[MAXX]=qpow(lad_and[MAXX],MOD-2);
	for (int i=MAXX-1;~i;i--)
		inv_and[i]=inv_and[i+1]*(i+1)%MOD;
	int chunk=(int)sqrt((double)MAXX);
	int cnt=1;
	 for (int i = 1; i <= MAXX; i += chunk, ++ cnt)
        for (int j = i; j < i + chunk && j <= MAXX; ++ j)
            in_chunk[j] = cnt;
    cnt --;
   int T;
   scan(T);
   for (int i=1;i<=T;i++)
   {
   	  dscan(query[i].n,query[i].m);
   	  query[i].id=i;
   	  left_area[in_chunk[query[i].m]].push_back(query[i]);
   	  
   }
 //  sort(query+1,query+T+1);
   for (int i=1;i<=cnt;i++)
   	if (left_area[i].size())
   	{
   		sort(left_area[i].begin(),left_area[i].end());
   		int val = 0, tmpn = left_area[i][0].n;
   		int tmpm = -1;
   		for (int j=0;j<left_area[i].size();j++)
   		{
   			while (tmpn<left_area[i][j].n)
   			{
   				val=(MOD-C(tmpn,tmpm)+val+val)%MOD;
   				tmpn++;
   			}
   			while (tmpm<left_area[i][j].m)
   			{
   				tmpm++;
   				val+=C(tmpn,tmpm);
   				val%=MOD;
   			}
   			while (tmpm>left_area[i][j].m)
   			{
   				val=(val+MOD-C(tmpn,tmpm))%MOD;
   				tmpm--;
   			}
   			ans_to_query[left_area[i][j].id]=val;
   		}
   	}
   	for (int i=1;i<=T;i++)
   		cout<<ans_to_query[i]<<endl;
	return 0;
}