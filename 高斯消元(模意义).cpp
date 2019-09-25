typedef long long ll;
const int maxn = 75;
ll abs[maxn][maxn];
const ll p = 1e9 + 7;
ll qpow(ll a, ll b) //  a^b%p
{
	a%=p;
	ll ans=1;
	for (;b;b>>=1)
	{
		if (b&1)
		{
			ans*=a;
			ans%=p;
		}
		a*=a;	
		a%=p;
	}
	return ans;
}

void solve(int _n)
{
    for( int i=1;i<=_n;i++)
	{
		if(!a[i][i])//主元不能为0
		{
			maxi=0;
			for(int j=i+1;j<=_n&&!maxi;j++)
				if(a[j][i]) maxi=j;
			if(!maxi) continue;//如果一整列都为0，不需要消元
			for(int j=i;j<=_n+1;j++)
				tmp=a[maxi][j],a[maxi][j]=a[i][j],a[i][j]=tmp;
		}
		for(int j=i+1;j<=_n;j++)
		{
			tmp=a[j][i];
			if(!tmp) continue;//已经为0，不需要消元
			for(int k=i;k<=_n+1;k++)
				a[j][k]=((a[j][k]*a[i][i]-a[i][k]*tmp)%p+p)%p;
		}
    }
    for(int i=_n;i;i--)
	{
		for(int j=i+1;j<=_n;j++)
			a[i][_n+1]=((a[i][_n+1]-ans[j]*a[i][j])%p+p)%p;
		ans[i]=a[i][_n+1]*qpow(a[i][i],p-2)%p;
	}
}