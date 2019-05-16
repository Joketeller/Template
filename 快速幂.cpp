#define ll long long 

ll qpow(ll a,ll b,ll p)   //  a^b%p   
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


