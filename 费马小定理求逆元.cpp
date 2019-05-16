//费马小定理求逆元，要求p为素数 ,O(logp)
ll qpow(ll a,ll b,ll p)
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

ll inv(ll a,ll p) // 返回a在膜p情况下得逆元 
{
	return qpow(a,p-2);
 } 



