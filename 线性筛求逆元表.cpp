//线性筛求逆元表
const ll p=mod;
ll inv[p-1]; 

void solve()
{
   inv[1]=1;
   for (int i=2;i<=p-1;i++)
   {
   	   ll a=p/i;
	   ll b=p%i;
	   inv[i]=(-a*inv[b])%p;
   }
   return;
}

