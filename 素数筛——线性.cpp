const int MAXN=;

bool vis[MAXN];
int prime[MAXN],cnt;
void get_prime(int n)
{
	memset(vis,0,sizeof(vis));
	cnt=0;
	for (int i=2;i<=n;i++)
	{
	 if(!vis[i])
	 {
       prime[++cnt]=i;
       vis[i]=1;
      }
 	 for (int j=1;j<=cnt && (prime[j]*i<=n);j++)
	 	{
	 		vis[i*prime[j]]=1;
	 		if (i%prime[j]==0)
	 		 break;
		 }	
	 }
}




