#define MX 70 
int mu[MX];
int prime[MX];
bool vis[MX];
void getmu()
{
	memset(vis,0,sizeof(vis));
	int cnt = 0;
	mu[1] = 1;
	for (int i = 2; i < MX ; ++i)
	{
		if(!vis[i])
		{
			prime[cnt++] = i;
			mu[i] = -1;
		}
		for (int j = 0; j < cnt&&i*prime[j]<MX ; ++j)
		{
			vis[i*prime[j]] = true;
			if(i%prime[j]) mu[i*prime[j]] = -mu[i];
			else
			{
				mu[i*prime[j]] = 0;
				break;
			}
		}
	}
}