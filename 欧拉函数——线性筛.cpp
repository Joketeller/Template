const int maxn=;

int prime[maxn],cnt,euler[maxn];

void get_euler(int n)
{
	memset(euler,0,sizeof(euler));
	euler[1]=1;
    cnt=0;
    for (int i=2;i<=n;i++)
    {
    	if (!euler[i])
    	{
    		euler[i]=i-1;
    		prime[++cnt]=i;
		}
		for (int j=1;j<=cnt && i*prime[j]<=n;j++)
		{
			if (i%prime[j]) euler[prime[j]*i]=euler[i]*(prime[j]-1);
			else 
			{
				euler[prime[j]*i]=euler[i]*prime[j];
				break;
			}
		}
	}
	  
}




