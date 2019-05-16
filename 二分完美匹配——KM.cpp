
struct KM
{
	int linker[MAX],lx[MAX],ly[MAX];
	bool visx[MAX],visy[MAX];
	int map[MAX][MAX];
	int n; 
	bool hungary(int u)
	{
			visx[u]=true;
		for(int i=0;i<n;i++)
			if(!visy[i]&&lx[u]+ly[i]==map[u][i])
			{//满足最大匹配 
				visy[i]=true;
				if(linker[i]==-1||hungary(linker[i])){
					linker[i]=u;
					return true;
				}	
			}
		return false;
	}
	void K_M()
	{
		int temp;
		memset(lx,0,sizeof(lx));
		memset(ly,0,sizeof(ly));
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				lx[i]=max(lx[i],map[i][j]);
			}
		}
		for(int i=0;i<n;i++)
		{
			while(1)
			{
				memset(visx,0,sizeof(visx));
				memset(visy,0,sizeof(visy));
				if(hungary(i))
					break;
				else{
					temp=INT_MAX;
						for(int j=0;j<n;j++) 
						if(visx[j])for(int k=0;k<n;k++){//x在交叉树内
							if(!visy[k]&&temp>lx[j]+ly[k]-map[j][k]){//y在外 
								temp=lx[j]+ly[k]-map[j][k];
							}
						}
						for(int j=0;j<n;j++){
							if(visx[j])
							lx[j]-=temp;
							if(visy[j])
							ly[j]+=temp;
						}
					}
			}
		}
	}
	void solve()
	{
		while(scanf("%d",&n)!=EOF)
		{
			memset(linker,-1,sizeof(linker));
			for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
			scanf("%d",&map[i][j]);
			K_M();
			int ans=0;
			for(int i=0;i<n;i++)
				ans+=map[linker[i]][i];
			printf("%d\n",ans);
		}
	}
}km;


