/**********************************编号从一开始，注意二分图对称结果要/2*****************/
#define maxm                                                               //边的个数
#define maxn                                                               //点的个数

struct Edge
{
    int to,next;
};

struct Hun
{
	Edge edge[maxm];
	int head[maxn],tot;
	int uN,N;
	void init()
	{
	    tot =0;
	    memset(head,-1,sizeof(head));
	    return ;
	}

	void addedge(int u,int v)
	{
	    edge[tot].to = v;
	    edge[tot].next = head[u];
	    head[u]=tot++;
	    return ;
	}

	int linker[maxn];
	bool used[maxn];


	bool dfs(int u)
	{
	    for(int i=head[u];i!=-1;i=edge[i].next)
	    {
	        int v =  edge[i].to;
	        if(!used[v])
	        {
	            used[v]=true;
	            if(linker[v]==-1||dfs(linker[v]))
	            {
	                linker[v]=u;
	                return true;
	            }
	        }
	    }
	    return false;
	}

	int hungary()
	{
	    int res = 0;
	    memset(linker,-1,sizeof(linker));
	    for(int u = 1;u<=uN;u++)                                          //编号从1开始
	    {
	        memset(used,false,sizeof(used));
	        if(dfs(u))
	        {
	            res++;
	        }
	    }
	    return res;                                                         //对称时要/2
	}
} hun;

