#define MAX                                //点的个数
#define INF 0x3f3f3f3f

struct Edge
{
	int u,v,c,f;
} ;

int m,n,sum;
struct Dinic
{
	int s,t;
	vector<Edge> E;
	vector<int> G[MAX]; 
	bool vis[MAX];
	int lev[MAX];
	int cur[MAX];
	void init(int l,int r)
	{
		E.clear();
		for (int i=l;i<=r;i++)
			G[i].clear();
	}
	void addedge(int from,int to,int cap)
	{
		E.push_back((Edge){from,to,cap,0});
		E.push_back((Edge){to,from,0,0});
		int m=E.size();
		G[from].push_back(m-2);
		G[to].push_back(m-1);
	}
	bool bfs()
	{
		me0(vis);
		queue<int>  q;
		q.push(s);
		lev[s]=0;
		vis[s]=1;
		while (!q.empty())
		{
			int now=q.front();
			q.pop();
			for (int i=0,_size=G[now].size();i<_size;i++)
			{
				Edge edge=E[G[now][i]];
				int nex=edge.v;
				if (!vis[nex] && edge.c>edge.f)
				{
					lev[nex]=lev[now]+1;
					q.push(nex);
					vis[nex]=1;
				}
			}
		}
		return vis[t];
	}
	int dfs(int now,int aug)
	{
		if (now==t || aug==0) return aug;
		int flow=0,f;
		for (int &i=cur[now],_size=G[now].size();i<_size;i++)
		{
			Edge& edge=E[G[now][i]];
			int nex=edge.v;
			if (lev[now]+1 ==lev[nex] && (f=dfs(nex,min(aug,edge.c-edge.f)))>0)
			{
				edge.f+=f;
				E[G[now][i]^1].f-=f;
				flow+=f;
				aug-=f;
				if (!aug) break;
			}
		}
		return flow;
	}
	int maxflow()
	{
		int  flow=0;
		while (bfs())
		{
			me0(cur);
			flow+=dfs(s,INF);
		}
		return flow;
	}
}dinic;
