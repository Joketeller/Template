#define MAXN                                      //点的个数
#define MAXM                                      //边的条数
#define INF 0x7fffffff
struct Edge
{
    int to, next, cap, flow, cost;
    int x, y;
};

struct Dinic_mincost
{
	Edge edge[MAXM];
	int head[MAXN],tol;
	int pre[MAXN],dis[MAXN];
	bool vis[MAXN];
	int N;
//	char map[MAXN][MAXN];
	void init()
	{
	    N = MAXN;
	    tol = 0;
	    memset(head, -1, sizeof(head));
	}
	void addedge(int u, int v, int cap, int cost)//左端点，右端点，容量，花费
	{
	    edge[tol]. to = v;
	    edge[tol]. cap = cap;
	    edge[tol]. cost = cost;
	    edge[tol]. flow = 0;
	    edge[tol]. next = head[u];
	    head[u] = tol++;
	    edge[tol]. to = u;
	    edge[tol]. cap = 0;
	    edge[tol]. cost = -cost;
	    edge[tol]. flow = 0;
	    edge[tol]. next = head[v];
	    head[v] = tol++;
	}
	bool spfa(int s, int t)
	{
	    queue<int> q;
	    for(int i = 0; i < N; i++)
	    {
	        dis[i] = INF;
	        vis[i] = false;
	        pre[i] = -1;
	    }
	    dis[s] = 0;
	    vis[s] = true;
	    q.push(s);
	    while(!q.empty())
	    {
	        int u = q.front();
	        q.pop();
	        vis[u] = false;
	        for(int i = head[u]; i != -1; i = edge[i]. next)
	        {
	            int v = edge[i]. to;
	            if(edge[i]. cap > edge[i]. flow &&
	                    dis[v] > dis[u] + edge[i]. cost )
	            {
	                dis[v] = dis[u] + edge[i]. cost;
	                pre[v] = i;
	                if(!vis[v])
	                {
	                    vis[v] = true;
	                    q.push(v);
	                }
	            }
	        }
	    }
	    if(pre[t] == -1) return false;
	    else return true;
	}
	//返回的是最大流， cost存的是最小费用
	int minCostMaxflow(int s, int t, int &cost)
	{
	    int flow = 0;
	    cost = 0;
	    while(spfa(s,t))
	    {
	        int Min = INF;
	        for(int i = pre[t]; i != -1; i = pre[edge[i^1]. to])
	        {
	            if(Min > edge[i]. cap - edge[i]. flow)
	                Min = edge[i]. cap - edge[i]. flow;
	        }
	        for(int i = pre[t]; i != -1; i = pre[edge[i^1]. to])
	        {
	            edge[i]. flow += Min;
	            edge[i^1]. flow -= Min;
	            cost += edge[i]. cost * Min;
	        }
	        flow += Min;
	    }
	    return flow;
	}
} dinic_mincost;
