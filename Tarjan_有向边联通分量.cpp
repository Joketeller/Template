stack<int> s; //栈
int dfn[MAXN], low[MAXN], scc[MAXN];  //scc标明x节点属于哪个联通分量 dfn为dfs序 
bool ins[MAXN];   //标明x是不是在栈中  
vector<int> g[MAXN];
int num = 0, tot = 0;

int addedge(int u,int v)
{
    g[u].push_back(v);
}

int dfs(int x)
{
	if (dfn[x])
	{
		if (ins[x] == 0)
			return INF;
		return low[x];
	}
	
	low[x] = dfn[x] = ++num;   //为节点x设定次序编号和low初值
	ins[x] = 1;
	s.push(x);
	for (int i = 0; i < g[x].size(); i ++)
		low[x] = min(low[x], dfs(g[x][i]));
	if (low[x] == dfn[x])
	{
		while (dfn[s.top()] != low[s.top()])
		{
			ins[s.top()] = 0;
			scc[s.top()] = tot;
			s.pop();
		}
		ins[s.top()] = 0;
		scc[s.top()] = tot++;
		s.pop();
	}
	return low[x];
}

主函数内：
	for (int i = 1; i <= n; i ++)
		if (!dfn[i])
			dfs(i);

//算法复杂度O(N+M)



/*****非递归版本*****/
stack<int> ss;
int nx[Maxn]; //下一个儿子
int lx[Maxn]; //上一个儿子
void tarjan(int u){
    while(!ss.empty()) ss.pop();
    ss.push(u);
    for(int i=1;i<=n;i++) nx[i]=head[i],lx[i]=-1;
    while(!ss.empty()){
        int v=ss.top();
        if(!dfn[v]){ //第一次访问
            dfn[v]=low[v]=++tmpdfn;
            st[++top]=v;
            in[v]=1;
        }
        if(lx[v]!=-1) low[v]=min(low[v],low[lx[v]]); //访问完儿子
        if(nx[v]!=-1){ //有儿子
            while(nx[v]!=-1){ //寻找下一个还未访问的儿子
                if(!dfn[p[nx[v]].to]){ //树边
                    lx[v]=p[nx[v]].to;
                    nx[v]=p[nx[v]].next;
                    ss.push(lx[v]);
                    break;
                }
                else if(in[p[nx[v]].to]) //回边
                    low[v]=min(low[v],dfn[p[nx[v]].to]);
                nx[v]=p[nx[v]].next;
            }
        }
        else{ //全部儿子访问完毕
            if(low[v]==dfn[v]){
                scc++;
                do{
                    belong[st[top]]=scc;
                    in[st[top]]=0;
                    cnt[scc]+=val[st[top]];
                }while(st[top--]!=v);
            }
            ss.pop();
        }
    }
}
