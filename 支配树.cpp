/********************** DAG ************************/
struct Graph  //有一个虚拟点0,可以在树中DFS来获得答案
{
static const int MAXN = 100050;
static const int MAXM = 200100;
static const int DC =   18; //log2(MAXN)
    int to[MAXM],nxt[MAXM],head[MAXN],tot,st[MAXN],N,in[MAXN],toplist[MAXN],listnum,fa[MAXN][DC],depth[MAXN];
    vector<int> pre[MAXN],tr[MAXN];
    void init(int _N)
    {
        for (int i = 0; i <= _N;i++)
            head[i] = 0, pre[i].clear(),tr[i].clear();
        tot = 0;
        N = _N;
        depth[0]=0;
        for (int i=0;i<DC;i++)
            fa[0][i]=0;
    }

    void addedge(int u,int v)
    {
        to[++tot] = v;
        nxt[tot] = head[u];
        head[u] = tot;
        pre[v].push_back(u);
        ++in[v];
    }

    bool topsort()
    {
        int top = 1;
        listnum=0;
        for (int i = 1; i <= N;i++)
            if (!in[i])
            {
                addedge(0,i);
            }
        st[top]=0;
        while (top)
        {
            int u=st[top];
            --top;
            for (int ed=head[u];ed;ed=nxt[ed])
            {
                --in[to[ed]];
                if (!in[to[ed]])
                    st[++top]=toplist[++listnum]=to[ed];
            }
        }
        return listnum==N;
    }
    int LCA(int x,int y)
    {
        if (depth[x]<depth[y]) 
            swap(x,y);
        for (int i=DC-1;i>=0;i--)
            if (depth[fa[x][i]]>=depth[y])
                x=fa[x][i];
        if (x==y)
            return x;
        for (int i=DC-1;i>=0;i--)
            if (fa[x][i]!=fa[y][i])
                x=fa[x][i],y=fa[y][i];
        return fa[x][0];
    }

    inline void update(int father,int son,int distance=1)
    {
        tr[father].push_back(son);
        depth[son]=depth[father]+distance;
        fa[son][0]=father;
        for (int i=1;i<DC;i++) 
            fa[son][i]=fa[fa[son][i-1]][i-1];   
    }

    void build_dominator_tree()
    {
        for (int i=1;i<=listnum;i++)
        {
            int u=toplist[i];
            int fa=pre[u][0];
            for (int i=1;i<pre[u].size();i++)
                fa=LCA(fa,pre[u][i]);
            update(fa,u);
        }
    }

    void auto_build()
    {
        topsort();
        build_dominator_tree();
    }
}DT;


/************************ 一般有向图 **************************/
#include<bits/stdc++.h>
using namespace std;
#define RI register int
int read() {
	int q=0;char ch=' ';
	while(ch<'0'||ch>'9') ch=getchar();
	while(ch>='0'&&ch<='9') q=q*10+ch-'0',ch=getchar();
	return q;
}
typedef long long LL;
const int N=50005,M=100005;
int n,m,tim;
int dfn[N],repos[N],mi[N],fa[N],f[N],semi[N],idom[N],ans[N];
struct graph{
	int tot,h[N],ne[M],to[M];
	void clear() {tot=0;for(RI i=0;i<=n;++i) h[i]=0;}
	void add(int x,int y) {to[++tot]=y,ne[tot]=h[x],h[x]=tot;}
}g,rg,ng,tr;

void init() {
	tim=0;g.clear(),rg.clear(),ng.clear(),tr.clear();
	for(RI i=1;i<=n;++i)
		repos[i]=dfn[i]=idom[i]=fa[i]=ans[i]=0,mi[i]=semi[i]=f[i]=i;
}
void tarjan(int x) {
	dfn[x]=++tim,repos[tim]=x;
	for(RI i=g.h[x];i;i=g.ne[i])
		if(!dfn[g.to[i]]) fa[g.to[i]]=x,tarjan(g.to[i]);
}
int find(int x) {
	if(x==f[x]) return x;
	int tmp=f[x];f[x]=find(f[x]);
	if(dfn[semi[mi[tmp]]]<dfn[semi[mi[x]]]) mi[x]=mi[tmp];
	return f[x];
}
void dfs(int x,LL num) {
	ans[x]=num+x;
	for(RI i=tr.h[x];i;i=tr.ne[i]) dfs(tr.to[i],num+x);
}
void work() {
	for(RI i=n;i>=2;--i) {
		int x=repos[i],tmp=n;
		for(RI j=rg.h[x];j;j=rg.ne[j]) {
			if(!dfn[rg.to[j]]) continue;//此题数据有误
			if(dfn[rg.to[j]]<dfn[x]) tmp=min(tmp,dfn[rg.to[j]]);
			else find(rg.to[j]),tmp=min(tmp,dfn[semi[mi[rg.to[j]]]]);
		}
		semi[x]=repos[tmp],f[x]=fa[x],ng.add(semi[x],x);
		
		x=repos[i-1];
		for(RI j=ng.h[x];j;j=ng.ne[j]) {
			int y=ng.to[j];find(y);
			if(semi[mi[y]]==semi[y]) idom[y]=semi[y];
			else idom[y]=mi[y];//此时idom[mi[y]]可能并未找到
		}
	}
	for(RI i=2;i<=n;++i) {
		int x=repos[i];
		if(idom[x]!=semi[x]) idom[x]=idom[idom[x]];
		tr.add(idom[x],x);
	}
	dfs(n,0);
}
int main()
{
	int x,y;
	while(~scanf("%d%d",&n,&m)) {
		init();
		for(RI i=1;i<=m;++i)
			x=read(),y=read(),g.add(x,y),rg.add(y,x);
		tarjan(n);work();
		for(RI i=1;i<n;++i) printf("%d ",ans[i]);
		printf("%d\n",ans[n]);
	}
    return 0;
}