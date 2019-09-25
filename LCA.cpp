/*   tarjan 离线 O((n+q)) 实际上还有并查集饿的复杂度 */
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <cmath>
using namespace std;
const int N=40000+5;
struct Edge{
    int cnt,x[N],y[N],z[N],nxt[N],fst[N];
    void set(){
        cnt=0;
        memset(x,0,sizeof x);
        memset(y,0,sizeof y);
        memset(z,0,sizeof z);
        memset(nxt,0,sizeof nxt);
        memset(fst,0,sizeof fst);
    }
    void add(int a,int b,int c){
        x[++cnt]=a;
        y[cnt]=b;
        z[cnt]=c;
        nxt[cnt]=fst[a];
        fst[a]=cnt;
    }
}e,q;
int T,n,m,from,to,dist,in[N],rt,dis[N],fa[N],ans[N];
bool vis[N];
void dfs(int rt){
    for (int i=e.fst[rt];i;i=e.nxt[i]){
        dis[e.y[i]]=dis[rt]+e.z[i];
        dfs(e.y[i]);
    }
}
int getf(int k){
    return fa[k]==k?k:fa[k]=getf(fa[k]);
}
void LCA(int rt){
    for (int i=e.fst[rt];i;i=e.nxt[i]){
        LCA(e.y[i]);
        fa[getf(e.y[i])]=rt;
    }
    vis[rt]=1;
    for (int i=q.fst[rt];i;i=q.nxt[i])
        if (vis[q.y[i]]&&!ans[q.z[i]])
            ans[q.z[i]]=dis[q.y[i]]+dis[rt]-2*dis[getf(q.y[i])];
}
int main(){
    scanf("%d",&T);
    while (T--){
        q.set(),e.set();
        memset(in,0,sizeof in);
        memset(vis,0,sizeof vis);
        memset(ans,0,sizeof ans);
        scanf("%d%d",&n,&m);
        for (int i=1;i<n;i++)
            scanf("%d%d%d",&from,&to,&dist),e.add(from,to,dist),in[to]++;
        for (int i=1;i<=m;i++)
            scanf("%d%d",&from,&to),q.add(from,to,i),q.add(to,from,i);
        rt=0;
        for (int i=1;i<=n&&rt==0;i++)
            if (in[i]==0)
                rt=i;
        dis[rt]=0;
        dfs(rt);
        for (int i=1;i<=n;i++)
            fa[i]=i;
        LCA(rt);
        for (int i=1;i<=m;i++)
            printf("%d\n",ans[i]);
    }
    return 0;
}






/*  倍增 O((n+q)logn) 离线*/
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;
const int N=10000+5;
vector <int> son[N];
int T,n,depth[N],fa[N][20],in[N],a,b;
void dfs(int prev,int rt){
    depth[rt]=depth[prev]+1;
    fa[rt][0]=prev;
    for (int i=1;i<20;i++)
        fa[rt][i]=fa[fa[rt][i-1]][i-1];
    for (int i=0;i<son[rt].size();i++)
        dfs(rt,son[rt][i]);
}
int LCA(int x,int y){
    if (depth[x]<depth[y])
        swap(x,y);
    for (int i=19;i>=0;i--)
        if (depth[x]-(1<<i)>=depth[y])
            x=fa[x][i];
    if (x==y)
        return x;
    for (int i=19;i>=0;i--)
        if (fa[x][i]!=fa[y][i])
            x=fa[x][i],y=fa[y][i];
    return fa[x][0];
}
int main(){
    scanf("%d",&T);
    while (T--){
        scanf("%d",&n);
        for (int i=1;i<=n;i++)
            son[i].clear();
        memset(in,0,sizeof in);
        for (int i=1;i<n;i++){
            scanf("%d%d",&a,&b);
            son[a].push_back(b);
            in[b]++;
        }
        depth[0]=-1;
        int rt=0;
        for (int i=1;i<=n&&rt==0;i++)
            if (in[i]==0)
                rt=i;
        dfs(0,rt);
        scanf("%d%d",&a,&b);
        printf("%d\n",LCA(a,b));
    }
    return 0;
}




/* RMQ 构造O(nlogn) 查询O(1) */
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;
const int N=10000+5;
vector <int> son[N];
int T,n,depth[N],fa[N][20],in[N],a,b;
void dfs(int prev,int rt){
    depth[rt]=depth[prev]+1;
    fa[rt][0]=prev;
    for (int i=1;i<20;i++)
        fa[rt][i]=fa[fa[rt][i-1]][i-1];
    for (int i=0;i<son[rt].size();i++)
        dfs(rt,son[rt][i]);
}
int LCA(int x,int y){
    if (depth[x]<depth[y])
        swap(x,y);
    for (int i=19;i>=0;i--)
        if (depth[x]-(1<<i)>=depth[y])
            x=fa[x][i];
    if (x==y)
        return x;
    for (int i=19;i>=0;i--)
        if (fa[x][i]!=fa[y][i])
            x=fa[x][i],y=fa[y][i];
    return fa[x][0];
}
int main(){
    scanf("%d",&T);
    while (T--){
        scanf("%d",&n);
        for (int i=1;i<=n;i++)
            son[i].clear();
        memset(in,0,sizeof in);
        for (int i=1;i<n;i++){
            scanf("%d%d",&a,&b);
            son[a].push_back(b);
            in[b]++;
        }
        depth[0]=-1;
        int rt=0;
        for (int i=1;i<=n&&rt==0;i++)
            if (in[i]==0)
                rt=i;
        dfs(0,rt);
        scanf("%d%d",&a,&b);
        printf("%d\n",LCA(a,b));
    }
    return 0;
}