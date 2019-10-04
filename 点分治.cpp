/******************  树上距离为k的点对是否存在  *******************/
#include<iostream>
#include<vector>
#include<algorithm>
#include<queue>
#include<cstring>
#include<cstdio>
using namespace std;

int read()
{
    int f=1,x=0;
    char ss=getchar();
    while(ss<'0'||ss>'9'){if(ss=='-')f=-1;ss=getchar();}
    while(ss>='0'&&ss<='9'){x=x*10+ss-'0';ss=getchar();}
    return f*x;
}

const int inf=10000000;
const int maxn=100010;
int n,m;
struct node{int v,dis,nxt;}E[maxn<<1];
int tot,head[maxn];
int maxp[maxn],size[maxn],dis[maxn],rem[maxn];
int vis[maxn],test[inf],judge[inf],q[maxn];
int query[1010];
int sum,rt;
int ans;

void add(int u,int v,int dis)
{
    E[++tot].nxt=head[u];
    E[tot].v=v;
    E[tot].dis=dis;
    head[u]=tot;
}

void getrt(int u,int pa)
{
    size[u]=1; maxp[u]=0;
    for(int i=head[u];i;i=E[i].nxt) 
    {
        int v=E[i].v;
        if(v==pa||vis[v]) continue;
        getrt(v,u);
        size[u]+=size[v];
        maxp[u]=max(maxp[u],size[v]);
    }
    maxp[u]=max(maxp[u],sum-size[u]);
    if(maxp[u]<maxp[rt]) rt=u;
}

void getdis(int u,int fa)
{
    rem[++rem[0]]=dis[u];
    for(int i=head[u];i;i=E[i].nxt)
    {
        int v=E[i].v;
        if(v==fa||vis[v])continue;
        dis[v]=dis[u]+E[i].dis;
        getdis(v,u);
    }
}

void calc(int u)
{
    int p=0;
    for(int i=head[u];i;i=E[i].nxt)
    {
        int v=E[i].v;
        if(vis[v])continue;
        rem[0]=0; dis[v]=E[i].dis;
        getdis(v,u);//处理u的每个子树的dis

        for(int j=rem[0];j;--j)//遍历当前子树的dis
        for(int k=1;k<=m;++k)//遍历每个询问
        if(query[k]>=rem[j])
        test[k]|=judge[query[k]-rem[j]];
        //如果query[k]-rem[j]d的路径存在就标记第k个询问

        for(int j=rem[0];j;--j)//保存出现过的dis于judge
        q[++p]=rem[j],judge[rem[j]]=1;
    }
    for(int i=1;i<=p;++i)//处理完这个子树就清空judge
    judge[q[i]]=0;//特别注意一定不要用memset，会T

}

void solve(int u)
{   
    //judge[i]表示到根距离为i的路径是否存在
    vis[u]=judge[0]=1; calc(u);//处理以u为根的子树
    for(int i=head[u];i;i=E[i].nxt)//对每个子树进行分治
    {
        int v=E[i].v;
        if(vis[v])continue;
        sum=size[v]; maxp[rt=0]=inf;
        getrt(v,0); solve(rt);//在子树中找重心并递归处理
    }
}

int main()
{
    n=read();m=read();
    for(int i=1;i<n;++i)
    {
        int u=read(),v=read(),dis=read();
        add(u,v,dis);add(v,u,dis);
    }
    for(int i=1;i<=m;++i)
    query[i]=read();//先记录每个询问以离线处理

    maxp[rt]=sum=n;//第一次先找整棵树的重心
    getrt(1,0); 
    solve(rt);//对树进行点分治

    for(int i=1;i<=m;++i)
    {
        if(test[i]) printf("AYE\n");
        else printf("NAY\n");
    }
    return 0;
}

/***************************************************/

/*****************树上小于等于k的点对个数********************/
#include<bits/stdc++.h>
#define ll long long
using namespace std;
inline int read(){
    int x=0;char ch=' ';int f=1;
    while(ch!='-'&&(ch<'0'||ch>'9'))ch=getchar();
    if(ch=='-')f=-1,ch=getchar();
    while(ch>='0'&&ch<='9')x=x*10+ch-'0',ch=getchar();
    return x*f;
}
struct edge{
    int to,next,w;
}e[80001];
int n,tot,root;
ll k;
int head[40001];
inline void addedge(int x,int y,int l){
    e[++tot].to=y;e[tot].next=head[x];e[tot].w=l;head[x]=tot;
}
int size[40001],vis[40001],mx,sz;
ll dis[40001],q[40001],l,r;
void getroot(int x,int fa){
    size[x]=1;int num=0;
    for(int i=head[x];i;i=e[i].next){
        int u=e[i].to;
        if(u==fa||vis[u])continue;
        getroot(u,x);
        size[x]+=size[u];
        num=max(num,size[u]);
    }
    num=max(num,sz-size[x]);
    if(num<mx){
        mx=num;root=x;
    }
}
void getdis(int x,int fa){
    q[++r]=dis[x];
    for(int i=head[x];i;i=e[i].next){
        int u=e[i].to;
        if(u==fa||vis[u])continue;
        dis[u]=dis[x]+e[i].w;
        getdis(u,x);
    }
}
ll calc(int x,int v){
    r=0;
    dis[x]=v;
    getdis(x,0);
    ll sum=0;
    l=1;
    sort(q+1,q+r+1);
    while(l<r){
        if(q[l]+q[r]<=k)sum+=r-l,l++;
        else r--;
    }
    return sum;
}
ll ans;
void dfs(int x){
    ans+=calc(x,0);
    vis[x]=1;
    for(int i=head[x];i;i=e[i].next){
        int u=e[i].to;
        if(vis[u])continue;
        ans-=calc(u,e[i].w);
        sz=size[u];
        mx=0x3f3f3f3f;
        getroot(u,0);
        dfs(root);
    }
}
int main(){
    n=read();
    for(int i=1;i<n;i++){
        int x=read(),y=read(),l=read();
        addedge(x,y,l);addedge(y,x,l);
    }
    k=read();
    sz=n;
    mx=0x3f3f3f3f;
    getroot(1,0);
    dfs(root);
    printf("%lld",ans);
    return 0;
}

/**********************************************************/


#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<cstring>
#include<cmath>
#include<algorithm>
#define maxn 10002
#define INF 2147483646
#define ll long long
using namespace std;

inline int gi()
{
  int date = 0,m = 1; char ch = 0;
  while(ch!='-'&&(ch<'0'||ch>'9'))ch = getchar();
  if(ch=='-'){m = -1; ch = getchar();}
  while(ch>='0' && ch<='9')
    {
      date = date*10+ch-'0';
      ch = getchar();
    }return date*m;
}
inline void write(ll qw)
{
  if(qw<0){putchar('-'); qw = -qw;}
  if(qw>9)write(qw/10);
  putchar(qw%10+'0');
}

struct node{int to,next,len;}t[2*maxn+5];
int head[maxn+5];
int sim[maxn+5],mxson[maxn+5];
bool vis[maxn+5];
int MX,root,dis[maxn+5],summar;
ll ans;
int n,k,cnt,S;

void addedge(int u,int v,int l)
{
  cnt ++;
  t[cnt].to = v;
  t[cnt].next = head[u];
  t[cnt].len = l;
  head[u] = cnt;
  return;
}

void getroot(int u,int fa)
{
  sim[u] = 1; mxson[u] = 0;
  for(int i = head[u];i;i = t[i].next)
    {
      int v = t[i].to;
      if(v == fa || vis[v])continue;
      getroot(v,u);
      sim[u] = sim[u] + sim[v];
      mxson[u] = max(sim[v],mxson[u]);
    }
  mxson[u] = max(mxson[u],S-sim[u]);
  if(mxson[u]<MX){root = u;MX = mxson[u];}
}

void getdis(int u,int fa,int dist)
{
  dis[++summar] = dist;
  for(int i = head[u];i;i=t[i].next)
    {
      int v = t[i].to;
      if(vis[v]||v == fa)continue;
      getdis(v,u,dist+t[i].len);
    }
  return;
}

int consolate(int sta,int len1)
{
  summar = 0;
  memset(dis,0,sizeof(dis));
  getdis(sta,0,len1);
  sort(dis+1,dis+summar+1);
  int L = 1,R = summar,tep=0;
  while(L<=R)
    {
      if(dis[R] + dis[L]<=k){tep = tep + R - L; L++;}
      else R--;
    }
  return tep;
}

void Divide(int tr)
{
  ans = ans + consolate(tr,0);
  vis[tr] = true;
  for(int i = head[tr];i;i = t[i].next)
    {
      int v = t[i].to;
      if(vis[v])continue;
      ans = ans - consolate(v,t[i].len);
      S = sim[v]; root = 0;
      MX = INF; getroot(v,0);
      Divide(root);
    }
  return;
}

int main()
{
  freopen("Tree.in","r",stdin);
  int u,v,l;
  while(1)
    {
      n = gi(); k = gi();
      if(n == 0 && k == 0)break;
      cnt = 0; memset(head,0,sizeof(head));
      for(int i = 1; i <= n - 1; i ++)
    {
      u = gi(); v = gi(); l = gi();
      addedge(u,v,l); addedge(v,u,l);
    }
      ans = 0;
      memset(vis,0,sizeof(vis));
      MX = INF ; S = n; getroot(1,0);
      Divide(root);
      write(ans);putchar('\n');
    }
  return 0;
}

