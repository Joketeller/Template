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


/*********************距离为k点对个数****************************/


