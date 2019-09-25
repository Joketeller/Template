#include <bits/stdc++.h>
using namespace std;

typedef long long LL;

const LL maxn = 110000;
const LL INF = 1e18;
const LL dimension = 2;
LL D;

struct Node
{
    LL d[dimension], maxpos[dimension], minpos[dimension], v, sum, lazy, cnt;
    //以中心点d作为空间的代表，max和min分别是空间的边界
    LL l, r;
    bool operator<(const Node &b) const
    {
        return d[D] < b.d[D];
    }
    bool operator==(const Node &b) const
    {
        bool ans = 1;
        for (int i = 0; i < dimension; i++)
        {
            ans &= d[i] == b.d[i];
        }
        return ans;
    }
} p[maxn];

bool in(int x1, int y1, int x2, int y2, int X1, int Y1, int X2, int Y2)
{
    return x1 <= X1 && X2 <= x2 && y1 <= Y1 && Y2 <= y2;
}

bool out(int x1, int y1, int x2, int y2, int X1, int Y1, int X2, int Y2)
{
    return x1 > X2 || x2 < X1 || y1 > Y2 || y2 < Y1;
}

struct KDT
{
    LL root, cnt, block;
    Node tr[maxn], now;
    //now,用来单点插入
    void pushup(int rt)
    {
        tr[rt].cnt = 1;
        int l = tr[rt].l, r = tr[rt].r;
        if (l)
            tr[rt].cnt += tr[l].cnt;
        if (r)
            tr[rt].cnt += tr[r].cnt;
        for (int i = 0; i < dimension; i++)
        {
            tr[rt].maxpos[i] = tr[rt].minpos[i] = tr[rt].d[i];
            if (l)
            {
                tr[rt].minpos[i] = min(tr[rt].minpos[i], tr[l].minpos[i]);
                tr[rt].maxpos[i] = max(tr[rt].maxpos[i], tr[l].maxpos[i]);
            }
            if (r)
            {
                tr[rt].minpos[i] = min(tr[rt].minpos[i], tr[r].minpos[i]);
                tr[rt].maxpos[i] = max(tr[rt].maxpos[i], tr[r].maxpos[i]);
            }
        }
        tr[rt].sum = tr[l].sum + tr[r].sum + tr[rt].v;
    }

    void pushdown(int rt)
    {
        if (tr[rt].lazy)
        {
            int l = tr[rt].l, r = tr[rt].r;
            if (l)
            {
                tr[l].lazy += tr[rt].lazy;
                tr[l].v += tr[rt].lazy;
                tr[l].sum += tr[l].cnt * tr[rt].lazy;
            }
            if (r)
            {
                tr[r].lazy += tr[rt].lazy;
                tr[r].v += tr[rt].lazy;
                tr[r].sum += tr[r].cnt * tr[rt].lazy;
            }
            tr[rt].lazy = 0;
        }
    }

    LL rebuild(LL l, LL r, LL dep)
    { //重构
        if (l > r)
            return 0;
        D = dep;
        LL mid = (l + r) >> 1;
        nth_element(p + l, p + mid, p + r + 1);
        tr[mid] = p[mid];
        tr[mid].l = rebuild(l, mid - 1, (dep + 1) % dimension);
        tr[mid].r = rebuild(mid + 1, r, (dep + 1) % dimension);
        pushup(mid);
        return mid;
    }

    void checkSize()
    {
        if (cnt == block)
        {
            for (int i = 1; i <= cnt; i++)
                p[i] = tr[i];
            root = rebuild(1, cnt, 0);
            block += 10000;
        }
    }

    void ins(LL &rt, bool D)
    { //单点插入，如果没有就新开点
        if (!rt)
        {
            rt = ++cnt;
            for (int i = 0; i < dimension; i++)
                tr[rt].d[i] = tr[rt].maxpos[i] = tr[rt].minpos[i] = now.d[i];
            tr[rt].v = tr[rt].sum = now.v;
            return;
        }
        if (now == tr[rt])
        {
            tr[rt].v += now.v, tr[rt].sum += now.v;
            return;
        }
        pushdown(rt);
        if (now.d[D] < tr[rt].d[D])
            ins(tr[rt].l, D ^ 1);
        else
            ins(tr[rt].r, D ^ 1);
        pushup(rt);
    }
    
    LL query(int rt, int x1, int y1, int x2, int y2)
    {
        if (!rt)
            return 0;
        LL res = 0;
        if (out(x1, y1, x2, y2, tr[rt].minpos[0], tr[rt].minpos[1], tr[rt].maxpos[0], tr[rt].maxpos[1]))
            return 0;
        if (in(x1, y1, x2, y2, tr[rt].minpos[0], tr[rt].minpos[1], tr[rt].maxpos[0], tr[rt].maxpos[1]))
            return tr[rt].sum;
        pushdown(rt);
        if (in(x1, y1, x2, y2, tr[rt].d[0], tr[rt].d[1], tr[rt].d[0], tr[rt].d[1]))
            res += tr[rt].v;
        res += query(tr[rt].l, x1, y1, x2, y2) + query(tr[rt].r, x1, y1, x2, y2);
        pushup(rt);
        return res;
    }

    void update(int rt, int x1, int y1, int x2, int y2, LL add)
    {
        if (!rt)
            return;
        if (out(x1, y1, x2, y2, tr[rt].minpos[0], tr[rt].minpos[1], tr[rt].maxpos[0], tr[rt].maxpos[1]))
            return;
        if (in(x1, y1, x2, y2, tr[rt].minpos[0], tr[rt].minpos[1], tr[rt].maxpos[0], tr[rt].maxpos[1]))
        {
            tr[rt].lazy += add;
            tr[rt].sum += add * tr[rt].cnt;
            tr[rt].v += add;
            return;
        }
        pushdown(rt);
        if (in(x1, y1, x2, y2, tr[rt].d[0], tr[rt].d[1], tr[rt].d[0], tr[rt].d[1]))
        {
            tr[rt].v += add;
        }
        update(tr[rt].l, x1, y1, x2, y2, add);
        update(tr[rt].r, x1, y1, x2, y2, add);
        pushup(rt);
    }

    void init()
    {
        root = cnt = 0;
        block = 10000;
    }

    LL getdis(LL val[dimension], LL rt)
    { //估价函数，用来寻找区间
        LL res = 0;
        for (LL i = 0; i < dimension; i++)
        {
            if (val[i] < tr[rt].minpos[i])
                res += (tr[rt].minpos[i] - val[i]) * (tr[rt].minpos[i] - val[i]);
            if (val[i] > tr[rt].maxpos[i])
                res += (val[i] - tr[rt].maxpos[i]) * (val[i] - tr[rt].maxpos[i]);
        }
        return res;
    }

    LL ans;

    void ask(LL val[dimension], LL rt)
    { //询问最近点 欧式距离的平方
        LL dis = 0;
        for (LL i = 0; i < dimension; i++)
            dis += ((tr[rt].d[i] - val[i]) * (tr[rt].d[i] - val[i]));
        if (dis == 0)
            dis = INF;
        if (dis < ans)
            ans = dis;
        LL dl = tr[rt].l ? getdis(val, tr[rt].l) : INF;
        LL dr = tr[rt].r ? getdis(val, tr[rt].r) : INF;
        if (dl < dr)
        {
            if (dl < ans)
                ask(val, tr[rt].l);
            if (dr < ans)
                ask(val, tr[rt].r);
        }
        else
        {
            if (dr < ans)
                ask(val, tr[rt].r);
            if (dl < ans)
                ask(val, tr[rt].l);
        }
    }
} Tree;

int n, m;
int id;
int L[maxn], R[maxn];
int dep[maxn];
vector<int> G[maxn];

void dfs(int u, int fa)
{
    dep[u] = dep[fa] + 1;
    L[u] = ++id;
    for (auto v : G[u])
    {
        if (v == fa)
            continue;
        dfs(v, u);
    }
    R[u] = id;
}

int main()
{
    while (~scanf("%d%d", &n, &m))
    {
        Tree.init();
        id = 0;
        for (int i = 0; i < n - 1; i++)
        {
            int u, v;
            scanf("%d%d", &u, &v);
            G[u].push_back(v);
            G[v].push_back(u);
        }
        dep[0] = -1;
        dfs(1, 0);
        for (int i = 1; i <= n; i++)
        {
            Node &now = p[i];
            now.d[0] = L[i];
            now.d[1] = dep[i];
        }
        Tree.root = Tree.rebuild(1, n, 0);
        while (m--)
        {
            int tp, l, x;
            scanf("%d", &tp);
            if (tp == 1)
            {
                scanf("%d%d", &l, &x);
                Tree.update(Tree.root, 1, l, n, l, x);
            }
            else
            {
                scanf("%d", &x);
                printf("%lld\n", Tree.query(Tree.root, L[x], 1, R[x], n));
            }
        }
    }
    return 0;
}

/******************mine*************************/
#include <bits/stdc++.h>
using namespace std;
#define ll int
#define dimension 2
const int maxn = 5e5 + 10;
#define INF 1000000000
int D;
// int root;

struct Node
{
    ll d[dimension], l, r, mn[dimension], mx[dimension];
   
    ll& operator [](int n)
    {
        return d[n];
    }
     friend bool operator < (Node a, Node b) 
    {
        return a[D] < b[D];
    }
    bool operator==(const Node &b) const
    {
        bool tmp = 1;
        for (int i = 0; i < dimension; i++)
            tmp &= (d[i] == b.d[i]);
        return tmp;
    }
    friend ll dis(Node x, Node y) //曼哈顿距离
    {
        ll ret=0;
        for (int i=0;i<dimension;++i)
        {
            ret+=abs(x[i]-y[i]);
        }
        return ret;
    }
} p[maxn];

// namespace KDTree
// {
struct Tree
{
    int root;
    Node tr[maxn],T;
    ll ans;
    void pushup(int rt)
    {
        for (int i = 0; i < dimension; ++i)
            tr[rt].mn[i] = tr[rt].mx[i] = tr[rt].d[i];
        if (tr[rt].l)
        {
            for (int i = 0; i < dimension; i++)
            {
                tr[rt].mn[i]=min(tr[rt].mn[i], tr[tr[rt].l].mn[i]);
                tr[rt].mx[i]=max(tr[rt].mx[i], tr[tr[rt].l].mx[i]);
            }
        }
        if (tr[rt].r)
        {
            for (int i = 0; i < dimension; i++)
            {
                tr[rt].mn[i]=min(tr[rt].mn[i], tr[tr[rt].r].mn[i]);
                tr[rt].mx[i]=max(tr[rt].mx[i], tr[tr[rt].r].mx[i]);
            }
        }
    }

    int build(int l, int r, int dim)
    {
        int mid = (l + r) >> 1;
        D=dim;
        nth_element(p + l, p + mid, p + r + 1);
        tr[mid] = p[mid];
        if (l < mid)
            tr[mid].l = build(l, mid - 1, dim^1);
        if (mid < r)
            tr[mid].r = build(mid + 1, r, dim^1 );
        pushup(mid);
        return mid;
    }

    ll getmn(Node x) //估值
    {
        ll ans=0;
        for (int i=0;i<dimension;i++)
        {
            ans+=max(T[i] -x.mx[i],0);
            ans+=max(x.mn[i] - T[i],0);
        }
        return ans;
    }

    ll getmx(Node x) //估值 
    {
        ll ans = 0;
        for (int i = 0; i < dimension;i++)
        {
            ans += max(abs(T[i] - x.mx[i]), abs(T[i] - x.mn[i]));
        }
        return ans;
    }

    inline void querymaxdistance(int k)
    {
        ans=max(ans,dis(T,tr[k]));
        ll dl=-INF,dr=-INF;
        int l=tr[k].l,r=tr[k].r;
        if (l)
            dl=getmx(tr[l]);
        if (r)
            dr=getmx(tr[r]);
        if (dl>dr)
        {
            if (dl>ans)
                querymaxdistance(l);
            if (dr>ans) 
                querymaxdistance(r);
        }
        else
        {
            if (dr>ans)
                querymaxdistance(r);
            if (dl>ans)
                querymaxdistance(l);
        }
    }

    inline void querymindistance(int k)
    {
        if (dis(T,tr[k]))
            ans=min(ans,dis(T,tr[k]));
        ll dl=INF,dr=INF;
        int l=tr[k].l,r=tr[k].r;
        if (tr[k].l)
            dl=getmn(tr[l]);
        if (tr[k].r)
            dr=getmn(tr[r]);
        if (dl<dr)
        {
            if (dl<ans)
                querymindistance(l);
            if (dr<ans) 
                querymindistance(r);
        }
        else
        {
            if (dr<ans)
                querymindistance(r);
            if (dl<ans)
                querymindistance(l);
        }
    }

    ll query(int op,ll x,ll y)
    {
        T[0]=x,T[1]=y;
        if (op==0)
            ans=-INF,querymaxdistance(root);
        else
            ans=INF,querymindistance(root);
        return ans;
    }

}kd;
ll x[maxn],y[maxn];

int n;
int main()
{
    #ifdef _IRONHEAD_
        assert(freopen("/Users/ironhead/algorithm/in.in", "r", stdin));
       // assert(freopen("/Users/ironhead/algorithm/out.out", "w", stdout));
    #endif
    scanf("%d",&n);
    for (int i=1;i<=n;i++)
    {
        scanf("%d%d",&x[i],&y[i]);
        p[i][0]=x[i],p[i][1]=y[i];
    }
    kd.root=kd.build(1,n,0);
    ll ans=INF;
    int cnt=0;
    for (int i=1;i<=n;i++)
    {
        ll mx=kd.query(0,x[i],y[i]);
        ll mn=kd.query(1,x[i],y[i]);
        ans=min(ans,mx-mn);
    }
    printf("%d\n",ans);
    return 0;
}

/*可能有点问题的板子，支持插入和查询曼哈顿距离*/
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#define N 100005
#define inf (1<<30)
#define scan(x) scanf("%d",&x)
using namespace std;
const int k=3;
int n,m,dim,rt,ans,num=0;
struct node{int p[k],minn[k],maxx[k];}a[N];
bool cmp(node x,node y){ return x.p[dim]<y.p[dim]; }
struct kd_tree{
    int c[N][2];
    node s[N],q;
    void update(int o)
    {//管辖范围
        int l=c[o][0],r=c[o][1];
        for(int i=0;i<k;i++){
            if(l){ s[o].minn[i]=min(s[o].minn[i],s[l].minn[i]); s[o].maxx[i]=max(s[o].maxx[i],s[l].maxx[i]); }
            if(r){ s[o].minn[i]=min(s[o].minn[i],s[r].minn[i]); s[o].maxx[i]=max(s[o].maxx[i],s[r].maxx[i]); }
        }
    }
    void add(int o,node t){ for(int i=0;i<k;i++)s[o].minn[i]=s[o].maxx[i]=s[o].p[i]=t.p[i]; }
    int dist(node t,int o)
    {
        int tmp=0;
        for(int i=0;i<k;i++) tmp+=max(0,s[o].minn[i]-t.p[i]);
        for(int i=0;i<k;i++) tmp+=max(0,t.p[i]-s[o].maxx[i]);
        return tmp;
    }//?
    void build(int &o,int l,int r,int now)
    {
        o=(l+r)>>1; dim=now%k;
        nth_element(a+l,a+o,a+r+1,cmp);
        add(o,a[o]);
        if(l<o) build(c[o][0],l,o-1,now+1);
        if(o<r) build(c[o][1],o+1,r,now+1);
        update(o);
    }

    void ins(int o,int now){
        now%=k;
        if(q.p[now]<s[o].p[now]){
            if(c[o][0]) ins(c[o][0],now+1);
            else c[o][0]=++n,add(n,q);
        }
        else{
            if(c[o][1]) ins(c[o][1],now+1);
            else c[o][1]=++n,add(n,q);
        }
        update(o);
    }
    void qry(int o){//曼哈顿距离,且只求最短，dis是最短距离
        int tmp=0;
        for(int i=0;i<k;i++) tmp+=abs(s[o].p[i]-q.p[i]);
        ans=min(ans,tmp);
        int dl=c[o][0]?dist(q,c[o][0]):inf,dr=c[o][1]?dist(q,c[o][1]):inf;
        if(dl<dr)
        {
            if(dl<ans) qry(c[o][0]);
            if(dr<ans) qry(c[o][1]);
        }else{
            if(dr<ans) qry(c[o][1]);
            if(dl<ans) qry(c[o][0]);
        }
    }
}kd;
int tmp;
int main()
{
    // k=3;
    scan(tmp),scan(tmp),scan(tmp),scan(m);
    scan(tmp),scan(a[1].p[0]),scan(a[1].p[1]),scan(a[1].p[2]);
    ++num;
    kd.build(rt,1,1,0);
    --m;
    while(m--){
        scan(tmp);
        scanf("%d%d%d",&kd.q.p[0],&kd.q.p[1],&kd.q.p[2]);
        if(tmp==1)
        { 
            a[++num]=kd.q;
            if (num%5000==0)
                kd.build(rt,1,num,0);
            else
            {
                kd.ins(rt,0);
            }
            
        }
        else{
            ans=inf; kd.qry(rt); printf("%d\n",ans);
        }
    }
    return 0;
}