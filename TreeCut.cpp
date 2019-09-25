#include <bits/stdc++.h>
using namespace std;
const int MOD=1e9+7;
/******************************树链剖分******************************/
struct treenode
{
    int val;
    int father, treesize, index, depth, topfather, maxson;
};

struct segmenttreenode
{
    int l, r,size;
    int val, lazy;
};

struct Edge
{
    int u, v;
};

class TreeCut
{
#define maxn 100005
#define maxm 1000005
#define lson rt<<1
#define rson rt<<1|1
private:
    int TreeSize,EdgeSize,cnt;
    treenode t[maxn];
    int head[maxn],nxt[maxm];
    int to[maxm];
    int val[maxn];  //中间值
    segmenttreenode segt[maxn << 2];

public:
    void init(int _size=0)
    {
        memset(head, 0, sizeof(head));
        TreeSize = EdgeSize = cnt  = 0;
    }

    void addedge(int _u,int _v)
    {
        ++EdgeSize;
        to[EdgeSize] = _v;
        nxt[EdgeSize] = head[_u];
        head[_u] = EdgeSize;
        ++EdgeSize;
        to[EdgeSize] = _u;
        nxt[EdgeSize] = head[_v];
        head[_v] = EdgeSize;
    }

    int relationshipDFS(int now, int father, int depth)
    {
        t[now].index = 0;
        t[now].maxson = 0;
        t[now].depth = depth;
        t[now].father = father;
        t[now].treesize = 1;
        int maxsonsize = -1;
        for (int i = head[now]; i;i=nxt[i])
        {
            if (to[i]==father)
                continue;
            t[now].treesize += relationshipDFS(to[i], now, depth + 1);
            if (t[to[i]].treesize>maxsonsize)
            {
                t[now].maxson = to[i];
                maxsonsize = t[to[i]].treesize;
            }
        }
        return t[now].treesize;
    }

    void ReIndexDFS(int now,int topfather) //topfather该链最高的儿子链顶
    {
        t[now].index = ++cnt;
        val[cnt] = t[now].val; 
        t[now].topfather = topfather;
        if (!t[now].maxson)
            return;
        ReIndexDFS(t[now].maxson, topfather);
        for (int i = head[now]; i;i=nxt[i])
        {
            if (!t[to[i]].index)
                ReIndexDFS(to[i], to[i]);
        }
    }

    void UpdateNode(int rt)
    {
        segt[rt].val=(segt[lson].val+segt[rson].val)%MOD;
    }

    void SegmentTreeBuild(int now,int _l,int _r)
    {
        segt[now].l = _l;
        segt[now].r = _r;
        segt[now].size=_r-_l+1;
        if (_l == _r)
        {
            segt[now].val = val[_l];
            segt[now].lazy = 0;
            return;
        }
        int mid = (_l + _r) >> 1;
        SegmentTreeBuild(now << 1, _l, mid);
        SegmentTreeBuild(now << 1 | 1, mid + 1, _r);
        UpdateNode(now);
    }

    void PushDown(int rt)   //下传标记
    {
        if (!segt[rt].lazy)
            return;
        segt[lson].val=(segt[lson].val+ segt[rt].val*segt[lson].size)%MOD;
        segt[rson].val=(segt[rson].val+ segt[rt].val*segt[rson].size)%MOD;
        segt[lson].lazy=(segt[lson].lazy+segt[rt].lazy)%MOD;
        segt[rson].lazy=(segt[rson].lazy+segt[rt].lazy)%MOD;
        segt[rt].lazy=0;
    }

    void IntervalAdd(int rt,int l,int r,int v)
    {
        if (l<=segt[rt].l && segt[rt].r <= r)
        {
            segt[rt].val+=v*segt[rt].size;
            segt[rt].lazy+=v;
            return;
        }
        PushDown(rt);
        int mid=(segt[rt].l+segt[rt].r)>>1;
        if (l<=mid)
        {
            IntervalAdd(lson,l,r,v);
        }
        if (r>mid)
        {
            IntervalAdd(rson,l,r,v);
        }
        UpdateNode(rt);
    }

    int IntervalSum(int rt,int l,int r)
    {
        int ans=0;
        if (l<= segt[rt].l && segt[rt].r<=r)
            return segt[rt].val;
        PushDown(rt);
        int mid=(segt[rt].l+segt[rt].r)>>1;
        if (l<=mid)
            ans=(ans+IntervalSum(lson,l,r))%MOD;
        if (r>mid)
            ans=(ans+IntervalSum(rson,l,r))%MOD;
        return ans;
    }

    int TreeSum(int x,int y)   //x,y路径求和
    {
        int ans=0;
        while (t[x].topfather!=t[y].topfather)
        {
            if (t[t[x].topfather].depth < t[t[y].topfather].depth)
                swap(x,y);
            ans=(ans+IntervalSum(1,t[t[x].topfather].index,t[x].index))%MOD;
            x=t[t[x].topfather].father;
        }
        if (t[x].depth>t[y].depth)
            swap(x,y);
        ans=(ans+IntervalSum(1,t[x].index,t[y].index))%MOD;
        return ans;
    }

    void TreeAdd(int x,int y,int v)
    {
         while (t[x].topfather!=t[y].topfather)
        {
            if (t[t[x].topfather].depth < t[t[y].topfather].depth)
                swap(x,y);
            IntervalAdd(1,t[t[x].topfather].index,t[x].index,v);
            x=t[t[x].topfather].father;
        }
        if (t[x].depth>t[y].depth)
            swap(x,y);
         IntervalAdd(1,t[x].index,t[y].index,v);
    }

    
#undef maxn
#undef maxm
#undef lson
#undef rson
};