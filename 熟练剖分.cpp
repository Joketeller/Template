#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define scan(x) scanf("%d",&x)

struct node 
{
    int idx, fa, val,depth,son,siz;
};

struct tree
{
    #define maxn 50000
    #define maxm 2500000
    node t[maxn];
    int head[maxn], v[maxn], nxt[maxn],cnt;

    void init()
    {
        cnt = 0;
        memset(head, 0, sizeof(head));
    }

    void addedge(int x,int y)
    {
        v[++cnt] = x;
        nxt[cnt] = head[y];
        head[y] = cnt;
        v[++cnt] = y;
        nxt[cnt] = head[x];
        head[x] = cnt;
    }

    void SizeDFS(int root, int pre)
    {
        
    }
};

int main()
{
#ifdef _IRONHEAD_
    assert(freopen("/Users/ironhead/algorithm/in.in", "r", stdin));
    assert(freopen("/Users/ironhead/algorithm/out.out", "w", stdout));
#endif
    return 0;
}