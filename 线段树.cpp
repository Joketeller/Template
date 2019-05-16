
#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <queue>
#include <cstdio>
#include <stack>
#include <map>
#include <string>
#include <set>
#include <unordered_set>
#include <iomanip>
#include <bitset>
#pragma GCC optimize(2)
#define eps 1e-5
#define mod 1000000007
#define pi acos(-1)
#define MAXN 100005
#define ee 2.71828182845904523536
using namespace std;
int a[MAXN][30];
int tmp;
template<size_t N, typename T>
struct Segment_tree
{
    struct Node
    {
        int l[2],r[2];
        long long s[2];
        int d,ls,rs;
    } seg[N * 2];
    int tot=1;
    void maintain(int x)
    {
        int dl=seg[seg[x].ls].d,dr=seg[seg[x].rs].d;
        seg[x].d=seg[seg[x].ls].d^seg[seg[x].rs].d;
        seg[x].l[0]=seg[seg[x].ls].l[0]+seg[seg[x].rs].l[(dl==0) ? 0:1];
        seg[x].l[1]=seg[seg[x].ls].l[1]+seg[seg[x].rs].l[(dl==0) ? 1:0];
        seg[x].r[0]=seg[seg[x].rs].r[0]+seg[seg[x].ls].r[(dr==0) ? 0:1];
        seg[x].r[1]=seg[seg[x].rs].r[1]+seg[seg[x].ls].r[(dr==0) ? 1:0];
        seg[x].s[0]=seg[seg[x].ls].s[0]+seg[seg[x].rs].s[0]+seg[seg[x].ls].r[0]*seg[seg[x].rs].l[0]+seg[seg[x].ls].r[1]*seg[seg[x].rs].l[1];
        seg[x].s[1]=seg[seg[x].ls].s[1]+seg[seg[x].rs].s[1]+seg[seg[x].ls].r[0]*seg[seg[x].rs].l[1]+seg[seg[x].ls].r[1]*seg[seg[x].rs].l[0];
    }
    void init()
    {
        memset(seg,0,sizeof(seg));
        tot=1;
    }
    void build(int left, int right, int x, int idx)
    {
        if (left == right)
        {
            seg[x].d=a[left][idx];
            seg[x].l[a[left][idx]]=seg[x].r[a[left][idx]]=seg[x].s[a[left][idx]]=1;
            return;
        }
        int mid = (left + right) / 2;
        seg[x].ls=++tot;
        build(left, mid, tot, idx);
        seg[x].rs=++tot;
        build(mid + 1, right, tot, idx);
        maintain(x);
    }
    T query(int l, int r, int left, int right, int x)
    {
        if(left==l && right==r) return seg[x].s[1];
        int mid = (left + right) / 2;
        if (r <= mid) return query(l, r, left, mid, seg[x].ls);
        else if (l > mid) return query(l, r, mid + 1, right, seg[x].rs);
    }
};