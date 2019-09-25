
# include <bits/stdc++.h>
# define lson l,mid,id<<1
# define rson mid+1,r,id<<1|1
using namespace std;
typedef long long LL;
const int maxn = 7e5+30;
struct node{
    int bit[31], len;
}a[3*maxn];
void fun(node &x, int val){
    if(x.len == 31) return;
    for(int i=30; ~i; --i){
        if(val&(1<<i)){
            if(x.bit[i] == 0){
                x.bit[i] = val;
                ++x.len;
                break;
            }
            else val ^= x.bit[i];
        }
    }
}
void update(int pos, int val, int l, int r, int id){
    fun(a[id], val);
    if(l == r) return;
    int mid = l+r>>1;
    if(pos <= mid) update(pos, val, lson);
    else update(pos, val, rson);
}
node Merge(node x, node y){
    if(x.len == 31) return x;
    if(y.len == 31) return y;
    node tmp = x;
    for(int i=30; ~i; --i)
        if(y.bit[i])
            fun(tmp, y.bit[i]);
    return tmp;
}
node query(int L, int R, int l, int r, int id){
    if(L<=l && R>=r) return a[id];
    int mid = l+r>>1;
    if(R <= mid) return query(L, R, lson);
    else if(L > mid) return query(L, R, rson);
    else return Merge(query(L, R, lson), query(L, R, rson));
}
int main()
{
      #ifdef _IRONHEAD_
            assert(freopen("/Users/ironhead/algorithm/in.in", "r", stdin));
            // assert(freopen("/Users/ironhead/algorithm/out.out", "w", stdout));
        #endif
    int T;
    int n, m, op, x, y;
    int temp;
    scanf("%d",&T);
    while (T--)
    {
        scanf("%d%d",&n,&m);
        int lastans=0;
        for (int i=1;i<=n;i++)
        {
            scanf("%d",&temp);
            update(i,temp,1,n+m,1);
        }
        int tot=n;
        for (int i=1;i<=m;i++)
        {
            scanf("%d",&op);
            if (op)
            {
                scanf("%d",&temp);
                temp^=lastans;
                update(++tot,temp,1,n+m,1);
            }
            else
            {
                scanf("%d%d",&x,&y);
                x^=lastans;
                y^=lastans;
                x=x%tot+1;
                y=y%tot+1;
                if (x>y) swap(x,y);
                node tmp=query(x,y,1,n+m,1);
                int ans = 0;
                for(int i=30; ~i; --i)
                    if((ans^tmp.bit[i]) > ans) ans ^= tmp.bit[i];
                printf("%d\n",ans);
                lastans=ans;
            }
        }
        memset(a,0,sizeof(a));
    }
    
}