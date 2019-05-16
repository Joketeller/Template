const int maxn = 100 + 10;
int par[maxn];     //父亲,  当par[x] = x时,x是所在的树的根
int Rank[maxn];    //树的高度

//初始化n个元素
void init(int n)
{
    for (int i = 0; i < n; i++) {
        par[i] = i;
        Rank[i] = 0;
    }
}

//查询树的根
int find(int x) {
    if (par[x] == x) {
        return x;
    }
    else {
        return par[x] = find(par[x]);
    }
}

//合并x和y所属集合
void unite(int x, int y) {
    x = find(x);
    y = find(y);
    if (x == y) return;
    
    if (Rank[x] < Rank[y]) {
        par[x] = y;
    } else {
        par[y] = x;
        if (Rank[x] == Rank[y]) Rank[x]++;    //如果x,y的树高相同,就让x的树高+1
    }
}

//判断x和y是否属于同一个集合
bool same(int x, int y) {
    return find(x) == find(y);
}