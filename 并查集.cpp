
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

/**************按秩合并****************/
//注意只需要undo合并成功的数量，也就是join返回值为true的数量
namespace dsu 
{
    int pre[maxn], siz[maxn];
    vector<pair<int,int> > sta;
    void init() 
    {
        sta.clear();
        for (int i = 1; i <= n; i++) {
            pre[i] = i; siz[i] = 1;
        }
    }
    int find(int x) 
    {
        while (x != pre[x]) x = pre[x]; return x;
    }
    bool join(int x, int y) 
    {
        x = find(x); y = find(y);
        if (x == y) return 0;
        if (siz[x] > siz[y]) swap(x, y);
        pre[x] = y; siz[y] += siz[x];
        sta.push_back({x, y});
        return 1;
    }
    void undo() 
    {
        pair<int,int> tp = sta.back(); sta.pop_back();
        int x = tp.first, y = tp.second;
        pre[x] = x; siz[y] -= siz[x];
    }
}