#include <bits/stdc++.h>
using namespace std;
const int maxn = 1e5 + 10;
int n, m;
int cnt;
struct node
{
    int L, R; //分别指向左右子树
    int sum;  //该节点所管辖区间范围内数的个数
    node()
    {
        sum = 0;
    }
} Tree[maxn * 20];

struct value
{
    int x;  //值的大小
    int id; //离散之前在原数组中的位置
    bool operator<(const value &b)
        const
    {
        return x < b.x;
    }
} Value[maxn];

int root[maxn]; //多颗线段树的根节点
int rk[maxn];   //原数组离散之后的数组
void init()
{
    cnt = 1;
    root[0] = 0;
    Tree[0].L = Tree[0].R = Tree[0].sum = 0;
}

void update(int num, int &rt, int l, int r)
{
    Tree[cnt++] = Tree[rt];
    rt = cnt - 1;
    Tree[rt].sum++;
    if (l == r)
        return;
    int mid = (l + r) >> 1;
    if (num <= mid)
        update(num, Tree[rt].L, l, mid);
    else
        update(num, Tree[rt].R, mid + 1, r);
}

int query(int i, int j, int k, int l, int r)
{
    int d = Tree[Tree[j].L].sum - Tree[Tree[i].L].sum;
    if (l == r)
        return l;
    int mid = (l + r) >> 1;
    if (k <= d)
        return query(Tree[i].L, Tree[j].L, k, l, mid);
    else
        return query(Tree[i].R, Tree[j].R, k - d, mid + 1, r);
}

int main()
{
    scanf("%d%d", &n, &m);
    for (int i = 1; i <= n; i++)
    {
        scanf("%d", &Value[i].x);
        Value[i].id = i;
    }
    //进行离散化
    sort(Value + 1, Value + n + 1);
    for (int i = 1; i <= n; i++)
    {
        rk[Value[i].id] = i;
    }
    init();
    for (int i = 1; i <= n; i++)
    {
        root[i] = root[i - 1];
        update(rk[i], root[i], 1, n);
    }
    int left, right, k;
    for (int i = 1; i <= m; i++)
    {
        scanf("%d%d%d", &left, &right, &k);
        printf("%d\n", Value[query(root[left - 1], root[right], k, 1, n)].x);
    }
    return 0;
}
