//****************************后缀自动机*******************************//
#include <bits/stdc++.h>
typedef long long ll;
using namespace std;
namespace SAM
{
const int Alphabet = 26;
const int MAXN = 5e5 + 10000;
struct SAMTree
{
    int size;
    struct node
    {
        int len;    //能够识别的最长子串长度
        int cnt;    //当前节点的right集合大小
        int rmax;   //right集中最大的位置
        node *fail; //失配指针
                    //  int val;    //information
        node *nxt[Alphabet];
    } t[MAXN], *last, *root;

    node *newnode(int _len)
    {
        node *k = t + (size++);
        k->len = _len;
        memset(k->nxt, 0, sizeof(k->nxt));
        k->fail = 0;
        k->rmax = 0;
        k->cnt = 0;
        return k;
    }

    void extend(int x, int id)
    {
        node *p = last, *np = newnode(last->len + 1);
        np->rmax = id;
        np->cnt = 1;
        while (p && !p->nxt[x])
            p->nxt[x] = np, p = p->fail;
        if (!p)
            np->fail = root;
        else
        {
            node *q = p->nxt[x];
            if (q->len == p->len + 1)
                np->fail = q;
            else
            {
                node *nq = newnode(p->len + 1);
                memcpy(nq->nxt, q->nxt, sizeof(q->nxt));
                nq->fail = q->fail;
                q->fail = np->fail = nq;
                while (p && p->nxt[x] == q)
                    p->nxt[x] = nq, p = p->fail;
            }
        }
        last = np;
    }

    void init()
    {
        size = 0;
        root = last = newnode(0);
    }

    int in[MAXN];
    queue<node *> q;
    void calculate_right_union()
    {
        memset(in, 0, sizeof(in));
        for (int i = 1; i < size; i++)
        {
            node *k = t + i;
            in[k->fail - t]++;
        }
        for (int i = 0; i < size; i++)
            if (!in[i])
                q.push(t + i);
        while (q.size())
        {
            node *f = q.front();
            q.pop();
            if (!f->fail)
                continue;
            f->fail->cnt += f->cnt;
            f->fail->rmax = max(f->fail->rmax, f->rmax);
            if (!(--in[f->fail - t]))
                q.push(f->fail);
        }
    }
} A;
    long long number(int len)
    {
        long long num = 0;
        for (int i = 0; i < A.size;i++)
            if (A.t[i].len>=len)
                num++;
        return num;
    }
}; // namespace SAM

char s[500005];
int main()
{
    #ifdef _IRONHEAD_
        assert(freopen("/Users/ironhead/algorithm/in.in", "r", stdin));
        assert(freopen("/Users/ironhead/algorithm/out.out", "w", stdout));
    #endif
    SAM::A.init();
    scanf("%s", s + 1);
    int len = strlen(s + 1);
    for (int i = 1; i <= len; i++)
    {
        SAM::A.extend(s[i] - 'a', i);
    }

    return 0;
}