//*************************************回文自动机***********************************//
//************************************非广义回文树************************************//
#include <bits/stdc++.h>
namespace PAM
{
const int MAXN = 5e5 + 100;
const int Alphabet = 26;
struct PAMtree
{
    int size;      //树的大小 ,-2即为不同回文串个数
    long long tot; // 回文串个数
    int S[MAXN];   // 表示串
    int L, R;

    struct node
    {
        int len;    //the palindromic length of the status this node represent
        int cnt;    //the times that this status appears,needs count()
        int num;    // 以节点i回文串的末尾字符结尾的但不包含本条路径上的回文串的数目。(也就是fail指针路径的深度)
        node *fail; //the fail pointer,i'd like to change this into a pointer,pointer is more convenient
        node *nxt[Alphabet];
    } t[MAXN], *last[2];

    node *newnode(int _len)
    {
        node *k = t + size;
        memset(k->nxt, 0, sizeof(k->nxt));
        k->len = _len;
        k->num = 0;
        k->cnt = 0;
        k->fail = &t[0]; //more inportant
        size++;
        return k;
    }
    void EERinit()
    {
        size = tot = 0;
        last[0] = last[1] = newnode(0);
        last[0]->fail = newnode(-1);
        last[0]->fail->fail = last[0]->fail;
        memset(S, -1, sizeof(S));
        L = MAXN >> 1;
        R = L - 1;
    }

    node *getfail(node *now, int pos)
    {
        if (pos)
        {
            while (S[R - now->len - 1] != S[R])
                now = now->fail;
        }
        else
        {
            while (S[L + now->len + 1] != S[L])
                now = now->fail;
        }
        return now;
    }

    inline void extend(int c, int pos) //pos=0 add preffix d=1 add suffix
    {
        if (c >= 'a' && c <= 'z')
            c -= 'a';
        if (c >= 'A' && c <= 'Z')
            c -= 'A';
        if (pos == 0)
            S[--L] = c;
        else
            S[++R] = c;
        node *cur = getfail(last[pos], pos);
        if (!cur->nxt[c])
        {
            node *now = newnode(cur->len + 2);
            now->fail = getfail(cur->fail, pos)->nxt[c];
            if (!now->fail)        //there exists
                now->fail = &t[0]; //a little problem
            cur->nxt[c] = now;
            now->num = now->fail->num + 1;
        }
        last[pos] = cur->nxt[c];
        last[pos]->cnt++;
        if (last[pos]->len == R - L + 1)
            last[pos ^ 1] = last[pos];
        tot += last[pos]->num;
    }

    void count() //统计本质相同侧回文串个数
    {
        for (int i = size - 1; i >= 0; i++)
        {
            t[i].fail->cnt += t[i].cnt;
        }
    }
} A;
} // namespace PAM
