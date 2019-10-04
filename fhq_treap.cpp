/***********************按值划分***********************************/
const int maxn = 1e5+5;
struct Node
{
    int l,r;
    int val,key;
    int size;
}fhq[maxn];
int cnt,root;
#include <random>
std::mt19937 rnd(233);
inline int newnode(int val)
{
    fhq[++cnt].val=val;
    fhq[cnt].key=rnd();
    fhq[cnt].size=1;
    return cnt;
}
inline void update(int now)
{
    fhq[now].size=fhq[fhq[now].l].size+fhq[fhq[now].r].size+1;
}
void split(int now,int val,int &x,int &y)
{
    if(!now) x=y=0;
    else 
    {
        if(fhq[now].val<=val)
        {
            x=now;
            split(fhq[now].r,val,fhq[now].r,y);
        }
        else 
        {
            y=now;
            split(fhq[now].l,val,x,fhq[now].l);
        }
        update(now);
    }
}
int merge(int x,int y)
{
    if(!x||!y) return x+y;
    if(fhq[x].key>fhq[y].key)           // > >= < <=
    {
        fhq[x].r=merge(fhq[x].r,y);
        update(x);
        return x;
    }
    else 
    {
        fhq[y].l=merge(x,fhq[y].l);
        update(y);
        return y;
    }
}
int x,y,z;
inline void ins(int val)
{
    split(root,val,x,y);
    root=merge(merge(x,newnode(val)),y);
}
inline void del(int val)
{
    split(root,val,x,z);
    split(x,val-1,x,y);
    y=merge(fhq[y].l,fhq[y].r);
    root=merge(merge(x,y),z);
}
inline int getrank(int val)
{
    split(root,val-1,x,y);
    // print(fhq[x].size+1);
    int ret = fhq[x].size + 1;
    root = merge(x, y);
    return ret;
}
inline int getnum(int rank)
{
    int now=root;
    while(now)
    {
        if(fhq[fhq[now].l].size+1==rank)
            break;
        else if(fhq[fhq[now].l].size>=rank)
            now=fhq[now].l;
        else 
        {
            rank-=fhq[fhq[now].l].size+1;
            now=fhq[now].r;
        }
    }
    return fhq[now].val;
}
inline int pre(int val)
{
    split(root,val-1,x,y);
    int now = x;
    while(fhq[now].r)
        now = fhq[now].r;
    // print(fhq[now].val);
    int ret = fhq[now].val;
    root = merge(x, y);
    return ret;
}
inline int nxt(int val)
{
    split(root,val,x,y);
    int now = y;
    while(fhq[now].l)
        now = fhq[now].l;
    // print(fhq[now].val);
    int ret = fhq[now].val;
    root=merge(x,y);
    return ret;
}



/****************************按大小划分*****************************/
#include <bits/stdc++.h>
using namespace std;
//======================================
const int maxn = 1e5+5;
struct Node
{
    int l,r;
    int val,key;
    int size;
    bool reverse;
}fhq[maxn];
int cnt,root;
#include <random>
std::mt19937 rnd(233);
inline int newnode(int val)
{
    fhq[++cnt].val=val;
    fhq[cnt].key=rnd();
    fhq[cnt].size=1;
    return cnt;
}
inline void update(int now)
{
    fhq[now].size=fhq[fhq[now].l].size+fhq[fhq[now].r].size+1;
}
inline void pushdown(int now)
{
    std::swap(fhq[now].l,fhq[now].r);
    fhq[fhq[now].l].reverse^=1;
    fhq[fhq[now].r].reverse^=1;
    fhq[now].reverse=false;
}
void split(int now,int siz,int &x,int &y)
{
    if(!now) x=y=0;
    else 
    {
        if(fhq[now].reverse) pushdown(now);
        if(fhq[fhq[now].l].size<siz)
        {
            x=now;
            split(fhq[now].r,siz-fhq[fhq[now].l].size-1,fhq[now].r,y);
        }
        else 
        {
            y=now;
            split(fhq[now].l,siz,x,fhq[now].l);
        }
        update(now);
    }
}
int merge(int x,int y)
{
    if(!x||!y) return x+y;
    if(fhq[x].key<fhq[y].key)
    {
        if(fhq[x].reverse) pushdown(x);
        fhq[x].r=merge(fhq[x].r,y);
        update(x);
        return x;
    }
    else 
    {
        if(fhq[y].reverse) pushdown(y);
        fhq[y].l=merge(x,fhq[y].l);
        update(y);
        return y;
    }
}
void reverse(int l,int r)
{
    int x,y,z;
    split(root,l-1,x,y);
    split(y,r-l+1,y,z);
    fhq[y].reverse^=1;
    root=merge(merge(x,y),z);
}
void ldr(int now)
{
    if(!now) return;
    if(fhq[now].reverse) pushdown(now);
    ldr(fhq[now].l);
    std::cout<<fhq[now].val<<" ";
    ldr(fhq[now].r);
}