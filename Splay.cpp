#include <bits/stdc++.h>
using namespace std;
//======================================
const int maxn = 1e5+5;
struct Node
{
    int l,r;
    int val,size;
    int cnt;        //当前结点重复次数，默认为1
}spl[maxn];         //内存池
int cnt,root;       //内存池计数器与根节点编号
inline void newnode(int &now,int &val)
{
    spl[now=++cnt].val=val;
    spl[cnt].size++;
    spl[cnt].cnt++;
}
inline void update(int now) //更新size
{
    spl[now].size=spl[spl[now].l].size+spl[spl[now].r].size+spl[now].cnt;
}
inline void zig(int &now)
{
    int l = spl[now].l;
    spl[now].l = spl[l].r;
    spl[l].r = now;
    now = l;
    update(spl[now].r),update(now);
}
inline void zag(int &now)
{
    int r = spl[now].r;
    spl[now].r = spl[r].l;
    spl[r].l = now;
    now = r;
    update(spl[now].l),update(now);
}
void splaying(int x,int &y) //我要把x伸展到y那个位置！
{
    if(x==y) return;        //如果到了终点，return
    int &l = spl[y].l, &r = spl[y].r;   //temp
    if(x==l) zig(y);        //如果左儿子是终点，那就单旋
    else if(x==r) zag(y);   //右儿子是终点也是单旋
    else        //否则就一定是双旋
    {
        if(spl[x].val<spl[y].val)
        {
            if(spl[x].val<spl[l].val)
                splaying(x,spl[l].l),zig(y),zig(y);     //zigzig情况
            else splaying(x,spl[l].r),zag(l),zig(y);    //zagzig情况
        }
        else 
        {
            if(spl[x].val>spl[r].val)
                splaying(x,spl[r].r),zag(y),zag(y);     //zagzag情况
            else splaying(x,spl[r].l),zig(r),zag(y);    //zigzag情况
        }
    }
}
inline void delnode(int now)
{
    splaying(now,root);     //将要删除的结点伸展至根结点
    if(spl[now].cnt>1) spl[now].size--,spl[now].cnt--;  //如果有重复，令重复次数--
    else if(spl[root].r)    //否则如果当前结点有后继
    {
        int p = spl[root].r;
        while(spl[p].l) p=spl[p].l;     //找到后继
        splaying(p,spl[root].r);        //将其伸展至根结点右儿子
        spl[spl[root].r].l=spl[root].l; //右儿子左儿子变为根结点
        root=spl[root].r;               //根结点变为根结点右儿子
        update(root);                   //更新一下size信息
    }
    else root = spl[root].l;    //伸展之后没有后继，说明它是最大的了，那就直接删除
}
void ins(int &now,int &val)
{
    if(!now) newnode(now,val),splaying(now,root);
    else if(val<spl[now].val) ins(spl[now].l,val);
    else if(val>spl[now].val) ins(spl[now].r,val);
    else spl[now].size++,spl[now].cnt++,splaying(now,root); //特判相同的情况
}
void del(int now,int &val)
{
    if(spl[now].val==val) delnode(now);
    else if(val<spl[now].val) del(spl[now].l,val);
    else del(spl[now].r,val);
}
//以下大致与以前的代码相同，有大变动的地方给出了注释
int getrank(int val)
{
    int now = root, rank = 1;
    while(now)
    {
        if(spl[now].val==val)   //找到了要的结点，这个之前的没有
        {
            rank+=spl[spl[now].l].size;
            splaying(now,root);
            break;
        }
        if(val<=spl[now].val)
            now=spl[now].l;
        else 
        {
            rank+=spl[spl[now].l].size+spl[now].cnt;
            now=spl[now].r;
        }
    }
    return rank;
}
int getnum(int rank)
{
    int now = root;
    while(now)
    {
        int lsize = spl[spl[now].l].size;
        if(lsize+1<=rank&&rank<=lsize+spl[now].cnt) //如果在这个范围内，那就是当前结点
        {
            splaying(now,root);
            break;
        }
        else if(lsize>=rank)
            now=spl[now].l;
        else 
        {
            rank-=lsize+spl[now].cnt;
            now=spl[now].r;
        }
    }
    return spl[now].val;
}