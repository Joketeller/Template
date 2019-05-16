struct segtree{
    int lson,rson,tag;ll sum;
}t[200010<<5];//数组要开的够大，比一般线段树大
int root[200010],treesize;
void push_down(int k,int lenl,int lenr)//还要传左右儿子的区间长度
{
    int tmp=t[k].tag;t[k].tag=0;
	//检查左右儿子是否存在，有可能修改打了lazytag没递归下去建点
    if (!t[k].lson)t[t[k].lson=++treesize]=(segtree){0,0,0,0ll};
	if (!t[k].rson)t[t[k].rson=++treesize]=(segtree){0,0,0,0ll};
    t[t[k].lson].sum+=lenl*tmp;t[t[k].lson].tag+=tmp;
    t[t[k].rson].sum+=lenr*tmp;t[t[k].rson].tag+=tmp;
}
void modify(int &k,int old,int l,int r,int x,int y,int d)//在旧节点上修改并新建节点
{
    t[k=++treesize]=(segtree){0,0,0,0ll};
    //新建一个节点
    int mid=(l+r)>>1;
	if (t[old].tag&&l!=r)push_down(old,mid-l+1,r-mid);
    //如果旧节点有tag而且不是叶节点就要下传
    if (l==x&&r==y){t[k].sum=(ll)(r-l+1)*d;t[k].tag=d;return;}
    if (y<=mid)modify(t[k].lson,t[old].lson,l,mid,x,y,d),t[k].rson=t[old].rson;
    //左边递归处理，右边照抄上一版本
    else if (x>mid)t[k].lson=t[old].lson,modify(t[k].rson,t[old].rson,mid+1,r,x,y,d);
    //左边照抄上一版本，右边递归处理
    else modify(t[k].lson,t[old].lson,l,mid,x,mid,d),modify(t[k].rson,t[old].rson,mid+1,r,mid+1,y,d);
    //左右都递归处理
    push_up(k);
}
//modify(root[i],root[i-1],1,n,x,y,d);这样只在上一个版本上多加log个点