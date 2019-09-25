void build(int l, int r, int &rt)
{
    rt = ++tot;
    sum[rt] = num[rt] = 0;
    if (l==r)
        return;
    int mid = (l + r) >> 1;
    build(l, mid, ls[rt]);
    build(mid + 1, r, rs[rt]);
};
 
void update(int last,ll p,int l,int r,int &rt)
{
    rt = ++tot;
    ls[rt] = ls[last];
    rs[rt] = rs[last];
    sum[rt] = sum[last] + lsh[p];
    num[rt] = num[last] + 1;
    if (l==r)
        return;
    int mid = (l + r) >> 1;
    if (p<=mid)
        update(ls[last], p, l, mid, ls[rt]);
    else
    {
        update(rs[last],p,mid+1,r,rs[rt]);
    }
}
 