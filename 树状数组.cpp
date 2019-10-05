
struct BIT
{
    int n;
    int value[MAXN];
    inline int low_bit(int x)
    {   
        return x & -x;
    }
    void init(int n)
    {
        this->n = n;
        memset(value + 1, 0, sizeof(int) * n);
    }
    void modify(int pos, int val)
    {
        while (pos <= n)
        {
            value[pos] = max(value[pos], val);
            pos += low_bit(pos);
        }
    }
    int query(int pos)
    {
        int ret = 0;
        while (pos)
        {
            ret = max(ret, value[pos]);
            pos -= low_bit(pos);
        }
        return ret;
    }
};
const int maxn = 2000, maxm = 2000;
struct DDBIT
{
    
    inline int lowbit(int x)
    {
        return x & -x;
    }
    int val[maxn][maxm];
    void change(int x, int y,int inc_val)
    {
        int j;
        while(x<=maxn)
        {
           j=y;
           while(j<=maxm)
            {
               val[x][j]+=inc_val;
                j+=lowbit(j);
            }
            x+=lowbit(x);
        }
    }   
    int ask(int i,int y)
    {
        int ans=0;
        int j;
        while(i!=0)
        {
            j=y;
            while(j!=0)
            {
                ans+=val[i][j];
                j-=lowbit(j);
            }
            i-=lowbit(i);
        }
        return ans;
    }
};


void ADDPRE(int x, int c)
{
     for (int i=x; i<=n; i+=i&(-i)) a[i] += c;
}
int SUMPRE(int x)
{
    int s = 0;
    for (int i=x; i>0; i-=i&(-i)) s += a[i];
    return s;
}