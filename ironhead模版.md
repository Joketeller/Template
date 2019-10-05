[TOC]



### 并查集

```c++
#include <bits/stdc++.h>
using namespace std;
/*********************按秩合并***********************/
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

/**************可回退并查集****************/
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


/*******************可持久化并查集*********************/
namespace FastIO
{
char buf[1 << 21], buf2[1 << 21], a[20], *p1 = buf, *p2 = buf, hh = '\n';
int p, p3 = -1;
void read() {}
void print() {}
inline int getc()
{
    return p1 == p2 && (p2 = (p1 = buf) + fread(buf, 1, 1 << 21, stdin), p1 == p2) ? EOF : *p1++;
}
inline void flush()
{
    fwrite(buf2, 1, p3 + 1, stdout), p3 = -1;
}
template <typename T, typename... T2>
inline void read(T &x, T2 &... oth)
{
    int f = 0;
    x = 0;
    char ch = getc();
    while (!isdigit(ch))
    {
        if (ch == '-')
            f = 1;
        ch = getc();
    }
    while (isdigit(ch))
    {
        x = x * 10 + ch - 48;
        ch = getc();
    }
    x = f ? -x : x;
    read(oth...);
}
template <typename T, typename... T2>
inline void print(T x, T2... oth)
{
    if (p3 > 1 << 20)
        flush();
    if (x < 0)
        buf2[++p3] = 45, x = -x;
    do
    {
        a[++p] = x % 10 + 48;
    } while (x /= 10);
    do
    {
        buf2[++p3] = a[p];
    } while (--p);
    buf2[++p3] = hh;
    print(oth...);
}
} // namespace FastIO
#define read FastIO::read
#define print FastIO::print
//======================================
const int maxn = 2e5+5;
int n;
struct Node
{
    int l,r,val;
}hjt[maxn*40*2];
int cnt,rootfa[maxn],rootdep[maxn],tot;
void build(int l,int r,int &now)
{
    now = ++cnt;
    if(l==r)
    {
        hjt[now].val=++tot;
        return;
    }
    int m = (l+r)>>1;
    build(l,m,hjt[now].l);
    build(m+1,r,hjt[now].r);
}
void modify(int l,int r,int ver,int &now,int pos,int val)
{
    hjt[now=++cnt]=hjt[ver];
    if(l==r)
    {
        hjt[now].val=val;
        return;
    }
    int m = (l+r)>>1;
    if(pos<=m) modify(l,m,hjt[ver].l,hjt[now].l,pos,val);
    else modify(m+1,r,hjt[ver].r,hjt[now].r,pos,val);
}
int query(int l,int r,int now,int pos)
{
    if(l==r) return hjt[now].val;
    int m = (l+r)>>1;
    if(pos<=m) return query(l,m,hjt[now].l,pos);
    else return query(m+1,r,hjt[now].r,pos);
}
int find(int ver,int x)
{
    int fx = query(1,n,rootfa[ver],x);
    return fx==x?x:find(ver,fx);
}
void merge(int ver,int x,int y)
{
    x = find(ver-1,x);          //ver-1
    y = find(ver-1,y);
    if(x==y)
    {
        rootfa[ver]=rootfa[ver-1];
        rootdep[ver]=rootdep[ver-1];
    }
    else
    {
        int depx = query(1,n,rootdep[ver-1],x);
        int depy = query(1,n,rootdep[ver-1],y);
        if(depx<depy)
        {
            modify(1,n,rootfa[ver-1],rootfa[ver],x,y);
            rootdep[ver]=rootdep[ver-1];
        }
        else if(depx>depy)
        {
            modify(1,n,rootfa[ver-1],rootfa[ver],y,x);
            rootdep[ver]=rootdep[ver-1];
        }
        else
        {
            modify(1,n,rootfa[ver-1],rootfa[ver],x,y);
            modify(1,n,rootdep[ver-1],rootdep[ver],y,depy+1);
        }
    }
}
int main(int argc, char const *argv[])
{
#ifndef ONLINE_JUDGE
    freopen("in.in", "r", stdin);
    freopen("out.out", "w", stdout);
#endif
    clock_t c1 = clock();
    //======================================
    int m;
    read(n,m);
    build(1,n,rootfa[0]);
    for(int ver=1;ver<=m;ver++)
    {
        int opt,x,y;
        read(opt);
        switch(opt)
        {
        case 1:
            read(x,y);
            merge(ver,x,y);
            break;
        case 2:
            read(x);
            rootfa[ver]=rootfa[x];
            rootdep[ver]=rootdep[x];
            break;
        case 3:
            read(x,y);
            rootfa[ver]=rootfa[ver-1];
            rootdep[ver]=rootdep[ver-1];
            int fx = find(ver,x);
            int fy = find(ver,y);
            print(fx==fy?1:0);
            break;
        }
    }
    //======================================
    FastIO::flush();
    std::cerr << "Time:" << clock() - c1 << "ms" << std::endl;
    return 0;
}
```



### 错排

```c++
long long d[26]={1,0};//要用long long 
 
void D()//求i个数错排有多少种
{
    for(int i=2;i<=14;i++)
        d[i]=(i-1)*(d[i-1]+d[i-2]);
}
```



### 大数模版

```c++
#include <bits/stdc++.h>

using namespace std; 

#define MAXN 9999
#define MAXSIZE 10
#define DLEN 4

class BigNum
{ 
private: 
	int a[50000];    //可以控制大数的位数 
	int len;       //大数长度
public: 
	BigNum(){ len = 1;memset(a,0,sizeof(a)); }   //构造函数
	BigNum(const int);       //将一个int类型的变量转化为大数
	BigNum(const char*);     //将一个字符串类型的变量转化为大数
	BigNum(const BigNum &);  //拷贝构造函数
	BigNum &operator=(const BigNum &);   //重载赋值运算符，大数之间进行赋值运算

	friend istream& operator>>(istream&,  BigNum&);   //重载输入运算符
	friend ostream& operator<<(ostream&,  BigNum&);   //重载输出运算符

	BigNum operator+(const BigNum &) const;   //重载加法运算符，两个大数之间的相加运算 
	BigNum operator-(const BigNum &) const;   //重载减法运算符，两个大数之间的相减运算 
	BigNum operator*(const BigNum &) const;   //重载乘法运算符，两个大数之间的相乘运算 
	BigNum operator/(const int   &) const;    //重载除法运算符，大数对一个整数进行相除运算

	BigNum operator^(const int  &) const;    //大数的n次方运算
	int    operator%(const int  &) const;    //大数对一个int类型的变量进行取模运算    
	bool   operator>(const BigNum & T)const;   //大数和另一个大数的大小比较
	bool   operator>(const int & t)const;      //大数和一个int类型的变量的大小比较

	void print();       //输出大数
}; 
BigNum::BigNum(const int b)     //将一个int类型的变量转化为大数
{ 
	int c,d = b;
	len = 0;
	memset(a,0,sizeof(a));
	while(d > MAXN)
	{
		c = d - (d / (MAXN + 1)) * (MAXN + 1); 
		d = d / (MAXN + 1);
		a[len++] = c;
	}
	a[len++] = d;
}
BigNum::BigNum(const char*s)     //将一个字符串类型的变量转化为大数
{
	int t,k,index,l,i;
	memset(a,0,sizeof(a));
	l=strlen(s);   
	len=l/DLEN;
	if(l%DLEN)
		len++;
	index=0;
	for(i=l-1;i>=0;i-=DLEN)
	{
		t=0;
		k=i-DLEN+1;
		if(k<0)
			k=0;
		for(int j=k;j<=i;j++)
			t=t*10+s[j]-'0';
		a[index++]=t;
	}
}
BigNum::BigNum(const BigNum & T) : len(T.len)  //拷贝构造函数
{ 
	int i; 
	memset(a,0,sizeof(a)); 
	for(i = 0 ; i < len ; i++)
		a[i] = T.a[i]; 
} 
BigNum & BigNum::operator=(const BigNum & n)   //重载赋值运算符，大数之间进行赋值运算
{
	int i;
	len = n.len;
	memset(a,0,sizeof(a)); 
	for(i = 0 ; i < len ; i++) 
		a[i] = n.a[i]; 
	return *this; 
}
istream& operator>>(istream & in,  BigNum & b)   //重载输入运算符
{
	char ch[MAXSIZE*4];
	int i = -1;
	in>>ch;
	int l=strlen(ch);
	int count=0,sum=0;
	for(i=l-1;i>=0;)
	{
		sum = 0;
		int t=1;
		for(int j=0;j<4&&i>=0;j++,i--,t*=10)
		{
			sum+=(ch[i]-'0')*t;
		}
		b.a[count]=sum;
		count++;
	}
	b.len =count++;
	return in;

}
ostream& operator<<(ostream& out,  BigNum& b)   //重载输出运算符
{
	int i;  
	cout << b.a[b.len - 1]; 
	for(i = b.len - 2 ; i >= 0 ; i--)
	{ 
		cout.width(DLEN); 
		cout.fill('0'); 
		cout << b.a[i]; 
	} 
	return out;
}

BigNum BigNum::operator+(const BigNum & T) const   //两个大数之间的相加运算
{
	BigNum t(*this);
	int i,big;      //位数   
	big = T.len > len ? T.len : len; 
	for(i = 0 ; i < big ; i++) 
	{ 
		t.a[i] +=T.a[i]; 
		if(t.a[i] > MAXN) 
		{ 
			t.a[i + 1]++; 
			t.a[i] -=MAXN+1; 
		} 
	} 
	if(t.a[big] != 0)
		t.len = big + 1; 
	else
		t.len = big;   
	return t;
}
BigNum BigNum::operator-(const BigNum & T) const   //两个大数之间的相减运算 
{  
	int i,j,big;
	bool flag;
	BigNum t1,t2;
	if(*this>T)
	{
		t1=*this;
		t2=T;
		flag=0;
	}
	else
	{
		t1=T;
		t2=*this;
		flag=1;
	}
	big=t1.len;
	for(i = 0 ; i < big ; i++)
	{
		if(t1.a[i] < t2.a[i])
		{ 
			j = i + 1; 
			while(t1.a[j] == 0)
				j++; 
			t1.a[j--]--; 
			while(j > i)
				t1.a[j--] += MAXN;
			t1.a[i] += MAXN + 1 - t2.a[i]; 
		} 
		else
			t1.a[i] -= t2.a[i];
	}
	t1.len = big;
	while(t1.a[len - 1] == 0 && t1.len > 1)
	{
		t1.len--; 
		big--;
	}
	if(flag)
		t1.a[big-1]=0-t1.a[big-1];
	return t1; 
} 

BigNum BigNum::operator*(const BigNum & T) const   //两个大数之间的相乘运算 
{ 
	BigNum ret; 
	int i,j,up; 
	int temp,temp1;   
	for(i = 0 ; i < len ; i++)
	{ 
		up = 0; 
		for(j = 0 ; j < T.len ; j++)
		{ 
			temp = a[i] * T.a[j] + ret.a[i + j] + up; 
			if(temp > MAXN)
			{ 
				temp1 = temp - temp / (MAXN + 1) * (MAXN + 1); 
				up = temp / (MAXN + 1); 
				ret.a[i + j] = temp1; 
			} 
			else
			{ 
				up = 0; 
				ret.a[i + j] = temp; 
			} 
		} 
		if(up != 0) 
			ret.a[i + j] = up; 
	} 
	ret.len = i + j; 
	while(ret.a[ret.len - 1] == 0 && ret.len > 1)
		ret.len--; 
	return ret; 
} 
BigNum BigNum::operator/(const int & b) const   //大数对一个整数进行相除运算
{ 
	BigNum ret; 
	int i,down = 0;   
	for(i = len - 1 ; i >= 0 ; i--)
	{ 
		ret.a[i] = (a[i] + down * (MAXN + 1)) / b; 
		down = a[i] + down * (MAXN + 1) - ret.a[i] * b; 
	} 
	ret.len = len; 
	while(ret.a[ret.len - 1] == 0 && ret.len > 1)
		ret.len--; 
	return ret; 
}
int BigNum::operator %(const int & b) const    //大数对一个int类型的变量进行取模运算    
{
	int i,d=0;
	for (i = len-1; i>=0; i--)
	{
		d = ((d * (MAXN+1))% b + a[i])% b;  
	}
	return d;
}
BigNum BigNum::operator^(const int & n) const    //大数的n次方运算
{
	BigNum t,ret(1);
	int i;
	if(n<0)
		exit(-1);
	if(n==0)
		return 1;
	if(n==1)
		return *this;
	int m=n;
	while(m>1)
	{
		t=*this;
		for( i=1;i<<1<=m;i<<=1)
		{
			t=t*t;
		}
		m-=i;
		ret=ret*t;
		if(m==1)
			ret=ret*(*this);
	}
	return ret;
}
bool BigNum::operator>(const BigNum & T) const   //大数和另一个大数的大小比较
{ 
	int ln; 
	if(len > T.len)
		return true; 
	else if(len == T.len)
	{ 
		ln = len - 1; 
		while(a[ln] == T.a[ln] && ln >= 0)
			ln--; 
		if(ln >= 0 && a[ln] > T.a[ln])
			return true; 
		else
			return false; 
	} 
	else
		return false; 
}
bool BigNum::operator >(const int & t) const    //大数和一个int类型的变量的大小比较
{
	BigNum b(t);
	return *this>b;
}

void BigNum::print()    //输出大数
{ 
	int i;   
	cout << a[len - 1]; 
	for(i = len - 2 ; i >= 0 ; i--)
	{ 
		cout.width(DLEN); 
		cout.fill('0'); 
		cout << a[i]; 
	} 
	cout << endl;
}

BigNum ans;
int main()
{
    ios::sync_with_stdio(false);
    int n;
    while (cin>>n)
    {
        ans = 1;
        for (int i = 1; i <= n;i++)
            ans = ans * i;
        cout << ans << endl;
    }
}

```



### 大素数判断及质因子分解

```c++
//质因子分解，小数据用筛法直接判，大数据用pollard_rho
//map[i]是含有多少个质因子i
//map<LL,int>::iterator c,c->first表示质因子，c->second表示次方
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;


typedef long long LL;
#define maxn_for_division 1000001
bool is_prime[maxn_for_division];
vector<int>prime;
map<LL,int>factor;
inline void get_prime() //素数筛
{
    for(int i=0;i<maxn_for_division;i++)is_prime[i]=1;
    is_prime[0]=is_prime[1]=0;
    for(int i=2;i<maxn_for_division;i++)
        if (is_prime[i])
        {
            prime.push_back(i);
            for (int j=i;j<maxn_for_division;j+=i)is_prime[j]=0;
        }
}
inline LL gcd(LL a,LL b){if (!b)return a;return gcd(b,a%b);}

inline LL mod_mul(LL a, LL b, LL p) //大数乘
{
    LL ans=0ll;
    a%=p,b%=p;
    if (a>b)    swap(a,b);
    while (b)
    {
        if (b&1)ans=(ans+a)%p;
        a=(a+a)%p;
        b>>=1;
    }
    return ans;
}
LL mod_pow(LL a,LL b,LL p)
{
    LL ans=1ll;
    a%=p;
    while (b)
    {
        if (b&1)ans=mod_mul(ans,a,p);
        a=mod_mul(a,a,p);
        b>>=1;
    }
    return ans;
}
bool witness(LL a,LL n)
{
    LL m=n-1;
    int j=0;
    while(!(m&1))j++,m>>=1;
    LL x=mod_pow(a,m,n);
    if (x==1||x==n-1)return 0;
    while(j--)
    {
        x=mod_mul(x,x,n);
        if(x==n-1)return 0;
    }
    return 1;
}
#define rep_times 20
bool Miller_Rabin(LL n)//判断n是否为素数
{
    srand(time(0));
    if(n<2)return 0;
    if(n==2)return 1;
    if (!(n&1))return 0;
    for(int i=0;i<rep_times;i++)
    {
        LL a=rand()%(n-1)+1;
        if (witness(a,n))return 0;
    }
    return 1;
}
#undef rep_times
LL Pollard_Rho(LL n,int c)
{
    LL x=2,y=2,d=1;
    while (d==1)
    {
        x=mod_mul(x,x,n)+c;
        y=mod_mul(y,y,n)+c;
        y=mod_mul(y,y,n)+c;
        d=gcd((x-y>=0?x-y:y-x),n);
    }
    if (d==n)return Pollard_Rho(n,c+1);
    return d;
}
bool Is_Prime(LL n)
{
    return n<maxn_for_division&&is_prime[n]||n>=maxn_for_division&&Miller_Rabin(n);
}
void Find_Factor(LL n)
{
    if (Is_Prime(n)){factor[n]++;return;}
    for (int i=0;i<prime.size()&&prime[i]<=n;i++)
        if (n%prime[i]==0)
        {
            while (n%prime[i]==0)
            {
                factor[prime[i]]++;
                n/=prime[i];
            }
        }
    if (n!=1)
    {
        if (Is_Prime(n))factor[n]++;
        else
        {
            LL p=Pollard_Rho(n,1);
            Find_Factor(p);
            Find_Factor(n/p);
        }
    }
}
int main()
{
    LL n;
	get_prime();
    while(0)
    {
        factor.clear();
        Find_Factor(n);
	//----------------------------output-------------------------------
        for(map<LL,int>::iterator c=factor.begin();c!=factor.end();)
        {
            printf("%lld^%d",c->first,c->second);
            if((++c)!=factor.end())printf("*");
        }
	//-----------------------------------------------------------------
        puts("");
    }
}
```



### 带修莫队

```c++
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cstring>
#include<algorithm>
#include<iostream>
#define fo(i,a,b) for(i=a;i<=b;i++)
#define fod(i,a,b) for(i=a;i>=b;i--)
using namespace std;
struct note 
{
    int x,y,x1,y1,z,wz;
}cz[50005];
bool cmp(note x,note y)
{
    return (x.x1<y.x1||(x.x1==y.x1&&x.y1<y.y1)||(x.x1==y.x1&&x.y1==y.y1&&x.z<y.z));
}
int n,m,a[50005],b[1000001],cg[50005][3],s[50005],size;
void up(int q,int k)
{
    int i;
    fo(i,cz[q].z+1,cz[k].z) 
    {
        if(k>1&&cg[i][0]>=cz[q].x&&cg[i][0]<=cz[q].y)
        {
            if(b[cg[i][1]]==0) s[cz[k].wz]++;
            b[cg[i][1]]++;  
        }
        cg[i][2]=a[cg[i][0]];
        a[cg[i][0]]=cg[i][1];
        if(k>1&&cg[i][0]>=cz[q].x&&cg[i][0]<=cz[q].y)
        {
            b[cg[i][2]]--; 
            if(b[cg[i][2]]==0) s[cz[k].wz]--;
        }
    }
}
void down(int q,int k)
{
    int i;
    fod(i,cz[q].z,cz[k].z+1) 
    {
        if(k>1&&cg[i][0]>=cz[q].x&&cg[i][0]<=cz[q].y)
        {
            if(b[cg[i][2]]==0) s[cz[k].wz]++;
            b[cg[i][2]]++; 
        }
        a[cg[i][0]]=cg[i][2]; 
        if(k>1&&cg[i][0]>=cz[q].x&&cg[i][0]<=cz[q].y)
        {
            b[cg[i][1]]--; 
            if(b[cg[i][1]]==0) s[cz[k].wz]--;
        }
        cg[i][2]=0;
    }
}
void make(int k,int l,int r,int v)
{
    int i;
    fo(i,l,r)
    {
        if (b[a[i]]==0&&v==1) s[cz[k].wz]++;
        b[a[i]]+=v;
        if (b[a[i]]==0&&v==-1) s[cz[k].wz]--; 
    }
}
int main()
{
    int i,j;
    cin>>n>>m;
    memset(b,0,sizeof(b));
    fo(i,1,n) scanf("%d",&a[i]); 
    scanf("\n");
    size=(int)pow(n,2.0/3);
    int time=0,num=0;
    fo(i,1,m)
    {
        char ch;
        int x,y;
        scanf("%c%d%d",&ch,&x,&y);
        scanf("\n");
        x++;
        if (ch=='Q')
        {
            cz[++num].x=x;
            cz[num].y=y;
            cz[num].x1=(x-1)/size+1;
            cz[num].y1=(y-1)/size+1;
            cz[num].z=time;
            cz[num].wz=num;
        }
        else 
        {
            cg[++time][0]=x;
            cg[time][1]=y;
        }
    }
    sort(cz+1,cz+num+1,cmp);
    cz[0].z=0;
    up(0,1);
    fo(i,cz[1].x,cz[1].y) 
    {
        if (b[a[i]]==0) s[cz[1].wz]++;
        b[a[i]]++;
    }
    fo(i,2,num)
    {
        s[cz[i].wz]=s[cz[i-1].wz]; 
        if(cz[i].z<cz[i-1].z) down(i-1,i);
        else up(i-1,i);
        int x=cz[i-1].x,x1=cz[i].x,y=cz[i-1].y,y1=cz[i].y;
        if(x<x1) make(i,x,x1-1,-1);
        else make(i,x1,x-1,1);
        if(y<y1) make(i,y+1,y1,1);
        else make(i,y1+1,y,-1); 
    } 
    fo(i,1,num) printf("%d\n",s[i]);
} 

```



### 点

```c++
// point Plus
int sgn(const double &x){ return x < -eps? -1 : (x > eps);}
inline double sqr(const double &x){ return x * x;}
struct Point
{
    double x, y;
    Point(const double &x = 0, const double &y = 0):x(x), y(y){}
    Point operator -(const Point &a)const{ return Point(x - a.x, y - a.y);}
    Point operator +(const Point &a)const{ return Point(x + a.x, y + a.y);}
    Point operator *(const double &a)const{ return Point(x * a, y * a);}
    Point operator /(const double &a)const{ return Point(x / a, y / a);}
    bool operator < (const Point &a)const{ return sgn(x - a.x) < 0 || (sgn(x - a.x) == 0 && sgn(y - a.y) < 0);}
    bool operator == (const Point &a)const{ return sgn(sgn(x - a.x) == 0 && sgn(y - a.y) == 0);}
    friend double det(const Point &a, const Point &b){ return a.x * b.y - a.y * b.x;}
    friend double dot(const Point &a, const Point &b){ return a.x * b.x + a.y * b.y;}
    friend double dist(const Point &a, const Point &b){ return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));}
    void in(){ scanf("%lf %lf", &x, &y); }
    void out()const{ printf("%lf %lf\n", x, y); }
};
```



### 点分治

```c++
/******************  树上距离为k的点对是否存在  *******************/
#include<iostream>
#include<vector>
#include<algorithm>
#include<queue>
#include<cstring>
#include<cstdio>
using namespace std;

int read()
{
    int f=1,x=0;
    char ss=getchar();
    while(ss<'0'||ss>'9'){if(ss=='-')f=-1;ss=getchar();}
    while(ss>='0'&&ss<='9'){x=x*10+ss-'0';ss=getchar();}
    return f*x;
}

const int inf=10000000;
const int maxn=100010;
int n,m;
struct node{int v,dis,nxt;}E[maxn<<1];
int tot,head[maxn];
int maxp[maxn],size[maxn],dis[maxn],rem[maxn];
int vis[maxn],test[inf],judge[inf],q[maxn];
int query[1010];
int sum,rt;
int ans;

void add(int u,int v,int dis)
{
    E[++tot].nxt=head[u];
    E[tot].v=v;
    E[tot].dis=dis;
    head[u]=tot;
}

void getrt(int u,int pa)
{
    size[u]=1; maxp[u]=0;
    for(int i=head[u];i;i=E[i].nxt) 
    {
        int v=E[i].v;
        if(v==pa||vis[v]) continue;
        getrt(v,u);
        size[u]+=size[v];
        maxp[u]=max(maxp[u],size[v]);
    }
    maxp[u]=max(maxp[u],sum-size[u]);
    if(maxp[u]<maxp[rt]) rt=u;
}

void getdis(int u,int fa)
{
    rem[++rem[0]]=dis[u];
    for(int i=head[u];i;i=E[i].nxt)
    {
        int v=E[i].v;
        if(v==fa||vis[v])continue;
        dis[v]=dis[u]+E[i].dis;
        getdis(v,u);
    }
}

void calc(int u)
{
    int p=0;
    for(int i=head[u];i;i=E[i].nxt)
    {
        int v=E[i].v;
        if(vis[v])continue;
        rem[0]=0; dis[v]=E[i].dis;
        getdis(v,u);//处理u的每个子树的dis

        for(int j=rem[0];j;--j)//遍历当前子树的dis
        for(int k=1;k<=m;++k)//遍历每个询问
        if(query[k]>=rem[j])
        test[k]|=judge[query[k]-rem[j]];
        //如果query[k]-rem[j]d的路径存在就标记第k个询问

        for(int j=rem[0];j;--j)//保存出现过的dis于judge
        q[++p]=rem[j],judge[rem[j]]=1;
    }
    for(int i=1;i<=p;++i)//处理完这个子树就清空judge
    judge[q[i]]=0;//特别注意一定不要用memset，会T

}

void solve(int u)
{   
    //judge[i]表示到根距离为i的路径是否存在
    vis[u]=judge[0]=1; calc(u);//处理以u为根的子树
    for(int i=head[u];i;i=E[i].nxt)//对每个子树进行分治
    {
        int v=E[i].v;
        if(vis[v])continue;
        sum=size[v]; maxp[rt=0]=inf;
        getrt(v,0); solve(rt);//在子树中找重心并递归处理
    }
}

int main()
{
    n=read();m=read();
    for(int i=1;i<n;++i)
    {
        int u=read(),v=read(),dis=read();
        add(u,v,dis);add(v,u,dis);
    }
    for(int i=1;i<=m;++i)
    query[i]=read();//先记录每个询问以离线处理

    maxp[rt]=sum=n;//第一次先找整棵树的重心
    getrt(1,0); 
    solve(rt);//对树进行点分治

    for(int i=1;i<=m;++i)
    {
        if(test[i]) printf("AYE\n");
        else printf("NAY\n");
    }
    return 0;
}

/***************************************************/

/*****************树上小于等于k的点对个数********************/
#include<bits/stdc++.h>
#define ll long long
using namespace std;
inline int read(){
    int x=0;char ch=' ';int f=1;
    while(ch!='-'&&(ch<'0'||ch>'9'))ch=getchar();
    if(ch=='-')f=-1,ch=getchar();
    while(ch>='0'&&ch<='9')x=x*10+ch-'0',ch=getchar();
    return x*f;
}
struct edge{
    int to,next,w;
}e[80001];
int n,tot,root;
ll k;
int head[40001];
inline void addedge(int x,int y,int l){
    e[++tot].to=y;e[tot].next=head[x];e[tot].w=l;head[x]=tot;
}
int size[40001],vis[40001],mx,sz;
ll dis[40001],q[40001],l,r;
void getroot(int x,int fa){
    size[x]=1;int num=0;
    for(int i=head[x];i;i=e[i].next){
        int u=e[i].to;
        if(u==fa||vis[u])continue;
        getroot(u,x);
        size[x]+=size[u];
        num=max(num,size[u]);
    }
    num=max(num,sz-size[x]);
    if(num<mx){
        mx=num;root=x;
    }
}
void getdis(int x,int fa){
    q[++r]=dis[x];
    for(int i=head[x];i;i=e[i].next){
        int u=e[i].to;
        if(u==fa||vis[u])continue;
        dis[u]=dis[x]+e[i].w;
        getdis(u,x);
    }
}
ll calc(int x,int v){
    r=0;
    dis[x]=v;
    getdis(x,0);
    ll sum=0;
    l=1;
    sort(q+1,q+r+1);
    while(l<r){
        if(q[l]+q[r]<=k)sum+=r-l,l++;
        else r--;
    }
    return sum;
}
ll ans;
void dfs(int x){
    ans+=calc(x,0);
    vis[x]=1;
    for(int i=head[x];i;i=e[i].next){
        int u=e[i].to;
        if(vis[u])continue;
        ans-=calc(u,e[i].w);
        sz=size[u];
        mx=0x3f3f3f3f;
        getroot(u,0);
        dfs(root);
    }
}
int main(){
    n=read();
    for(int i=1;i<n;i++){
        int x=read(),y=read(),l=read();
        addedge(x,y,l);addedge(y,x,l);
    }
    k=read();
    sz=n;
    mx=0x3f3f3f3f;
    getroot(1,0);
    dfs(root);
    printf("%lld",ans);
    return 0;
}

/**********************************************************/


#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<cstring>
#include<cmath>
#include<algorithm>
#define maxn 10002
#define INF 2147483646
#define ll long long
using namespace std;

inline int gi()
{
  int date = 0,m = 1; char ch = 0;
  while(ch!='-'&&(ch<'0'||ch>'9'))ch = getchar();
  if(ch=='-'){m = -1; ch = getchar();}
  while(ch>='0' && ch<='9')
    {
      date = date*10+ch-'0';
      ch = getchar();
    }return date*m;
}
inline void write(ll qw)
{
  if(qw<0){putchar('-'); qw = -qw;}
  if(qw>9)write(qw/10);
  putchar(qw%10+'0');
}

struct node{int to,next,len;}t[2*maxn+5];
int head[maxn+5];
int sim[maxn+5],mxson[maxn+5];
bool vis[maxn+5];
int MX,root,dis[maxn+5],summar;
ll ans;
int n,k,cnt,S;

void addedge(int u,int v,int l)
{
  cnt ++;
  t[cnt].to = v;
  t[cnt].next = head[u];
  t[cnt].len = l;
  head[u] = cnt;
  return;
}

void getroot(int u,int fa)
{
  sim[u] = 1; mxson[u] = 0;
  for(int i = head[u];i;i = t[i].next)
    {
      int v = t[i].to;
      if(v == fa || vis[v])continue;
      getroot(v,u);
      sim[u] = sim[u] + sim[v];
      mxson[u] = max(sim[v],mxson[u]);
    }
  mxson[u] = max(mxson[u],S-sim[u]);
  if(mxson[u]<MX){root = u;MX = mxson[u];}
}

void getdis(int u,int fa,int dist)
{
  dis[++summar] = dist;
  for(int i = head[u];i;i=t[i].next)
    {
      int v = t[i].to;
      if(vis[v]||v == fa)continue;
      getdis(v,u,dist+t[i].len);
    }
  return;
}

int consolate(int sta,int len1)
{
  summar = 0;
  memset(dis,0,sizeof(dis));
  getdis(sta,0,len1);
  sort(dis+1,dis+summar+1);
  int L = 1,R = summar,tep=0;
  while(L<=R)
    {
      if(dis[R] + dis[L]<=k){tep = tep + R - L; L++;}
      else R--;
    }
  return tep;
}

void Divide(int tr)
{
  ans = ans + consolate(tr,0);
  vis[tr] = true;
  for(int i = head[tr];i;i = t[i].next)
    {
      int v = t[i].to;
      if(vis[v])continue;
      ans = ans - consolate(v,t[i].len);
      S = sim[v]; root = 0;
      MX = INF; getroot(v,0);
      Divide(root);
    }
  return;
}

int main()
{
  freopen("Tree.in","r",stdin);
  int u,v,l;
  while(1)
    {
      n = gi(); k = gi();
      if(n == 0 && k == 0)break;
      cnt = 0; memset(head,0,sizeof(head));
      for(int i = 1; i <= n - 1; i ++)
    {
      u = gi(); v = gi(); l = gi();
      addedge(u,v,l); addedge(v,u,l);
    }
      ans = 0;
      memset(vis,0,sizeof(vis));
      MX = INF ; S = n; getroot(1,0);
      Divide(root);
      write(ans);putchar('\n');
    }
  return 0;
}
```



### 动态DP

```c++
#include <bits/stdc++.h>
using namespace std;

const int maxn = 1e5+5;
const int INF = 0x3f3f3f3f;     //极大值，-INF就是极小值
//! 这里注意极小值不能使用INT_MIN,因为修改操作有负数，会导致小爆INT
//TODO 链式前向星
struct E
{
    int to,next;
}Edge[maxn<<1];
int tot,Head[maxn];
inline void AddEdge(int u,int v)
{
    Edge[tot].to=v;
    Edge[tot].next=Head[u];
    Head[u]=tot++;
}
//TODO 为了学这个又写了个简化版的矩阵模板，方便快捷好写好用，可以学习一下，Github上有
template <int row,int col>
struct Matrix
{
    int r,c;
    int ele[row][col];          //元素
    Matrix():r(row),c(col) {}
    int& operator()(int a,int b) { return ele[a][b]; }  //用括号运算符获取矩阵中元素
};
template <int m,int n,int p>
auto operator*(Matrix<m,n> m1,Matrix<n,p> m2)   //如果不用C++14就要把前面的auto换成Matrix<m,p>
{
    Matrix<m,p> ret;            //用于存储结果
    //? 因为新的矩乘法则，初始化为极小值
    memset(ret.ele,0xcf,sizeof(ret.ele));
    for(int i=0;i<m;i++)        //同样可以使用ikj优化
        for(int k=0;k<n;k++)
            for(int j=0;j<p;j++)
                ret(i,j)=std::max(ret(i,j),m1(i,k)+m2(k,j));
    return ret;
}
Matrix<2,2> ident,g[maxn];  //ident为单位矩阵，在主函数中初始化，g为矩阵池，存放x结点的矩阵是什么
int val[maxn],f[2][maxn];   //val存放各结点权值，f用于初始化矩阵池g
//TODO 树链剖分第一遍dfs，与普通树链剖分无区别
int fa[maxn],dep[maxn],siz[maxn],son[maxn];
void dfs1(int u,int f)
{
    fa[u]=f,siz[u]=1,dep[u]=dep[f]+1;
    for(int i=Head[u];~i;i=Edge[i].next)
    {
        int v = Edge[i].to;
        if(v==f) continue;
        dfs1(v,u);
        siz[u]+=siz[v];
        if(siz[v]>siz[son[u]])
            son[u]=v;
    }
}
//TODO 树链剖分第二遍dfs，在dfs过程中初始化矩阵池
int tim,dfn[maxn],nfd[maxn],top[maxn],end[maxn];
void dfs2(int u,int t)
{
    top[u]=t,dfn[u]=++tim,nfd[tim]=u;   //nfd[x]存放dfn为x的是哪个结点
    if(!son[u])                         //如果没有重儿子（没有儿子）即为叶子结点
    {
        f[1][u]=val[u];                 //叶子结点的f
        g[u]=ident;                     //叶子结点的矩阵为单位矩阵
        end[u]=u;                       //叶子结点所在重链尾端即为此叶子结点
        return;
    }
    g[u](1,0)=val[u],g[u](1,1)=-INF;    //初始化矩阵，左下角放上点权，右下角一直是INF
    dfs2(son[u],t);
    end[u]=end[son[u]];         //整条重链的end值都等于最下方叶结点，这样一层一层传上来
    for(int i=Head[u];~i;i=Edge[i].next)
    {
        int v = Edge[i].to;
        if(v==fa[u]||v==son[u]) continue;
        dfs2(v,v);                      //这都是轻儿子
        g[u](0,0)=g[u](0,1)+=std::max(f[0][v],f[1][v]);//每个轻儿子处理好后，将值加到g(i,0)和g(i,1)里
        g[u](1,0)+=f[0][v];
    }
    //? g处理好后即可处理g，f再用于处理其父亲结点的g
    f[0][u]=g[u](0,0)+std::max(f[0][son[u]],f[1][son[u]]);
    f[1][u]=g[u](1,0)+f[0][son[u]];
}
//TODO 线段树
struct Node
{
    int l,r;
    Matrix<2,2> val;
}sgt[maxn<<2];
#define ls(x) (x<<1)
#define rs(x) (x<<1|1)
#define pushup(x) sgt[x].val=sgt[ls(x)].val*sgt[rs(x)].val
void build(int l,int r,int k=1)
{
    sgt[k].l=l,sgt[k].r=r;
    if(l==r)
    {
        sgt[k].val=g[nfd[l]];   //在矩阵池里获取应该赋值的矩阵
        return;
    }
    int m = (l+r)>>1;
    build(l,m,ls(k));
    build(m+1,r,rs(k));
    pushup(k);
}
auto query(int x,int y,int k=1) //区间查询乘积
{
    int l=sgt[k].l,r=sgt[k].r;
    if(x<=l&&y>=r) return sgt[k].val;
    int m = (l+r)>>1;
    Matrix<2,2> ret = ident;    //应初始化为单位矩阵
    if(x<=m) ret = query(x,y,ls(k));
    if(y>m) ret = ret * query(x,y,rs(k));   //注意右乘
    return ret;
}
void modify(int x,int p,int k=1) //单点修改x结点的矩阵为矩阵池中p位置上的矩阵
{
    int l=sgt[k].l,r=sgt[k].r;
    if(l==r)
    {
        sgt[k].val=g[p];
        return;
    }
    int m = (l+r)>>1;
    if(x<=m) modify(x,p,ls(k));
    else modify(x,p,rs(k));
    pushup(k);
}
auto queryt(int x)          //查询x结点的f(x,0),f(x,1)矩阵
{
    Matrix<2,1> tmp;
    tmp(0,0)=0,tmp(1,0)=val[end[x]];        //初始矩阵，PPT里有说
    return query(dfn[x],dfn[end[x]]) * tmp; //注意左乘
}
void modifyt(int x,int z)   //单点修改x结点的权值为z
{
    Matrix<2,1> od,nw;
    g[x](1,0)+=z-val[x];    //相当于先减val[x]再加z，加法交换律
    od = queryt(top[x]);    //修改之前先查询一个旧的矩阵
    val[x]=z;               //然后修改，因为queryt函数里初始矩阵是用val[x]算出来的
    while(x)                //只要x没跳到头，这样写的话主函数里dfs1调用时第二个参数就必须为0
    {
        modify(dfn[x],x);       //先修改
        nw = queryt(top[x]);    //再查一个新的矩阵出来
        x = fa[top[x]];
        //? 在矩阵池里修改，用于下次线段树上修改
        g[x](0,0)=g[x](0,1)+=std::max(nw(0,0),nw(1,0))-std::max(od(0,0),od(1,0));
        g[x](1,0)+=nw(0,0)-od(0,0);
        od = queryt(top[x]);    //下次修改之前先查一个老的出来
    }
}
int main(int argc, char const *argv[])
{
#ifndef ONLINE_JUDGE
    freopen("in.in", "r", stdin);
    freopen("out.out", "w", stdout);
#endif
    clock_t c1 = clock();
    //======================================
    memset(Head,-1,sizeof(Head));   //链式前向星-1写法
    ident(1,0)=ident(0,1)=-INF;     //单位矩阵赋值
    int n,m;
    read(n,m);
    for(int i=1;i<=n;i++) read(val[i]);
    for(int i=1;i<n;i++)
    {
        int u,v;
        read(u,v);
        AddEdge(u,v),AddEdge(v,u);
    }
    dfs1(1,0);                      //这里注意(1,0)
    dfs2(1,1);
    build(1,n);
    while(m--)
    {
        int x,z;
        read(x,z);
        modifyt(x,z);
        auto ans = queryt(1);               //查询根结点矩阵，结果为一个2×1矩阵
        print(std::max(ans(0,0),ans(1,0))); //查出的矩阵中较大的元素即为答案
    }
    //======================================
    FastIO::flush();
    std::cerr << "Time:" << clock() - c1 << "ms" << std::endl;
    return 0;
}
```



### 读入挂

```c++

/**********************read**************************/
template<typename T>void read(T &x)
{
    x=0; int w=0; char c=0;
    while (c<'0'||c>'9') w|=c=='-',c=getchar();
    while (c>='0'&&c<='9') x=(x<<1)+(x<<3)+(c^48),c=getchar();
    x=w?-x:x;
}
template <typename T,typename... Args> inline void read(T& t, Args&... args)
{
    read(t);read(args...);
}
```



### 杜教筛

```c++
#include<bits/stdc++.h>
#include<unordered_map>
#define N 6000010
using namespace std;
template<typename T>inline void read(T &x)
{
    x=0;
    static int p;p=1;
    static char c;c=getchar();
    while(!isdigit(c)){if(c=='-')p=-1;c=getchar();}
    while(isdigit(c)) {x=(x<<1)+(x<<3)+(c-48);c=getchar();}
    x*=p;
}
bool vis[N];
int mu[N],sum1[N],phi[N];
long long sum2[N];
int cnt,prim[N];
unordered_map<long long,long long>w1;
unordered_map<int,int>w;
void get(int maxn)
{
    phi[1]=mu[1]=1;
    for(int i=2;i<=maxn;i++)
    {
        if(!vis[i])
        {
            prim[++cnt]=i;
            mu[i]=-1;phi[i]=i-1;
        }
        for(int j=1;j<=cnt&&prim[j]*i<=maxn;j++)
        {
            vis[i*prim[j]]=1;
            if(i%prim[j]==0)
            {
                phi[i*prim[j]]=phi[i]*prim[j];
                break;
            }
            else mu[i*prim[j]]=-mu[i],phi[i*prim[j]]=phi[i]*(prim[j]-1);
        }
    }
    for(int i=1;i<=maxn;i++)sum1[i]=sum1[i-1]+mu[i],sum2[i]=sum2[i-1]+phi[i];
}
int djsmu(int x)
{
    if(x<=6000000)return sum1[x];
    if(w[x])return w[x];
    int ans=1;
    for(int l=2,r;l>=0&&l<=x;l=r+1)
    {
        r=x/(x/l);
        ans-=(r-l+1)*djsmu(x/l);
    }
    return w[x]=ans;
}
long long djsphi(long long x)
{
    if(x<=6000000)return sum2[x];
    if(w1[x])return w1[x];
    long long ans=x*(x+1)/2;
    for(long long l=2,r;l<=x;l=r+1)
    {
        r=x/(x/l);
        ans-=(r-l+1)*djsphi(x/l);
    }
    return w1[x]=ans;
}
int main()
{
    int t,n;
    read(t);
    get(6000000);
    while(t--)
    {
        read(n);
        printf("%lld %d\n",djsphi(n),djsmu(n));
    }
    return 0;
}
```



### 凸包和旋转卡壳

```c++
// point Plus
int sgn(const double &x){ return x < -eps? -1 : (x > eps);}
inline double sqr(const double &x){ return x * x;}
struct Point
{
    double x, y;
    Point(const double &x = 0, const double &y = 0):x(x), y(y){}
    Point operator -(const Point &a)const{ return Point(x - a.x, y - a.y);}
    Point operator +(const Point &a)const{ return Point(x + a.x, y + a.y);}
    Point operator *(const double &a)const{ return Point(x * a, y * a);}
    Point operator /(const double &a)const{ return Point(x / a, y / a);}
    bool operator < (const Point &a)const{ return sgn(x - a.x) < 0 || (sgn(x - a.x) == 0 && sgn(y - a.y) < 0);}
    bool operator == (const Point &a)const{ return sgn(sgn(x - a.x) == 0 && sgn(y - a.y) == 0);}
    friend double det(const Point &a, const Point &b){ return a.x * b.y - a.y * b.x;}
    friend double dot(const Point &a, const Point &b){ return a.x * b.x + a.y * b.y;}
    friend double dist(const Point &a, const Point &b){ return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));}
    void in(){ scanf("%lf %lf", &x, &y); }
    void out()const{ printf("%lf %lf\n", x, y); }
};

struct Poly //多边形类
{
    vector<Point>p; //顺时针凸包
    vector<Point>tb;// 逆时针凸包
    void in(const int &r)
    {
        p.resize(r);  //不早凸包的时候可以把p改为a
        for(int i = 0; i < r; i++) p[i].in();
    }
    //判断点集是否为凸包(返回m-1==n),或者用凸包点算出凸包顶点tb(本题即是)
    void isCanHull()
    {
        sort(p.begin(), p.end());
        p.erase(unique(p.begin(), p.end()), p.end());
        int n = p.size();
        tb.resize(n * 2 + 5);
        int m = 0;
        for(int i = 0; i < n; i++)
        {
            while(m > 1 && sgn(det(tb[m - 1] - tb[m - 2], p[i] - tb[m - 2])) <= 0)m--;
            tb[m++] = p[i];
        }
        int k = m;
        for(int i = n - 2; i >= 0; i--)
        {
            while(m > k && sgn(det(tb[m - 1] - tb[m -2], p[i] - tb[m - 2])) <= 0)m--;
            tb[m++] = p[i];
        }
        tb.resize(m);
        if(m > 1)tb.resize(m - 1);
        //for(int i = 0; i < m - 1; i++) tb[i].out();
    }  
    //旋转卡壳算法：输入凸包，输入最远两个点的坐标及最远距离
    int maxdist()//int &first,int &second,若要输出坐标可放入
    {
        int n=tb.size(),first,second;
        int Max=-INF;
        if(n==1) return Max;//first=second=0;
        for(int i=0,j=1;i<n;i++)
        {
            while(sgn(det(tb[(i+1)%n]-tb[i],tb[j]-tb[i])-det(tb[(i+1)%n]-tb[i],tb[(j+1)%n]-tb[i]))<0)
                j=(j+1)%n;
            int d=dist(tb[i],tb[j]);
            if(d>Max) Max=d;//first=i,second=j;
            d=dist(tb[(i+1)%n],tb[(j+1)%n]);
            if(d>Max) Max=d;//first=i,second=j;
        }
        return Max;
    }
}poly;
```



### 二分图最大匹配——匈牙利算法

```c++
/**********************************编号从一开始，注意二分图对称结果要/2*****************/
#define maxm                                                               //边的个数
#define maxn                                                               //点的个数

struct Edge
{
    int to,next;
};

struct Hun
{
	Edge edge[maxm];
	int head[maxn],tot;
	int uN,N;
	void init()
	{
	    tot =0;
	    memset(head,-1,sizeof(head));
	    return ;
	}

	void addedge(int u,int v)
	{
	    edge[tot].to = v;
	    edge[tot].next = head[u];
	    head[u]=tot++;
	    return ;
	}

	int linker[maxn];
	bool used[maxn];


	bool dfs(int u)
	{
	    for(int i=head[u];i!=-1;i=edge[i].next)
	    {
	        int v =  edge[i].to;
	        if(!used[v])
	        {
	            used[v]=true;
	            if(linker[v]==-1||dfs(linker[v]))
	            {
	                linker[v]=u;
	                return true;
	            }
	        }
	    }
	    return false;
	}

	int hungary()
	{
	    int res = 0;
	    memset(linker,-1,sizeof(linker));
	    for(int u = 1;u<=uN;u++)                                          //编号从1开始
	    {
	        memset(used,false,sizeof(used));
	        if(dfs(u))
	        {
	            res++;
	        }
	    }
	    return res;                                                         //对称时要/2
	}
} hun;
```



### 二分完美匹配——KM

```c++

struct KM
{
	int linker[MAX],lx[MAX],ly[MAX];
	bool visx[MAX],visy[MAX];
	int map[MAX][MAX];
	int n; 
	bool hungary(int u)
	{
			visx[u]=true;
		for(int i=0;i<n;i++)
			if(!visy[i]&&lx[u]+ly[i]==map[u][i])
			{//满足最大匹配 
				visy[i]=true;
				if(linker[i]==-1||hungary(linker[i])){
					linker[i]=u;
					return true;
				}	
			}
		return false;
	}
	void K_M()
	{
		int temp;
		memset(lx,0,sizeof(lx));
		memset(ly,0,sizeof(ly));
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				lx[i]=max(lx[i],map[i][j]);
			}
		}
		for(int i=0;i<n;i++)
		{
			while(1)
			{
				memset(visx,0,sizeof(visx));
				memset(visy,0,sizeof(visy));
				if(hungary(i))
					break;
				else{
					temp=INT_MAX;
						for(int j=0;j<n;j++) 
						if(visx[j])for(int k=0;k<n;k++){//x在交叉树内
							if(!visy[k]&&temp>lx[j]+ly[k]-map[j][k]){//y在外 
								temp=lx[j]+ly[k]-map[j][k];
							}
						}
						for(int j=0;j<n;j++){
							if(visx[j])
							lx[j]-=temp;
							if(visy[j])
							ly[j]+=temp;
						}
					}
			}
		}
	}
	void solve()
	{
		while(scanf("%d",&n)!=EOF)
		{
			memset(linker,-1,sizeof(linker));
			for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
			scanf("%d",&map[i][j]);
			K_M();
			int ans=0;
			for(int i=0;i<n;i++)
				ans+=map[linker[i]][i];
			printf("%d\n",ans);
		}
	}
}km;
```



### 费马小定理求逆元

```c++
//O(logp)
ll qpow(ll a,ll b,ll p)
{
	a%=p;
	ll ans=1;
	for (;b;b>>=1)
	{
		if (b&1)
		{
			ans*=a;
			ans%=p;
		}
		a*=a;	
		a%=p;
	}
	return ans;
}

ll inv(ll a,ll p)
{
	return qpow(a,p-2);
 } 
```



### 费用流

```c++
struct edge {
	int to, capacity, cost, rev;
	edge() {}
	edge(int to, int _capacity, int _cost, int _rev) :to(to), capacity(_capacity), cost(_cost), rev(_rev) {}
};
struct Min_Cost_Max_Flow {
	int V, H[maxn + 5], dis[maxn + 5], PreV[maxn + 5], PreE[maxn + 5];
	vector<edge> G[maxn + 5];
	//调用前初始化
	void Init(int n) {
		V = n;
		for (int i = 0; i <= V; ++i)G[i].clear();
	}
	//加边
	void Add_Edge(int from, int to, int cap, int cost) {
		G[from].push_back(edge(to, cap, cost, G[to].size()));
		G[to].push_back(edge(from, 0, -cost, G[from].size() - 1));
	}
	//flow是自己传进去的变量，就是最后的最大流，返回的是最小费用
	int Min_cost_max_flow(int s, int t, int f, int& flow) {
		int res = 0; fill(H, H + 1 + V, 0);
		while (f) {
			priority_queue <pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>> > q;
			fill(dis, dis + 1 + V, INF);
			dis[s] = 0; q.push(pair<int, int>(0, s));
			while (!q.empty()) {
				pair<int, int> now = q.top(); q.pop();
				int v = now.second;
				if (dis[v] < now.first)continue;
				for (int i = 0; i < G[v].size(); ++i) {
					edge& e = G[v][i];
					if (e.capacity > 0 && dis[e.to] > dis[v] + e.cost + H[v] - H[e.to]) {
						dis[e.to] = dis[v] + e.cost + H[v] - H[e.to];
						PreV[e.to] = v;
						PreE[e.to] = i;
						q.push(pair<int, int>(dis[e.to], e.to));
					}
				}
			}
			if (dis[t] == INF)break;
			for (int i = 0; i <= V; ++i)H[i] += dis[i];
			int d = f;
			for (int v = t; v != s; v = PreV[v])d = min(d, G[PreV[v]][PreE[v]].capacity);
			f -= d; flow += d; res += d*H[t];
			for (int v = t; v != s; v = PreV[v]) {
				edge& e = G[PreV[v]][PreE[v]];
				e.capacity -= d;
				G[v][e.rev].capacity += d;
			}
		}
		return res;
	}
	int Max_cost_max_flow(int s, int t, int f, int& flow) {
		int res = 0;
		fill(H, H + 1 + V, 0);
		while (f) {
			priority_queue <pair<int, int>> q;
			fill(dis, dis + 1 + V, -INF);
			dis[s] = 0;
			q.push(pair<int, int>(0, s));
			while (!q.empty()) {
				pair<int, int> now = q.top(); q.pop();
				int v = now.second;
				if (dis[v] > now.first)continue;
				for (int i = 0; i < G[v].size(); ++i) {
					edge& e = G[v][i];
					if (e.capacity > 0 && dis[e.to] < dis[v] + e.cost + H[v] - H[e.to]) {
						dis[e.to] = dis[v] + e.cost + H[v] - H[e.to];
						PreV[e.to] = v;
						PreE[e.to] = i;
						q.push(pair<int, int>(dis[e.to], e.to));
					}
				}
			}
			if (dis[t] == -INF)break;
			for (int i = 0; i <= V; ++i)H[i] += dis[i];
			int d = f;
			for (int v = t; v != s; v = PreV[v])d = min(d, G[PreV[v]][PreE[v]].capacity);
			f -= d; flow += d;
			res += d*H[t];
			for (int v = t; v != s; v = PreV[v]) {
				edge& e = G[PreV[v]][PreE[v]];
				e.capacity -= d;
				G[v][e.rev].capacity += d;
			}
		}
		return res;
	}
};
```



### 高斯消元(模意义)

```c++
typedef long long ll;
const int maxn = 75;
ll abs[maxn][maxn];
const ll p = 1e9 + 7;
ll qpow(ll a, ll b) //  a^b%p
{
	a%=p;
	ll ans=1;
	for (;b;b>>=1)
	{
		if (b&1)
		{
			ans*=a;
			ans%=p;
		}
		a*=a;	
		a%=p;
	}
	return ans;
}

void solve(int _n)
{
    for( int i=1;i<=_n;i++)
	{
		if(!a[i][i])//主元不能为0
		{
			maxi=0;
			for(int j=i+1;j<=_n&&!maxi;j++)
				if(a[j][i]) maxi=j;
			if(!maxi) continue;//如果一整列都为0，不需要消元
			for(int j=i;j<=_n+1;j++)
				tmp=a[maxi][j],a[maxi][j]=a[i][j],a[i][j]=tmp;
		}
		for(int j=i+1;j<=_n;j++)
		{
			tmp=a[j][i];
			if(!tmp) continue;//已经为0，不需要消元
			for(int k=i;k<=_n+1;k++)
				a[j][k]=((a[j][k]*a[i][i]-a[i][k]*tmp)%p+p)%p;
		}
    }
    for(int i=_n;i;i--)
	{
		for(int j=i+1;j<=_n;j++)
			a[i][_n+1]=((a[i][_n+1]-ans[j]*a[i][j])%p+p)%p;
		ans[i]=a[i][_n+1]*qpow(a[i][i],p-2)%p;
	}
}
```



### 后缀数组

```c++
/*
*suffix array
*倍增算法 O(n*logn)
*待排序数组长度为 n, 放在 0 n-1 中，在最后面补一个 0
*da(str ,sa,rk,height, n , );//注意是 n;
*例如：
*n = 8;
* num[] = { 1, 1, 2, 1, 1, 1, 1, 2, $ }; 注意 num 最后一位为 0，其他
大于 0
*rk[] = 4, 6, 8, 1, 2, 3, 5, 7, 0 ;rk[0 n-1] 为有效值，rk[n]
必定为 0 无效值
*sa[] = 8, 3, 4, 5, 0, 6, 1, 7, 2 ;sa[1 n] 为有效值，sa[0] 必定为 n 是
无效值
*height[]= 0, 0, 3, 2, 3, 1, 2, 0, 1 ;height[2 n] 为有效值
kuangbin
ACM Template of kuangbin
 *
*/
#include <bits/stdc++.h>
using namespace std;
const int MAXN = 20010;
int t1[MAXN], t2[MAXN], c[MAXN]; //求 SA 数组需要的中间变量，不需要赋值
//待排序的字符串放在 s 数组中，从 s[0] 到 s[n-1], 长度为 n, 且最大值小于 m,
//除 s[n-1] 外的所有 s[i] 都大于 0，r[n-1]=0
//函数结束以后结果放在 sa 数组中
bool cmp(int *r, int a, int b, int l) {
	return r[a] == r[b] && r[a + l] == r[b + l];
}
void da(int str[], int sa[], int rk[], int height[], int n, int m) {
	str[n]=0;
	n++;
	int i, j, p, *x = t1, *y = t2;
//第一轮基数排序，如果 s 的最大值很大，可改为快速排序
	for (i = 0; i < m; i++)c[i] = 0;
	for (i = 0; i < n; i++)c[x[i] = str[i]]++;
	for (i = 1; i < m; i++)c[i] += c[i - 1];
	for (i = n - 1; i >= 0; i--)  sa[--c[x[i]]] = i;
	for (j = 1; j <= n; j <<= 1) {
		p = 0;
//直接利用 sa 数组排序第二关键字
		for (i = n - j; i < n; i++)y[p++] = i; //后面的 j 个数第二关键字为空的最小
		for (i = 0; i < n; i++)if (sa[i] >= j) y[p++] = sa[i] - j;
//这样数组 y 保存的就是按照第二关键字排序的结果
//基数排序第一关键字
		for (i = 0; i < m; i++)c[i] = 0;
		for (i = 0; i < n; i++)c[x[y[i]]]++;
		for (i = 1; i < m; i++)c[i] += c[i - 1];
		for (i = n - 1; i >= 0; i--)sa[--c[x[y[i]]]] = y[i];
//根据 sa 和 x 数组计算新的 x 数组
		swap(x, y);
		p = 1; x[sa[0]] = 0;
		for (i = 1; i < n; i++)
			x[sa[i]] = cmp(y, sa[i - 1], sa[i], j) ? p - 1 : p++;
		if (p >= n)break;
		m = p;//下次基数排序的最大值
	}
	int k = 0;
	n--;
	for (i = 0; i <= n; i++)rk[sa[i]] = i;
	for (i = 0; i < n; i++) {
		if (k) k--;
		j = sa[rk[i] - 1];
		while (str[i + k] == str[j + k])k++;
		height[rk[i]] = k;
	}
}
int rk[MAXN], height[MAXN];
int RMQ[MAXN];
int mm[MAXN];

int best[20][MAXN];
void initRMQ(int n) {
	mm[0] = -1;
	for (int i = 1; i <= n; i++)
		mm[i] = ((i & (i - 1)) == 0) ? mm[i - 1] + 1 : mm[i - 1];
	for (int i = 1; i <= n; i++)best[0][i] = i;
	for (int i = 1; i <= mm[n]; i++)
		for (int j = 1; j + (1 << i) - 1 <= n; j++) {
			int a = best[i - 1][j];
			int b = best[i - 1][j + (1 << (i - 1))];
			if (RMQ[a] < RMQ[b])best[i][j] = a;
			else best[i][j] = b;
		}
}
int askRMQ(int a, int b) {
	int t;
	t = mm[b - a + 1];
	b -= (1 << t) - 1;
	a = best[t][a]; b = best[t][b];
	return RMQ[a] < RMQ[b] ? a : b;
}
int lcp(int a, int b) {
	a = rk[a]; b = rk[b];
	if (a > b)swap(a, b);
	return height[askRMQ(a + 1, b)];
}
char str[MAXN];
int r[MAXN];
int sa[MAXN];
int main()
{
	while (scanf("%s", str) == 1) {
		int len = strlen(str);
		int n = 2 * len + 1;
		for (int i = 0; i < len; i++)r[i] = str[i];
		for (int i = 0; i < len; i++)r[len + 1 + i] = str[len - 1 - i];
		r[len] = 1;
		r[n] = 0;
		da(r, sa, rk, height, n, 128);
		for (int i = 1; i <= n; i++)RMQ[i] = height[i];
		initRMQ(n);
		int ans = 0, st;
		int tmp;
		for (int i = 0; i < len; i++) {
			tmp = lcp(i, n - i); //偶对称
			if (2 * tmp > ans) {
				ans = 2 * tmp;
				st = i - tmp;
			}
			tmp = lcp(i, n - i - 1); //奇数对称
			if (2 * tmp - 1 > ans) {
				ans = 2 * tmp - 1;
				st = i - tmp + 1;
			}
		}
		str[st + ans] = 0;
		printf("%s\n", str + st);
	}
	return 0;
}
```



### 后缀数组DC3

```c++
/*
* 后缀数组
* DC3 算法，复杂度 O(n)
* 所有的相关数组都要开三倍
*/

#define F(x) ((x)/3+((x)%3==1?0:tb))
#define G(x) ((x)<tb?(x)*3+1:((x)-tb)*3+2)

int wa[MAXN*3],wb[MAXN*3],wv[MAXN*3],wss[MAXN*3];
int c0(int *r,int a,int b){
    return r[a] == r[b] && r[a+1] == r[b+1] && r[a+2] == r[b+2];
}
int c12(int k,int *r,int a,int b){
    if(k == 2)
        return r[a] < r[b] || ( r[a] == r[b] && c12(1,r,a+1,b+1) );
    else return r[a] < r[b] || ( r[a] == r[b] && wv[a+1] < wv[b+1] );
}
void sort(int *r,int *a,int *b,int n,int m)
{
    int i;
    for(i=0;i<n;i++)
        wv[i]=r[a[i]];
    for(i=0;i<m;i++)
        wss[i]=0;
    for(i=0;i<n;i++)
        wss[wv[i]]++;
    for(i=1;i<m;i++)
        wss[i]+=wss[i-1];
    for(i=n-1;i>=0;i--)
        b[--wss[wv[i]]]=a[i];
}
void dc3(int *r,int *sa,int n,int m){
    
    int i,j,*rn = r+n;
    int *san=sa+n,ta=0,tb=(n+1)/3,tbc=0,p;
    r[n]=r[n+1]=0;
    for(i=0;i<n;i++)
        if(i%3!=0)
            wa[tbc++]=i;
    sort(r+2,wa,wb,tbc,m);
    sort(r+1,wb,wa,tbc,m);
    sort(r,wa,wb,tbc,m);
    for(p=1,rn[F(wb[0])]=0,i=1;i<tbc;i++)
        rn[F(wb[i])]=c0(r,wb[i-1],wb[i])?p-1:p++;
    if(p<tbc)dc3(rn,san,tbc,p);
    else for(i=0;i<tbc;i++) san[rn[i]]=i;
    for(i=0;i<tbc;i++) if(san[i]<tb)wb[ta++]=san[i]*3;
    if(n%3==1) wb[ta++]=n-1;
    sort(r,wb,wa,ta,m);
    for(i=0;i<tbc;i++)wv[wb[i]=G(san[i])] =i;
    for(i=0,j=0,p=0;i<ta&&j<tbc;p++)
        sa[p]=c12(wb[j]%3,r,wa[i],wb[j])?wa[i++]:wb[j++];
    for(;i<ta;p++)sa[p]=wa[i++];
    for(;j<tbc;p++)sa[p]=wb[j++];
}

void da(int str[],int sa[],int rank[],int height[],int n,int m){
    for(int i=n;i<n*3;i++)
        str[i]=0;
    dc3(str,sa,n+1,m);
    int i,j,k=0;
    for(i=0;i<=n;i++)rank[sa[i]]=i;
    for(i=0;i<n;i++){
        if(k)k--;
        j=sa[rank[i]-1];
        while(str[i+k]==str[j+k])k++;
        height[rank[i]]=k;
    }
}

```



### 矩阵快速幂

```c++

const int MAX = 4;
const long long mod = 1e9 + 7;
struct Matrix {
    long long mp[MAX][MAX];
    int n, m;

    Matrix(int _n, int _m)
    {
        n = _n, m = _m;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                mp[i][j] = 0;
    }

    Matrix operator+(const Matrix &b)const {
        Matrix tmp(n, m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                tmp.mp[i][j] = mp[i][j] + b.mp[i][j];
                tmp.mp[i][j] %= mod;
            }
        return tmp;
    }

    Matrix operator*(const Matrix &b)const {
        Matrix ret(n, b.m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                for (int k = 0; k < m; k++) {
                    ret.mp[i][j] += mp[i][k] * b.mp[k][j];
                    ret.mp[i][j] %= mod;
                }
        return ret;
    }

    Matrix operator^(long long k)const {
        Matrix ret(n, m), b(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                b.mp[i][j] = mp[i][j];
            ret.mp[i][i] = 1;
        }
        while (k) {
            if (k & 1)
                ret = ret * b;
            b = b * b;
            k >>= 1;
        }
        return ret;
    }
};
```



### 快速幂

```c++
#define ll long long 

ll qpow(ll a,ll b,ll p)   //  a^b%p   
{
	a%=p;
	ll ans=1;
	for (;b;b>>=1)
	{
		if (b&1)
		{
			ans*=a;
			ans%=p;
		}
		a*=a;	
		a%=p;
	}
	return ans;
}
```



### 模拟退火

```c++
#include<bits/stdc++.h>
using namespace std;
inline void read(int &x)
{
    x=0;
    static int p;p=1;
    static char c;c=getchar();
    while(!isdigit(c)){if(c=='-')p=-1;c=getchar();}
    while(isdigit(c)) {x=(x<<1)+(x<<3)+(c-48);c=getchar();}
    x*=p;
}
double a,b,c,d,e,f;
double ansx,ansy,ansz;
const double eps=1e-10;
double dis(double x,double y,double z)
{
    return sqrt(x*x+y*y+z*z);
}
double calc(double x,double y)
{
    double A=c;
    double B=d*y+e*x;
    double C=a*x*x+b*y*y+f*x*y-1.0;
    double delta=B*B-4.0*A*C;
    if(delta<0)return 210000000.0;
    double x1=(-B+sqrt(delta))/(2.0*A);
    double x2=(-B-sqrt(delta))/(2.0*A);
    if(dis(x,y,x1)>dis(x,y,x2))return x2;
    return x1;
}
void MNTH()
{
    for(int times=1;times<=1;times++)
    {
        double T=10000;
        while(T>eps)
        {
            double nowx=ansx+(rand()*2-RAND_MAX)*T;
            double nowy=ansy+(rand()*2-RAND_MAX)*T;
            double nowz=calc(nowx,nowy);
            if(nowz==210000000.0){T*=0.99;continue;}
            double delta=dis(nowx,nowy,nowz)-dis(ansx,ansy,ansz);
            if(delta<0)ansx=nowx,ansy=nowy,ansz=nowz;
            else if(exp(delta/T)*RAND_MAX<rand())ansx=nowx,ansy=nowy,ansz=nowz;
            T*=0.99;
        }
    }
    printf("%.6lf\n",dis(ansx,ansy,ansz));
}
int main()
{
    srand(19890604);
    while(scanf("%lf%lf%lf%lf%lf%lf",&a,&b,&c,&d,&e,&f)!=EOF)
    {
        ansx=0;
        ansy=0;
        ansz=sqrt(1.0/c);
        MNTH();
    }
    return 0;
}
```



### 欧拉函数——单个

```c++
// ŷ��������������,O(n^(1/2))
#include <bits/stdc++.h>
#define me0(x) memset(x,0,sizeof(x))
#define ll long long
#define ull unsigned long long
#define scan(x) scanf("%d",&x)
#define scanll(x) scanf("%lld",&x)
#define dscan(x,y) scanf("%d%d",&x,&y)
#define rep(x,be,en) for (x=be;x<=en;x++)
#define fr1(n) for (int i=1;i<=n;i++)
#define fr(n) for (int i=0;i<n;i++)
using namespace std;
ll euler(ll n)
{
	ll ans=n;
	for (int i=2;i*i<=n;i++)
	{
		if (n%i==0)
		{
			ans=ans/i*(i-1);
			while (n%i==0)
			 n/=i;
		}
	}
	if (n>1) ans=ans/n*(n-1);
	return ans;
 } 
const ll mod=1e9+7;
ll num[100];
int main()
{
	int i,j;
	int top=1;
	num[top]=mod;
	for (i=2;i<=  ;  i++)
	 num[i]=euler(num[i-1]);
    for (i=1;i<=100;i++) 
    cout<<num[i]<<" ";
}
```



### 欧拉函数——线性筛

```c++
const int maxn=;

int prime[maxn],cnt,euler[maxn];

void get_euler(int n)
{
	memset(euler,0,sizeof(euler));
	euler[1]=1;
    cnt=0;
    for (int i=2;i<=n;i++)
    {
    	if (!euler[i])
    	{
    		euler[i]=i-1;
    		prime[++cnt]=i;
		}
		for (int j=1;j<=cnt && i*prime[j]<=n;j++)
		{
			if (i%prime[j]) euler[prime[j]*i]=euler[i]*(prime[j]-1);
			else 
			{
				euler[prime[j]*i]=euler[i]*prime[j];
				break;
			}
		}
	}	  
}
```



### 树链剖分

```c++
#include <bits/stdc++.h>
using namespace std;
const int MOD=1e9+7;
/******************************树链剖分******************************/
/*  
 */
struct treenode
{
    int val;
    int father, treesize, index, depth, topfather, maxson;
};

struct segmenttreenode
{
    int l, r,size;
    int val, lazy;
};

struct Edge
{
    int u, v;
};

class TreeCut
{
#define maxn 100005
#define maxm maxn<<1
#define lson rt<<1
#define rson rt<<1|1
private:
    int TreeSize,EdgeSize,cnt;
    treenode t[maxn];
    int head[maxn],nxt[maxm];
    int to[maxm];
    int val[maxn];  //中间值
    segmenttreenode segt[maxn << 2];

public:
    void init(int _size=0)
    {
        memset(head, 0, sizeof(head));
        EdgeSize = cnt  = 0;
        TreeSize=_size;
    }

    void SetValues(vector<int> &vals)
    {
        int id=0;
        for (auto &tmp:vals)
        {
            t[++id].val=tmp;
        }
        if (!TreeSize)
            TreeSize=id;
    }

    void ReadValues()
    {
        for (int i = 1; i <= TreeSize;i++)
        {
            cin >> t[i].val;
        }
    }

    void AddEdge(int _u,int _v)
    {
        ++EdgeSize;
        to[EdgeSize] = _v;
        nxt[EdgeSize] = head[_u];
        head[_u] = EdgeSize;
        ++EdgeSize;
        to[EdgeSize] = _u;
        nxt[EdgeSize] = head[_v];
        head[_v] = EdgeSize;
    }

    int relationshipDFS(int now, int father, int depth)
    {
        t[now].index = 0;
        t[now].maxson = 0;
        t[now].depth = depth;
        t[now].father = father;
        t[now].treesize = 1;
        int maxsonsize = -1;
        for (int i = head[now]; i;i=nxt[i])
        {
            if (to[i]==father)
                continue;
            t[now].treesize += relationshipDFS(to[i], now, depth + 1);
            if (t[to[i]].treesize>maxsonsize)
            {
                t[now].maxson = to[i];
                maxsonsize = t[to[i]].treesize;
            }
        }
        return t[now].treesize;
    }

    void ReIndexDFS(int now,int topfather) //topfather该链最高的儿子链顶
    {
        t[now].index = ++cnt;
        val[cnt] = t[now].val; 
        t[now].topfather = topfather;
        if (!t[now].maxson)
            return;
        ReIndexDFS(t[now].maxson, topfather);
        for (int i = head[now]; i;i=nxt[i])
        {
            if (!t[to[i]].index)
                ReIndexDFS(to[i], to[i]);
        }
    }

    void UpdateNode(int rt)
    {
        segt[rt].val=(segt[lson].val+segt[rson].val)%MOD;
    }

    void SegmentTreeBuild(int now,int _l,int _r)
    {
        segt[now].l = _l;
        segt[now].r = _r;
        segt[now].size=_r-_l+1;
        if (_l == _r)
        {
            segt[now].val = val[_l];
            segt[now].lazy = 0;
            return;
        }
        int mid = (_l + _r) >> 1;
        SegmentTreeBuild(now << 1, _l, mid);
        SegmentTreeBuild(now << 1 | 1, mid + 1, _r);
        UpdateNode(now);
    }

    void PushDown(int rt)   //下传标记
    {
        if (!segt[rt].lazy)
            return;
        segt[lson].val=(segt[lson].val+ segt[rt].lazy*segt[lson].size)%MOD;
        segt[rson].val=(segt[rson].val+ segt[rt].lazy*segt[rson].size)%MOD;
        segt[lson].lazy=(segt[lson].lazy+segt[rt].lazy)%MOD;
        segt[rson].lazy=(segt[rson].lazy+segt[rt].lazy)%MOD;
        segt[rt].lazy=0;
    }

    void IntervalAdd(int rt,int l,int r,int v)
    {
        if (l<=segt[rt].l && segt[rt].r <= r)
        {
            segt[rt].val+=v*segt[rt].size;
            segt[rt].lazy+=v;
            return;
        }
        PushDown(rt);
        int mid=(segt[rt].l+segt[rt].r)>>1;
        if (l<=mid)
        {
            IntervalAdd(lson,l,r,v);
        }
        if (r>mid)
        {
            IntervalAdd(rson,l,r,v);
        }
        UpdateNode(rt);
    }

    int IntervalSum(int rt,int l,int r)
    {
        int ans=0;
        if (l<= segt[rt].l && segt[rt].r<=r)
            return segt[rt].val;
        PushDown(rt);
        int mid=(segt[rt].l+segt[rt].r)>>1;
        if (l<=mid)
            ans=(ans+IntervalSum(lson,l,r))%MOD;
        if (r>mid)
            ans=(ans+IntervalSum(rson,l,r))%MOD;
        return ans;
    }

    //x,y路径求和
    int TreeSum(int x,int y)   
    {
        int ans=0;
        while (t[x].topfather!=t[y].topfather)
        {
            if (t[t[x].topfather].depth < t[t[y].topfather].depth)
                swap(x,y);
            ans = (ans + IntervalSum(1, t[t[x].topfather].index, t[x].index)) % MOD;
            x=t[t[x].topfather].father;
        }
        if (t[x].depth>t[y].depth)
            swap(x,y);
        ans=(ans+IntervalSum(1,t[x].index,t[y].index))%MOD;
        return ans;
    }

    //x节点到y节点路径上的点权值+v
    void TreeAdd(int x,int y,int v)
    {
         while (t[x].topfather!=t[y].topfather)
        {
            if (t[t[x].topfather].depth < t[t[y].topfather].depth)
                swap(x,y);
            IntervalAdd(1,t[t[x].topfather].index,t[x].index,v);
            x=t[t[x].topfather].father;
        }
        if (t[x].depth>t[y].depth)
            swap(x,y);
         IntervalAdd(1,t[x].index,t[y].index,v);
    }

    //x的子树权值+v
    void SubTreeAdd(int x,int v)
    {
        IntervalAdd(1, t[x].index, t[x].index + t[x].treesize - 1, v);
    }

    //单点x 权值+v
    void NodeAdd(int x,int v)
    {
        IntervalAdd(1, t[x].index, t[x].index, v);
    }

    //询问x到根上权值和
    int QueryNodeToRoot(int x)
    {
        return (TreeSum(1,x)+MOD)%MOD;
    }

    //完成树链剖分
    void CutTree(int _root=1)
    {
        relationshipDFS(_root, 0, 1);
        ReIndexDFS(_root, _root);
        SegmentTreeBuild(1, 1, cnt);
    }

    int LCA(int x,int y)
    {
        while (t[x].topfather!=t[y].topfather)
        {
            if (t[t[x].topfather].depth<t[t[y].topfather].depth)
                swap(x,y);
            x=t[t[x].topfather].father;
        }
        if (t[x].depth>t[y].depth)
            swap(x,y);
        return x;
    }

    void debug()
    {
        cout << "vals" << endl;
        for (int i = 1; i <= cnt;i++)
            cout << val[i] << " ";
        cout << endl
             << "-------------------------------" << endl;
    }

#undef maxn
#undef maxm
#undef lson
#undef rson
};
```



### 素数筛——线性

```c++
const int MAXN=;

bool vis[MAXN];
int prime[MAXN],cnt;
void get_prime(int n)
{
	memset(vis,0,sizeof(vis));
	cnt=0;
	for (int i=2;i<=n;i++)
	{
	 if(!vis[i])
	 {
       prime[++cnt]=i;
       vis[i]=1;
      }
 	 for (int j=1;j<=cnt && (prime[j]*i<=n);j++)
	 	{
	 		vis[i*prime[j]]=1;
	 		if (i%prime[j]==0)
	 		 break;
		 }	
	 }
}
```



### 拓展欧几里得求逆元

```c++
typedef long long ll;
void extgcd(ll a, ll b, ll &d, ll &x, ll &y)
{
	if (!b)
	{
		d = a;
		x = 1;
		y = 0;
	}
	else
	{
		extgcd(b, a % b, d, y, x);
		y -= x * (a / b);
	}
}
ll inverse(ll a, ll n)
{
	ll d, x, y;
	extgcd(a, n, d, x, y);
	return d == 1 ? (x + n) % n : -1;
}
```



### 替罪羊树 

```c++
#include <vector>
using namespace std;

namespace Scapegoat_Tree
{
#define MAXN (100000 + 10)
const double alpha = 0.75;
struct Node
{
    Node *ch[2];
    int key, size, cover; // size为有效节点的数量，cover为节点总数量
    bool exist;           // 是否存在（即是否被删除）
    void PushUp(void)
    {
        size = ch[0]->size + ch[1]->size + (int)exist;
        cover = ch[0]->cover + ch[1]->cover + 1;
    }
    bool isBad(void)
    { // 判断是否需要重构
        return ((ch[0]->cover > cover * alpha + 5) ||
                (ch[1]->cover > cover * alpha + 5));
    }
};
struct STree
{
  protected:
    Node mem_poor[MAXN];      //内存池，直接分配好避免动态分配内存占用时间
    Node *tail, *root, *null; // 用null表示NULL的指针更方便，tail为内存分配指针，root为根
    Node *bc[MAXN];
    int bc_top; // 储存被删除的节点的内存地址，分配时可以再利用这些地址

    Node *NewNode(int key)
    {
        Node *p = bc_top ? bc[--bc_top] : tail++;
        p->ch[0] = p->ch[1] = null;
        p->size = p->cover = 1;
        p->exist = true;
        p->key = key;
        return p;
    }
    void Travel(Node *p, vector<Node *> &v)
    {
        if (p == null)
            return;
        Travel(p->ch[0], v);
        if (p->exist)
            v.push_back(p); // 构建序列
        else
            bc[bc_top++] = p; // 回收
        Travel(p->ch[1], v);
    }
    Node *Divide(vector<Node *> &v, int l, int r)
    {
        if (l >= r)
            return null;
        int mid = (l + r) >> 1;
        Node *p = v[mid];
        p->ch[0] = Divide(v, l, mid);
        p->ch[1] = Divide(v, mid + 1, r);
        p->PushUp(); // 自底向上维护，先维护子树
        return p;
    }
    void Rebuild(Node *&p)
    {
        static vector<Node *> v;
        v.clear();
        Travel(p, v);
        p = Divide(v, 0, v.size());
    }
    Node **Insert(Node *&p, int val)
    {
        if (p == null)
        {
            p = NewNode(val);
            return &null;
        }
        else
        {
            p->size++;
            p->cover++;

            // 返回值储存需要重构的位置，若子树也需要重构，本节点开始也需要重构，以本节点为根重构
            Node **res = Insert(p->ch[val >= p->key], val);
            if (p->isBad())
                res = &p;
            return res;
        }
    }
    void Erase(Node *p, int id)
    {
        p->size--;
        int offset = p->ch[0]->size + p->exist;
        if (p->exist && id == offset)
        {
            p->exist = false;
            return;
        }
        else
        {
            if (id <= offset)
                Erase(p->ch[0], id);
            else
                Erase(p->ch[1], id - offset);
        }
    }

  public:
    void Init(void)
    {
        tail = mem_poor;
        null = tail++;
        null->ch[0] = null->ch[1] = null;
        null->cover = null->size = null->key = 0;
        root = null;
        bc_top = 0;
    }
    STree(void) { Init(); }

    void Insert(int val)
    {
        Node **p = Insert(root, val);
        if (*p != null)
            Rebuild(*p);
    }
    int Rank(int val)
    {
        Node *now = root;
        int ans = 1;
        while (now != null)
        { // 非递归求排名
            if (now->key >= val)
                now = now->ch[0];
            else
            {
                ans += now->ch[0]->size + now->exist;
                now = now->ch[1];
            }
        }
        return ans;
    }
    int Kth(int k)
    {
        Node *now = root;
        while (now != null)
        { // 非递归求第K大
            if (now->ch[0]->size + 1 == k && now->exist)
                return now->key;
            else if (now->ch[0]->size >= k)
                now = now->ch[0];
            else
                k -= now->ch[0]->size + now->exist, now = now->ch[1];
        }
    }
    void Erase(int k)
    {
        Erase(root, Rank(k));
        if (root->size < alpha * root->cover)
            Rebuild(root);
    }
    void Erase_kth(int k)
    {
        Erase(root, k);
        if (root->size < alpha * root->cover)
            Rebuild(root);
    }
};
#undef MAXN

} // namespace Scapegoat_Tree
```



### 线段树

```c++

#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <queue>
#include <cstdio>
#include <stack>
#include <map>
#include <string>
#include <set>
#include <unordered_set>
#include <iomanip>
#include <bitset>
#pragma GCC optimize(2)
#define eps 1e-5
#define mod 1000000007
#define pi acos(-1)
#define MAXN 100005
#define ee 2.71828182845904523536
using namespace std;
int a[MAXN][30];
int tmp;
template<size_t N, typename T>
struct Segment_tree
{
    struct Node
    {
        int l[2],r[2];
        long long s[2];
        int d,ls,rs;
    } seg[N * 2];
    int tot=1;
    void maintain(int x)
    {
        int dl=seg[seg[x].ls].d,dr=seg[seg[x].rs].d;
        seg[x].d=seg[seg[x].ls].d^seg[seg[x].rs].d;
        seg[x].l[0]=seg[seg[x].ls].l[0]+seg[seg[x].rs].l[(dl==0) ? 0:1];
        seg[x].l[1]=seg[seg[x].ls].l[1]+seg[seg[x].rs].l[(dl==0) ? 1:0];
        seg[x].r[0]=seg[seg[x].rs].r[0]+seg[seg[x].ls].r[(dr==0) ? 0:1];
        seg[x].r[1]=seg[seg[x].rs].r[1]+seg[seg[x].ls].r[(dr==0) ? 1:0];
        seg[x].s[0]=seg[seg[x].ls].s[0]+seg[seg[x].rs].s[0]+seg[seg[x].ls].r[0]*seg[seg[x].rs].l[0]+seg[seg[x].ls].r[1]*seg[seg[x].rs].l[1];
        seg[x].s[1]=seg[seg[x].ls].s[1]+seg[seg[x].rs].s[1]+seg[seg[x].ls].r[0]*seg[seg[x].rs].l[1]+seg[seg[x].ls].r[1]*seg[seg[x].rs].l[0];
    }
    void init()
    {
        memset(seg,0,sizeof(seg));
        tot=1;
    }
    void build(int left, int right, int x, int idx)
    {
        if (left == right)
        {
            seg[x].d=a[left][idx];
            seg[x].l[a[left][idx]]=seg[x].r[a[left][idx]]=seg[x].s[a[left][idx]]=1;
            return;
        }
        int mid = (left + right) / 2;
        seg[x].ls=++tot;
        build(left, mid, tot, idx);
        seg[x].rs=++tot;
        build(mid + 1, right, tot, idx);
        maintain(x);
    }
    T query(int l, int r, int left, int right, int x)
    {
        if(left==l && right==r) return seg[x].s[1];
        int mid = (left + right) / 2;
        if (r <= mid) return query(l, r, left, mid, seg[x].ls);
        else if (l > mid) return query(l, r, mid + 1, right, seg[x].rs);
    }
};
```



### 线性基+线段树

```c++

# include <bits/stdc++.h>
# define lson l,mid,id<<1
# define rson mid+1,r,id<<1|1
using namespace std;
typedef long long LL;
const int maxn = 7e5+30;
struct node{
    int bit[31], len;
}a[3*maxn];
void fun(node &x, int val){
    if(x.len == 31) return;
    for(int i=30; ~i; --i){
        if(val&(1<<i)){
            if(x.bit[i] == 0){
                x.bit[i] = val;
                ++x.len;
                break;
            }
            else val ^= x.bit[i];
        }
    }
}
void update(int pos, int val, int l, int r, int id){
    fun(a[id], val);
    if(l == r) return;
    int mid = l+r>>1;
    if(pos <= mid) update(pos, val, lson);
    else update(pos, val, rson);
}
node Merge(node x, node y){
    if(x.len == 31) return x;
    if(y.len == 31) return y;
    node tmp = x;
    for(int i=30; ~i; --i)
        if(y.bit[i])
            fun(tmp, y.bit[i]);
    return tmp;
}
node query(int L, int R, int l, int r, int id){
    if(L<=l && R>=r) return a[id];
    int mid = l+r>>1;
    if(R <= mid) return query(L, R, lson);
    else if(L > mid) return query(L, R, rson);
    else return Merge(query(L, R, lson), query(L, R, rson));
}
int main()
{
      #ifdef _IRONHEAD_
            assert(freopen("/Users/ironhead/algorithm/in.in", "r", stdin));
            // assert(freopen("/Users/ironhead/algorithm/out.out", "w", stdout));
        #endif
    int T;
    int n, m, op, x, y;
    int temp;
    scanf("%d",&T);
    while (T--)
    {
        scanf("%d%d",&n,&m);
        int lastans=0;
        for (int i=1;i<=n;i++)
        {
            scanf("%d",&temp);
            update(i,temp,1,n+m,1);
        }
        int tot=n;
        for (int i=1;i<=m;i++)
        {
            scanf("%d",&op);
            if (op)
            {
                scanf("%d",&temp);
                temp^=lastans;
                update(++tot,temp,1,n+m,1);
            }
            else
            {
                scanf("%d%d",&x,&y);
                x^=lastans;
                y^=lastans;
                x=x%tot+1;
                y=y%tot+1;
                if (x>y) swap(x,y);
                node tmp=query(x,y,1,n+m,1);
                int ans = 0;
                for(int i=30; ~i; --i)
                    if((ans^tmp.bit[i]) > ans) ans ^= tmp.bit[i];
                printf("%d\n",ans);
                lastans=ans;
            }
        }
        memset(a,0,sizeof(a));
    }
}
```



### 线性筛求逆元表

```c++
const ll p=mod;
ll inv[p-1]; 

void solve()
{
   inv[1]=1;
   for (int i=2;i<=p-1;i++)
   {
   	ll a=p/i;
	   ll b=p%i;
	   inv[i]=(-a*inv[b])%p;
   }
   return;
}
```



### 匈牙利指派

```c++
#include <bits/stdc++.h>
using namespace std;
const double inf = 1e100;
const int ms = 105;
double u[ms], v[ms];
int p[ms], way[ms];
double minv[ms];
bool used[ms];
pair<vector<int>, double> solve(const vector<vector<double>> &matrix)
{
    int n = matrix.size();
    if (n == 0)
        return {vector<int>(), 0};
    for (int i = 1; i <= n; i++)
    {
        for (int i = 0; i <= n; i++)
            minv[i] = inf;
        memset(way, 0, (n + 1) * sizeof(int));
        for (int j = 0; j <= n; j++)
            used[j] = false;
        p[0] = i;
        int k0 = 0;
        do
        {
            used[k0] = true;
            int i0 = p[k0], k1;
            double delta = inf;
            for (int j = 1; j <= n; j++)
            {
                if (!used[j])
                {
                    double cur = matrix[i0 - 1][j - 1] - u[i0] - v[j];
                    if (cur < minv[j])
                    {
                        minv[j] = cur;
                        way[j] = k0;
                    }
                    if (minv[j] < delta)
                    {
                        delta = minv[j];
                        k1 = j;
                    }
                }
            }
            for (int j = 0; j <= n; j++)
            {
                if (used[j])
                {
                    u[p[j]] += delta;
                    v[j] -= delta;
                }
                else
                {
                    minv[j] -= delta;
                }
            }
            k0 = k1;
        } while (p[k0] != 0);
        do
        {
            int k1 = way[k0];
            p[k0] = p[k1];
            k0 = k1;
        } while (k0 != 0);
    }
    vector<int> ans(n, -1);
    for (int j = 1; j <= n; j++)
    {
        if (p[j] == 0)
            continue;
        ans[p[j] - 1] = j - 1;
    }
    return {ans, -v[0]};
}
/**** ans 是每个行取的列号 0-indexed ****/
```



### 圆与简单三角形面积交

```c++
	/**************************************

Function : Direct area of a circle and triangle
***********/
const double eps = 1e-8;            //浮点数精度控制

struct point                        //点或者向量结构
{
    double x,y;
    point(double _x=0.0,double _y=0.0)
        : x(_x),y(_y) {}
    point operator - (const point & v)
    {
        return point(x-v.x,y-v.y);
    }
    double sqrx()                    //向量的模
    {
        return sqrt(x*x+y*y);
    }
};
double xmult(point & p1,point & p2,point & p0)        //叉乘
{
    return (p1.x-p0.x)*(p2.y-p0.y)-(p1.y-p0.y)*(p2.x-p0.x);
}
double distancex(point & p1,point & p2)
{
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}
point intersection(point u1,point u2,point v1,point v2)        //两直线交点
{
    point ret=u1;
    double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
            /((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
    ret.x+=(u2.x-u1.x)*t;
    ret.y+=(u2.y-u1.y)*t;
    return ret;
}
void intersection_line_circle(point c,double r,point l1,point l2,point& p1,point& p2){
    point p=c;
    double t;
    p.x+=l1.y-l2.y;
    p.y+=l2.x-l1.x;
    p=intersection(p,c,l1,l2);
    t=sqrt(r*r-distancex(p,c)*distancex(p,c))/distancex(l1,l2);
    p1.x=p.x+(l2.x-l1.x)*t;
    p1.y=p.y+(l2.y-l1.y)*t;
    p2.x=p.x-(l2.x-l1.x)*t;
    p2.y=p.y-(l2.y-l1.y)*t;
}

point ptoseg(point p,point l1,point l2)            //点到线段的最近距离
{
    point t=p;
    t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
    if (xmult(l1,t,p)*xmult(l2,t,p)>eps)
    return distancex(p,l1)<distancex(p,l2)?l1:l2;
    return intersection(p,t,l1,l2);
}
double distp(point & a,point & b)
{
    return (a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y);
}
double Direct_Triangle_Circle_Area(point a,point b,point o,double r)
{
    double sign=1.0;
    a=a-o;
    b=b-o;
    o=point(0.0,0.0);
    if(fabs(xmult(a,b,o))<eps) return 0.0;
    if(distp(a,o)>distp(b,o))
    {
        swap(a,b);
        sign=-1.0;
    }
    if(distp(a,o)<r*r+eps)
    {
        if(distp(b,o)<r*r+eps) return xmult(a,b,o)/2.0*sign;
        point p1,p2;
        intersection_line_circle(o,r,a,b,p1,p2);
        if(distancex(p1,b)>distancex(p2,b)) swap(p1,p2);
        double ret1=fabs(xmult(a,p1,o));
        double ret2=acos( p1*b/p1.sqrx()/b.sqrx() )*r*r;
        double ret=(ret1+ret2)/2.0;
        if(xmult(a,b,o)<eps && sign>0.0 || xmult(a,b,o)>eps && sign<0.0) ret=-ret;
        return ret;
    }
    point ins=ptoseg(o,a,b);
    if(distp(o,ins)>r*r-eps)
    {
        double ret=acos( a*b/a.sqrx()/b.sqrx() )*r*r/2.0;
        if(xmult(a,b,o)<eps && sign>0.0 || xmult(a,b,o)>eps && sign<0.0) ret=-ret;
        return ret;
    }
    point p1,p2;
    intersection_line_circle(o,r,a,b,p1,p2);
    double cm=r/(distancex(o,a)-r);
    point m=point( (o.x+cm*a.x)/(1+cm) , (o.y+cm*a.y)/(1+cm) );
    double cn=r/(distancex(o,b)-r);
    point n=point( (o.x+cn*b.x)/(1+cn) , (o.y+cn*b.y)/(1+cn) );
    double ret1 = acos( m*n/m.sqrx()/n.sqrx() )*r*r;
    double ret2 = acos( p1*p2/p1.sqrx()/p2.sqrx() )*r*r-fabs(xmult(p1,p2,o));
    double ret=(ret1-ret2)/2.0;
    if(xmult(a,b,o)<eps && sign>0.0 || xmult(a,b,o)>eps && sign<0.0) ret=-ret;
    return ret;
}
```



### 支配树

```c++
/********************** DAG ************************/
struct Graph  //有一个虚拟点0,可以在树中DFS来获得答案
{
static const int MAXN = 100050;
static const int MAXM = 200100;
static const int DC =   18; //log2(MAXN)
    int to[MAXM],nxt[MAXM],head[MAXN],tot,st[MAXN],N,in[MAXN],toplist[MAXN],listnum,fa[MAXN][DC],depth[MAXN];
    vector<int> pre[MAXN],tr[MAXN];
    void init(int _N)
    {
        for (int i = 0; i <= _N;i++)
            head[i] = 0, pre[i].clear(),tr[i].clear();
        tot = 0;
        N = _N;
        depth[0]=0;
        for (int i=0;i<DC;i++)
            fa[0][i]=0;
    }

    void addedge(int u,int v)
    {
        to[++tot] = v;
        nxt[tot] = head[u];
        head[u] = tot;
        pre[v].push_back(u);
        ++in[v];
    }

    bool topsort()
    {
        int top = 1;
        listnum=0;
        for (int i = 1; i <= N;i++)
            if (!in[i])
            {
                addedge(0,i);
            }
        st[top]=0;
        while (top)
        {
            int u=st[top];
            --top;
            for (int ed=head[u];ed;ed=nxt[ed])
            {
                --in[to[ed]];
                if (!in[to[ed]])
                    st[++top]=toplist[++listnum]=to[ed];
            }
        }
        return listnum==N;
    }
    int LCA(int x,int y)
    {
        if (depth[x]<depth[y]) 
            swap(x,y);
        for (int i=DC-1;i>=0;i--)
            if (depth[fa[x][i]]>=depth[y])
                x=fa[x][i];
        if (x==y)
            return x;
        for (int i=DC-1;i>=0;i--)
            if (fa[x][i]!=fa[y][i])
                x=fa[x][i],y=fa[y][i];
        return fa[x][0];
    }

    inline void update(int father,int son,int distance=1)
    {
        tr[father].push_back(son);
        depth[son]=depth[father]+distance;
        fa[son][0]=father;
        for (int i=1;i<DC;i++) 
            fa[son][i]=fa[fa[son][i-1]][i-1];   
    }

    void build_dominator_tree()
    {
        for (int i=1;i<=listnum;i++)
        {
            int u=toplist[i];
            int fa=pre[u][0];
            for (int i=1;i<pre[u].size();i++)
                fa=LCA(fa,pre[u][i]);
            update(fa,u);
        }
    }

    void auto_build()
    {
        topsort();
        build_dominator_tree();
    }
}DT;


/************************ 一般有向图 **************************/
#include<bits/stdc++.h>
using namespace std;
#define RI register int
int read() {
	int q=0;char ch=' ';
	while(ch<'0'||ch>'9') ch=getchar();
	while(ch>='0'&&ch<='9') q=q*10+ch-'0',ch=getchar();
	return q;
}
typedef long long LL;
const int N=50005,M=100005;
int n,m,tim;
int dfn[N],repos[N],mi[N],fa[N],f[N],semi[N],idom[N],ans[N];
struct graph{
	int tot,h[N],ne[M],to[M];
	void clear() {tot=0;for(RI i=0;i<=n;++i) h[i]=0;}
	void add(int x,int y) {to[++tot]=y,ne[tot]=h[x],h[x]=tot;}
}g,rg,ng,tr;

void init() {
	tim=0;g.clear(),rg.clear(),ng.clear(),tr.clear();
	for(RI i=1;i<=n;++i)
		repos[i]=dfn[i]=idom[i]=fa[i]=ans[i]=0,mi[i]=semi[i]=f[i]=i;
}
void tarjan(int x) {
	dfn[x]=++tim,repos[tim]=x;
	for(RI i=g.h[x];i;i=g.ne[i])
		if(!dfn[g.to[i]]) fa[g.to[i]]=x,tarjan(g.to[i]);
}
int find(int x) {
	if(x==f[x]) return x;
	int tmp=f[x];f[x]=find(f[x]);
	if(dfn[semi[mi[tmp]]]<dfn[semi[mi[x]]]) mi[x]=mi[tmp];
	return f[x];
}
void dfs(int x,LL num) {
	ans[x]=num+x;
	for(RI i=tr.h[x];i;i=tr.ne[i]) dfs(tr.to[i],num+x);
}
void work() {
	for(RI i=n;i>=2;--i) {
		int x=repos[i],tmp=n;
		for(RI j=rg.h[x];j;j=rg.ne[j]) {
			if(!dfn[rg.to[j]]) continue;//此题数据有误
			if(dfn[rg.to[j]]<dfn[x]) tmp=min(tmp,dfn[rg.to[j]]);
			else find(rg.to[j]),tmp=min(tmp,dfn[semi[mi[rg.to[j]]]]);
		}
		semi[x]=repos[tmp],f[x]=fa[x],ng.add(semi[x],x);
		
		x=repos[i-1];
		for(RI j=ng.h[x];j;j=ng.ne[j]) {
			int y=ng.to[j];find(y);
			if(semi[mi[y]]==semi[y]) idom[y]=semi[y];
			else idom[y]=mi[y];//此时idom[mi[y]]可能并未找到
		}
	}
	for(RI i=2;i<=n;++i) {
		int x=repos[i];
		if(idom[x]!=semi[x]) idom[x]=idom[idom[x]];
		tr.add(idom[x],x);
	}
	dfs(n,0);
}
int main()
{
	int x,y;
	while(~scanf("%d%d",&n,&m)) {
		init();
		for(RI i=1;i<=m;++i)
			x=read(),y=read(),g.add(x,y),rg.add(y,x);
		tarjan(n);work();
		for(RI i=1;i<n;++i) printf("%d ",ans[i]);
		printf("%d\n",ans[n]);
	}
    return 0;
}
```



### 中国剩余定理

```c++
//a 为 余数， m为模数
#include <cstdio>
#include <algorithm>
#define maxn 100000
#define ll long long
using namespace std;
ll a[maxn], m[maxn], N;
void exgcd(ll a, ll b, ll &x, ll &y)
{
	if(!b){x=1,y=0;return;}
	ll xx, yy;
	exgcd(b,a%b,xx,yy);
	x=yy, y=xx-a/b*yy;
}
ll gcd(ll a, ll b){return !b?a:gcd(b,a%b);}
ll lcm(ll a, ll b){return a/gcd(a,b)*b;}
ll CRT(ll *a, ll *m, ll n)
{
	ll a1, m1, a2, m2, d, x1, x2;
	a1=a[1], m1=m[1];
	for(ll i=2;i<=n;i++)
	{
		a2=a[i], m2=m[i];
		d=gcd(m1,m2);
		if((a2-a1)%d!=0)return -1;
		exgcd(m1,m2,x1,x2);
		x1=(a2-a1)/d*x1%m2;
		a1=(a1+x1*m1)%lcm(m1,m2);
		m1=lcm(m1,m2);
	}
	return (a1+m1)%m1;
}
int main()
{
	ll i, N;
	while(~scanf("%lld",&N))
	{
		for(i=1;i<=N;i++)scanf("%lld%lld",m+i,a+i);
		printf("%lld\n",CRT(a,m,N));
	}
	return 0;
}
```



### 主席树

```c++
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
```



### 字典树(lrj)

```c++
//  main.cpp
//  字典树(高性能)

#include <cstring>
#define TYPES_OF_CHAR 26
#define MAX_WORDS 4000
#define MAX_WORD_LENGTH 101
struct Trie{
    int ch[MAX_WORDS*MAX_WORD_LENGTH][TYPES_OF_CHAR];
    int is_last[MAX_WORDS*MAX_WORD_LENGTH];
    int tree_size;
    Trie()
    {
        tree_size=1;
        memset(ch[0],0,sizeof(ch[0]));
        memset(is_last,0,sizeof(is_last));
    }
    void clear()
    {
        memset(ch[0],0,sizeof(ch[0]));
        memset(is_last,0,sizeof(is_last));
    }
    //⚠️需要按照情况更改
    int index_get(char c)
    {
        return c-'a';
    }
    
    void insert_word(char *a,int v)
    {
        int u=0,n=(int)strlen(a);
        for(int i=0;i<n;i++)
        {
            int c=index_get(a[i]);
            if(!ch[u][c]){
                memset(ch[tree_size],0,sizeof(ch[tree_size]));
                is_last[tree_size]=0;
                ch[u][c] = tree_size++;
            }
            u=ch[u][c];
        }
        is_last[u]=v;
    }
};
```



### 字典树(STL)

```c++
//  字典树（STL）

#include <string>
#include <map>
using namespace std;
struct _node{
    map<char,_node*> m;
    int last=0;
}root;
void add(string s)
{
    _node *p=&root;
    for(int i=0;i<s.length();i++)
    {
        if(!p->m.count(s[i]))
        {
            p->m[s[i]]=new _node;
        }
        p=p->m[s[i]];
    }
    p->last++;
}
```



### 组合数前缀和——莫队

```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<queue>
#include<vector>
#include<iomanip>

#define ll long long
#define scan(x) scanf("%d",&x)
#define dscan(x,y) scanf("%d%d",&x,&y)
using namespace std;
const int  MOD=1e9+7;
const int MAXN=1e5+100;
struct node
{
	int m,n,id;
	bool operator < (const node& b)
	{
	   	return n<b.n;
	}
} query[MAXN];
vector<node> left_area[MAXN];
ll lad_and[MAXN];
ll inv_and[MAXN];
int ans_to_query[MAXN];
int in_chunk[MAXN];
const int MAXX=1e5;
ll qpow(ll base, ll po)
{
	ll ans=1;
	for (;po;base*=base,base%=MOD,po>>=1)
	{
		if (po&1) 
		{
			ans*=base;
			ans%=MOD;
		}
	}
	return ans;
}

ll C(int n,int m)
{
	return (lad_and[n] * inv_and[n-m] %MOD *inv_and[m]) %MOD;
}


int main(int argc, char const *argv[])
{
	lad_and[0]=1;
	for (int i=1;i<=MAXX;i++)
	{
		lad_and[i]=lad_and[i-1]*(ll)i % MOD;
	}
	inv_and[MAXX]=qpow(lad_and[MAXX],MOD-2);
	for (int i=MAXX-1;~i;i--)
		inv_and[i]=inv_and[i+1]*(i+1)%MOD;
	int chunk=(int)sqrt((double)MAXX);
	int cnt=1;
	 for (int i = 1; i <= MAXX; i += chunk, ++ cnt)
        for (int j = i; j < i + chunk && j <= MAXX; ++ j)
            in_chunk[j] = cnt;
    cnt --;
   int T;
   scan(T);
   for (int i=1;i<=T;i++)
   {
   	  dscan(query[i].n,query[i].m);
   	  query[i].id=i;
   	  left_area[in_chunk[query[i].m]].push_back(query[i]);
   	  
   }
 //  sort(query+1,query+T+1);
   for (int i=1;i<=cnt;i++)
   	if (left_area[i].size())
   	{
   		sort(left_area[i].begin(),left_area[i].end());
   		int val = 0, tmpn = left_area[i][0].n;
   		int tmpm = -1;
   		for (int j=0;j<left_area[i].size();j++)
   		{
   			while (tmpn<left_area[i][j].n)
   			{
   				val=(MOD-C(tmpn,tmpm)+val+val)%MOD;
   				tmpn++;
   			}
   			while (tmpm<left_area[i][j].m)
   			{
   				tmpm++;
   				val+=C(tmpn,tmpm);
   				val%=MOD;
   			}
   			while (tmpm>left_area[i][j].m)
   			{
   				val=(val+MOD-C(tmpn,tmpm))%MOD;
   				tmpm--;
   			}
   			ans_to_query[left_area[i][j].id]=val;
   		}
   	}
   	for (int i=1;i<=T;i++)
   		cout<<ans_to_query[i]<<endl;
	return 0;
}
```



### 最大流——Dinic

```c++
#define MAX                                //点的个数
#define INF 0x3f3f3f3f

struct Edge
{
	int u,v,c,f;
} ;

int m,n,sum;
struct Dinic
{
	int s,t;
	vector<Edge> E;
	vector<int> G[MAX]; 
	bool vis[MAX];
	int lev[MAX];
	int cur[MAX];
	void init(int l,int r)
	{
		E.clear();
		for (int i=l;i<=r;i++)
			G[i].clear();
	}
	void addedge(int from,int to,int cap)
	{
		E.push_back((Edge){from,to,cap,0});
		E.push_back((Edge){to,from,0,0});
		int m=E.size();
		G[from].push_back(m-2);
		G[to].push_back(m-1);
	}
	bool bfs()
	{
		me0(vis);
		queue<int>  q;
		q.push(s);
		lev[s]=0;
		vis[s]=1;
		while (!q.empty())
		{
			int now=q.front();
			q.pop();
			for (int i=0,_size=G[now].size();i<_size;i++)
			{
				Edge edge=E[G[now][i]];
				int nex=edge.v;
				if (!vis[nex] && edge.c>edge.f)
				{
					lev[nex]=lev[now]+1;
					q.push(nex);
					vis[nex]=1;
				}
			}
		}
		return vis[t];
	}
	int dfs(int now,int aug)
	{
		if (now==t || aug==0) return aug;
		int flow=0,f;
		for (int &i=cur[now],_size=G[now].size();i<_size;i++)
		{
			Edge& edge=E[G[now][i]];
			int nex=edge.v;
			if (lev[now]+1 ==lev[nex] && (f=dfs(nex,min(aug,edge.c-edge.f)))>0)
			{
				edge.f+=f;
				E[G[now][i]^1].f-=f;
				flow+=f;
				aug-=f;
				if (!aug) break;
			}
		}
		return flow;
	}
	int maxflow(int s,int t)
	{
		this->s = s;
		this->t = t;
		int  flow=0;
		while (bfs())
		{
			me0(cur);
			flow+=dfs(s,INF);
		}
		return flow;
	}
}dinic;
```



### 最小费用最大流

```c++

#define MAXN                                      //点的个数
#define MAXM                                      //边的条数
#define INF 0x7fffffff
struct Edge
{
    int to, next, cap, flow, cost;
    int x, y;
};

struct Dinic_mincost
{
	Edge edge[MAXM];
	int head[MAXN],tol;
	int pre[MAXN],dis[MAXN];
	bool vis[MAXN];
	int N;
//	char map[MAXN][MAXN];
	void init()
	{
	    N = MAXN;
	    tol = 0;
	    memset(head, -1, sizeof(head));
	}
	void addedge(int u, int v, int cap, int cost)//左端点，右端点，容量，花费
	{
	    edge[tol]. to = v;
	    edge[tol]. cap = cap;
	    edge[tol]. cost = cost;
	    edge[tol]. flow = 0;
	    edge[tol]. next = head[u];
	    head[u] = tol++;
	    edge[tol]. to = u;
	    edge[tol]. cap = 0;
	    edge[tol]. cost = -cost;
	    edge[tol]. flow = 0;
	    edge[tol]. next = head[v];
	    head[v] = tol++;
	}
	bool spfa(int s, int t) //源点，汇点
	{
	    queue<int> q;
	    for(int i = 0; i < N; i++)
	    {
	        dis[i] = INF;
	        vis[i] = false;
	        pre[i] = -1;
	    }
	    dis[s] = 0;
	    vis[s] = true;
	    q.push(s);
	    while(!q.empty())
	    {
	        int u = q.front();
	        q.pop();
	        vis[u] = false;
	        for(int i = head[u]; i != -1; i = edge[i]. next)
	        {
	            int v = edge[i]. to;
	            if(edge[i]. cap > edge[i]. flow &&
	                    dis[v] > dis[u] + edge[i]. cost )
	            {
	                dis[v] = dis[u] + edge[i]. cost;
	                pre[v] = i;
	                if(!vis[v])
	                {
	                    vis[v] = true;
	                    q.push(v);
	                }
	            }
	        }
	    }
	    if(pre[t] == -1) return false;
	    else return true;
	}
	//返回的是最大流容量， cost存的是最小费用
	int minCostMaxflow(int s, int t, int &cost)
	{
	    int flow = 0;
	    cost = 0;
	    while(spfa(s,t))
	    {
	        int Min = INF;
	        for(int i = pre[t]; i != -1; i = pre[edge[i^1]. to])
	        {
	            if(Min > edge[i]. cap - edge[i]. flow)
	                Min = edge[i]. cap - edge[i]. flow;
	        }
	        for(int i = pre[t]; i != -1; i = pre[edge[i^1]. to])
	        {
	            edge[i]. flow += Min;
	            edge[i^1]. flow -= Min;
	            cost += edge[i]. cost * Min;
	        }
	        flow += Min;
	    }
	    return flow;
	}
} dinic_mincost;
```



### 最小覆盖圆

```c++
/*最小圆覆盖*/
/*给定n个点, 让求半径最小的圆将n个点全部包围,可以在圆上*/
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#define EPS 1e-8
using namespace std;
const int maxn = 550;
struct point{
    double x, y;
};
int sgn(double x)
{
    if (fabs(x) < EPS)
        return 0;
    return x < 0 ? -1 : 1;
}
double get_distance(const point a, const point b)//两点之间的距离
{
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}
point get_circle_center(const point a, const point b, const point c)//得到三角形外接圆的圆心
{
    point center;
    double a1 = b.x - a.x;
    double b1 = b.y - a.y;
    double c1 = (a1 * a1 + b1 * b1) / 2.0;
    double a2 = c.x - a.x;
    double b2 = c.y - a.y;
    double c2 = (a2 * a2 + b2 * b2) / 2.0;
    double d = a1 * b2 - a2 * b1;
    center.x = a.x + (c1 * b2 - c2 * b1) / d;
    center.y = a.y + (a1 * c2 - a2 * c1) / d;
    return center;
}
//p表示定点, n表示顶点的个数, c代表最小覆盖圆圆心, r是半径
void min_cover_circle(point *p, int n, point &c, double &r)//找最小覆盖圆(这里没有用全局变量p[], 因为是为了封装一个函数便于调用)
{
    random_shuffle(p, p + n);//随机函数,使用了之后使程序更快点,也可以不用
    c = p[0];
    r = 0;
    for (int i = 1; i < n; i++)
    {
        if (sgn(get_distance(p[i], c) - r) > 0)//如果p[i]在当前圆的外面, 那么以当前点为圆心开始找
        {
            c = p[i];//圆心为当前点
            r = 0;//这时候这个圆只包括他自己.所以半径为0
            for (int j = 0; j < i; j++)//找它之前的所有点
            {
                if (sgn(get_distance(p[j], c) - r) > 0)//如果之前的点有不满足的, 那么就是以这两点为直径的圆
                {
                    c.x = (p[i].x + p[j].x) / 2.0;
                    c.y = (p[i].y + p[j].y) / 2.0;
                    r = get_distance(p[j], c);
                    for (int k = 0; k < j; k++)
                    {
                        if (sgn(get_distance(p[k], c) - r) > 0)//找新作出来的圆之前的点是否还有不满足的, 如果不满足一定就是三个点都在圆上了
                        {
                            c = get_circle_center(p[i], p[j], p[k]);
                            r = get_distance(p[i], c);
                        }
                    }
                }
            }
        }
    }
}
int main()
{
    int n;
    point p[maxn];
    point c; double r;
    while (~scanf("%d", &n) && n)
    {
        for (int i = 0; i < n; i++)
            scanf("%lf %lf", &p[i].x, &p[i].y);
        min_cover_circle(p, n, c, r);
        printf("%.2lf %.2lf %.2lf\n", c.x, c.y, r);
    }
    return 0;
}
```



### AC自动机

```c++
#include<bits/stdc++.h>
#define pb push_back
using namespace std;
const int MAXN=5e5+5;
struct Trie
{
	int next[MAXN][26],fail[MAXN],end[MAXN];
	int root,L;
	int newnode()
	{
		for(int i=0;i<26;i++)
			next[L][i]=-1;
		end[L++]=0;
		return L-1;
	}
	void init()
	{
		L=0;
		root=newnode();
	}
	void insert(char buf[])
	{
		int len=strlen(buf);
		int now=root;
		for(int i=0;i<len;i++)
		{
			if(next[now][buf[i]-'a']==-1) 
				next[now][buf[i]-'a']=newnode();
			now=next[now][buf[i]-'a'];
		}
		end[now]++;
	}
	void build()
	{
		queue<int> Q;
		fail[root]=root;
		for(int i=0;i<26;i++)
			if(next[root][i]==-1)
				next[root][i]=root;
			else
			{
				fail[next[root][i]]=root;
				Q.push(next[root][i]);
			}
		while(!Q.empty())
		{
			int now=Q.front();
			Q.pop();
			for(int i=0;i<26;i++)
				if(next[now][i]==-1)
					next[now][i]=next[fail[now]][i];
				else
				{
					fail[next[now][i]]=next[fail[now]][i];
					Q.push(next[now][i]);
				}
			end[now]+=end[fail[now]];
		}
	}
}ac;
```



### 树状数组

```c++

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
```



### Dijkstra

```c++
//单源最短路Dijkstra
/*
 * 使用优先队列优化Dijkstra算法
 * 复杂度O(ElogE)
 * 注意对vector<Edge>E[MAXN]进行初始化后加边
 */
#include <iostream>
#include <cstring>
#include <queue>
using namespace std;
const int INF = 0x3f3f3f3f;
const int MAXN=1000010;
struct qnode
{
    int v;
    int c;
    qnode(int _v=0,int _c=0):v(_v),c(_c){}
    bool operator <(const qnode &r)const
    {
        return c>r.c;
    }
};
struct Edge
{
    int v,cost;
    Edge(int _v=0,int _cost=0):v(_v),cost(_cost){}
};
vector<Edge>E[MAXN];
bool vis[MAXN];
int dist[MAXN];
int pre[MAXN];
void showpath(int from,int to)
{
    vector<int> path;
    cout << "the total length is " << dist[to] << endl;
    cout << "Path" << endl;
    cout << "-------------------------" << endl;
    while (to != from)
    {
        path.push_back(to);
        to = pre[to];
    }
    path.push_back(from);
    for (int i = path.size() - 1; i >= 0;i--)
        if (i!=0)
            cout << path[i] << " -> ";
        else
            cout << path[i] << endl;
    cout << "--------------------------" << endl;
}
void Dijkstra(int n, int start) //点的编号从1开始
{
    memset(vis,false,sizeof(vis));
    for(int i=1;i<=n;i++)dist[i]=INF;
    priority_queue<qnode>que;
    while(!que.empty())que.pop();
    dist[start]=0;
    que.push(qnode(start,0));
    qnode tmp;
    while(!que.empty())
    {
        tmp=que.top();
        que.pop();
        int u=tmp.v;
        if(vis[u])continue;
        vis[u]=true;
        for(int i=0;i<E[u].size();i++)
        {
           int v=E[tmp.v][i].v;
           int cost=E[u][i].cost;
           if(!vis[v]&&dist[v]>dist[u]+cost)
           {
               dist[v]=dist[u]+cost;
               pre[v] = u;
               que.push(qnode(v, dist[v]));
           }
        }
    }
}
void addedge(int u,int v,int w)
{
    E[u].push_back(Edge(v,w));
}
```



### FastIO

```c++
namespace fastIO//输入外挂
{
	#define BUF_SIZE 100000
	//fread -> read
	bool IOerror = 0;
	inline char nc() {
		static char buf[BUF_SIZE], *p1 = buf + BUF_SIZE, *pend = buf + BUF_SIZE;
		if(p1 == pend) {
			p1 = buf;
			pend = buf + fread(buf, 1, BUF_SIZE, stdin);
			if(pend == p1) {
				IOerror = 1;
				return -1;
			}
		}
		return *p1++;
	}
	inline bool blank(char ch) {
		return ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t';
	}
	inline void read(int &x) {
		char ch;
		while(blank(ch = nc()));
		if(IOerror)
			return;
		for(x = ch - '0'; (ch = nc()) >= '0' && ch <= '9'; x = x * 10 + ch - '0');
	}
	#undef BUF_SIZE
};
using namespace fastIO;

void Out(int a)//输出外挂
{
    if(a>9)
        Out(a/10);
    putchar(a%10+'0');
}
```



### FFT

```c++

namespace myFFT
{
	#define maxn ((1<<17)+1)
	const double pi = acos(-1.0);
	struct cp
	{
		double a, b;
		cp operator +(const cp &o)const {return (cp) {a + o.a, b + o.b};}
		cp operator -(const cp &o)const {return (cp) {a - o.a, b - o.b};}
		cp operator *(const cp &o)const {return (cp) {a*o.a - b*o.b, b*o.a + a*o.b};}
		cp operator *(const double &o)const {return (cp) {a*o, b*o};}
		cp operator !() const {return (cp) {a, -b};}
	} w[maxn];
	int pos[maxn];
	void fft_init(int len)
	{
		int j = 0;
		while ((1 << j) < len)j++;
		j--;
		for (int i = 0; i < len; i++)
			pos[i] = pos[i >> 1] >> 1 | ((i & 1) << j);
	}
	void dft(cp *x, int len, int sta)
	{
		for (int i = 0; i < len; i++)
			if (i < pos[i])swap(x[i], x[pos[i]]);
		w[0] = (cp) {1, 0};
		for (unsigned i = 2; i <= len; i <<= 1)
		{
			cp g = (cp) {cos(2 * pi / i), sin(2 * pi / i)*sta};
			for (int j = i >> 1; j >= 0; j -= 2)w[j] = w[j >> 1];
			for (int j = 1; j<i >> 1; j += 2)w[j] = w[j - 1] * g;
			for (int j = 0; j < len; j += i)
			{
				cp *a = x + j, *b = a + (i >> 1);
				for (int l = 0; l<i >> 1; l++)
				{
					cp o = b[l] * w[l];
					b[l] = a[l] - o;
					a[l] = a[l] + o;
				}
			}
		}
		if (sta == -1)for (int i = 0; i < len; i++)x[i].a /= len, x[i].b /= len;
	}
	cp x[maxn], y[maxn], z[maxn];
	void FFT1(int *a, int *b, int n, int m, int *c)     //for int
	{
		int len = 1;
		while (len <= (n + m) >> 1)len <<= 1;
		fft_init(len);
		for (int i = n / 2; i < len; i++)x[i].a = x[i].b = 0;
		for (int i = m / 2; i < len; i++)y[i].a = y[i].b = 0;
		for (int i = 0; i < n; i++)(i & 1 ? x[i >> 1].b : x[i >> 1].a) = a[i];
		for (int i = 0; i < m; i++)(i & 1 ? y[i >> 1].b : y[i >> 1].a) = b[i];
		dft(x, len, 1), dft(y, len, 1);
		for (int i = 0; i < len / 2; i++)
		{
			int j = len - 1 & len - i;
			z[i] = x[i] * y[i] - (x[i] - !x[j]) * (y[i] - !y[j]) * (w[i] + (cp) {1, 0}) * 0.25;
		}
		for (int i = len / 2; i < len; i++)
		{
			int j = len - 1 & len - i;
			z[i] = x[i] * y[i] - (x[i] - !x[j]) * (y[i] - !y[j]) * ((cp) {1, 0} -w[i ^ len >> 1]) * 0.25;
		}
		dft(z, len, -1);
		for (int i = 0; i < n + m; i++)
			if (i & 1)c[i] = (int)(z[i >> 1].b + 0.5);
			else c[i] = (int)(z[i >> 1].a + 0.5);
	}
	void FFT2(int *a, int *b, int n, int m, ll *c)	//for long long
	{
		int len = 1;
		while (len <= (n + m) >> 1)len <<= 1;
		fft_init(len);
		for (int i = n / 2; i < len; i++)x[i].a = x[i].b = 0;
		for (int i = m / 2; i < len; i++)y[i].a = y[i].b = 0;
		for (int i = 0; i < n; i++)(i & 1 ? x[i >> 1].b : x[i >> 1].a) = a[i];
		for (int i = 0; i < m; i++)(i & 1 ? y[i >> 1].b : y[i >> 1].a) = b[i];
		dft(x, len, 1), dft(y, len, 1);
		for (int i = 0; i < len / 2; i++)
		{
			int j = len - 1 & len - i;
			z[i] = x[i] * y[i] - (x[i] - !x[j]) * (y[i] - !y[j]) * (w[i] + (cp) {1, 0}) * 0.25;
		}
		for (int i = len / 2; i < len; i++)
		{
			int j = len - 1 & len - i;
			z[i] = x[i] * y[i] - (x[i] - !x[j]) * (y[i] - !y[j]) * ((cp) {1, 0} -w[i ^ len >> 1]) * 0.25;
		}
		dft(z, len, -1);
		for (int i = 0; i < n + m; i++)
			if (i & 1)c[i] = (ll)(z[i >> 1].b + 0.5);
			else c[i] = (ll)(z[i >> 1].a + 0.5);
	}
	#undef maxn
} // myFFT

/**********************************************************DFT***************************************************/
const double pi = acos(-1.0);
struct Complex
{
 
    double R, I;
 
    inline Complex(double real = 0.0, double image = 0.0)
    {
        R = real, I = image;
    }
 
    inline Complex operator + (Complex const &b) const
    {
        return Complex(R + b.R, I + b.I);
    }
 
    inline Complex operator - (Complex const &b) const
    {
        return Complex(R - b.R, I - b.I);
    }
 
    inline Complex operator * (Complex const &b) const
    {
        return Complex(R * b.R - I * b.I, I * b.R + R * b.I);
    }
};
 
inline int turn(int n)
{
    int i = 1;
    for(; i < n; i <<= 1);
    return i;
}
void Change(Complex y[], int len)
{
    int i, j, k;
    for(i = 1, j = len / 2; i < len - 1; i++)
    {
        if(i < j)
            swap(y[i], y[j]);
        k = len / 2;
        while(j >= k)
        {
            j -= k;
            k /= 2;
        }
        if(j < k)
            j += k;
    }
}
void FFT(Complex P[], int n, int op)
{
    Change(P, n);
    for(int len = 2; len <= n; len <<= 1)
    {
        int m = len >> 1;
        Complex unit = Complex(cos(pi / m * op), sin(pi / m * op));
        for(int i = 0; i < n; i += len)
        {
            Complex W = Complex(1, 0);
            for(int j = 0; j < m; j++, W = W * unit)
            {
                Complex p = P[i + j], q = P[i + j + m];
                P[i + j] = p + W * q;
                P[i + j + m] = p - W * q;
            }
        }
    }
}
```



### fhq_treap

```c++
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
```



### FWT

```c++
#include<bits/stdc++.h>
using namespace std;

const int maxn = 112345;
#define LL long long 
const LL INFF = 0x3f3f3f3f3f3f3f3fll;

void fwt (ll a[] , int n ,bool on) {
    for ( int d = 1 ; d < n ; d <<= 1 ) {
        for ( int k = d << 1 , i = 0 ; i < n ; i += k ) {
            for ( int j = 0 ; j < d ; ++ j ) {
                LL x = a[i + j] , y = a[i + j + d] ;
                if(on){
                    a[i + j] = ( x + y )  ;
                    a[i + j + d] = ( x - y  )  ;
                }
                else{
                    a[i + j] = ( x + y ) / 2 ;
                    a[i + j + d] = ( x - y ) / 2;
                }
            }
        }
    }
}

char inp[22][maxn];
int arr[maxn];
LL A[1<<22],B[1<<22];

int getone(int x,int n){
    int ret = 0 ;
    while(x){
        ret += x & 1;
        x >>= 1;
    }
    return min(ret,n - ret);
}

int main(){
    int n,m;
    scanf("%d %d",&n,&m);
    for(int i = 0 ; i < n; i ++) scanf("%s",inp[i]);
    for(int i = 0 ; i < m ; i++){
        arr[i] = 0;
        for(int j = 0 ; j < n ; j ++){
            arr[i] |= (inp[j][i] - '0') << j;
        }
    }
    memset(A,0,sizeof(A));
    for(int i = 0 ; i < m ; i++) A[arr[i]] ++;
    for(int i = 0 ; i < (1<<n) ; i ++){
        B[i] = getone(i,n);
    }
    fwt(A,1<<n,true);
    fwt(B,1<<n,true);
    for(int i = 0; i < (1<<n);i++) A[i] *= B[i];
    fwt(A,1<<n,false);
    LL ans = INFF;
    for(int i = 0 ; i < (1<<n);i++) ans = min(ans,A[i]);
    printf("%I64d\n",ans);
    return 0;
}
```



### K_LIS最大覆盖

```c++
#include<cstdio>
#include<map>
#include<algorithm>
using namespace std;
typedef long long ll;
const int K=5;
int Case, n, i, x;
ll ans;
map<int, ll> T[K];
void ins(int o,int x,ll p){
  if(o>=K)return;
  T[o][x]+=p;
  ans+=p;
  while(p){
    map<int,ll>::iterator it=T[o].lower_bound(x+1);
    if(it==T[o].end())return;
    ll t=min(p,it->second);
    ans-=t;
    p-=t;
    ins(o+1,it->first,t);
    if(t==it->second)T[o].erase(it);else it->second-=t;
  }
}
int main(){
  scanf("%d",&Case);
  while(Case--){
    scanf("%d",&n);
    ans=0;
    for(i=0;i<K;i++)T[i].clear();
    for(i=1;i<=n;i++){
      scanf("%d",&x);
      ins(0,x,x);
      printf("%lld%c",ans,i<n?' ':'\n');
    }
  }
}
```



### KDTree

```c++
#include <bits/stdc++.h>
using namespace std;

typedef long long LL;

const LL maxn = 110000;
const LL INF = 1e18;
const LL dimension = 2;
LL D;

struct Node
{
    LL d[dimension], maxpos[dimension], minpos[dimension], v, sum, lazy, cnt;
    //以中心点d作为空间的代表，max和min分别是空间的边界
    LL l, r;
    bool operator<(const Node &b) const
    {
        return d[D] < b.d[D];
    }
    bool operator==(const Node &b) const
    {
        bool ans = 1;
        for (int i = 0; i < dimension; i++)
        {
            ans &= d[i] == b.d[i];
        }
        return ans;
    }
} p[maxn];

bool in(int x1, int y1, int x2, int y2, int X1, int Y1, int X2, int Y2)
{
    return x1 <= X1 && X2 <= x2 && y1 <= Y1 && Y2 <= y2;
}

bool out(int x1, int y1, int x2, int y2, int X1, int Y1, int X2, int Y2)
{
    return x1 > X2 || x2 < X1 || y1 > Y2 || y2 < Y1;
}

struct KDT
{
    LL root, cnt, block;
    Node tr[maxn], now;
    //now,用来单点插入
    void pushup(int rt)
    {
        tr[rt].cnt = 1;
        int l = tr[rt].l, r = tr[rt].r;
        if (l)
            tr[rt].cnt += tr[l].cnt;
        if (r)
            tr[rt].cnt += tr[r].cnt;
        for (int i = 0; i < dimension; i++)
        {
            tr[rt].maxpos[i] = tr[rt].minpos[i] = tr[rt].d[i];
            if (l)
            {
                tr[rt].minpos[i] = min(tr[rt].minpos[i], tr[l].minpos[i]);
                tr[rt].maxpos[i] = max(tr[rt].maxpos[i], tr[l].maxpos[i]);
            }
            if (r)
            {
                tr[rt].minpos[i] = min(tr[rt].minpos[i], tr[r].minpos[i]);
                tr[rt].maxpos[i] = max(tr[rt].maxpos[i], tr[r].maxpos[i]);
            }
        }
        tr[rt].sum = tr[l].sum + tr[r].sum + tr[rt].v;
    }

    void pushdown(int rt)
    {
        if (tr[rt].lazy)
        {
            int l = tr[rt].l, r = tr[rt].r;
            if (l)
            {
                tr[l].lazy += tr[rt].lazy;
                tr[l].v += tr[rt].lazy;
                tr[l].sum += tr[l].cnt * tr[rt].lazy;
            }
            if (r)
            {
                tr[r].lazy += tr[rt].lazy;
                tr[r].v += tr[rt].lazy;
                tr[r].sum += tr[r].cnt * tr[rt].lazy;
            }
            tr[rt].lazy = 0;
        }
    }

    LL rebuild(LL l, LL r, LL dep)
    { //重构
        if (l > r)
            return 0;
        D = dep;
        LL mid = (l + r) >> 1;
        nth_element(p + l, p + mid, p + r + 1);
        tr[mid] = p[mid];
        tr[mid].l = rebuild(l, mid - 1, (dep + 1) % dimension);
        tr[mid].r = rebuild(mid + 1, r, (dep + 1) % dimension);
        pushup(mid);
        return mid;
    }

    void checkSize()
    {
        if (cnt == block)
        {
            for (int i = 1; i <= cnt; i++)
                p[i] = tr[i];
            root = rebuild(1, cnt, 0);
            block += 10000;
        }
    }

    void ins(LL &rt, bool D)
    { //单点插入，如果没有就新开点
        if (!rt)
        {
            rt = ++cnt;
            for (int i = 0; i < dimension; i++)
                tr[rt].d[i] = tr[rt].maxpos[i] = tr[rt].minpos[i] = now.d[i];
            tr[rt].v = tr[rt].sum = now.v;
            return;
        }
        if (now == tr[rt])
        {
            tr[rt].v += now.v, tr[rt].sum += now.v;
            return;
        }
        pushdown(rt);
        if (now.d[D] < tr[rt].d[D])
            ins(tr[rt].l, D ^ 1);
        else
            ins(tr[rt].r, D ^ 1);
        pushup(rt);
    }
    
    LL query(int rt, int x1, int y1, int x2, int y2)
    {
        if (!rt)
            return 0;
        LL res = 0;
        if (out(x1, y1, x2, y2, tr[rt].minpos[0], tr[rt].minpos[1], tr[rt].maxpos[0], tr[rt].maxpos[1]))
            return 0;
        if (in(x1, y1, x2, y2, tr[rt].minpos[0], tr[rt].minpos[1], tr[rt].maxpos[0], tr[rt].maxpos[1]))
            return tr[rt].sum;
        pushdown(rt);
        if (in(x1, y1, x2, y2, tr[rt].d[0], tr[rt].d[1], tr[rt].d[0], tr[rt].d[1]))
            res += tr[rt].v;
        res += query(tr[rt].l, x1, y1, x2, y2) + query(tr[rt].r, x1, y1, x2, y2);
        pushup(rt);
        return res;
    }

    void update(int rt, int x1, int y1, int x2, int y2, LL add)
    {
        if (!rt)
            return;
        if (out(x1, y1, x2, y2, tr[rt].minpos[0], tr[rt].minpos[1], tr[rt].maxpos[0], tr[rt].maxpos[1]))
            return;
        if (in(x1, y1, x2, y2, tr[rt].minpos[0], tr[rt].minpos[1], tr[rt].maxpos[0], tr[rt].maxpos[1]))
        {
            tr[rt].lazy += add;
            tr[rt].sum += add * tr[rt].cnt;
            tr[rt].v += add;
            return;
        }
        pushdown(rt);
        if (in(x1, y1, x2, y2, tr[rt].d[0], tr[rt].d[1], tr[rt].d[0], tr[rt].d[1]))
        {
            tr[rt].v += add;
        }
        update(tr[rt].l, x1, y1, x2, y2, add);
        update(tr[rt].r, x1, y1, x2, y2, add);
        pushup(rt);
    }

    void init()
    {
        root = cnt = 0;
        block = 10000;
    }

    LL getdis(LL val[dimension], LL rt)
    { //估价函数，用来寻找区间
        LL res = 0;
        for (LL i = 0; i < dimension; i++)
        {
            if (val[i] < tr[rt].minpos[i])
                res += (tr[rt].minpos[i] - val[i]) * (tr[rt].minpos[i] - val[i]);
            if (val[i] > tr[rt].maxpos[i])
                res += (val[i] - tr[rt].maxpos[i]) * (val[i] - tr[rt].maxpos[i]);
        }
        return res;
    }

    LL ans;

    void ask(LL val[dimension], LL rt)
    { //询问最近点 欧式距离的平方
        LL dis = 0;
        for (LL i = 0; i < dimension; i++)
            dis += ((tr[rt].d[i] - val[i]) * (tr[rt].d[i] - val[i]));
        if (dis == 0)
            dis = INF;
        if (dis < ans)
            ans = dis;
        LL dl = tr[rt].l ? getdis(val, tr[rt].l) : INF;
        LL dr = tr[rt].r ? getdis(val, tr[rt].r) : INF;
        if (dl < dr)
        {
            if (dl < ans)
                ask(val, tr[rt].l);
            if (dr < ans)
                ask(val, tr[rt].r);
        }
        else
        {
            if (dr < ans)
                ask(val, tr[rt].r);
            if (dl < ans)
                ask(val, tr[rt].l);
        }
    }
} Tree;

int n, m;
int id;
int L[maxn], R[maxn];
int dep[maxn];
vector<int> G[maxn];

void dfs(int u, int fa)
{
    dep[u] = dep[fa] + 1;
    L[u] = ++id;
    for (auto v : G[u])
    {
        if (v == fa)
            continue;
        dfs(v, u);
    }
    R[u] = id;
}

int main()
{
    while (~scanf("%d%d", &n, &m))
    {
        Tree.init();
        id = 0;
        for (int i = 0; i < n - 1; i++)
        {
            int u, v;
            scanf("%d%d", &u, &v);
            G[u].push_back(v);
            G[v].push_back(u);
        }
        dep[0] = -1;
        dfs(1, 0);
        for (int i = 1; i <= n; i++)
        {
            Node &now = p[i];
            now.d[0] = L[i];
            now.d[1] = dep[i];
        }
        Tree.root = Tree.rebuild(1, n, 0);
        while (m--)
        {
            int tp, l, x;
            scanf("%d", &tp);
            if (tp == 1)
            {
                scanf("%d%d", &l, &x);
                Tree.update(Tree.root, 1, l, n, l, x);
            }
            else
            {
                scanf("%d", &x);
                printf("%lld\n", Tree.query(Tree.root, L[x], 1, R[x], n));
            }
        }
    }
    return 0;
}

/******************mine*************************/
#include <bits/stdc++.h>
using namespace std;
#define ll int
#define dimension 2
const int maxn = 5e5 + 10;
#define INF 1000000000
int D;
// int root;

struct Node
{
    ll d[dimension], l, r, mn[dimension], mx[dimension];
   
    ll& operator [](int n)
    {
        return d[n];
    }
     friend bool operator < (Node a, Node b) 
    {
        return a[D] < b[D];
    }
    bool operator==(const Node &b) const
    {
        bool tmp = 1;
        for (int i = 0; i < dimension; i++)
            tmp &= (d[i] == b.d[i]);
        return tmp;
    }
    friend ll dis(Node x, Node y) //曼哈顿距离
    {
        ll ret=0;
        for (int i=0;i<dimension;++i)
        {
            ret+=abs(x[i]-y[i]);
        }
        return ret;
    }
} p[maxn];

// namespace KDTree
// {
struct Tree
{
    int root;
    Node tr[maxn],T;
    ll ans;
    void pushup(int rt)
    {
        for (int i = 0; i < dimension; ++i)
            tr[rt].mn[i] = tr[rt].mx[i] = tr[rt].d[i];
        if (tr[rt].l)
        {
            for (int i = 0; i < dimension; i++)
            {
                tr[rt].mn[i]=min(tr[rt].mn[i], tr[tr[rt].l].mn[i]);
                tr[rt].mx[i]=max(tr[rt].mx[i], tr[tr[rt].l].mx[i]);
            }
        }
        if (tr[rt].r)
        {
            for (int i = 0; i < dimension; i++)
            {
                tr[rt].mn[i]=min(tr[rt].mn[i], tr[tr[rt].r].mn[i]);
                tr[rt].mx[i]=max(tr[rt].mx[i], tr[tr[rt].r].mx[i]);
            }
        }
    }

    int build(int l, int r, int dim)
    {
        int mid = (l + r) >> 1;
        D=dim;
        nth_element(p + l, p + mid, p + r + 1);
        tr[mid] = p[mid];
        if (l < mid)
            tr[mid].l = build(l, mid - 1, dim^1);
        if (mid < r)
            tr[mid].r = build(mid + 1, r, dim^1 );
        pushup(mid);
        return mid;
    }

    ll getmn(Node x) //估值
    {
        ll ans=0;
        for (int i=0;i<dimension;i++)
        {
            ans+=max(T[i] -x.mx[i],0);
            ans+=max(x.mn[i] - T[i],0);
        }
        return ans;
    }

    ll getmx(Node x) //估值 
    {
        ll ans = 0;
        for (int i = 0; i < dimension;i++)
        {
            ans += max(abs(T[i] - x.mx[i]), abs(T[i] - x.mn[i]));
        }
        return ans;
    }

    inline void querymaxdistance(int k)
    {
        ans=max(ans,dis(T,tr[k]));
        ll dl=-INF,dr=-INF;
        int l=tr[k].l,r=tr[k].r;
        if (l)
            dl=getmx(tr[l]);
        if (r)
            dr=getmx(tr[r]);
        if (dl>dr)
        {
            if (dl>ans)
                querymaxdistance(l);
            if (dr>ans) 
                querymaxdistance(r);
        }
        else
        {
            if (dr>ans)
                querymaxdistance(r);
            if (dl>ans)
                querymaxdistance(l);
        }
    }

    inline void querymindistance(int k)
    {
        if (dis(T,tr[k]))
            ans=min(ans,dis(T,tr[k]));
        ll dl=INF,dr=INF;
        int l=tr[k].l,r=tr[k].r;
        if (tr[k].l)
            dl=getmn(tr[l]);
        if (tr[k].r)
            dr=getmn(tr[r]);
        if (dl<dr)
        {
            if (dl<ans)
                querymindistance(l);
            if (dr<ans) 
                querymindistance(r);
        }
        else
        {
            if (dr<ans)
                querymindistance(r);
            if (dl<ans)
                querymindistance(l);
        }
    }

    ll query(int op,ll x,ll y)
    {
        T[0]=x,T[1]=y;
        if (op==0)
            ans=-INF,querymaxdistance(root);
        else
            ans=INF,querymindistance(root);
        return ans;
    }

}kd;
ll x[maxn],y[maxn];

int n;
int main()
{
    #ifdef _IRONHEAD_
        assert(freopen("/Users/ironhead/algorithm/in.in", "r", stdin));
       // assert(freopen("/Users/ironhead/algorithm/out.out", "w", stdout));
    #endif
    scanf("%d",&n);
    for (int i=1;i<=n;i++)
    {
        scanf("%d%d",&x[i],&y[i]);
        p[i][0]=x[i],p[i][1]=y[i];
    }
    kd.root=kd.build(1,n,0);
    ll ans=INF;
    int cnt=0;
    for (int i=1;i<=n;i++)
    {
        ll mx=kd.query(0,x[i],y[i]);
        ll mn=kd.query(1,x[i],y[i]);
        ans=min(ans,mx-mn);
    }
    printf("%d\n",ans);
    return 0;
}

/*可能有点问题的板子，支持插入和查询曼哈顿距离*/
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#define N 100005
#define inf (1<<30)
#define scan(x) scanf("%d",&x)
using namespace std;
const int k=3;
int n,m,dim,rt,ans,num=0;
struct node{int p[k],minn[k],maxx[k];}a[N];
bool cmp(node x,node y){ return x.p[dim]<y.p[dim]; }
struct kd_tree{
    int c[N][2];
    node s[N],q;
    void update(int o)
    {//管辖范围
        int l=c[o][0],r=c[o][1];
        for(int i=0;i<k;i++){
            if(l){ s[o].minn[i]=min(s[o].minn[i],s[l].minn[i]); s[o].maxx[i]=max(s[o].maxx[i],s[l].maxx[i]); }
            if(r){ s[o].minn[i]=min(s[o].minn[i],s[r].minn[i]); s[o].maxx[i]=max(s[o].maxx[i],s[r].maxx[i]); }
        }
    }
    void add(int o,node t){ for(int i=0;i<k;i++)s[o].minn[i]=s[o].maxx[i]=s[o].p[i]=t.p[i]; }
    int dist(node t,int o)
    {
        int tmp=0;
        for(int i=0;i<k;i++) tmp+=max(0,s[o].minn[i]-t.p[i]);
        for(int i=0;i<k;i++) tmp+=max(0,t.p[i]-s[o].maxx[i]);
        return tmp;
    }//?
    void build(int &o,int l,int r,int now)
    {
        o=(l+r)>>1; dim=now%k;
        nth_element(a+l,a+o,a+r+1,cmp);
        add(o,a[o]);
        if(l<o) build(c[o][0],l,o-1,now+1);
        if(o<r) build(c[o][1],o+1,r,now+1);
        update(o);
    }

    void ins(int o,int now){
        now%=k;
        if(q.p[now]<s[o].p[now]){
            if(c[o][0]) ins(c[o][0],now+1);
            else c[o][0]=++n,add(n,q);
        }
        else{
            if(c[o][1]) ins(c[o][1],now+1);
            else c[o][1]=++n,add(n,q);
        }
        update(o);
    }
    void qry(int o){//曼哈顿距离,且只求最短，dis是最短距离
        int tmp=0;
        for(int i=0;i<k;i++) tmp+=abs(s[o].p[i]-q.p[i]);
        ans=min(ans,tmp);
        int dl=c[o][0]?dist(q,c[o][0]):inf,dr=c[o][1]?dist(q,c[o][1]):inf;
        if(dl<dr)
        {
            if(dl<ans) qry(c[o][0]);
            if(dr<ans) qry(c[o][1]);
        }else{
            if(dr<ans) qry(c[o][1]);
            if(dl<ans) qry(c[o][0]);
        }
    }
}kd;
int tmp;
int main()
{
    // k=3;
    scan(tmp),scan(tmp),scan(tmp),scan(m);
    scan(tmp),scan(a[1].p[0]),scan(a[1].p[1]),scan(a[1].p[2]);
    ++num;
    kd.build(rt,1,1,0);
    --m;
    while(m--){
        scan(tmp);
        scanf("%d%d%d",&kd.q.p[0],&kd.q.p[1],&kd.q.p[2]);
        if(tmp==1)
        { 
            a[++num]=kd.q;
            if (num%5000==0)
                kd.build(rt,1,num,0);
            else
            {
                kd.ins(rt,0);
            }
            
        }
        else{
            ans=inf; kd.qry(rt); printf("%d\n",ans);
        }
    }
    return 0;
}
```



### KMP&exKMP

```c++
/*******************************EXKMP**************************************/

//计算next数组，保存到next参数中。
void getExtendNext (const char *t, int *next)
{
	int lt = strlen(t);
	for (int i = 1, j = -1, a, p; i < lt; i++, j--)
		if (j < 0 || i + next[i - a] >= p)
		{
			if (j < 0)	j = 0, p = i;
			while (p < lt && t[j] == t[p])	j++, p++;
			next[i] = j, a = i;
		}
		else	next[i] = next[i - a];
}

//扩展KMP算法求extend数组，s为文本串，t为模式串
void getExtend (const char *s, const char *t, int *extend)
{
	int ls = strlen(s), lt = strlen(t);
	int next[mx]； getExtendNext(t, next);	//计算next数组。mx是表示最大长度的常数。
	for (int i = 0, j = -1, a, p; i < ls; i++, j--)
		if (j < 0 || i + next[i - a] >= p)
		{
			if (j < 0)	j = 0, p = i;
			while (p < ls && j < lt && s[p] == t[j]) 	j++, p++;
			extend[i] = j, a = i;
		}
		else	extend[i] = next[i - a];
}
//********************************************************************//
//******************************KMP***********************************//
 void get_next(char *s,int *nxt)//数组为nxt,从0开始的序号
{
    int len=strlen(s);
	int i=0,k=-1;
	nxt[0]=-1;
	while (i<len)
	{ 
		if (k==-1||s[i]==s[k])
		{
			i++;
			k++;
			if (s[i]!=s[k]) nxt[i]=k;
			else	nxt[i]=nxt[k];
		}
		else 
		 k=nxt[k];
	}
}
	//*****************************************************************// 
```



### LCA

```c++
/*   tarjan 离线 O((n+q)) 实际上还有并查集饿的复杂度 */
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <cmath>
using namespace std;
const int N=40000+5;
struct Edge{
    int cnt,x[N],y[N],z[N],nxt[N],fst[N];
    void set(){
        cnt=0;
        memset(x,0,sizeof x);
        memset(y,0,sizeof y);
        memset(z,0,sizeof z);
        memset(nxt,0,sizeof nxt);
        memset(fst,0,sizeof fst);
    }
    void add(int a,int b,int c){
        x[++cnt]=a;
        y[cnt]=b;
        z[cnt]=c;
        nxt[cnt]=fst[a];
        fst[a]=cnt;
    }
}e,q;
int T,n,m,from,to,dist,in[N],rt,dis[N],fa[N],ans[N];
bool vis[N];
void dfs(int rt){
    for (int i=e.fst[rt];i;i=e.nxt[i]){
        dis[e.y[i]]=dis[rt]+e.z[i];
        dfs(e.y[i]);
    }
}
int getf(int k){
    return fa[k]==k?k:fa[k]=getf(fa[k]);
}
void LCA(int rt){
    for (int i=e.fst[rt];i;i=e.nxt[i]){
        LCA(e.y[i]);
        fa[getf(e.y[i])]=rt;
    }
    vis[rt]=1;
    for (int i=q.fst[rt];i;i=q.nxt[i])
        if (vis[q.y[i]]&&!ans[q.z[i]])
            ans[q.z[i]]=dis[q.y[i]]+dis[rt]-2*dis[getf(q.y[i])];
}
int main(){
    scanf("%d",&T);
    while (T--){
        q.set(),e.set();
        memset(in,0,sizeof in);
        memset(vis,0,sizeof vis);
        memset(ans,0,sizeof ans);
        scanf("%d%d",&n,&m);
        for (int i=1;i<n;i++)
            scanf("%d%d%d",&from,&to,&dist),e.add(from,to,dist),in[to]++;
        for (int i=1;i<=m;i++)
            scanf("%d%d",&from,&to),q.add(from,to,i),q.add(to,from,i);
        rt=0;
        for (int i=1;i<=n&&rt==0;i++)
            if (in[i]==0)
                rt=i;
        dis[rt]=0;
        dfs(rt);
        for (int i=1;i<=n;i++)
            fa[i]=i;
        LCA(rt);
        for (int i=1;i<=m;i++)
            printf("%d\n",ans[i]);
    }
    return 0;
}






/*  倍增 O((n+q)logn) 离线*/
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;
const int N=10000+5;
vector <int> son[N];
int T,n,depth[N],fa[N][20],in[N],a,b;
void dfs(int prev,int rt){
    depth[rt]=depth[prev]+1;
    fa[rt][0]=prev;
    for (int i=1;i<20;i++)
        fa[rt][i]=fa[fa[rt][i-1]][i-1];
    for (int i=0;i<son[rt].size();i++)
        dfs(rt,son[rt][i]);
}
int LCA(int x,int y){
    if (depth[x]<depth[y])
        swap(x,y);
    for (int i=19;i>=0;i--)
        if (depth[x]-(1<<i)>=depth[y])
            x=fa[x][i];
    if (x==y)
        return x;
    for (int i=19;i>=0;i--)
        if (fa[x][i]!=fa[y][i])
            x=fa[x][i],y=fa[y][i];
    return fa[x][0];
}
int main(){
    scanf("%d",&T);
    while (T--){
        scanf("%d",&n);
        for (int i=1;i<=n;i++)
            son[i].clear();
        memset(in,0,sizeof in);
        for (int i=1;i<n;i++){
            scanf("%d%d",&a,&b);
            son[a].push_back(b);
            in[b]++;
        }
        depth[0]=-1;
        int rt=0;
        for (int i=1;i<=n&&rt==0;i++)
            if (in[i]==0)
                rt=i;
        dfs(0,rt);
        scanf("%d%d",&a,&b);
        printf("%d\n",LCA(a,b));
    }
    return 0;
}




/* RMQ 构造O(nlogn) 查询O(1) */
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;
const int N=10000+5;
vector <int> son[N];
int T,n,depth[N],fa[N][20],in[N],a,b;
void dfs(int prev,int rt){
    depth[rt]=depth[prev]+1;
    fa[rt][0]=prev;
    for (int i=1;i<20;i++)
        fa[rt][i]=fa[fa[rt][i-1]][i-1];
    for (int i=0;i<son[rt].size();i++)
        dfs(rt,son[rt][i]);
}
int LCA(int x,int y){
    if (depth[x]<depth[y])
        swap(x,y);
    for (int i=19;i>=0;i--)
        if (depth[x]-(1<<i)>=depth[y])
            x=fa[x][i];
    if (x==y)
        return x;
    for (int i=19;i>=0;i--)
        if (fa[x][i]!=fa[y][i])
            x=fa[x][i],y=fa[y][i];
    return fa[x][0];
}
int main(){
    scanf("%d",&T);
    while (T--){
        scanf("%d",&n);
        for (int i=1;i<=n;i++)
            son[i].clear();
        memset(in,0,sizeof in);
        for (int i=1;i<n;i++){
            scanf("%d%d",&a,&b);
            son[a].push_back(b);
            in[b]++;
        }
        depth[0]=-1;
        int rt=0;
        for (int i=1;i<=n&&rt==0;i++)
            if (in[i]==0)
                rt=i;
        dfs(0,rt);
        scanf("%d%d",&a,&b);
        printf("%d\n",LCA(a,b));
    }
    return 0;
}
```



### LCT

```c++
#include<bits/stdc++.h>
#define N 300005
using namespace std;
int n,m,val[N];
struct Link_Cut_Tree{
    int top,c[N][2],fa[N],xr[N],q[N],rev[N];
    inline void pushup(int x){xr[x]=xr[c[x][0]]^xr[c[x][1]]^val[x];}
    inline void pushdown(int x){
        int l=c[x][0],r=c[x][1];
        if(rev[x]){
            rev[l]^=1;rev[r]^=1;rev[x]^=1;
            swap(c[x][0],c[x][1]);
        }
    }
    inline bool isroot(int x){return c[fa[x]][0]!=x&&c[fa[x]][1]!=x;}
    void rotate(int x){
        int y=fa[x],z=fa[y],l,r;
        if(c[y][0]==x)l=0;else l=1;r=l^1;
        if(!isroot(y)){if(c[z][0]==y)c[z][0]=x;else c[z][1]=x;}
        fa[x]=z;fa[y]=x;fa[c[x][r]]=y;
        c[y][l]=c[x][r];c[x][r]=y;
        pushup(y);pushup(x);
    }
    void splay(int x){
        top=1;q[top]=x;
        for(int i=x;!isroot(i);i=fa[i])q[++top]=fa[i];
        for(int i=top;i;i--)pushdown(q[i]);
        while(!isroot(x)){
            int y=fa[x],z=fa[y];
            if(!isroot(y)){
                if((c[y][0]==x)^(c[z][0]==y))rotate(x);
                else rotate(y);
            }rotate(x);
        }
    }
    void access(int x){for(int t=0;x;t=x,x=fa[x])splay(x),c[x][1]=t,pushup(x);}
    void makeroot(int x){access(x);splay(x);rev[x]^=1;}
    int find(int x){access(x);splay(x);while(c[x][0])x=c[x][0];return x;}
    void split(int x,int y){makeroot(x);access(y);splay(y);}
    void cut(int x,int y){split(x,y);if(c[y][0]==x)c[y][0]=0,fa[x]=0;}
    void link(int x,int y){makeroot(x);fa[x]=y;}
}T;
inline int read(){
    int f=1,x=0;char ch;
    do{ch=getchar();if(ch=='-')f=-1;}while(ch<'0'||ch>'9');
    do{x=x*10+ch-'0';ch=getchar();}while(ch>='0'&&ch<='9');
    return f*x;
}
int main(){
    n=read();m=read();
    for(int i=1;i<=n;i++)val[i]=read(),T.xr[i]=val[i];
    while(m--){
        int opt=read();
        if(opt==0){
            int x=read(),y=read();T.split(x,y);
            printf("%d\n",T.xr[y]);
        }
        if(opt==1){
            int x=read(),y=read(),xx=T.find(x),yy=T.find(y);
            if(xx!=yy)T.link(x,y);
        }
        if(opt==2){
            int x=read(),y=read(),xx=T.find(x),yy=T.find(y);
            if(xx==yy)T.cut(x,y);
        }
        if(opt==3){
            int x=read(),y=read();
            T.access(x);T.splay(x);val[x]=y;T.pushup(x);
        }
    }
    return 0;
}
```



### Manachar

```c++
//s+1原数组
//ss+1是添加了字符之后的数组
//r[]半径
//d[]差分后是以i开头回文串的个数
 
int Init()
{
    int len = strlen(s+1);
    ss[0] = '$';ss[1] = '#';
    int j = 2;
    for (int i = 1; i <= len; i++)ss[j++] = s[i],ss[j++] = '#';
    ss[j] = '\0';
    return j;
}

void Manacher()
{
    int len=Init();
    int p,mx=0;
    for (int i = 1; i < len; i++)
    {
        if (i<mx)
            r[i]=min(r[2*p-i],mx-i);
        else
                r[i] = 1;
        while (ss[i-r[i]]==ss[i+r[i]])
            r[i]++;
        if (mx<i+r[i])
        {
            p=i;
            mx=i+r[i];
        }
    }
    for (int i=2;i<len;i++)
    {
        if (ss[i]=='#' && r[i]==1) continue;
        int x=i/2-r[i]/2+1,y=i/2+r[i]/2-!(i&1);
        d[x]++;
        d[(x+y)/2+1]--;
    }
}
```



### mobius

```c++
#define MX 70 
int mu[MX];
int prime[MX];
bool vis[MX];
void getmu()
{
	memset(vis,0,sizeof(vis));
	int cnt = 0;
	mu[1] = 1;
	for (int i = 2; i < MX ; ++i)
	{
		if(!vis[i])
		{
			prime[cnt++] = i;
			mu[i] = -1;
		}
		for (int j = 0; j < cnt&&i*prime[j]<MX ; ++j)
		{
			vis[i*prime[j]] = true;
			if(i%prime[j]) mu[i*prime[j]] = -mu[i];
			else
			{
				mu[i*prime[j]] = 0;
				break;
			}
		}
	}
}
```



### NTT

```c++
const int XN = 2e6 + 11, P = 998244353;

int Pow(long long base, int v)
{
    long long res;
    for (res = 1; v; v >>= 1, (base *= base) %= P)
        if (v & 1)
            (res *= base) %= P;
    return res;
}

int D(int x)
{
    ((x >= P) && (x -= P)) || ((x < 0) && (x += P));
    return x;
}

void NTT(int a[], int n, int op)
{
    for (int i = 1, j = n >> 1; i < n - 1; ++i)
    {
        if (i < j)
            std::swap(a[i], a[j]);
        int k = n >> 1;
        while (k <= j)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    for (int len = 2; len <= n; len <<= 1)
    {
        int rt = Pow(3, (P - 1) / len);
        for (int i = 0; i < n; i += len)
        {
            int w = 1;
            for (int j = i; j < i + len / 2; ++j)
            {
                int u = a[j], t = 1LL * a[j + len / 2] * w % P;
                a[j] = D(u + t), a[j + len / 2] = D(u - t);
                w = 1LL * w * rt % P;
            }
        }
    }
    if (op == -1)
    {
        std::reverse(a + 1, a + n);
        int in = Pow(n, P - 2);
        for (int i = 0; i < n; ++i)
            a[i] = 1LL * a[i] * in % P;
    }
}

std::vector<int> Conv(std::vector<int> const &A, std::vector<int> const &B, int N)
{
    static int a[XN], b[XN];
    auto Make2 = [](int x) -> int {
        return 1 << ((32 - __builtin_clz(x)) + ((x & (-x)) != x));
    };
    int n = Make2(A.size() + B.size() - 1);
    for (int i = 0; i < n; ++i)
    {
        a[i] = i < A.size() ? A[i] : 0;
        b[i] = i < B.size() ? B[i] : 0;
    }
    NTT(a, n, 1);
    NTT(b, n, 1);
    for (int i = 0; i < n; ++i)
        a[i] = 1LL * a[i] * b[i] % P;
    NTT(a, n, -1);
    std::vector<int> C(N);
    for (int i = 0; i < N; i++)
        C[i] = a[i];
    return C;
}
```



### SPFA

```c++
#include <bits/stdc++.h>
/*
 * 单源最短路SPFA
 * 时间复杂度 0(kE)
 * 这个是队列实现，有时候改成栈实现会更加快，很容易修改
 * 这个复杂度是不定的
 */
const int MAXN=1010;
const int INF=0x3f3f3f3f;
struct Edge
{
    int v;
    int cost;
    Edge(int _v=0,int _cost=0):v(_v),cost(_cost){}
};
vector<Edge>E[MAXN];
void addedge(int u,int v,int w)
{
    E[u].push_back(Edge(v,w));
}
bool vis[MAXN];//在队列标志
int cnt[MAXN];//每个点的入队列次数
int dist[MAXN];
bool SPFA(int start,int n)
{
    memset(vis,false,sizeof(vis));
    for(int i=1;i<=n;i++)dist[i]=INF;
    vis[start]=true;
    dist[start]=0;
    queue<int>que;
    while(!que.empty())que.pop();
    que.push(start);
    memset(cnt,0,sizeof(cnt));
    cnt[start]=1;
    while(!que.empty())
    {
        int u=que.front();
        que.pop();
        vis[u]=false;
        for(int i=0;i<E[u].size();i++)
        {
           int v=E[u][i].v;
           if(dist[v]>dist[u]+E[u][i].cost)
           {
               dist[v]=dist[u]+E[u][i].cost;
               if(!vis[v])
               {
                   vis[v]=true;
                   que.push(v);
                   if(++cnt[v]>n)return false;
                   //cnt[i]为入队列次数，用来判定是否存在负环回路
               }
           }
        }
    }
    return true;
}
```



### Splay

```c++
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
```



### Suffix Automaton

```c++
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
```



### Tarjan_有向边联通分量

```c++
stack<int> s; //栈
int dfn[MAXN], low[MAXN], scc[MAXN];  //scc标明x节点属于哪个联通分量 dfn为dfs序 
bool ins[MAXN];   //标明x是不是在栈中  
vector<int> g[MAXN];
int num = 0, tot = 0;

int addedge(int u,int v)
{
    g[u].push_back(v);
}

int dfs(int x)
{
	if (dfn[x])
	{
		if (ins[x] == 0)
			return INF;
		return low[x];
	}
	
	low[x] = dfn[x] = ++num;   //为节点x设定次序编号和low初值
	ins[x] = 1;
	s.push(x);
	for (int i = 0; i < g[x].size(); i ++)
		low[x] = min(low[x], dfs(g[x][i]));
	if (low[x] == dfn[x])
	{
		while (dfn[s.top()] != low[s.top()])
		{
			ins[s.top()] = 0;
			scc[s.top()] = tot;
			s.pop();
		}
		ins[s.top()] = 0;
		scc[s.top()] = tot++;
		s.pop();
	}
	return low[x];
}

主函数内：
	for (int i = 1; i <= n; i ++)
		if (!dfn[i])
			dfs(i);

//算法复杂度O(N+M)



/*****非递归版本*****/
stack<int> ss;
int nx[Maxn]; //下一个儿子
int lx[Maxn]; //上一个儿子
void tarjan(int u){
    while(!ss.empty()) ss.pop();
    ss.push(u);
    for(int i=1;i<=n;i++) nx[i]=head[i],lx[i]=-1;
    while(!ss.empty()){
        int v=ss.top();
        if(!dfn[v]){ //第一次访问
            dfn[v]=low[v]=++tmpdfn;
            st[++top]=v;
            in[v]=1;
        }
        if(lx[v]!=-1) low[v]=min(low[v],low[lx[v]]); //访问完儿子
        if(nx[v]!=-1){ //有儿子
            while(nx[v]!=-1){ //寻找下一个还未访问的儿子
                if(!dfn[p[nx[v]].to]){ //树边
                    lx[v]=p[nx[v]].to;
                    nx[v]=p[nx[v]].next;
                    ss.push(lx[v]);
                    break;
                }
                else if(in[p[nx[v]].to]) //回边
                    low[v]=min(low[v],dfn[p[nx[v]].to]);
                nx[v]=p[nx[v]].next;
            }
        }
        else{ //全部儿子访问完毕
            if(low[v]==dfn[v]){
                scc++;
                do{
                    belong[st[top]]=scc;
                    in[st[top]]=0;
                    cnt[scc]+=val[st[top]];
                }while(st[top--]!=v);
            }
            ss.pop();
        }
    }
}
```



### Tarjan——割点割边点边双联通分量

```c++
注意如果有多组样例的初始化
/***************************割点********************************/
const int MAXN=;
int n,m,stamp=0,low[MAXN],dfn[MAXN],iscut[MAXN];
vector<int> vec[MAXN];
void tarjan(int index,int fa){
    int child=0;
    low[index]=dfn[index]=++stamp;
    for(int i=0;i<vec[index].size();i++)
    {
        int tmp=vec[index][i];
        if(!dfn[tmp])
        {
            child++;
            tarjan(tmp,index);
            low[index]=min(low[index],low[tmp]);
            if(low[tmp]>=dfn[index])
                iscut[index]=1;           //是割点
        }
        else if(dfn[tmp]<dfn[index] && tmp!=fa)
        {
            low[index]=min(low[index],dfn[tmp]);
        }
    }
    if(fa<0 && child==1)
        iscut[index]=0;             //不是割点
}

/***********************************************************/


/**************************割边******************************/
int n,stamp=0,dfn[MAXN],low[MAXN];
int cnt,ansu[10005],ansv[10005];  //割边的，
vector<int> vec[MAXN];


void addAns(int x,int y)
{
    if(x>y)
        swap(x,y);
    ansu[cnt]=x, ansu[cnt]=y;
    cnt++;
}

void tarjan(int index,int fa)
{
    int tmp;
    dfn[index]=low[index]=++stamp;
    for(int i=0;i<vec[index].size();i++)
    {
        tmp=vec[index][i];
        if(!dfn[tmp])
        {
            tarjan(tmp,index);
            low[index]=min(low[index],low[tmp]);
            if(low[tmp]>dfn[index])
                addAns(index,tmp);
        }

        else if(dfn[tmp]<dfn[index] && tmp!=fa)
        {
            low[index]=min(low[index],dfn[tmp]);
        }
    }
}

/**************************************************************/

/***************************点双联通分量**************************/
struct Edge{
    int u,v;
    Edge(int u=0,int v=0):u(u),v(v){}
}e[maxm];
int n,m,stamp=0,dfn[maxn],low[maxn],iscut[maxn],bccno[maxn];
int scnt,stack[maxm],bcc_cnt;
vector<int> vec[maxn],bcc[maxn];

void tarjan(int index,int fa)
{
    int child=0,tmp;
    dfn[index]=low[index]=++stamp;
    for(int i=0;i<vec[index].size();i++)
    {
        tmp=e[vec[index][i]].v;
        if(!dfn[tmp])
        {
            stack[++scnt]=vec[index][i],child++;
            tarjan(tmp,index);
            low[index]=min(low[index],low[tmp]);
            if(low[tmp]>=dfn[index])
            {
                iscut[index]=1;
                bcc[++bcc_cnt].clear();
                while(1)
                {
                    int num=stack[scnt--];
                    if(bccno[e[num].u]!=bcc_cnt)
                    {
                        bcc[bcc_cnt].push_back(e[num].u);
                        bccno[e[num].u]=bcc_cnt;
                    }
                    if(bccno[e[num].v]!=bcc_cnt)
                    {
                        bcc[bcc_cnt].push_back(e[num].v);
                        bccno[e[num].v]=bcc_cnt;
                    }
                    if(e[num].u==index && e[num].v==tmp)
                        break;
                }
            }
        }
        else if(dfn[tmp]<dfn[index] && tmp!=fa)
        {
            stack[++scnt]=vec[index][i];
            low[index]=min(low[index], dfn[tmp]);
        }
    }
    if(fa<0 && child==1)
        iscut[index]=0;
}

void find_bcc()
{
    // 割顶的bccno值毫无意义，因为它是属于多个连通分量中的点

    memset(dfn,0,sizeof(dfn));
    memset(low,0,sizeof(low));
    memset(iscut,0,sizeof(iscut));
    memset(bccno,0,sizeof(bccno));
    memset(bcc,0,sizeof(bcc));
    stamp=scnt=bcc_cnt=0;
    for(int i=1;i<=n;i++)
        if(!dfn[i])
            tarjan(i,-1);
}
/***********************************************************************/

/********************************边双连通分量************************************/
struct Edge{
    int u,v;
    Edge(int u=0,int v=0):u(u),v(v){}
}e[maxm];
int n,m,stamp,dfn[maxn],low[maxn],bccno[maxn],bcc_cnt;
vector<int> vec[maxn],bcc[maxn];
bool g[maxn][maxn],isbridge[maxm];

void tarjan(int index,int fa)
{
    int tmp;
    dfn[index]=low[index]=++stamp;
    for(int i=0;i<vec[index].size();i++)
    {
        tmp=e[vec[index][i]].v;
        if(!dfn[tmp])
        {
            tarjan(tmp,index);
            low[index]=min(low[index],low[tmp]);
            if(low[tmp]>dfn[index])
                isbridge[vec[index][i]]=isbridge[vec[index][i]^1]=1;
        }
        else if(dfn[tmp]<dfn[index] && tmp!=fa)
        {
            low[index]=min(low[index], dfn[tmp]);
        }
    }
}

void dfs(int index)
{
    dfn[index]=1;
    bccno[index]=bcc_cnt;
    for(int i=0;i<vec[index].size();i++)
    {
        int tmp=vec[index][i];
        if(isbridge[tmp])
            continue;
        if(!dfn[e[tmp].v])
        {
            dfs(e[tmp].v);
        }
    }
}

void find_ebcc(){
    bcc_cnt=stamp=0;
    memset(dfn,0,sizeof(dfn));
    memset(low,0,sizeof(low));
    memset(isbridge,0,sizeof(isbridge));
    memset(bccno,0,sizeof(bccno));
    memset(bcc,0,sizeof(bcc));
    for(int i=1;i<=n;i++)
        if(!dfn[i])
            tarjan(i, -1);
    memset(dfn,0,sizeof(dfn));
    for(int i=1;i<=n;i++)
    {
        if(!dfn[i])
        {
            bcc_cnt++;
            dfs(i);
        }
    }               
}

/*********************************************************************/
```

