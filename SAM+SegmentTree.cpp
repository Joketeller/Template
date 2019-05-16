#include<bits/stdc++.h>
using namespace std;

const int N=11e5+10,M=N*30;
int n,m,Q,ans,ps;
int rt[N],b[N],c[N],id[N],bz[N][21];
char s[N];

int read()
{
	int ret=0;char c=getchar();
	while(!isdigit(c)) c=getchar();
	while(isdigit(c)) ret=ret*10+(c^48),c=getchar();
	return ret;
}

struct SAM
{
	int sz,las,fa[N],mx[N],ch[N][26];
	int extend(int c)
	{
		int p,q,np,nq;
		if (ch[las][c])
	    {
	        p=las,q=ch[p][c];
	        if (mx[q]==mx[p]+1) return (las=q);
	        nq=++sz;mx[nq]=mx[p]+1;
	        memcpy(ch[nq],ch[q],sizeof(ch[q]));
	        fa[nq]=fa[q];fa[q]=nq;
	        while(p && ch[p][c]==q) ch[p][c]=nq,p=fa[p];
	        return (las=nq);
	    }
		p=las;las=np=++sz;mx[np]=mx[p]+1;
		while(p && !ch[p][c]) ch[p][c]=np,p=fa[p];
		if(!p) fa[np]=1;
		else
		{
			q=ch[p][c];
			if(mx[p]+1==mx[q]) fa[np]=q;
			else
			{
				nq=++sz;mx[nq]=mx[p]+1;
				memcpy(ch[nq],ch[q],sizeof(ch[q]));
				fa[nq]=fa[q];fa[np]=fa[q]=nq;
				while(p && ch[p][c]==q) ch[p][c]=nq,p=fa[p];
			}
		}
		return las;
	}
}S;

struct Segment
{
	int sz,ls[M],rs[M],mx[M],pos[M];
	void upd(int x)
	{
		if(mx[ls[x]]>=mx[rs[x]]) mx[x]=mx[ls[x]],pos[x]=pos[ls[x]];
		else mx[x]=mx[rs[x]],pos[x]=pos[rs[x]];
	}
	int merge(int x,int y,int l,int r)
	{
		if(!x || !y) return x+y;
		int d=++sz,mid=(l+r)>>1;
		if(l==r){mx[d]=mx[x]+mx[y];pos[d]=l;return d;}
		ls[d]=merge(ls[x],ls[y],l,mid);
		rs[d]=merge(rs[x],rs[y],mid+1,r);
		upd(d); return d;
	}
	void insert(int &x,int l,int r,int p)
	{
		if(!x) x=++sz;
		if(l==r){++mx[x];pos[x]=l;return;}
		int mid=(l+r)>>1;
		if(p<=mid) insert(ls[x],l,mid,p);
		else insert(rs[x],mid+1,r,p);
		upd(x);
	}
	void query(int x,int l,int r,int L,int R)
	{
		if(L<=l && r<=R) {if(mx[x]>ans)ans=mx[x],ps=pos[x];return;}
		int mid=(l+r)>>1;
		if(L<=mid) query(ls[x],l,mid,L,R);
		if(R>mid) query(rs[x],mid+1,r,L,R);
	}
	void init()
	{
		for(int i=1;i<=S.sz;++i) bz[i][0]=S.fa[i];
		for(int j=1;j<=20;++j) for(int i=1;i<=S.sz;++i) 
			bz[i][j]=bz[bz[i][j-1]][j-1];
		for(int i=1;i<=S.sz;++i) b[S.mx[i]]++;
		for(int i=1;i<=S.sz;++i) b[i]+=b[i-1];
		for(int i=S.sz;i;--i) c[b[S.mx[i]]--]=i;
		for(int i=S.sz;i>1;--i) rt[S.fa[c[i]]]=merge(rt[S.fa[c[i]]],rt[c[i]],1,m);
	}
}tr;

int main()
{
#ifndef ONLINE_JUDGE
	freopen("CF666E.in","r",stdin);
	freopen("CF666E.out","w",stdout);
#endif
	S.sz=S.las=1;
	scanf("%s",s+1);n=strlen(s+1);
	for(int i=1;i<=n;++i) id[i]=S.extend(s[i]-'a');
	
	m=read();
	for(int i=1,l;i<=m;++i)
	{
		scanf("%s",s+1);l=strlen(s+1);S.las=1;
		for(int j=1,x;j<=l;++j) x=S.extend(s[j]-'a'),tr.insert(rt[x],1,m,i);
	}
	//for(int i=1;i<=S.sz;++i) printf("%d\n",S.fa[i]);
	tr.init();
	
	Q=read();
	while(Q--)
	{
		int l=read(),r=read(),x=read(),y=read(),now=id[y];
		for(int i=20;~i;--i) if(S.mx[bz[now][i]]>=y-x+1) now=bz[now][i];
		ans=0;ps=l;tr.query(rt[now],1,m,l,r);
		printf("%d %d\n",ps,ans);
	}

	return 0;
}
