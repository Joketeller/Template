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


