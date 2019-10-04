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
