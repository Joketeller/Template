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