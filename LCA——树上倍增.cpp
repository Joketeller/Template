const int max0=(int)log2(maxn)+1;

int father(int i,int k)  
{  
    for(int x=0;x<=int(log2(k));x++)  
        if((1<<x)&k)    //(1<<x)&k可以判断k的二进制表示中，<strong>第(x-1)位上</strong>是否为1  
            i=fa[i][x];     //把i往上提  
    return i;  
}  


void dfs(int x)  
{  
    for(int i=1;i<=max0;i++)  
        if(fa[x][i-1])   //在dfs(x)之前，x的父辈们的fa数组都已经计算完毕，所以可以用来计算x  
            fa[x][i]=fa[fa[x][i-1]][i-1];  
        else break;    //如果x已经没有第2^(i-1)个父亲了，那么也不会有更远的父亲，直接break  
    for(/*每一个与x相连的节点i*/)  
        if(i!=fa[x][0])     //如果i不是x的父亲就是x的儿子  
        {  
            fa[i][0]=x;       //记录儿子的第一个父亲是x  
            dep[i]=dep[x]+1;      //处理深度  
            dfs(i);  
        }  
}  

LCA(int u,int v)  
{  
    if(dep[u]<dep[v])swap(u,v);  //我们默认u的深度一开始大于v，那么如果u的深度小就交换u和v  
    int delta=dep[u]-dep[v];    //计算深度差  
    for(int x=0;x<=max0;x++)    //此循环用于提到深度相同。  
        if((1<<x)&delta)  
            u=fa[u][x];  
    if(u==v)return u;  
    for(int x=max0;x>=0;x--)     //<strong>注意！此处循环必须是从大到小!</strong>因为我们应该越提越“精确”，  
        if(fa[u][x]!=fa[v][x])   //如果从小到大的话就有可能无法提到正确位置，自己可以多想一下  
        {  
            u=fa[u][x];  
            v=fa[v][x];  
        }  
    return fa[u][0];    //此时u、v的第一个父亲就是LCA。  
}  