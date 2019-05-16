struct SAM{
    int last,cnt;//last表示目前自动机的最终态[能接受最长串的态],cnt是正在打下的标号
    int a[maxn][26],fa[maxn],mx[maxn];
     
    SAM(){last=++cnt;} //初始化,先要建立一个初始态
     
    void extend(int c){
        int p=last,np=last=++cnt;
        mx[np]=mx[p]+1;
        while(!a[p][c] && p)//一直沿着parent指针往上找到第一个有c出边的终态
            a[p][c]=np,p=fa[p];
        if(!p) fa[np]=1;//如果所有的点之前都没有含c的出边，那么包含新后缀的Right集合不存在，就将np的parent指向初始态
        else{
            int q=a[p][c];
            if(mx[p]+1==mx[q])//如果这个后缀正好能包含前面所有的子串，即不会出现一个以后缀形式包含此串却不是上一个串的后缀的串也连接到q，这种情况就可以保证直接将np的Parent练接到q即可[q的Right集合包含了np]
                fa[np]=q;
            else{
                int nq=++cnt;mx[nq]=mx[p]+1;
                memcpy(a[nq],a[q],sizeof(a[q]));
                fa[nq]=fa[q];//新建一个nq节点，并且将q复制一份到nq节点上
                fa[np]=fa[q]=nq;//之前的q上连接的就是包含这个后缀的其它串，而这个后缀上的串因为末尾添加可以在Right集合中加入一个新的位置，从而连到一个新的状态
                while(a[p][c]==q) a[p][c]=nq,p=fa[p];//这个新的状态也是之前那些不能转移节点所到达状态的Parent
            }
        }
    }
};

