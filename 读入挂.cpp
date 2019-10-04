
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