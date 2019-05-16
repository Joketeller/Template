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
