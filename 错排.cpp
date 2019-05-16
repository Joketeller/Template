long long d[26]={1,0};//要用long long 
 
void D()//求i个数错排有多少种
{
    for(int i=2;i<=14;i++)
        d[i]=(i-1)*(d[i-1]+d[i-2]);
}
 