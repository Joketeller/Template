

 


int exgcd(int a,int p,int &x,int &y)  
{
	int d=a;
	if (p!=0)
	{
		d=exgcd(p,a%p,y,x);
		y-=(a/p)*x; 
	}
	else 
	{
		x=1;
		y=0;
	}
	return d;
}

