/*******************************EXKMP**************************************/

//计算next数组，保存到next参数中。
void getExtendNext (const char *t, int *next)
{
	int lt = strlen(t);
	for (int i = 1, j = -1, a, p; i < lt; i++, j--)
		if (j < 0 || i + next[i - a] >= p)
		{
			if (j < 0)	j = 0, p = i;
			while (p < lt && t[j] == t[p])	j++, p++;
			next[i] = j, a = i;
		}
		else	next[i] = next[i - a];
}

//扩展KMP算法求extend数组，s为文本串，t为模式串
void getExtend (const char *s, const char *t, int *extend)
{
	int ls = strlen(s), lt = strlen(t);
	int next[mx]； getExtendNext(t, next);	//计算next数组。mx是表示最大长度的常数。
	for (int i = 0, j = -1, a, p; i < ls; i++, j--)
		if (j < 0 || i + next[i - a] >= p)
		{
			if (j < 0)	j = 0, p = i;
			while (p < ls && j < lt && s[p] == t[j]) 	j++, p++;
			extend[i] = j, a = i;
		}
		else	extend[i] = next[i - a];
}
//********************************************************************//
//******************************KMP***********************************//
 void get_next(char *s,int *nxt)//数组为nxt,从0开始的序号
{
    int len=strlen(s);
	int i=0,k=-1;
	nxt[0]=-1;
	while (i<len)
	{ 
		if (k==-1||s[i]==s[k])
		{
			i++;
			k++;
			if (s[i]!=s[k]) nxt[i]=k;
			else	nxt[i]=nxt[k];
		}
		else 
		 k=nxt[k];
	}
}
	//*****************************************************************// 
