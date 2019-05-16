//  main.cpp
//  字典树(高性能)

#include <cstring>
#define TYPES_OF_CHAR 26
#define MAX_WORDS 4000
#define MAX_WORD_LENGTH 101
struct Trie{
    int ch[MAX_WORDS*MAX_WORD_LENGTH][TYPES_OF_CHAR];
    int is_last[MAX_WORDS*MAX_WORD_LENGTH];
    int tree_size;
    Trie()
    {
        tree_size=1;
        memset(ch[0],0,sizeof(ch[0]));
        memset(is_last,0,sizeof(is_last));
    }
    void clear()
    {
        memset(ch[0],0,sizeof(ch[0]));
        memset(is_last,0,sizeof(is_last));
    }
    //⚠️需要按照情况更改
    int index_get(char c)
    {
        return c-'a';
    }
    
    void insert_word(char *a,int v)
    {
        int u=0,n=(int)strlen(a);
        for(int i=0;i<n;i++)
        {
            int c=index_get(a[i]);
            if(!ch[u][c]){
                memset(ch[tree_size],0,sizeof(ch[tree_size]));
                is_last[tree_size]=0;
                ch[u][c] = tree_size++;
            }
            u=ch[u][c];
        }
        is_last[u]=v;
    }
};

