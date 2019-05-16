//  字典树（STL）

#include <string>
#include <map>
using namespace std;
struct _node{
    map<char,_node*> m;
    int last=0;
}root;
void add(string s)
{
    _node *p=&root;
    for(int i=0;i<s.length();i++)
    {
        if(!p->m.count(s[i]))
        {
            p->m[s[i]]=new _node;
        }
        p=p->m[s[i]];
    }
    p->last++;
}

