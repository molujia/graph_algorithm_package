# include <queue>
# include <cstring>
# include <iostream>
# include <queue>
# define LL long long
using namespace std;
// #define getc() (((__S==__T)&&((__T=(__S=__B)+fread(__B,1,1<<15,stdin)),__S==__T))?EOF:*__S++)
template<class TT>
inline void qread(TT &x){
     x=0; bool f=false; char c=getchar();
     for(;c<48||c>57;c=getchar())f|=(c=='-');
     for(;c>47&&c<58;c=getchar())x=(x<<1)+(x<<3)+(c^48);
     (f)&&(x=-x);
}
const int N = 2002, M = 120002, INF_INT = 0x3f3f3f3f;
const LL INF_LL = 9223372036854775807;
int head[N+2], to[(M+2)<<1], nxt[(M+2)<<1], cnt=1;
int val[(M+2)<<1], maxFlow=0;
int n,m,s,t;
int FLAG[N+2][N+2], Height[N+2], gap[N<<1], Extra[N+2];
bool inQ[N+2];
struct cmp
{
	inline bool operator()(int a,int b) const
	{
		return Height[a]<Height[b];//因为在优先队列中的节点高度不会改变，所以可以直接比较
	}
};
priority_queue<int, vector<int>, cmp> PQ;
void AddEdge(int u, int v, LL w) {
    ++cnt;
    nxt[cnt]=head[u];
    to[cnt]=v;
    val[cnt]=w;
    head[u] = cnt;
    
}
int AddFLowEdge(int u, int v, LL w) {
    AddEdge(u,v,w);
    AddEdge(v,u,0);
    return cnt-1;
}
bool bfs(int s, int t) { //一次寻找一条增广路，不走重复路径
    queue<int> Q;
    memset(Height,INF_INT,sizeof(Height));
    Height[t]=0;
    Q.push(t);
    while (!Q.empty()) {
        int now = Q.front();
        Q.pop();
        for (int i=head[now];i;i=nxt[i]) {
            int dest = to[i];
			if (val[i^1]==0) continue;
            if (Height[dest]<=Height[now]+1) continue;
            Height[dest]=Height[now] + 1;
            Q.push(dest);
            
        }
    }
    return Height[s]!=INF_INT ;
}
void Push(int now) {
    for (int i=head[now];i&&Extra[now];i=nxt[i]) {
        if (val[i]==0) continue;
        int dest = to[i];
        if (Height[dest]!=Height[now]-1) continue; // 只向更低一级的地方推流。
        int flow = min(Extra[now], val[i]);
        val[i]-=flow;
        val[i^1]+=flow;
        Extra[now]-=flow;
        Extra[dest]+=flow;
        if (inQ[dest]==0) {
            inQ[dest]=1;
            PQ.push(dest);
        }
    }
    return;
}
void relabel(int now) {
    Height[now] = INF_INT;
    for (int i=head[now];i;i=nxt[i]) {
        int dest = to[i];
        if (val[i] && Height[now]>Height[dest]+1) {
            Height[now] = Height[dest]+1;
        }
    }
}
inline int HLPP(int s, int t)
{
	if(!bfs(s,t))//s和t不连通
		return 0;
	Height[s]=n;
	memset(gap,0,sizeof(gap));
	for(int i=1;i<=n;i++)
		if(Height[i]<=n)
			++gap[Height[i]];
            
    Extra[s]=INF_INT;
    inQ[s]=1;
    inQ[t]=1;
    for (int i=head[s];i;i=nxt[i]) {
        int dest = to[i];
        int flow = val[i];
		if(flow)
		{
			val[i]-=flow;
            val[i^1]+=flow;
            // Extra[s]-=flow;
            Extra[dest]+=flow;
            if (inQ[dest]||Height[dest]>n) continue;
			PQ.push(dest);
            inQ[dest]=1;
		}
    }
	while(!PQ.empty())
	{
        int now = PQ.top();
        PQ.pop();
        Push(now);
        inQ[now]=0;
		if(Extra[now])
		{
			if(!--gap[Height[now]])//gap优化，因为当前节点是最高的所以修改的节点一定不在优先队列中，不必担心修改对优先队列会造成影响
				for(int i=1;i<=n;i++)
					if(inQ[i]==0&&Height[i]>Height[now]&&Height[i]<INF_INT) {
						Height[i]=n+1;
                    }
			relabel(now);
            if (Height[now]!=INF_INT) {
            ++gap[Height[now]];
			PQ.push(now);
            inQ[now]=1;}
		}
	}
	return Extra[t];
}

int main() {
    // ios::synv_with_stdio(false);
    // cin>>n>>m>>s>>t;
    qread(n);
    qread(m);
    qread(s);
    qread(t);
    int u, v;
    LL w;
    for (int i=1; i<=m; ++i) {
        // cin>>u>>v>>w;
        qread(u);
        qread(v);
        qread(w);
        // 此处建边时不考虑有重边的情况
        // AddFLowEdge(u,v,w);
        if (u==v) continue;
        // 此处建边时考虑有重边的情况
        if (FLAG[u][v]) {
            int pt = FLAG[u][v];
            val[pt] +=w;
        }
        else {
            FLAG[u][v] = AddFLowEdge(u,v,w);
        }

    }
    cout<<HLPP(s ,t)<<endl;
}
#undef LL
#undef getc