// ISAP算法
// 之前已经有结论，一个点u在增广之后，到s的最短路径不会更短。同理，到t的最短路径也不会更短。并且
// 将不能到达t的点的最短路径定义为INF，可以证明，这样的定义逻辑上上是没问题的。
// dinic算法的核心是利用每个点最短路径的不降性质，通过bfs分层的做法，实际上是在对每个点到s的最短路径进行从小到大枚举。
// 但是dinic记录的是到s的最短距离，而在dfs中被枚举到的点无法得知比自己距离s更近的点的完整信息，从而无法及时作出是否应当更新自己深度的判断。
// 如果考虑为记录到t的最短距离，那么回溯到u点时，u点可以根据回溯信息进行判断是否更新自己的深度。具体方式如下：
// ISAP算法对每个点u始终维持d[u]<=D[u]。D[u]为u到t真正的最短路长度。
// 我们可知，如果一次经过u的增广成功，那么d[u]==D[u]。如果经过u的增广路径走遍了还剩余流量，那么此时u的最短路径不可能是d[u]了，所以d[u]++;
// 在最短路分层设定下，u只走到下一层可走的点。对于下一层之后的点v，必定不存在经过(u,v)合理的增广路径。否则与整个图的层次结构矛盾。所以u走遍增广
// 路只需要走遍到d[u]-1层的边。
// 初始分层时，由于只要满足d[u]<=D[u]，从t无论边的正反随意跑一个bfs。
// 可以使用gap优化。
// ISAP在一遍dfs内一般多路增广。同dinic，ISAP也可以采取当前弧优化。其原理实际跟dinic一样，故时间复杂度也是O(V^2E)
// 由洛谷P3376 提供正确性检验

# include <iostream>
# include <queue>
# include <cstring>
# define LL long long
using namespace std;

const int N = 200, M = 5000, INF_INT = 2147483647;
const LL INF_LL = 9223372036854775807;
int head[N+2], to[(M+2)<<1], nxt[(M+2)<<1], cnt=1, cur[N+2];
LL val[(M+2)<<1], maxFlow=0;
int n,m,s,t;
int FLAG[N+2][N+2], Depth[N+2], gap[N+2];
// 网络流的建边由于本身算法的时间复杂度比较高可以建邻接矩阵也不会爆空间。但是若用邻接表，最好从1开始建边，因为这样反边和正边的关系就是异或1.
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

LL dfs(int now, LL flow) {
    if (now==t) {
		maxFlow+=flow;
        return flow;
    }
	LL used = 0;
    // 注意循环中的写法可以直接改变cur的值
    for (int &i=cur[now];i;i=nxt[i]) { 
        int dest = to[i];
        if (Depth[dest]==Depth[now]-1&&(val[i])){
            LL ret = dfs(dest, min(flow-used, val[i]));
            if (ret>0) { 
                val[i]-=ret;
                val[i^1]+=ret;
                used+=ret;
            }
			if (used==flow) {
				return used;
			}
        }
    }
	--gap[Depth[now]];
	if (gap[Depth[now]]==0) {
		Depth[s]=n;
	}
	++Depth[now];
	gap[Depth[now]]++;
    return used;
}
void bfs(int s, int t) { //一次寻找一条增广路，不走重复路径
    queue<int> Q;
    memset(Depth,0x7f,sizeof(Depth));
    Depth[t]=0;
	gap[0]=1;
    Q.push(t);
    while (!Q.empty()) {
        int now = Q.front();
        Q.pop();
        for (int i=head[now];i;i=nxt[i]) {
            int dest = to[i];
			// 考虑反边最大限度建立比较精确的分层图。
			if (val[i^1]==0) continue;
            if (Depth[dest]<n) continue;
            Depth[dest]=Depth[now] + 1;
			gap[Depth[dest]]++;
            Q.push(dest);
            
        }
    }
    return ;
}
LL ISAP(int s, int t)
{
	bfs(s,t);
    while (Depth[s]<n) {
        memcpy(cur, head, sizeof(head));
        dfs(s,INF_LL);
    }
    return maxFlow;
}

int main() {
    cin>>n>>m>>s>>t;
    int u, v;
    LL w;
    for (int i=1; i<=m; ++i) {
        cin>>u>>v>>w;
        // 此处建边时不考虑有重边的情况
        // AddFLowEdge(u,v,w);
        
        // 此处建边时考虑有重边的情况
        if (FLAG[u][v]) {
            int pt = FLAG[u][v];
            val[pt] +=w;
        }
        else {
            FLAG[u][v] = AddFLowEdge(u,v,w);
        }

    }
    cout<<ISAP(s ,t)<<endl;
}
#undef LL