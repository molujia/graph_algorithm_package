// Edmonds-Karp算法 
// 算法不采用dfs，主要是因为dfs有可能经过一条容量极小的边和他的反边而到达汇点的可能，这样每次更新的增广流量都非常的小，使得时间复杂度与边的容量相关。
// 这条容量小的边记作z一定不是必走不可的选择，否则就会一定程度锁死最大流的上限。
// 既然这样，一定有一条容量大的边实际上完成了经过这条容量小的边的功能，也就是沟通了z边的上游和下游。这样的边在bfs中会被先搜索到或者先从队列中出来
// 所以采用bfs比dfs优秀。

// 复杂度上限O(VE^2)，一般远不及上限。
// 复杂度证明：
// 证明的核心要点是每条增广路径都是一条最短路。
// 先证明每次bfs增广之后每个点到源点的最短距离是不降的。因为假设v点到源点的最短距离减少了，而且v点是第一个减少的，并且这条路径经过了边(u,v)。
// 由于v点是第一个距离减少的，所以u的没减少。那么这就说明，(u,v)这条边在此次增广之前不存在。因为如果存在了，v在此次增广前，他的最短距离至多是
// d[u]+1,增广之后还减少了，至多是d[u]。可是增广之后v的最短路径是由u到v，而且u的最短距离是不降的，所以增广之后v的最短距离至少是d[u]+1。
// 所以，(u,v)这条边在此次增广之前不存在。而这次增广之后存在了，说明这次增广把(v,u)走了。走了(v,u)意味着u在增广前的最短路径有一截(v,u)，所以增广前的
// d[u]等于d[v]+1。而增广之后D[v]=D[u]+1,D[u]>=d[u]。所以D[v]>=d[v]+2。跟假设矛盾。所以不存在这样的v点。

// 定义一条边的容量在一次增广路径中被走光称为关键边。假设(u,v)某次成为了关键边。注意到(u,v)再次回到网络中必定是走了(v,u)这条边。结合上面证明的结论
// 这一来一去，d[u]至少增加了2。一个点到源点的最短距离最大也就是V。所以某条边(u,v)成为关键边最多成为V/2次。每一次增广都至少有一条
// 关键边，所以总共最多增广VE/2次。增广一次本身是一次bfs，复杂度做多是E。所以总的复杂度是O(VE^2)
// 由洛谷P3376 提供正确性检验

# include <iostream>
# include <queue>
# include <cstring>
# define LL long long
using namespace std;

const int N = 200, M = 5000, INF = 2005020600;
int head[N+2], to[(M+2)<<1], nxt[(M+2)<<1], cnt=1, pre[N+2];
LL val[(M+2)<<1], maxFlow=0, dist[(M+2)<<1];
int n,m,s,t;
bool vis[N+2];
int FLAG[N+2][N+2];
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


int bfs(int s, int t) { //一次寻找一条增广路，不走重复路径
    queue<int> Q;
    Q.push(s);
    memset(vis,0,sizeof(vis));
    int i=0;
    vis[s]=1;
    dist[s]=INF;
    while (!Q.empty()) {
        int cur = Q.front();
        Q.pop();
        for (int i=head[cur];i;i=nxt[i]) {
            if (val[i]==0) continue;
            int dest = to[i];
            if (vis[dest]) continue;
            dist[dest] = min(dist[cur], val[i]);
            pre[dest] = i;
            Q.push(dest);
            vis[dest] = 1;
            if(dest==t) return 1;
        }
    }
    return 0;
}
void EK(int s,int t)
{
    int increase=0;
    int i=0;
    while (bfs(s,t))
    {
        int k=t;
        increase = dist[t];
        while (k!=s)
        {
            
            int last=pre[k];//从后往前找路径
            val[last]-=increase;
            val[last^1]+=increase;
            k=to[last^1];
        }
        maxFlow+=increase;
    }
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
    EK(s ,t);
    cout<<maxFlow<<endl;
}