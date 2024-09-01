// Dinic算法，使用了多路增广，复杂度没有变化

# include <iostream>
# include <queue>
# include <cstring>
# define LL long long
using namespace std;

const int N = 200, M = 5000, INF_INT = 2147483647;
const LL INF_LL = 9223372036854775807;
int head[N+2], to[(M+2)<<1], nxt[(M+2)<<1], cnt=1, cur[N+2];
LL val[(M+2)<<1];
int n,m,s,t;
int FLAG[N+2][N+2], Depth[N+2];
LL maxFlow=0;
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
    // 注意循环中的写法可以直接改变cur的值
    LL used = 0;
    for (int &i=cur[now];i;i=nxt[i]) { 
        int dest = to[i];
        if (Depth[dest]==Depth[now]+1&&(val[i])){
            LL ret = dfs(dest, min(flow-used, val[i]));
            if (ret==0) { // 剪枝优化，和当前弧的意思是一样的，dest这个点已经没有意义了。
                Depth[dest]=0;
            }
            else {
                val[i]-=ret;
                val[i^1]+=ret;
                used+=ret;
                if (used==flow) break;
            }
        }
    }
    return used;
}
int bfs(int s, int t) { //一次寻找一条增广路，不走重复路径
    queue<int> Q;
    Q.push(s);
    memset(Depth,0,sizeof(Depth));
    Depth[s]=1;
    Q.push(s);
    while (!Q.empty()) {
        int now = Q.front();
        Q.pop();
        for (int i=head[now];i;i=nxt[i]) {
            if (val[i]==0) continue;
            int dest = to[i];
            if (Depth[dest]!=0) continue;
            Depth[dest]=Depth[now] + 1;
            Q.push(dest);
            if (dest==t) {
                return 1;
            }
        }
    }
    return 0;
}
LL Dinic(int s, int t)
{
    while (bfs(s,t)) {
        memcpy(cur, head, sizeof(head));
        dfs(s, INF_LL); // 内层循环不断找寻增广路，复杂度O(VE)
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
    cout<<Dinic(s ,t)<<endl;
}
#undef LL