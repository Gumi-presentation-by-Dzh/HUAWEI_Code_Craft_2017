
#include <deque>
#include <queue>
#include <vector>
#include <algorithm>
#include <cstring>
#include <utility>
#include <cstdio>
#include <climits>
#include <iostream>
#include <cstdlib>
#include <random>
#include <cmath>
#include <ctime>
#include <cmath>
#include "lib_io.h"
#include "deploy.h"

using namespace std;

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    char * topo_file;
    Graph graph;                                  
    graph.init_g(topo, line_num);              //读数据
    graph.allpair_sp();                                ////利用SPFA将所有点的最短路径找出来/这个里面的队列有优化的方式
    PSO pso(graph);                            

    vector<int> server, best_server;            
    vector<vector<int> > node;
    vector<int> flow;
    int best_index = graph.customer_num;
    int kmean_times = 2;
    int block_size = (int)(sqrt(graph.customer_num) + 0.1);
    
    graph.clustering(1, server);                     
    pso.add_candidate(server);
    graph.clustering(best_index, best_server);   


    long long best_cost = best_index * graph.server_cost;
    for (int i = best_index - 1; i > 1; i -= block_size) {
        for (int j = 0; j < kmean_times ; ++j) {
            graph.clustering(i, server);
            graph.add_server(server);
            long long cost = graph.costflow() + i * graph.server_cost;
            if (cost < best_cost) {
                best_cost = cost;
                best_index = i;
                best_server.swap(server);
            }
        }
    }
    pso.add_candidate(best_server);
    block_size = (block_size >> 1);
    int min_index = max(best_index - block_size, 1);
    int max_index = min(best_index + block_size, graph.customer_num);
    int max_size = 10;
    for (int i = min_index; i <= max_index; ++i) {
        for (int j = 0; j < kmean_times; ++j) {
            graph.clustering(i, server);
            pso.add_candidate(server);
        }
    }
    psoprocess -= pso.init(max_size);
    preprocess = psoprocess * 0.4;
    while ((double)clock() / CLOCKS_PER_SEC < preprocess)     
        pso.phase(1);
    pso.reset_particle();
    while ((double)clock() / CLOCKS_PER_SEC < psoprocess)
        pso.phase(2);
    
    
    pso.get_best_particle(best_server);                    
    graph.add_server(best_server);
    graph.costflow();
    graph.print_flow(node, flow);

    /*输出结果*/
    int node_size = node.size();
    topo_file = new char[node_size * vbound * 5];
    topo_file[0] = '\0';
    char line[vbound * 5];
    char tmp[100];
    sprintf(line, "%d\n\n", node_size);
    strcat(topo_file, line);
    for (int i = 0; i < node_size; ++i) {
        line[0] = '\0';
        int node_size_1 = node[i].size() - 1;
        for (int j = 0; j < node_size_1; ++j) {
            sprintf(tmp, "%d ", node[i][j]);
            strcat(line, tmp);
        }
        sprintf(tmp, "%d ", node[i][node_size_1] - graph.node_num);
        strcat(line, tmp);
        sprintf(tmp, "%d\n", (int)flow[i]);
        strcat(line, tmp);
        strcat(topo_file, line);
    }
    write_result(topo_file, filename);
    delete []topo_file;
}

bool cmp(const Particle & p1, const Particle & p2) {
    return p1.cost == p2.cost ? p1.cost_best < p2.cost_best : p1.cost < p2.cost;
}

template <class T>
void flush(vector<T> & v) {
    int i = v.size() - 1, j = 0;
    while (i >= 0) {
        j = rand() % (i + 1);
        swap(v[i], v[j]);
        --i;
    }
}

void Graph::init_g(char * topo[MAX_EDGE_NUM], int line_num) { 
    int line = 0;
    int u, v, c, w;
    if (line < line_num)
        sscanf(topo[line], "%d %d %d", &node_num, &edge_num, &customer_num);
    s = node_num + customer_num;
    t = s + 1;
    epool_flag = 0;
    memset(e, 0, sizeof(e));
    G.resize(node_num, vector<nnEdge>());
    d.resize(node_num, vector<int>(node_num, INF));
    line += 2;
    sscanf(topo[line], "%lld", &server_cost);
    line += 2;
    for (int i = 0; i < edge_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d %d", &u, &v, &w, &c);//链路起始节点、终止节点、总带宽、单位租用费
        G[u].emplace_back(v, c);
        G[v].emplace_back(u, c);
        add_edge(u, v, w, c);
        add_edge(v, u, w, c);
    }
    ++line;
    need_flow = 0;
    for (int i = 0; i < customer_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d", &u, &v, &w);
        consumers.emplace_back(v, w);
        add_edge(v, u + node_num, w, 0);
        add_edge(u + node_num, t, w, 0);
        need_flow += w;
    }
    tmp_epool_flat = epool_flag;
}

void Graph::allpair_sp() {                                 //利用SPFA将所有点的最短路径找出来
    for(int s = 0; s < node_num; ++s) {             
        d[s][s] = 0;
        deque<int> q;                
        memset(vis,0,sizeof(vis));
        vis[s] = 1;
        q.push_back(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            vis[u] = 0;
            for (unsigned int i = 0; i < G[u].size(); ++i) {
                int v = G[u][i].v;
                int dis =  d[s][u] + G[u][i].c;
                if (dis < d[s][v])
                { 
                    d[s][v] = dis;
                    if (!vis[v]) {
                        vis[v] = 1;
                        if (q.size () && d[s][v] < d[s][q[0]])
                            q.push_front(v);
                        else
                            q.push_back(v);
                    }
                }
            }
        }
    }
}

void Graph::clustering(int k, vector<int> & clusters) {
    clusters.resize(k);                             //重新 指定 有效元素的个数
    memset(vis, -1, sizeof(vis));                   //初始化vis标记vector
    vector<vector<int> > kmean_node(k);             //定义新kmean_node vector
    int min_dist, min_index;                       

    // flush(consumers);                  //将消费节点的顺序打乱    
    // for (int i = 0; i < k; ++i) {
    //     clusters[i] = consumers[i].v;          //取前k个消费节点里面的相连网络节点作为cluster
    // }
    
	vector<double> dist;
	dist.resize(consumers.size());
    clusters[0] = consumers[rand() % consumers.size()].v;           //kmeans++ 的初始化
    for (int n_cluster = 1; n_cluster < k; ++n_cluster) {
        double sum = 0;
        for (int i = 0; i < consumers.size(); ++i) {
            int Mindist = INF;
            for(int j = 0; j < n_cluster; ++j) {         //所有点到每一个聚类中心的最近距离，这个距离存放在d这个数组中   
                if(Mindist > d[consumers[i].v][clusters[j]])   //对于每个点，我们都计算其和最近的一个“种子点”的距离D(x)并保存在一个数组里
                    Mindist = d[consumers[i].v][clusters[j]];
            }
            sum += Mindist;
            dist[i] = Mindist;
        }
        sum = sum * rand() / (RAND_MAX - 1);                    //先取一个能落在Sum中的随机值Random
		for (int i = 0; i < consumers.size(); ++i) {
            if((sum -= dist[i]) > 0)	continue;
            clusters[n_cluster] = consumers[i].v;
            break;
        }
	}


    while (1) {
        for (int i = 0; i < k; ++i)
            kmean_node[i].clear();                  //初始化
        bool update = 0;                            //循环终止条件
        for (int i = 0; i < customer_num; ++i) {    //循环消费节点次数
            min_dist = INF;                         
            min_index = 0;                          //初始化最小值
            for (int j = 0; j < k; ++j) {
                if(d[clusters[j]][consumers[i].v] < min_dist) {        //找当前消费节点链接节点(consumers[i].v)的最近的cluster节点在cluster里面的位置
                    min_dist = d[clusters[j]][consumers[i].v];         
                    min_index = j;                                          //位置给min_index 大小给min_dist
                }
            }
            if (vis[i] != min_index) {              //如果当前vis[i]标记的不是consumers[i].v距离最近的cluster节点就更新
                update = 1;
                vis[i] = min_index;
            }
            kmean_node[vis[i]].push_back(i);        //可以理解这里就是聚类，将距离这个服务器最近的消费节点绑在这个服务器节点（kmean_node）的邻接表上
            //将consumers[i]的i放到kmean_node[（刚刚的出来离consumers[i].v距离最近的cluster节点）]
        }                                               
        if (!update)                            //如果没有更新 就不找了
            break;
        for (int j = 0; j < k; ++j) {           //如果有更新，说明聚类情况有改变，聚类还未收敛，遍历kmean_node
            min_dist = INF;                     
            min_index = 0;                          
            for (int l = 0; l < node_num; ++l) {            //遍历所有网络节点
                int dist = 0;
                for (unsigned int i = 0; i < kmean_node[j].size(); ++i) {   
                    dist += d[l][consumers[kmean_node[j][i]].v];       //当前节点到当前服务器节点链上所有消费节点的距离和
                }
                if (dist < min_dist) {              //有优化解
                    min_index = l;
                    min_dist = dist;
                }
            }
            clusters[j] = min_index;            //调整服务器位置
        }
    }
}

void Graph::add_server(vector<int> & v) {
    if (epool_flag != tmp_epool_flat) {
        epool_flag = tmp_epool_flat;
        for (Edge *j = e[s]; j; j = j->next) {
            int x = j->t;
            e[x] = e[x]->next;
        }
        e[s] = 0;
        Edge *j = epool + epool_flag;
        for (Edge *i = epool; i < j; ++i) {
            i->u = i->U;
            i->c = i->C;
        }
    }
    for (unsigned int i = 0; i < v.size(); ++i) {
        add_edge(s, v[i], INF, 0);
    }
}

void Graph::add_edge(int u, int v, int w, int c) {
    Edge *e1 = epool + epool_flag++, *e2 = epool + epool_flag++;
    *e1 = (Edge){v, w, c, w, c, e[u], e2}, e[u] = e1;
    *e2 = (Edge){u, 0, -c, 0, -c, e[v], e1}, e[v] = e2;
}


void Graph::print_flow(vector<vector<int> > & node, vector<int> &flow) {
    node.clear();
    flow.clear();
    while (1) {
        vector<int> Tmp;
        int u = s;
        int S = INF;
        while (u != t) {
            bool flag=0;
            for (Edge *i = e[u]; i; i = i->next) {
                int v = i->t;
                if (i->U > i->u) {
                    S = min(S, i->U - i->u);
                    u = v;
                    flag = 1;
                    break;
                }
            }
            if (!flag) break;
        }
        if (u != t) break;
        u = s;
        flow.push_back(S);
        while (u != t) {
            for (Edge *i = e[u]; i; i = i->next) {
                int v = i->t;
                if (i->U > i->u) {
                    i->u += S;
                    u = v;
                    break;
                }
            }
            if (u != t) Tmp.push_back(u);
        }
        node.push_back(Tmp);
    }
}

int Graph::aug(int u, int m) {
    if (u == t)
        return cost += (long long)dist * m, m;
    int d = m;
    vis[u] = 1;
    for (Edge *i = e[u]; i; i = i->next) {
        if (i->u && !i->c && !vis[i->t]) {
            int f = aug(i->t, min(d, i->u));
            i->u -= f;
            i->pair->u += f;
            d -= f;
            if (!d)
                return m;
        }
    }
    return m - d;
}

bool Graph::modlabel() {
    deque <int> q;
    memset(vis , 0, sizeof(vis));
    memset(D, 0x33 , sizeof(D));
    q.push_back(s);
    D[s] = 0;
    vis[s] = 1;
    while (!q.empty ()) {
        int u = q.front ();
        q.pop_front ();
        for (Edge *i = e[u]; i; i = i->next) {
            int v = i->t;
            int dis = D[u] + i->c;
            if (i->u && dis < D[v]) {
                D[v] = dis;
                if (!vis[v]) {
                    vis[v] = 1;
                    if (q.size () && D[v] < D[q[0]])
                        q.push_front(v);
                    else
                        q.push_back(v);
                }
            }
        }
        vis[u] = 0;
    }
    for (Edge *i = epool; i < epool + epool_flag; ++i) {
        i->c -= D[i->t] - D[i->pair->t];
    }
    dist += D[t];
    return D[t] < INF;
}

long long Graph::costflow() {
    flow = dist = 0;
    cost = 0;
    while (modlabel()) {
        int tmpf;
        do {
            memset(vis , 0, sizeof(vis));
            tmpf = aug(s, INT_MAX);
            flow += tmpf;
        } while (tmpf);
    }
    if (flow != need_flow)
        cost = INFLL;
    return cost;
}

PSO::PSO(Graph & g) {
    graph = &g;
    l = graph->node_num;
}

void PSO::decode(vector<double> & vd, vector<int> & vi) {
    vi.clear();
    for (int i = 0; i < l; ++i) {
        if (vd[i] > 0.5)
            vi.push_back(i);
    }
}

void PSO::add_candidate(vector<int> & v) {
    candidates.emplace_back(l, v, graph);
}

void PSO::get_best_particle(vector<int> & server) {
    decode(gbest.localbest, server);
}

void PSO::GA_cross(Particle & s1, Particle & s2) {
    int r1 = rand() % l, r2 = rand() % l;
    if (r1 > r2)
        swap(r1, r2);
    while (r1 < r2) {
        swap(s1.v[r1], s2.v[r1]);
        ++r1;
    }

}


void PSO::shuffle(Particle & s) {
    int r1, r2;
    do {
        r1 = rand() % l;
    } while (s.v[r1] > 0.5);
    do {
        r2 = rand() % l;
    } while (s.v[r2] < 0.5);
    swap(s.v[r1], s.v[r2]);
}

inline void PSO::PSO_update(Particle & s) {
   
    int  temp = 0;              //标记位
    int lg = (int)sqrt(graph->node_num);  //步长
   
    w = 0.88;                    //改变粒子之后还原W
   
    for (int i = 0; i < l; ++i) {
        s.vp[i] = w * s.vp[i] + c1 * rand() / RAND_MAX * (s.localbest[i] - s.v[i]) + c2 * rand() / RAND_MAX * (gbest.localbest[i] - s.v[i]); 
        s.v[i] = (1 / (1 + exp(100*(0.5-(s.v[i] + s.vp[i])))));

        //*自适应算法注入*// 
        temp++;   
        if(temp == lg && 2*i>l){        //调整W       //如果算法停滞(改用步长调整)，就对W采取适应
            get_f();                //算全局的F
            if(s.cost < fg_avg)
                w = w - (w - w_min)*abs((s.cost-fg_avg)/(gbest.cost_best-fg_avg));
            else
                if(s.cost > f_avg)
                    w = 1.5 - 1/(1 + k1 * exp(-k2 * deta));
            temp = 0;
        }
    }
}

void
PSO::get_f(){                                           //APSO global fitness
    int k = candidates.size();
    double sum_f = 0;
    for (int i = 0; i < k; ++i){
        sum_f += candidates[i].cost;
    }
    f_avg = sum_f / k;
    sum_f = 0;
    int count = 0;
    for (int i = 0; i < k; ++i){
        if(candidates[i].cost < f_avg){                 //此处的fitness为cost也就是意味着cost越低 fitness越好
            sum_f += candidates[i].cost;
            count++;
        }
    }
    fg_avg = sum_f / count;
    deta = abs(gbest.cost_best-fg_avg);
}



inline void PSO::update_best_particle(Particle & s) {
    vector<int> v;
    decode(s.v, v);
    graph->add_server(v);
    s.cost = graph->costflow() + v.size() * graph->server_cost;
    if (s.cost < s.cost_best) {
        s.localbest = s.v;
        s.cost_best = s.cost;
        if (s.cost_best < gbest.cost_best) {
            gbest.localbest = s.localbest;
            gbest.cost_best = s.cost_best;
            cnt = 0;
        }
    }
}

inline void PSO::reset_particle() {
    int k = candidates.size();
    for (int i = 0; i < k; ++i)
        candidates.push_back(candidates[i]);
}
void 
PSO::phase(int i) {
    switch(i)
    {
        case 1:{
            for (int i = 0; i < max_size; ++i) {
                shuffle(candidates[i]);
                update_best_particle(candidates[i]);
                PSO_update(candidates[i]);
            }
            if (++cnt > 200) {
                preprocess *= 0.5;
                psoprocess *= 0.5;                 
                cnt = 0;
            }
            break;
        }
        case 2:{
                int i = 0;
                int j = max_size >> 1;
                sort(candidates.begin(),candidates.end(), cmp);
                for (; i < j; ++i)
                    shuffle(candidates[i]);
                for (; i < max_size; ++i)
                    PSO_update(candidates[i]);
                for (i = 0, j = max_size - 1; i < j; ++i, --j) {
                    GA_cross(candidates[i], candidates[j]);
                    update_best_particle(candidates[i]);
                    update_best_particle(candidates[j]);
                }
                if (++cnt > 200) {
                    preprocess *= 0.5;
                    psoprocess *= 0.5;
                    cnt = 0;
                }
            break;
        }
    }

}

double PSO::init(int size) {
    max_size = size;
    c1 = 1.0;
    c2 = 1.6;
    w = 0.9;
    cnt = 0;

    sort(candidates.begin(), candidates.end(), cmp);
    gbest = candidates[0];
    //decode(gbest.localbest, v);

    int best_size = 0;
    for (int i = 0; i < l; ++i)
        best_size += gbest.localbest[i] > 0.5 ? 1 : 0;      //decode 过程
    if (best_size > 150)
        max_size = 4;
    else if (best_size > 100)
        max_size = 14;
    else 
        max_size = 10;
    best_size *= 0.7;

    int limit_size = min(max_size >> 1, (int)candidates.size());
    candidates.resize(limit_size);

    vector<int> v;

    for (int i = limit_size; i < max_size; ++i) {
        graph->clustering(best_size, v);
        add_candidate(v);
    }

    clock_t t1 = clock();
    phase(1);
    clock_t t2 = clock();
    return double(t2 - t1) / CLOCKS_PER_SEC;
}
