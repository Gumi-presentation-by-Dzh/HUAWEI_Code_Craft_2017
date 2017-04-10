#ifndef __ROUTE_H__
#define __ROUTE_H__

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

const int vbound = 2050;
const int INF = 0x33333333;
const long long INFLL = 0x3333333333333333;
static double psoprocess = 87;
static double preprocess;

struct nnEdge {
    int v;
    int c;//与之相连的顶点，单位租用费用
    nnEdge(int _v, int _c): v(_v), c(_c) {}
};

struct Consumer {
    int v, w;//与之相连的网络节点，需求带宽
    Consumer(int _v, int _w): v(_v), w(_w) {}
};

struct Edge {
    int t;
    int u;
    int c;
    int U; 
    int C;
    Edge *next, *pair;
};

class Graph {
public:
    void init_g(char * topo[MAX_EDGE_NUM], int line_num);         //读取信息
    void allpair_sp();                                            //allpair_sp最短路径
    void clustering(int k, vector<int> & clusters);             //聚类

    void add_server(vector<int> & v);
    void add_edge(int u, int v, int w, int c);
    long long costflow();
    int aug(int u, int m);
    bool modlabel();
    void print_flow(vector<vector<int> > &node, vector<int> &flow);  

    int need_flow, node_num, edge_num, customer_num;        //视频带宽消耗需求总数，网络节点、边、消费节点总数
    long long server_cost;  //服务器部署成本
    vector<vector<nnEdge> > G;                 //图信息
    vector<Consumer> consumers;         //消费节点情况
    vector<vector<int> > d;                          //放的最短路径
    int vis[vbound];                                  //一个bool数组用来一般做节点是不是经过的处理
    Edge epool[MAX_EDGE_NUM * 5];
    Edge *e[vbound];
    int epool_flag, s, t, tmp_epool_flat;//s为网络节点+消费节点  
    int flow, dist, D[vbound];
    long long cost;                                 
};

class Particle {
public:
    Particle(int length=0): v(length, 0), localbest(length, 0), vp(length, 0), cost_best(INFLL), cost(INFLL) {}

    Particle(int length, vector<int> & vi, Graph * & graph): v(length, 0), localbest(length, 0), vp(length, 0) {
        int size = vi.size();
        for (int i = 0; i < size; ++i) {
            localbest[vi[i]] = v[vi[i]] = 1;
        }
        graph->add_server(vi);
        cost = cost_best = graph->costflow() + size * graph->server_cost;
    }
    vector<double> v;
    vector<double> localbest;
    vector<double> vp;
    long long cost_best;
    long long cost;
};

// used for pso algorithm
class PSO {
public:
    PSO(Graph & graph);
    double init(int size);
    inline void reset_particle();
    void add_candidate(vector<int> & v);
    void phase(int i);
    void get_best_particle(vector<int> & server);
    void decode(vector<double> & vd, vector<int> & vi);
    void GA_cross(Particle & s1, Particle & s2);
    void shuffle(Particle & s);
    inline void PSO_update(Particle & s);
    inline void update_best_particle(Particle & s);
    vector<Particle> candidates;
    Particle gbest;
    int l;
    int max_size;
    int cnt;
    double c1, c2, w;
    Graph *graph;


    //apso优化 所需要的结构体
    double w_min = 0.5;

    int k1 = 1.5;       //APSO w自适应算法参数
    int k2 = 0.3;

    double f_avg;           //粒子群的平均适应值
    double fg_avg;          //适应值优于f_avg的适应值求平均得到 
    double deta;            //评价粒子群的早熟收敛程度,fm-fg_avg    fm = gbest.fi
    
    void get_f();               //获得当前粒子群平均适应值
    
};

template <class T>
void flush(vector<T> & v);

bool cmp(const Particle & p1, const Particle & p2);

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename);

#endif
