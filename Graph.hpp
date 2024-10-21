#ifndef Graph_hpp
#define Graph_hpp

#include <cstdint>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <queue>
#include "time.h"
#include <ctime>
#include <algorithm>
#include <cstring>
#include <chrono>
#include "MurmurHash3.h"
#include <sys/resource.h>


#include <limits>
using namespace std;


struct PS_edge {
    uint32_t edge_l;
    uint32_t edge_r;
    double priority;
    
    PS_edge(uint32_t edge_l, uint32_t edge_r, double priority) : edge_l(edge_l), edge_r(edge_r), priority(priority) {}
};

struct PS_edge2 {
    int edge_l;
    int edge_r;
    
    PS_edge2(int edge_l, int edge_r) : edge_l(edge_l), edge_r(edge_r) {}
};

struct ComparePriority {
    bool operator()(const PS_edge& a, const PS_edge& b) {
        return a.priority > b.priority;
    }
};

struct ComparePriority_min {
    bool operator()(const PS_edge& a, const PS_edge& b) {
        return a.priority < b.priority;
    }
};


// random seed
extern int random_seed;

// Parameters for Hash function
extern const unsigned long long arr[8];
extern unsigned long long Prime1;
extern unsigned long long Prime3;
extern unsigned long long Prime2;
extern unsigned long long Prime4;


class Graph{
private:
    string str;
    uint32_t l_n, r_n, m;
    uint32_t *edge_l;
    uint32_t *edge_r;
    vector<uint32_t> RS_sample_edge_l;
    vector<uint32_t> RS_sample_edge_r;
    uint32_t *degree_l;
    uint32_t *degree_r;
    uint32_t Reservoir_size;
    uint32_t *adj_l;
    uint32_t *adj_r;
    uint32_t *start_l;
    uint32_t *start_r;
    vector<vector<uint32_t> > RS_adj_l, RS_adj_r;
    vector<unordered_set<uint32_t> > RS_adj_l2, RS_adj_r2;
    vector<vector<uint32_t> > PS_adj_l, PS_adj_r;
    vector<vector<uint32_t> > DEABC_adj_l, DEABC_adj_r;
    vector<vector<uint32_t> > DEABC_PLUS_adj_l, DEABC_PLUS_adj_r;
    double zStar;
    int time_;
    double *count_l, *count_r;
    bool exactcnt;
    int TM;
    double delta; //decaying factor for past estimations (real value in [0,1))
    uint32_t J; //period for estimation update
    unordered_map<int, double> DEABC_count_l, DEABC_count_r; //estimations for basic method DEABC-0
    unordered_map<int, double> DEABC_estimations_l, DEABC_estimations_r; // estimations for main method DEABC
    
    vector<PS_edge> S_DEABC_PLUS;
    vector<PS_edge2> S_DEABC_PLUS2;
    uint32_t *hashmap;

    uint32_t Time_point;

    uint32_t Edge_num_point;

    vector<uint64_t> Time_point_butterflies;
    vector<uint64_t> Edge_num_point_butterflies;
    
public:
    Graph(const char *_dir, uint32_t BUCKET_BITS);
    ~Graph();
    void readGraph();
    uint64_t Butterfly_counting();

    vector<uint64_t> Butterfly_counting_time_point();
    vector<uint64_t> Butterfly_counting_edge_num_point();

    void init();

    
    uint64_t Reservoir_sampling(uint32_t BUCKET_BITS);
    int intersectionSize(const vector<uint32_t>& vec1, const vector<uint32_t>& vec2);
    
    uint64_t Priority_sampling(uint32_t BUCKET_BITS, int is_hash);
    double hash(uint32_t a, uint32_t b, uint32_t c, uint32_t d, mt19937 &generator);
    
    uint64_t Priority_sampling_DEABC0(uint32_t BUCKET_BITS);
    double EdgeHash_p(uint32_t u, uint32_t v, unsigned long long a1, unsigned long long b1);
    uint32_t hash_func(uint32_t a);
    
    uint64_t Priority_sampling_DEABC_PLUS(uint32_t BUCKET_BITS);
    int hash_bucket(uint32_t a, uint32_t b, uint32_t c, uint32_t d, mt19937 generator);
    double GetCounts_DEABC_PLUS(uint32_t v_l, uint32_t v_r, uint64_t &sum_degree_l, uint64_t &sum_degree_r, uint32_t i);
    double getConstant(uint32_t buckets);
    double hash_to_double(uint32_t v_l, uint32_t v_r);
    double hash(uint32_t a, uint32_t b, uint32_t c, uint32_t d);
    uint64_t murmur_hash_64(uint32_t a, uint32_t b, uint32_t c, uint32_t d);
    uint32_t hash_bucket(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t range);
    
};
#endif /* Graph_hpp */
