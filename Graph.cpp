#include "Graph.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

// random seed
int random_seed = 0;

// Parameters for Hash function
const unsigned long long arr[8] = {0xFC13C8E7,
            0xA2A9FFD4,
            0x597ECDDC,
            0x8AF8DA7E,
            0xAF531D42,
            0x842A21DD,
            0x1DEE299F,
            0xBFEC63E9};

unsigned long long Prime1 = 3584999771;
unsigned long long Prime3 = 67532401;
unsigned long long Prime2 = 4294967291;
unsigned long long Prime4 = 8532401;

Graph::Graph(const char *_dir, uint32_t BUCKET_BITS) {
    this->str = string(_dir);
    l_n = m = r_n = Reservoir_size = 0;
    edge_l = NULL;
    edge_r = NULL;
    degree_l = NULL;
    degree_r = NULL;
    adj_l = NULL;
    adj_r = NULL;
    start_l = NULL;
    start_r = NULL;
    
    count_l = NULL;
    count_r = NULL;
    hashmap = NULL;
    TM = 0;
    time_ = 0;
    Time_point = 100;
    Edge_num_point = 10000;

    this->Reservoir_size = pow(2, BUCKET_BITS);
}

Graph::~Graph(){
    if(edge_l!=NULL){
        delete [] edge_l; edge_l = NULL;
    }
    if(edge_r!=NULL){
        delete [] edge_r; edge_r = NULL;
    }
    if(degree_l!=NULL){
        delete [] degree_l; degree_l = NULL;
    }
    if(degree_r!=NULL){
        delete [] degree_r; degree_r = NULL;
    }
    if(adj_l!=NULL){
        delete [] adj_l; adj_l = NULL;
    }
    if(adj_r!=NULL){
        delete [] adj_r; adj_r = NULL;
    }
    if(start_l!=NULL){
        delete [] start_l; start_l = NULL;
    }
    if(start_r!=NULL){
        delete [] start_r; start_r = NULL;
    }
    if(count_l!=NULL){
        delete [] count_l; count_l = NULL;
    }
    if(count_r!=NULL){
        delete [] count_r; count_r = NULL;
    }

    if(hashmap!=NULL){
        delete [] hashmap; hashmap = NULL;
    }

    
}

void Graph::readGraph(){
    ifstream infile;   //输入流
    
    infile.open(str, ios::in);
    if (!infile.is_open()){
        cout<<"Open file failure"<<endl;
        exit(0);
    }
    infile>>l_n>>r_n>>m;
    
    // cout<<"l_n:"<<l_n<<" r_n:"<<r_n<<" m:"<<m<<endl;
    
//    m = m/1000;
    
    if(edge_l == NULL) edge_l = new uint32_t[m];
    if(edge_r == NULL) edge_r = new uint32_t[m];
    if(degree_l == NULL) degree_l = new uint32_t[l_n];
    if(degree_r == NULL) degree_r = new uint32_t[r_n];
    if(adj_l == NULL) adj_l = new uint32_t[m];
    if(adj_r == NULL) adj_r = new uint32_t[m];
    if(start_l == NULL) start_l = new uint32_t[l_n + 1];
    if(start_r == NULL) start_r = new uint32_t[r_n + 1];
    if(count_l == NULL) count_l = new double[l_n];
    if(count_r == NULL) count_r = new double[r_n];
    if(hashmap == NULL) hashmap = new uint32_t[std::max(l_n, r_n)];
    RS_adj_l.resize(l_n);
    RS_adj_r.resize(r_n);
    RS_adj_l2.resize(l_n);
    RS_adj_r2.resize(r_n);
    PS_adj_l.resize(l_n);
    PS_adj_r.resize(r_n);
    DEABC_adj_l.resize(l_n);
    DEABC_adj_r.resize(r_n);
    DEABC_PLUS_adj_l.resize(l_n);
    DEABC_PLUS_adj_r.resize(r_n);

    RS_sample_edge_l.resize(this->Reservoir_size);
    RS_sample_edge_r.resize(this->Reservoir_size);

    // S_DEABC_PLUS2.resize(this->Reservoir_size, PS_edge2(-1, -1));

    Time_point_butterflies.resize(Time_point,0);
    Edge_num_point_butterflies.resize(m/Edge_num_point, 0);

    memset(degree_l, 0, sizeof(uint32_t)*l_n);
    memset(degree_r, 0, sizeof(uint32_t)*r_n);
    
    memset(count_l, 0, sizeof(double)*l_n);
    memset(count_r, 0, sizeof(double)*r_n);
    
    uint32_t m_i = 0;
    uint32_t left_ = 0, right_ = 0;
    
    while (infile>>left_>>right_){
        edge_l[m_i] = left_;
        edge_r[m_i++] = right_;
        degree_l[left_] ++;
        degree_r[right_] ++;
        if(m_i>=m) break;
    }
    
    start_l[0] = 0;
    start_r[0] = 0;
    for(uint32_t i = 0;i<l_n;i++){
        start_l[i+1] = start_l[i] + degree_l[i];
        RS_adj_l[i].reserve(degree_l[i]);
        RS_adj_l2[i].reserve(degree_l[i]*1.2);
        PS_adj_l[i].reserve(degree_l[i]);
        DEABC_adj_l[i].reserve(degree_l[i]);
        DEABC_PLUS_adj_l[i].reserve(degree_l[i]);
    }
    for(uint32_t i = 0;i<r_n;i++){
        start_r[i+1] = start_r[i] + degree_r[i];
        RS_adj_r[i].reserve(degree_r[i]);
        RS_adj_r2[i].reserve(degree_r[i]*1.2);
        PS_adj_r[i].reserve(degree_r[i]);
        DEABC_adj_r[i].reserve(degree_r[i]);
        // DEABC_PLUS_adj_r[i].resize(0);
        DEABC_PLUS_adj_r[i].reserve(degree_r[i]);
    }
    for(uint32_t i = 0;i <m;i++){
        adj_l[start_l[edge_l[i]]++] = edge_r[i];
        adj_r[start_r[edge_r[i]]++] = edge_l[i];
    }
    start_l[0] = 0;
    start_r[0] = 0;
    for(uint32_t i = 0;i<l_n;i++){
        start_l[i+1] = start_l[i] + degree_l[i];
        sort(adj_l+start_l[i], adj_l+start_l[i+1]);
    }
    for(uint32_t i = 0;i<r_n;i++){
        start_r[i+1] = start_r[i] + degree_r[i];
        sort(adj_r+start_r[i], adj_r+start_r[i+1]);
    }

}

uint64_t Graph::Butterfly_counting(){
    

    uint64_t sun_degree_l = 0;
    uint64_t sun_degree_r = 0;
    uint32_t *adj_u;
    uint32_t *adj_v;
    uint32_t *start_u;
    uint32_t *start_v;
    uint32_t *hashmap;
    uint32_t *visit_vertex;
    uint32_t n_, visit_vertex_i;
    uint64_t butterfly_num = 0;
    
    for(uint32_t i = 0; i<l_n; i++){
        sun_degree_l += degree_l[i]*degree_l[i];
    }
    for(uint32_t i = 0; i<r_n; i++){
        sun_degree_r += degree_r[i]*degree_r[i];
    }
    
    
    if(sun_degree_l < sun_degree_r){
        adj_v = adj_r;
        adj_u = adj_l;
        start_v = start_r;
        start_u = start_l;
        hashmap = new uint32_t[l_n];
        memset(hashmap, 0, sizeof(uint32_t)*l_n);
        visit_vertex = new uint32_t[r_n];
        n_ = r_n;
    }else{
        adj_v = adj_l;
        adj_u = adj_r;
        start_v = start_l;
        start_u = start_r;
        hashmap = new uint32_t[r_n];
        memset(hashmap, 0, sizeof(uint32_t)*r_n);
        visit_vertex = new uint32_t[l_n];
        n_ = l_n;
    }
    
    for(uint32_t v = 0; v < n_; v++){
        visit_vertex_i = 0;
//        if(v%1000 == 0) cout<<v<<endl;
        
        for(uint32_t u = start_v[v]; u<start_v[v+1]; u++){
//            cout<<adj_v[u]<<" "<< start_u[adj_v[u]]<<endl;
            for(uint32_t w = start_u[adj_v[u]]; w < start_u[adj_v[u] + 1]; w++){
                if(adj_u[w] >= v) continue;
                if(hashmap[adj_u[w]] == 0) visit_vertex[visit_vertex_i++] = adj_u[w];
                hashmap[adj_u[w]]++;
                
            }
            
        }
        
        for(uint32_t i = 0;i<visit_vertex_i;i++){
            butterfly_num += hashmap[visit_vertex[i]]*(hashmap[visit_vertex[i]]-1)/2;
            hashmap[visit_vertex[i]] = 0;
        }
        
    }
    
    return butterfly_num;
    
}

vector<uint64_t> Graph::Butterfly_counting_time_point(){

    uint32_t time_num = m/Time_point;

    uint32_t time_i = 0;

    vector<vector<uint32_t> > adj_l_(l_n);
    vector<vector<uint32_t> > adj_r_(r_n);

    uint64_t sum_degree_l = 0;
    uint64_t sum_degree_r = 0;
    uint32_t n_, visit_vertex_i;
    uint32_t *visit_vertex = new uint32_t[max(l_n, r_n)];
    uint32_t *hashmap = new uint32_t[max(l_n, r_n)];

    for(uint32_t i = 0;i<l_n;i++){
        adj_l_.reserve(degree_l[i]);
    }
    for(uint32_t i = 0;i<r_n;i++){
        adj_r_.reserve(degree_r[i]);
    }
    

    for(uint32_t i = 0;i < m;i++){

        adj_l_[edge_l[i]].push_back(edge_r[i]);
        adj_r_[edge_r[i]].push_back(edge_l[i]);

        if(i%time_num == 0){
            // cout<<i<<endl;
            for(uint32_t j = 0;j<l_n;j++){
                sort(adj_l_[j].begin(), adj_l_[j].end());
                sum_degree_l += adj_l_[j].size()*adj_l_[j].size();
            }
            for(uint32_t j = 0;j<r_n;j++){
                sort(adj_r_[j].begin(), adj_r_[j].end());
                sum_degree_r += adj_r_[j].size()*adj_r_[j].size();
            }

            bool use_left_side = (sum_degree_l < sum_degree_r);
            auto& adj_v = use_left_side ? adj_r_ : adj_l_;
            auto& adj_u = use_left_side ? adj_l_ : adj_r_;
            n_ = use_left_side ? r_n : l_n;

            for(uint32_t v = 0; v < n_; v++){
                visit_vertex_i = 0;
        
                for(uint32_t u = 0; u<adj_v[v].size(); u++){
                    
                    for(uint32_t w = 0; w < adj_u[adj_v[v][u]].size(); w++){
                        if(adj_u[adj_v[v][u]][w] >= v) continue;
                        if(hashmap[adj_u[adj_v[v][u]][w]] == 0) visit_vertex[visit_vertex_i++] = adj_u[adj_v[v][u]][w];
                        hashmap[adj_u[adj_v[v][u]][w]]++;
                
                    }
            
                }
        
                for(uint32_t i = 0;i<visit_vertex_i;i++){
                    Time_point_butterflies[time_i] += hashmap[visit_vertex[i]]*(hashmap[visit_vertex[i]]-1)/2;
                    hashmap[visit_vertex[i]] = 0;
                }
        
            }

            time_i++;
        }

    }

    return Time_point_butterflies;
}

vector<uint64_t> Graph::Butterfly_counting_edge_num_point(){

    uint32_t time_num = Edge_num_point;

    uint32_t time_i = 0;

    vector<vector<uint32_t> > adj_l_(l_n);
    vector<vector<uint32_t> > adj_r_(r_n);

    uint64_t sum_degree_l = 0;
    uint64_t sum_degree_r = 0;
    uint32_t n_, visit_vertex_i;
    uint32_t *visit_vertex = new uint32_t[max(l_n, r_n)];
    uint32_t *hashmap = new uint32_t[max(l_n, r_n)];

    for(uint32_t i = 0;i<l_n;i++){
        adj_l_.reserve(degree_l[i]);
    }
    for(uint32_t i = 0;i<r_n;i++){
        adj_r_.reserve(degree_r[i]);
    }
    

    for(uint32_t i = 0;i < m;i++){

        adj_l_[edge_l[i]].push_back(edge_r[i]);
        adj_r_[edge_r[i]].push_back(edge_l[i]);

        if(i%time_num == 0){
            for(uint32_t j = 0;j<l_n;j++){
                sort(adj_l_[j].begin(), adj_l_[j].end());
                sum_degree_l += adj_l_[j].size()*adj_l_[j].size();
            }
            for(uint32_t j = 0;j<r_n;j++){
                sort(adj_r_[j].begin(), adj_r_[j].end());
                sum_degree_r += adj_r_[j].size()*adj_r_[j].size();
            }

            bool use_left_side = (sum_degree_l < sum_degree_r);
            auto& adj_v = use_left_side ? adj_r_ : adj_l_;
            auto& adj_u = use_left_side ? adj_l_ : adj_r_;
            n_ = use_left_side ? r_n : l_n;

            for(uint32_t v = 0; v < n_; v++){
                visit_vertex_i = 0;
        
                for(uint32_t u = 0; u<adj_v[v].size(); u++){
                    
                    for(uint32_t w = 0; w < adj_u[adj_v[v][u]].size(); w++){
                        if(adj_u[adj_v[v][u]][w] >= v) continue;
                        if(hashmap[adj_u[adj_v[v][u]][w]] == 0) visit_vertex[visit_vertex_i++] = adj_u[adj_v[v][u]][w];
                        hashmap[adj_u[adj_v[v][u]][w]]++;
                
                    }
            
                }
        
                for(uint32_t i = 0;i<visit_vertex_i;i++){
                    Edge_num_point_butterflies[time_i] += hashmap[visit_vertex[i]]*(hashmap[visit_vertex[i]]-1)/2;
                    hashmap[visit_vertex[i]] = 0;
                }
        
            }

            time_i++;
        }

    }
    return Edge_num_point_butterflies;
}


uint64_t Graph::Reservoir_sampling(uint32_t BUCKET_BITS) {
    double butterfly_num = 0;
    this->Reservoir_size = pow(2, BUCKET_BITS);
    double currentEdges = 0;
    double increment = 0;
    double y = 0;
    double Pr = 0;
    uint64_t sum_degree_l = 0, sum_degree_r = 0;

    uint32_t time_num = m/Time_point;
    
    srand(static_cast<unsigned int>(time(0)));
    std::random_device rd;
    std::mt19937 gen(rd());

    int RS_sample_edge_l_i = 0;

    std::uniform_real_distribution<> dis(0.0, 1.0);

    // memset(hashmap, 0, sizeof(uint32_t) * std::max(l_n, r_n));

    for (uint32_t i = 0; i < m; i++) {
        //if (i % time_num == 0) cout << (uint64_t)butterfly_num << endl;
        if(find(RS_adj_l[edge_l[i]].begin(),RS_adj_l[edge_l[i]].end(),edge_r[i]) != RS_adj_l[edge_l[i]].end()) continue;
        currentEdges++;

        if(currentEdges < Reservoir_size) y = currentEdges;
        else y = (double)Reservoir_size;
        Pr = (y/currentEdges) * ((y-1)/(currentEdges-1)) * ((y-2)/(currentEdges-2));
        increment = 1/Pr;

        bool use_left_side = (sum_degree_l < sum_degree_r);
        auto& RS_adj_v = use_left_side ? RS_adj_r : RS_adj_l;
        auto& RS_adj_u = use_left_side ? RS_adj_l : RS_adj_r;
        uint32_t v = use_left_side ? edge_r[i] : edge_l[i];
        uint32_t u = use_left_side ? edge_l[i] : edge_r[i];

        for (uint32_t w : RS_adj_v[v]) {
            hashmap[w] = i+1;
        }
        for (uint32_t w : RS_adj_u[u]) {
            if (w == v) continue;
            for (uint32_t w2 : RS_adj_v[w]) {
                if (w2 == u) continue;
                
                if (hashmap[w2] == i+1) {
                    butterfly_num += increment;
                }
            }
        } 

        if (currentEdges <= Reservoir_size) {
            
            RS_sample_edge_l[RS_sample_edge_l_i] = edge_l[i];
            RS_sample_edge_r[RS_sample_edge_l_i++] = edge_r[i];
            RS_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + RS_adj_l[edge_l[i]].size()*RS_adj_l[edge_l[i]].size() - (RS_adj_l[edge_l[i]].size() - 1)*(RS_adj_l[edge_l[i]].size() - 1);
            RS_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + RS_adj_r[edge_r[i]].size()*RS_adj_r[edge_r[i]].size() - (RS_adj_r[edge_r[i]].size() - 1)*(RS_adj_r[edge_r[i]].size() - 1);
            
        } else if (dis(gen) < (double)Reservoir_size/currentEdges) {
            int del_pos = rand() % Reservoir_size;
            uint32_t del_l = RS_sample_edge_l[del_pos];
            uint32_t del_r = RS_sample_edge_r[del_pos];
            RS_sample_edge_l[del_pos] = edge_l[i];
            RS_sample_edge_r[del_pos] = edge_r[i];

            // uint32_t del_l = rand() % RS_adj_l.size();
            
            // while (RS_adj_l[del_l].size() <= 0) {
            //     del_l = rand() % RS_adj_l.size();
            // }
            
            // uint32_t del_r = RS_adj_l[del_l][rand() % RS_adj_l[del_l].size()];

            auto& del_adj_l = RS_adj_l[del_l];
            auto& del_adj_r = RS_adj_r[del_r];

            std::swap(del_adj_l.back(), *std::find(del_adj_l.begin(), del_adj_l.end(), del_r));
            del_adj_l.pop_back();
            sum_degree_l = sum_degree_l + RS_adj_l[del_l].size()*RS_adj_l[del_l].size() - (RS_adj_l[del_l].size() + 1)*(RS_adj_l[del_l].size() + 1);

            std::swap(del_adj_r.back(), *std::find(del_adj_r.begin(), del_adj_r.end(), del_l));
            del_adj_r.pop_back();
            sum_degree_r = sum_degree_r + RS_adj_r[del_r].size()*RS_adj_r[del_r].size() - (RS_adj_r[del_r].size() + 1)*(RS_adj_r[del_r].size() + 1);

            RS_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + RS_adj_l[edge_l[i]].size()*RS_adj_l[edge_l[i]].size() - (RS_adj_l[edge_l[i]].size() - 1)*(RS_adj_l[edge_l[i]].size() - 1);
            RS_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + RS_adj_r[edge_r[i]].size()*RS_adj_r[edge_r[i]].size() - (RS_adj_r[edge_r[i]].size() - 1)*(RS_adj_r[edge_r[i]].size() - 1);
        }
    }

    return (uint64_t)butterfly_num;
}

// uint64_t Graph::Reservoir_sampling(uint32_t BUCKET_BITS) {
//     uint64_t butterfly_num = 0;
//     this->Reservoir_size = pow(2, BUCKET_BITS);
//     double currentEdges = 0;
//     uint32_t increment = 0;
//     double y = 0;
//     double Pr = 0;
//    uint64_t sum_degree_l = 0, sum_degree_r = 0;

//    srand(static_cast<unsigned int>(time(0)));
//    std::random_device rd;
//    std::mt19937 gen(rd());

//    RS_sample_edge_l.reserve(Reservoir_size);
//    RS_sample_edge_r.reserve(Reservoir_size);

//    std::uniform_real_distribution<> dis(0.0, 1.0);

//    for (uint32_t i = 0; i < m; i++) {
//        if (i % 10000 == 0) cout << i << endl;
//        currentEdges++;

//        if(currentEdges < Reservoir_size) y = currentEdges;
//        else y = (double)Reservoir_size;
//        Pr = (y/currentEdges) * ((y-1)/(currentEdges-1)) * ((y-2)/(currentEdges-2));
//        increment = round(1/Pr);

//        bool use_left_side = (sum_degree_l < sum_degree_r);
//        auto& RS_adj_v = use_left_side ? RS_adj_r : RS_adj_l;
//        auto& RS_adj_u = use_left_side ? RS_adj_l : RS_adj_r;
//        uint32_t v = use_left_side ? edge_r[i] : edge_l[i];
//        uint32_t u = use_left_side ? edge_l[i] : edge_r[i];

//        for (uint32_t w : RS_adj_u[u]) {
//            if (w == v) continue;
//            for (uint32_t w2 : RS_adj_v[w]) {
//                if (w2 == u) continue;
//                if (std::find(RS_adj_v[v].begin(), RS_adj_v[v].end(), w2) != RS_adj_v[v].end()) {
//                    butterfly_num += increment;
//                }
//            }
//        }
       
// //        for (uint32_t w : RS_adj_u[u]) {
// //                if (w == v) continue;  // 跳过与自己相连的边
// //                for (uint32_t w2 : RS_adj_v[w]) {
// //                    if (w2 == u || w2 <= v) continue;  // 避免重复计数
// //                    if (std::find(RS_adj_v[v].begin(), RS_adj_v[v].end(), w2) != RS_adj_v[v].end()) {
// //                        butterfly_num += increment;
// //                    }
// //                }
// //            }

//        if (currentEdges <= Reservoir_size) {
//            RS_sample_edge_l.push_back(edge_l[i]);
//            RS_sample_edge_r.push_back(edge_r[i]);
//            RS_adj_l[edge_l[i]].push_back(edge_r[i]);
//            sum_degree_l = sum_degree_l + RS_adj_l[edge_l[i]].size()*RS_adj_l[edge_l[i]].size() - (RS_adj_l[edge_l[i]].size() - 1)*(RS_adj_l[edge_l[i]].size() - 1);
//            RS_adj_r[edge_r[i]].push_back(edge_l[i]);
//            sum_degree_r = sum_degree_r + RS_adj_r[edge_r[i]].size()*RS_adj_r[edge_r[i]].size() - (RS_adj_r[edge_r[i]].size() - 1)*(RS_adj_r[edge_r[i]].size() - 1);
//        } else if (dis(gen) < (double)Reservoir_size/currentEdges) {
//            int del_pos = rand() % Reservoir_size;
//            uint32_t del_l = RS_sample_edge_l[del_pos];
//            uint32_t del_r = RS_sample_edge_r[del_pos];
//            RS_sample_edge_l[del_pos] = edge_l[i];
//            RS_sample_edge_r[del_pos] = edge_r[i];

//            auto& del_adj_l = RS_adj_l[del_l];
//            auto& del_adj_r = RS_adj_r[del_r];

//            std::swap(del_adj_l.back(), *std::find(del_adj_l.begin(), del_adj_l.end(), del_r));
//            del_adj_l.pop_back();
//            sum_degree_l = sum_degree_l + RS_adj_l[del_l].size()*RS_adj_l[del_l].size() - (RS_adj_l[del_l].size() + 1)*(RS_adj_l[del_l].size() + 1);

//            std::swap(del_adj_r.back(), *std::find(del_adj_r.begin(), del_adj_r.end(), del_l));
//            del_adj_r.pop_back();
//            sum_degree_r = sum_degree_r + RS_adj_r[del_r].size()*RS_adj_r[del_r].size() - (RS_adj_r[del_r].size() + 1)*(RS_adj_r[del_r].size() + 1);

//            RS_adj_l[edge_l[i]].push_back(edge_r[i]);
//            sum_degree_l = sum_degree_l + RS_adj_l[edge_l[i]].size()*RS_adj_l[edge_l[i]].size() - (RS_adj_l[edge_l[i]].size() - 1)*(RS_adj_l[edge_l[i]].size() - 1);
//            RS_adj_r[edge_r[i]].push_back(edge_l[i]);
//            sum_degree_r = sum_degree_r + RS_adj_r[edge_r[i]].size()*RS_adj_r[edge_r[i]].size() - (RS_adj_r[edge_r[i]].size() - 1)*(RS_adj_r[edge_r[i]].size() - 1);
//        }
//    }

//    return butterfly_num;
// }

int Graph::intersectionSize(const vector<uint32_t>& vec1, const vector<uint32_t>& vec2) {
    unordered_set<uint32_t> set1(vec1.begin(), vec1.end());
    int count = 0;

    for (const auto& elem : vec2) {
        if (set1.find(elem) != set1.end()) {
            count++;
            set1.erase(elem);
        }
    }

    return count;
}

uint64_t Graph::Priority_sampling(uint32_t BUCKET_BITS, int is_hash){
    this->Reservoir_size = pow(2,BUCKET_BITS);
    double butterfly_num = 0;
    uint32_t currentEdges = 0;
    std::priority_queue<PS_edge, std::vector<PS_edge>, ComparePriority> pq;
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    uint32_t u, v;
    
    zStar = 0;
    
    for(uint32_t i = 0;i<l_n;i++){
        PS_adj_l[i].clear();
    }
    for(uint32_t i = 0;i<r_n;i++){
        PS_adj_r[i].clear();
        
    }
    
    uint32_t time_ = m/Time_point;
    
    /**
    生成随机数
     */
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    
    /**
    hash相关
     */
    random_seed = time(NULL);
    mt19937_64 generator_DEABC(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    unsigned long long c = dist2(generator_DEABC);
    unsigned long long d = dist2(generator_DEABC);
    
        
    for(uint32_t i = 0; i< m;i++){
        currentEdges++;
        
        //if(i%100000 == 0) cout<<i<<endl;
        
        
        bool use_left_side = (sum_degree_l < sum_degree_r);
        auto& PS_adj_v = use_left_side ? PS_adj_r : PS_adj_l;
        auto& PS_adj_u = use_left_side ? PS_adj_l : PS_adj_r;
        v = use_left_side ? edge_r[i] : edge_l[i];
        u = use_left_side ? edge_l[i] : edge_r[i];
        
        double q, q2;
        if(zStar == 0) {
            q = 1;
            q2 = 1;
        }
        else{
            q = 1-zStar;
            q2 = (double)(Reservoir_size - 3)/(double)Reservoir_size;
        }
        
        uint32_t increment = q2/(q*q*q);
        
        for(uint32_t w : PS_adj_u[u]){
            if(w == v) continue;
            for(uint32_t w2 : PS_adj_v[w]){
                if(w2 == u) continue;
                
                if(find(PS_adj_v[v].begin(),PS_adj_v[v].end(),w2) != PS_adj_v[v].end()){
                    butterfly_num += increment;
                }
                
            }
        }
        
        
        double random_value;
        
        if(is_hash == 0){
            random_value = EdgeHash_p(edge_l[i],edge_r[i],c,d);
            
            while (random_value == 0.0) {
                random_value = EdgeHash_p(edge_l[i],edge_r[i],c,d);
            }
        }else if(is_hash == 1){
            random_value = hash(edge_l[i], edge_r[i], c, d, gen);
            while (random_value == 0.0) {
                random_value = hash(edge_l[i], edge_r[i], c, d, gen);
            }
        }
        else{
            random_value = dis(gen);
            while (random_value == 0.0) {
                random_value = dis(gen);
            }
        }

        
                
        if(currentEdges < Reservoir_size){
            if(find(PS_adj_l[edge_l[i]].begin(),PS_adj_l[edge_l[i]].end(),edge_r[i]) == PS_adj_l[edge_l[i]].end()){
                PS_adj_l[edge_l[i]].push_back(edge_r[i]);
                sum_degree_l = sum_degree_l + PS_adj_l[edge_l[i]].size()*PS_adj_l[edge_l[i]].size() - (PS_adj_l[edge_l[i]].size() - 1)*(PS_adj_l[edge_l[i]].size() - 1);
                PS_adj_r[edge_r[i]].push_back(edge_l[i]);
                sum_degree_r = sum_degree_r + PS_adj_r[edge_r[i]].size()*PS_adj_r[edge_r[i]].size() - (PS_adj_r[edge_r[i]].size() - 1)*(PS_adj_r[edge_r[i]].size() - 1);
                pq.push(PS_edge(edge_l[i], edge_r[i], random_value));
            }
        }else{
            
            PS_edge min_p_edge = pq.top();
            if(random_value > min_p_edge.priority){
                if(find(PS_adj_l[edge_l[i]].begin(),PS_adj_l[edge_l[i]].end(),edge_r[i]) == PS_adj_l[edge_l[i]].end()){
                    PS_adj_l[edge_l[i]].push_back(edge_r[i]);
                    sum_degree_l = sum_degree_l + PS_adj_l[edge_l[i]].size()*PS_adj_l[edge_l[i]].size() - (PS_adj_l[edge_l[i]].size() - 1)*(PS_adj_l[edge_l[i]].size() - 1);
                    PS_adj_r[edge_r[i]].push_back(edge_l[i]);
                    sum_degree_r = sum_degree_r + PS_adj_r[edge_r[i]].size()*PS_adj_r[edge_r[i]].size() - (PS_adj_r[edge_r[i]].size() - 1)*(PS_adj_r[edge_r[i]].size() - 1);
                    pq.push(PS_edge(edge_l[i], edge_r[i], random_value));
                
                    pq.pop();
                
                    auto& del_adj_l = PS_adj_l[min_p_edge.edge_l];
                    auto& del_adj_r = PS_adj_r[min_p_edge.edge_r];
                    std::swap(del_adj_l.back(), *std::find(del_adj_l.begin(), del_adj_l.end(), min_p_edge.edge_r));
                    del_adj_l.pop_back();
                    sum_degree_l = sum_degree_l + PS_adj_l[min_p_edge.edge_l].size()*PS_adj_l[min_p_edge.edge_l].size() - (PS_adj_l[edge_l[i]].size() + 1)*(PS_adj_l[edge_l[i]].size() + 1);
                
                    std::swap(del_adj_r.back(), *std::find(del_adj_r.begin(), del_adj_r.end(), min_p_edge.edge_l));
                    del_adj_r.pop_back();
                    sum_degree_r = sum_degree_r + PS_adj_r[min_p_edge.edge_r].size()*PS_adj_r[min_p_edge.edge_r].size() - (PS_adj_r[min_p_edge.edge_r].size() + 1)*(PS_adj_r[min_p_edge.edge_r].size() + 1);
                
                    zStar = max(zStar, min_p_edge.priority);
                }
                
            }
            
        }
        
                
        
    }
    
    return butterfly_num;
}

// uint64_t Graph::Priority_sampling_DEABC0(uint32_t BUCKET_BITS){
//     this->Reservoir_size = pow(2,BUCKET_BITS);
//     double butterfly_num = 0;
//     uint32_t currentEdges = 0;
//     random_seed = time(NULL);
    
//     uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
//     mt19937_64 generator_DEABC(random_seed);
//     mt19937 gen(random_seed);
//     uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    
//     priority_queue<PS_edge, std::vector<PS_edge>, ComparePriority_min> buffer;
    
//     unsigned long long c = dist2(generator_DEABC);
//     unsigned long long d = dist2(generator_DEABC);
    
// //    vector<vector<uint32_t> > DEABC_adj_v, DEABC_adj_u;
//     uint32_t u, v;
//     double *count_u, *count_v;
    
//     exactcnt = 1;
    
//     for(uint32_t  i = 0; i< m; i++){
//         currentEdges++;
        
//         //if(i%100000 == 0) cout<<i<<endl;
        
//         double edge_rank = hash(edge_l[i],edge_r[i],c,d, gen);
// //        double edge_rank = EdgeHash_p(edge_l[i],edge_r[i],c,d);
        
//         bool sampled = 0;
        
//         if(find(DEABC_adj_l[edge_l[i]].begin(),DEABC_adj_l[edge_l[i]].end(),edge_r[i]) == DEABC_adj_l[edge_l[i]].end()){
//             if(currentEdges < Reservoir_size){
//                 buffer.push(PS_edge(edge_l[i],edge_r[i],edge_rank));
                
//                 DEABC_adj_l[edge_l[i]].push_back(edge_r[i]);
//                 sum_degree_l = sum_degree_l + DEABC_adj_l[edge_l[i]].size()*DEABC_adj_l[edge_l[i]].size() - (DEABC_adj_l[edge_l[i]].size() - 1)*(DEABC_adj_l[edge_l[i]].size() - 1);
//                 DEABC_adj_r[edge_r[i]].push_back(edge_l[i]);
//                 sum_degree_r = sum_degree_r + DEABC_adj_r[edge_r[i]].size()*DEABC_adj_r[edge_r[i]].size() - (DEABC_adj_r[edge_r[i]].size() - 1)*(DEABC_adj_r[edge_r[i]].size() - 1);
//                 sampled = 1;
                
//             }else{
//                 if(exactcnt){
//                     exactcnt = 0;
//                 }
//                 PS_edge max_p_edge = buffer.top();
//                 if(edge_rank < max_p_edge.priority){
                    
//                     auto& del_adj_l = DEABC_adj_l[max_p_edge.edge_l];
//                     auto& del_adj_r = DEABC_adj_r[max_p_edge.edge_r];
//                     std::swap(del_adj_l.back(), *std::find(del_adj_l.begin(), del_adj_l.end(), max_p_edge.edge_r));
//                     del_adj_l.pop_back();
//                     sum_degree_l = sum_degree_l + DEABC_adj_l[max_p_edge.edge_l].size()*DEABC_adj_l[max_p_edge.edge_l].size() - (DEABC_adj_l[max_p_edge.edge_l].size() + 1)*(DEABC_adj_l[max_p_edge.edge_l].size() + 1);
                    
//                     std::swap(del_adj_r.back(), *std::find(del_adj_r.begin(), del_adj_r.end(), max_p_edge.edge_l));
//                     del_adj_r.pop_back();
//                     sum_degree_r = sum_degree_r + DEABC_adj_r[max_p_edge.edge_r].size()*DEABC_adj_r[max_p_edge.edge_r].size() - (DEABC_adj_r[max_p_edge.edge_r].size() - 1)*(DEABC_adj_r[max_p_edge.edge_r].size() - 1);
                    
//                     buffer.pop();
                    
//                     buffer.push(PS_edge(edge_l[i],edge_r[i],edge_rank));
                    
//                     DEABC_adj_l[edge_l[i]].push_back(edge_r[i]);
//                     sum_degree_l = sum_degree_l + DEABC_adj_l[edge_l[i]].size()*DEABC_adj_l[edge_l[i]].size() - (DEABC_adj_l[edge_l[i]].size() - 1)*(DEABC_adj_l[edge_l[i]].size() - 1);
//                     DEABC_adj_r[edge_r[i]].push_back(edge_l[i]);
//                     sum_degree_r = sum_degree_r + DEABC_adj_r[edge_r[i]].size()*DEABC_adj_r[edge_r[i]].size() - (DEABC_adj_r[edge_r[i]].size() - 1)*(DEABC_adj_r[edge_r[i]].size() - 1);
//                     sampled = 1;
                    
//                 }
//             }
            
//             if(exactcnt){
                
//                 bool use_left_side = (sum_degree_l < sum_degree_r);
//                 auto& DEABC_adj_v = use_left_side ? DEABC_adj_r : DEABC_adj_l;
//                 auto& DEABC_adj_u = use_left_side ? DEABC_adj_l : DEABC_adj_r;
//                 v = use_left_side ? edge_r[i] : edge_l[i];
//                 u = use_left_side ? edge_l[i] : edge_r[i];
                
//                 for(uint32_t w : DEABC_adj_u[u]){
//                     if(w == v) continue;
//                     double weightSum_w = 0;
//                     for(uint32_t w2 : DEABC_adj_v[w]){
//                         if(w2 == u) continue;
                        
//                         if(find(DEABC_adj_v[v].begin(),DEABC_adj_v[v].end(),w2) != DEABC_adj_v[v].end()){
//                             butterfly_num++;
//                         }
                        
//                     }
//                 }
                
                    
//             }else{
//                 double qT = 0;
//                 if(sampled){
//                     qT = (((double)Reservoir_size - 4.0) / (double)Reservoir_size) / pow(buffer.top().priority, 4.0);
                    
//                     bool use_left_side = (sum_degree_l < sum_degree_r);
//                     auto& DEABC_adj_v = use_left_side ? DEABC_adj_r : DEABC_adj_l;
//                     auto& DEABC_adj_u = use_left_side ? DEABC_adj_l : DEABC_adj_r;
//                     v = use_left_side ? edge_r[i] : edge_l[i];
//                     u = use_left_side ? edge_l[i] : edge_r[i];

//                     for(uint32_t w : DEABC_adj_u[u]){
//                         if(w == v) continue;
//                         double weightSum_w = 0;
//                         for(uint32_t w2 : DEABC_adj_v[w]){
//                             if(w2 == u) continue;
                            
//                             if(find(DEABC_adj_v[v].begin(),DEABC_adj_v[v].end(),w2) != DEABC_adj_v[v].end()){
//                                 butterfly_num+=qT;
//                             }
                            
//                         }
                        
//                     }
                    
//                 }
//             }
//         }
        
        
//     }
    
//     return (uint64_t)butterfly_num;
// }

uint64_t Graph::Priority_sampling_DEABC0(uint32_t BUCKET_BITS){
    this->Reservoir_size = pow(2,BUCKET_BITS);
    double butterfly_num = 0;
    uint32_t currentEdges = 0;
    random_seed = time(NULL);
    uint32_t time_num = m/Time_point;
    
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
    mt19937_64 generator_DEABC(random_seed);
    mt19937 gen(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    
    priority_queue<PS_edge, std::vector<PS_edge>, ComparePriority_min> buffer;
    
    unsigned long long c = dist2(generator_DEABC);
    unsigned long long d = dist2(generator_DEABC);
    
    uint32_t u, v;
    
    // memset(hashmap, 0, sizeof(uint32_t) * std::max(l_n, r_n));
    
    for(uint32_t  i = 0; i< m; i++){
        currentEdges++;

        // if (i % time_num == 0) cout << (uint64_t)butterfly_num << endl;
        
        // if(i%10000 == 0) cout<<i<<endl;
        
        // double edge_rank = hash(edge_l[i],edge_r[i],c,d, gen);
        // double edge_rank = EdgeHash_p(edge_l[i],edge_r[i],c,d);
        double edge_rank = hash(edge_l[i],edge_r[i],c,d);
        
        
        
        if(currentEdges < Reservoir_size){
            if(find(DEABC_adj_l[edge_l[i]].begin(),DEABC_adj_l[edge_l[i]].end(),edge_r[i]) == DEABC_adj_l[edge_l[i]].end()){
                buffer.push(PS_edge(edge_l[i],edge_r[i],edge_rank));
                
                DEABC_adj_l[edge_l[i]].push_back(edge_r[i]);
                sum_degree_l = sum_degree_l + DEABC_adj_l[edge_l[i]].size()*DEABC_adj_l[edge_l[i]].size() - (DEABC_adj_l[edge_l[i]].size() - 1)*(DEABC_adj_l[edge_l[i]].size() - 1);
                DEABC_adj_r[edge_r[i]].push_back(edge_l[i]);
                sum_degree_r = sum_degree_r + DEABC_adj_r[edge_r[i]].size()*DEABC_adj_r[edge_r[i]].size() - (DEABC_adj_r[edge_r[i]].size() - 1)*(DEABC_adj_r[edge_r[i]].size() - 1);
                    

                bool use_left_side = (sum_degree_l < sum_degree_r);
                auto& DEABC_adj_v = use_left_side ? DEABC_adj_r : DEABC_adj_l;
                auto& DEABC_adj_u = use_left_side ? DEABC_adj_l : DEABC_adj_r;
                v = use_left_side ? edge_r[i] : edge_l[i];
                u = use_left_side ? edge_l[i] : edge_r[i];

                for (uint32_t w : DEABC_adj_v[v]) {
                    hashmap[w] = i+1;
                }
                for (uint32_t w : DEABC_adj_u[u]) {
                    if (w == v) continue;
                    for (uint32_t w2 : DEABC_adj_v[w]) {
                        if (w2 == u) continue;
                
                        if (hashmap[w2] == i+1) {
                            butterfly_num += 1;
                        }
                    }
                }
                
                // for(uint32_t w : DEABC_adj_u[u]){
                //     if(w == v) continue;
                //     double weightSum_w = 0;
                //     for(uint32_t w2 : DEABC_adj_v[w]){
                //         if(w2 == u) continue;
                        
                //         if(find(DEABC_adj_v[v].begin(),DEABC_adj_v[v].end(),w2) != DEABC_adj_v[v].end()){
                //             butterfly_num++;
                //         }
                        
                //     }
                // }
                
            }
                
        }else{
            
            PS_edge max_p_edge = buffer.top();
            // PS_edge max_p_edge = PS_edge(edge_l[i],edge_r[i],1);
            if(edge_rank < max_p_edge.priority){
                if(find(DEABC_adj_l[edge_l[i]].begin(),DEABC_adj_l[edge_l[i]].end(),edge_r[i]) == DEABC_adj_l[edge_l[i]].end()){
                    
                    auto& del_adj_l = DEABC_adj_l[max_p_edge.edge_l];
                    auto& del_adj_r = DEABC_adj_r[max_p_edge.edge_r];
                    std::swap(del_adj_l.back(), *std::find(del_adj_l.begin(), del_adj_l.end(), max_p_edge.edge_r));
                    del_adj_l.pop_back();
                    sum_degree_l = sum_degree_l + DEABC_adj_l[max_p_edge.edge_l].size()*DEABC_adj_l[max_p_edge.edge_l].size() - (DEABC_adj_l[max_p_edge.edge_l].size() + 1)*(DEABC_adj_l[max_p_edge.edge_l].size() + 1);
                    
                    std::swap(del_adj_r.back(), *std::find(del_adj_r.begin(), del_adj_r.end(), max_p_edge.edge_l));
                    del_adj_r.pop_back();
                    sum_degree_r = sum_degree_r + DEABC_adj_r[max_p_edge.edge_r].size()*DEABC_adj_r[max_p_edge.edge_r].size() - (DEABC_adj_r[max_p_edge.edge_r].size() - 1)*(DEABC_adj_r[max_p_edge.edge_r].size() - 1);
                
                    buffer.pop();
                    
                    buffer.push(PS_edge(edge_l[i],edge_r[i],edge_rank));
                    
                    DEABC_adj_l[edge_l[i]].push_back(edge_r[i]);
                    sum_degree_l = sum_degree_l + DEABC_adj_l[edge_l[i]].size()*DEABC_adj_l[edge_l[i]].size() - (DEABC_adj_l[edge_l[i]].size() - 1)*(DEABC_adj_l[edge_l[i]].size() - 1);
                    DEABC_adj_r[edge_r[i]].push_back(edge_l[i]);
                    sum_degree_r = sum_degree_r + DEABC_adj_r[edge_r[i]].size()*DEABC_adj_r[edge_r[i]].size() - (DEABC_adj_r[edge_r[i]].size() - 1)*(DEABC_adj_r[edge_r[i]].size() - 1);
                    
                    double qT = (((double)Reservoir_size - 4.0) / (double)Reservoir_size) / pow(buffer.top().priority, 4.0);
                    
                    bool use_left_side = (sum_degree_l < sum_degree_r);
                    auto& DEABC_adj_v = use_left_side ? DEABC_adj_r : DEABC_adj_l;
                    auto& DEABC_adj_u = use_left_side ? DEABC_adj_l : DEABC_adj_r;
                    v = use_left_side ? edge_r[i] : edge_l[i];
                    u = use_left_side ? edge_l[i] : edge_r[i];

                    for (uint32_t w : DEABC_adj_v[v]) {
                        hashmap[w] = i+1;
                    }
                    for (uint32_t w : DEABC_adj_u[u]) {
                        if (w == v) continue;
                        for (uint32_t w2 : DEABC_adj_v[w]) {
                            if (w2 == u) continue;
                
                            if (hashmap[w2] == i+1) {
                                butterfly_num += qT;
                            }
                        }
                    }

                    // for(uint32_t w : DEABC_adj_u[u]){
                    //     if(w == v) continue;
                    //     double weightSum_w = 0;
                    //     for(uint32_t w2 : DEABC_adj_v[w]){
                    //         if(w2 == u) continue;
                            
                    //         if(find(DEABC_adj_v[v].begin(),DEABC_adj_v[v].end(),w2) != DEABC_adj_v[v].end()){
                    //             butterfly_num+=qT;
                    //         }
                            
                    //     }
                        
                    // }
                }
            }
        }
    }
    return (uint64_t)butterfly_num;
}


double Graph::EdgeHash_p(uint32_t u, uint32_t v, unsigned long long a1, unsigned long long b1) {
    uint64_t edge = (static_cast<uint64_t>(hash_func(u ^ arr[random_seed % 8])) << 32) + 
                static_cast<uint64_t>(hash_func(v ^ arr[(random_seed + 3) % 8]));
    return double((((a1*(edge%Prime1)+b1)%Prime1)%Prime4))/Prime4;
}

uint32_t Graph::hash_func(uint32_t a)
{
    a = (a+0x479ab41d) + (a<<8);
    a = (a^0xe4aa10ce) ^ (a>>5);
    a = (a+0x9942f0a6) - (a<<14);
    a = (a^0x5aedd67d) ^ (a>>3);
    a = (a+0x17bea992) + (a<<7);
    return a;
}

double Graph::hash(uint32_t a, uint32_t b, uint32_t c, uint32_t d, mt19937 &generator) {
    // 使用std::hash函数生成哈希值
    size_t hashValue = a ^ (b << 1) * c ^ (d >> 1);
    
    // 创建一个以哈希值为种子的均匀分布的随机数生成器
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    generator.seed(hashValue);
    return distribution(generator);

}

int Graph::hash_bucket(uint32_t a, uint32_t b, uint32_t c, uint32_t d, mt19937 generator){
    size_t hashValue = a ^ (b << 1) * c ^ (d >> 1);
    
    // 创建一个以哈希值为种子的均匀分布的随机数生成器
    std::uniform_int_distribution<int> distribution(0, Reservoir_size-1);
    generator.seed(hashValue);
    return distribution(generator);
}


double Graph::getConstant(uint32_t buckets) {
    switch (buckets) {
        case 16:
            return 0.673;
        case 32:
            return 0.697;
        case 64:
            return 0.709;
        default:
            return (0.7213 / (1 + 1.079 / buckets));
//                return 0.85;
    }
}

double Graph::GetCounts_DEABC_PLUS(uint32_t v_l, uint32_t v_r, uint64_t &sum_degree_l, uint64_t &sum_degree_r, uint32_t i){
    
    
    bool use_left_side = (sum_degree_l < sum_degree_r);
    auto& DEABC_PLUS_adj_v = use_left_side ? DEABC_PLUS_adj_r : DEABC_PLUS_adj_l;
    auto& DEABC_PLUS_adj_u = use_left_side ? DEABC_PLUS_adj_l : DEABC_PLUS_adj_r;
    int v = use_left_side ? v_r : v_l;
    int u = use_left_side ? v_l : v_r;
    
    double count_num = 0;

    for (uint32_t w : DEABC_PLUS_adj_v[v]) {
        hashmap[w] = i+1;
    }
    for (uint32_t w : DEABC_PLUS_adj_u[u]) {
        if (w == v) continue;
        for (uint32_t w2 : DEABC_PLUS_adj_v[w]) {
            if (w2 == u) continue;
                
            if (hashmap[w2] == i+1) {
                count_num ++;
            }
        }
    }
    
    // for(uint32_t w : DEABC_PLUS_adj_u[u]){
    //     if(w == v) continue;
    //     for(uint32_t w2 : DEABC_PLUS_adj_v[w]){
    //         if(w2 == u) continue;
            
    //         if(find(DEABC_PLUS_adj_v[v].begin(),DEABC_PLUS_adj_v[v].end(),w2) != DEABC_PLUS_adj_v[v].end()){
    //             count_num++;
    //         }
    //     }
    // }
    return count_num;
}


uint64_t Graph::Priority_sampling_DEABC_PLUS(uint32_t BUCKET_BITS){
    this->Reservoir_size = pow(2,BUCKET_BITS);
    double butterfly_num = 0;
    uint32_t currentEdges = 0;
    random_seed = time(NULL);

    uint32_t time_num = m/Time_point;
    
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
    mt19937_64 generator_DEABC(random_seed);
    mt19937 gen(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    unsigned long long c = dist2(generator_DEABC);
    unsigned long long d = dist2(generator_DEABC);
    
    uint32_t num_in_S = 0;
    uint32_t u, v;
    
    double q = 1;
    // S_DEABC_PLUS.clear();
    // S_DEABC_PLUS.resize(Reservoir_size, PS_edge(0, 0, -1.0));
    
    uint32_t g_max;
    uint32_t y;
    double n_heat = 0, uu;

    double min_p = 1;

    uint32_t non_empty_bucket_num = 0;

    // memset(hashmap, 0, sizeof(uint32_t) * std::max(l_n, r_n));

    auto start = std::chrono::high_resolution_clock::now();
    
    
    for(uint32_t  i = 0; i< m; i++){
        currentEdges++;
        g_max=0;
        // cout<<i<<endl;

        // if(i%100000 == 0){
            // cout<<i<<endl;
        //     std::cout<<i<<" "<<(uint64_t)n_heat<<std::endl;

        // } 


        // double edge_rank = hash(edge_l[i],edge_r[i], c, d, gen);
        // uint32_t edge_bucket = hash_bucket(edge_l[i],edge_r[i], c, d, gen);

        // if (i % time_num == 0) cout << (uint64_t)butterfly_num << endl;

        double edge_rank = hash(edge_l[i],edge_r[i], c, d);
        uint32_t edge_bucket = hash_bucket(edge_l[i],edge_r[i], c, d, this->Reservoir_size);
        
        double p;

        
        // if(n_heat <= 0 || non_empty_bucket_num<=2){
        //     p = 1;
        // }else{
        //     p = ((double)non_empty_bucket_num/(double)n_heat)*((double)(non_empty_bucket_num-1)/((double)n_heat-1)) *((double)(non_empty_bucket_num-2)/((double)n_heat-2));
        // }
        
        // if(find(DEABC_PLUS_adj_l[edge_l[i]].begin(),DEABC_PLUS_adj_l[edge_l[i]].end(),edge_r[i]) == DEABC_PLUS_adj_l[edge_l[i]].end()){
        //     if(n_heat > 0 && dis(gen)*(double)(n_heat - non_empty_bucket_num)/(double)n_heat < ((double)n_heat/(double)currentEdges)){
        //         butterfly_num += (GetCounts_DEABC_PLUS(edge_l[i], edge_r[i], sum_degree_l, sum_degree_r, i))/p;
        //     }
        // }

        
        if(S_DEABC_PLUS[edge_bucket].priority == -1 || edge_rank < S_DEABC_PLUS[edge_bucket].priority){

            if(n_heat <= 0 || non_empty_bucket_num<=3){
                p = 1;
            }else{
                p = ((double)non_empty_bucket_num/(double)n_heat)*((double)(non_empty_bucket_num-1)/((double)n_heat-1)) *((double)(non_empty_bucket_num-2)/((double)n_heat-2)) *((double)(non_empty_bucket_num-3)/((double)n_heat-3));
            }

            // cout<<p<<endl;
        
            // if(find(DEABC_PLUS_adj_l[edge_l[i]].begin(),DEABC_PLUS_adj_l[edge_l[i]].end(),edge_r[i]) == DEABC_PLUS_adj_l[edge_l[i]].end()){
                // if(n_heat > 0 && dis(gen)*(double)(n_heat - non_empty_bucket_num)/(double)n_heat < ((double)n_heat/(double)currentEdges)){
                butterfly_num += (GetCounts_DEABC_PLUS(edge_l[i], edge_r[i], sum_degree_l, sum_degree_r, i))/p;
            // }


            

            if(S_DEABC_PLUS[edge_bucket].priority != -1){
                g_max = static_cast<int>(floor(-log2(S_DEABC_PLUS[edge_bucket].priority)));

                
                
                auto& del_adj_l = DEABC_PLUS_adj_l[S_DEABC_PLUS[edge_bucket].edge_l];
                auto& del_adj_r = DEABC_PLUS_adj_r[S_DEABC_PLUS[edge_bucket].edge_r];
                std::swap(del_adj_l.back(), *std::find(del_adj_l.begin(), del_adj_l.end(), S_DEABC_PLUS[edge_bucket].edge_r));
                del_adj_l.pop_back();
                sum_degree_l = sum_degree_l + DEABC_PLUS_adj_l[S_DEABC_PLUS[edge_bucket].edge_l].size()*DEABC_PLUS_adj_l[S_DEABC_PLUS[edge_bucket].edge_l].size() - (DEABC_PLUS_adj_l[S_DEABC_PLUS[edge_bucket].edge_l].size() + 1)*(DEABC_PLUS_adj_l[S_DEABC_PLUS[edge_bucket].edge_l].size() + 1);
                
                
                
                std::swap(del_adj_r.back(), *std::find(del_adj_r.begin(), del_adj_r.end(), S_DEABC_PLUS[edge_bucket].edge_l));
                del_adj_r.pop_back();
                sum_degree_r = sum_degree_r + DEABC_PLUS_adj_r[S_DEABC_PLUS[edge_bucket].edge_r].size()*DEABC_PLUS_adj_r[S_DEABC_PLUS[edge_bucket].edge_r].size() - (DEABC_PLUS_adj_r[S_DEABC_PLUS[edge_bucket].edge_r].size() - 1)*(DEABC_PLUS_adj_r[S_DEABC_PLUS[edge_bucket].edge_r].size() - 1);
            }else non_empty_bucket_num++;

            

            y = static_cast<int>(floor(-log2(edge_rank)));

            
            
            if(y > g_max){
                n_heat += 2.0/q;
                double uu1 = 1.0/pow(2.0, y);
                double uu2 = 1.0/pow(2.0, g_max);
                uu = uu1 - uu2;
                q += (1.0 / Reservoir_size) * uu;
            }

            

            S_DEABC_PLUS[edge_bucket] = PS_edge(edge_l[i], edge_r[i], edge_rank);

            DEABC_PLUS_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + DEABC_PLUS_adj_l[edge_l[i]].size()*DEABC_PLUS_adj_l[edge_l[i]].size() - (DEABC_PLUS_adj_l[edge_l[i]].size() - 1)*(DEABC_PLUS_adj_l[edge_l[i]].size() - 1);
            DEABC_PLUS_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + DEABC_PLUS_adj_r[edge_r[i]].size()*DEABC_PLUS_adj_r[edge_r[i]].size() - (DEABC_PLUS_adj_r[edge_r[i]].size() - 1)*(DEABC_PLUS_adj_r[edge_r[i]].size() - 1);

            // cout<<3<<endl;
        }
        
        
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Function execution time: " << duration.count() << std::endl;

    // std::cout<<m<<" "<<(uint64_t)n_heat<<std::endl;
    
    return (uint64_t)butterfly_num;
}



double Graph::hash_to_double(uint32_t a, uint32_t b){
    double da = static_cast<double>(a);
    double db = static_cast<double>(b);
        
        // 使用 sin 和 cos 函数生成一个近似均匀分布的值
    double result = fabs(sin(da) * cos(db));
        
        // 确保结果在 0 到 1 范围内
    return fmod(result, 1.0);

}

// MurmurHash3 生成64位哈希值
uint64_t Graph::murmur_hash_64(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    uint64_t hash[2];  // MurmurHash3_x64_128 的输出是128位，使用前64位
    uint32_t data[4] = {a, b, c, d};

    // 调用MurmurHash3_x64_128生成哈希
    MurmurHash3_x64_128(data, sizeof(data), 0, hash);

    return hash[0];  // 只使用前64位哈希值
}

// 将MurmurHash的输出映射到[0, 1)
double Graph::hash(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    uint64_t hash_value = murmur_hash_64(a, b, c, d);

    // 将哈希值映射到[0, 1)
    return static_cast<double>(hash_value) / static_cast<double>(UINT64_MAX);
}

// 将哈希结果映射到 [1, range] 范围的整数
uint32_t Graph::hash_bucket(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t range) {
    uint64_t hash_value = murmur_hash_64(a, b, c, d);

    // 将哈希值映射到 [1, range]，通过取模运算并加1
    return (hash_value % range);
}

void Graph::init(){
    memset(hashmap, 0, sizeof(uint32_t) * std::max(this->l_n, this->r_n));
    for(uint32_t i = 0;i<l_n;i++){
        DEABC_adj_l[i].clear();
        DEABC_PLUS_adj_l[i].clear();
        DEABC_adj_l[i].reserve(degree_l[i]);
        DEABC_PLUS_adj_l[i].reserve(degree_l[i]);
    }
    for(uint32_t i = 0;i<r_n;i++){
        DEABC_adj_r[i].clear();
        DEABC_PLUS_adj_r[i].clear();
        DEABC_adj_r[i].reserve(degree_r[i]);
        DEABC_PLUS_adj_r[i].reserve(degree_r[i]);
    }
    S_DEABC_PLUS.clear();
    S_DEABC_PLUS.resize(Reservoir_size, PS_edge(0, 0, -1.0));
    
}