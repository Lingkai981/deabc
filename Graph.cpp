//
//  Graph.cpp
//  bc_ABACUS
//
//  Created by LingkaiMeng on 2024/5/18.
//

#include "Graph.hpp"

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

Graph::Graph(const char *_dir){
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
    exactcnt = 1;
    TM = 0;
    time_ = 0;
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

    
}



void Graph::readGraph(){
    ifstream infile;   //输入流
    
    infile.open(str, ios::in);
    if (!infile.is_open()){
        cout<<"Open degree file failure"<<endl;
        exit(0);
    }
    infile>>l_n>>r_n>>m;
    
    cout<<"l_n:"<<l_n<<" r_n:"<<r_n<<" m:"<<m<<endl;
    
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
    RS_adj_l.resize(l_n);
    RS_adj_r.resize(r_n);
    PS_adj_l.resize(l_n);
    PS_adj_r.resize(r_n);
    FURL_adj_l.resize(l_n);
    FURL_adj_r.resize(r_n);
    PartitionCT_adj_l.resize(l_n);
    PartitionCT_adj_r.resize(r_n);

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
        PS_adj_l[i].reserve(degree_l[i]);
        FURL_adj_l[i].reserve(degree_l[i]);
        PartitionCT_adj_l[i].reserve(degree_l[i]);
    }
    for(uint32_t i = 0;i<r_n;i++){
        start_r[i+1] = start_r[i] + degree_r[i];
        RS_adj_r[i].reserve(degree_r[i]);
        PS_adj_r[i].reserve(degree_r[i]);
        FURL_adj_r[i].reserve(degree_r[i]);
        PartitionCT_adj_r[i].reserve(degree_r[i]);
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

//uint64_t Graph::Reservoir_sampling(uint32_t BUCKET_BITS){
//    uint64_t butterfly_num = 0;
//    this->Reservoir_size = pow(2,BUCKET_BITS);
//    double currentEdges = 0;
//    uint32_t increment = 0;
//    double y = 0;
//    double Pr = 0;
//    uint64_t sum_degree_l = 0, sum_degree_r = 0;
//    vector<vector<uint32_t> > RS_adj_v, RS_adj_u;
//    uint32_t u, v, del_l, del_r;
//    
//    srand(static_cast<unsigned int>(time(0)));
//    
//    std::random_device rd; // 获取随机设备
//    std::mt19937 gen(rd()); // 以随机设备生成随机种子
//    
//    RS_sample_edge_l.reserve(Reservoir_size);
//    RS_sample_edge_r.reserve(Reservoir_size);
//    
//    for(uint32_t i = 0; i < m; i++){
//        if(i%100000 == 0) cout<<i<<endl;
//        currentEdges++;
//        if(currentEdges < Reservoir_size) y = currentEdges;
//        else y = (double)Reservoir_size;
//        Pr = (y/currentEdges) * ((y-1)/(currentEdges-1)) * ((y-2)/(currentEdges-2));
//        increment = round(1/Pr);
//        
//        if(sum_degree_l < sum_degree_r){
//            RS_adj_v = RS_adj_r;
//            RS_adj_u = RS_adj_l;
//            v = edge_r[i];
//            u = edge_l[i];
//        }
//            
//        else{
//            RS_adj_v = RS_adj_l;
//            RS_adj_u = RS_adj_r;
//            v = edge_l[i];
//            u = edge_r[i];
//        }
//        
//        for(uint32_t w : RS_adj_u[u]){
//            if(w == v) continue;
//            for(uint32_t w2 : RS_adj_v[w]){
//                if(w2 == u) continue;
//                
//                if(find(RS_adj_v[v].begin(),RS_adj_v[v].end(),w2) != RS_adj_v[v].end()){
//                    butterfly_num += increment;
//                }
//                
//            }
//        }
//        
//        bernoulli_distribution bd((double)Reservoir_size/currentEdges);
//        
//        if((uint32_t)currentEdges <= Reservoir_size){
//            RS_sample_edge_l.push_back(edge_l[i]);
//            RS_sample_edge_r.push_back(edge_r[i]);
//            
//            RS_adj_l[edge_l[i]].push_back(edge_r[i]);
//            sum_degree_l = sum_degree_l + RS_adj_l[edge_l[i]].size()*RS_adj_l[edge_l[i]].size() - (RS_adj_l[edge_l[i]].size() - 1)*(RS_adj_l[edge_l[i]].size() - 1);
//            RS_adj_r[edge_r[i]].push_back(edge_l[i]);
//            sum_degree_r = sum_degree_r + RS_adj_r[edge_r[i]].size()*RS_adj_r[edge_r[i]].size() - (RS_adj_r[edge_r[i]].size() - 1)*(RS_adj_r[edge_r[i]].size() - 1);
//        }
//        else if(bd(gen)){
//            
//            int del_pos = rand() % Reservoir_size;
//            del_l = RS_sample_edge_l[del_pos];
//            del_r = RS_sample_edge_r[del_pos];
//            RS_sample_edge_l[del_pos] = edge_l[i];
//            RS_sample_edge_r[del_pos] = edge_r[i];
//            
//            *std::find(RS_adj_l[del_l].begin(), RS_adj_l[del_l].end(),
//                       del_r) = RS_adj_l[del_l][RS_adj_l[del_l].size() - 1];
//            RS_adj_l[del_l].resize(RS_adj_l[del_l].size() - 1);
////            RS_adj_l[del_l].pop_back();
//            
//            sum_degree_l = sum_degree_l + RS_adj_l[del_l].size()*RS_adj_l[del_l].size() - (RS_adj_l[del_l].size() + 1)*(RS_adj_l[del_l].size() + 1);
//            
//            *std::find(RS_adj_r[del_r].begin(), RS_adj_r[del_r].end(),
//                       del_l) = RS_adj_r[del_r][RS_adj_r[del_r].size() - 1];
//            RS_adj_r[del_r].resize(RS_adj_r[del_r].size() - 1);
////            RS_adj_l[del_r].pop_back();
//            
//            sum_degree_r = sum_degree_r + RS_adj_r[del_r].size()*RS_adj_r[del_r].size() - (RS_adj_r[del_r].size() + 1)*(RS_adj_r[del_r].size() + 1);
//            
//            
//            RS_adj_l[edge_l[i]].push_back(edge_r[i]);
//            sum_degree_l = sum_degree_l + RS_adj_l[edge_l[i]].size()*RS_adj_l[edge_l[i]].size() - (RS_adj_l[edge_l[i]].size() - 1)*(RS_adj_l[edge_l[i]].size() - 1);
//            RS_adj_r[edge_r[i]].push_back(edge_l[i]);
//            sum_degree_r = sum_degree_r + RS_adj_r[edge_r[i]].size()*RS_adj_r[edge_r[i]].size() - (RS_adj_r[edge_r[i]].size() - 1)*(RS_adj_r[edge_r[i]].size() - 1);
//        }
//            
//    }
//    
//    
//    return butterfly_num;
//}

uint64_t Graph::Reservoir_sampling(uint32_t BUCKET_BITS) {
    uint64_t butterfly_num = 0;
    this->Reservoir_size = pow(2, BUCKET_BITS);
    double currentEdges = 0;
    uint32_t increment = 0;
    double y = 0;
    double Pr = 0;
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
    srand(static_cast<unsigned int>(time(0)));
    std::random_device rd;
    std::mt19937 gen(rd());

    RS_sample_edge_l.reserve(Reservoir_size);
    RS_sample_edge_r.reserve(Reservoir_size);

    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (uint32_t i = 0; i < m; i++) {
        if (i % 100000 == 0) cout << i << endl;
        currentEdges++;

        if(currentEdges < Reservoir_size) y = currentEdges;
        else y = (double)Reservoir_size;
        Pr = (y/currentEdges) * ((y-1)/(currentEdges-1)) * ((y-2)/(currentEdges-2));
        increment = round(1/Pr);

        bool use_left_side = (sum_degree_l < sum_degree_r);
        auto& RS_adj_v = use_left_side ? RS_adj_r : RS_adj_l;
        auto& RS_adj_u = use_left_side ? RS_adj_l : RS_adj_r;
        uint32_t v = use_left_side ? edge_r[i] : edge_l[i];
        uint32_t u = use_left_side ? edge_l[i] : edge_r[i];

        for (uint32_t w : RS_adj_u[u]) {
            if (w == v) continue;
            for (uint32_t w2 : RS_adj_v[w]) {
                if (w2 == u) continue;
                if (std::find(RS_adj_v[v].begin(), RS_adj_v[v].end(), w2) != RS_adj_v[v].end()) {
                    butterfly_num += increment;
                }
            }
        }

        if (currentEdges <= Reservoir_size) {
            RS_sample_edge_l.push_back(edge_l[i]);
            RS_sample_edge_r.push_back(edge_r[i]);
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

    return butterfly_num;
}

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
    
    vector<vector<uint32_t> > PS_adj_v, PS_adj_u;
    uint32_t u, v;
    
    zStar = 0;
    
    for(uint32_t i = 0;i<l_n;i++){
        PS_adj_l[i].clear();
    }
    for(uint32_t i = 0;i<r_n;i++){
        PS_adj_r[i].clear();
        
    }
    
    
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
    mt19937_64 generator_FURL(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    unsigned long long c = dist2(generator_FURL);
    unsigned long long d = dist2(generator_FURL);
    
        
    for(uint32_t i = 0; i< m;i++){
        currentEdges++;
        
//        if(i%10000 == 0) cout<<i<<endl;
        
        if(sum_degree_l < sum_degree_r){
            PS_adj_v = PS_adj_r;
            PS_adj_u = PS_adj_l;
            v = edge_r[i];
            u = edge_l[i];
        }
            
        else{
            PS_adj_v = PS_adj_l;
            PS_adj_u = PS_adj_r;
            v = edge_l[i];
            u = edge_r[i];
        }
        
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
//            butterfly_num += (uint32_t)(increment)*intersectionSize(PS_adj_v[v], PS_adj_v[w]);
//            butterfly_num += (uint32_t)(increment)*intersectionSize(PS_adj_v[v], PS_adj_v[w]);
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
            PS_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + PS_adj_l[edge_l[i]].size()*PS_adj_l[edge_l[i]].size() - (PS_adj_l[edge_l[i]].size() - 1)*(PS_adj_l[edge_l[i]].size() - 1);
            PS_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + PS_adj_r[edge_r[i]].size()*PS_adj_r[edge_r[i]].size() - (PS_adj_r[edge_r[i]].size() - 1)*(PS_adj_r[edge_r[i]].size() - 1);
            pq.push(PS_edge(edge_l[i], edge_r[i], random_value));
        }else{
            
            PS_edge min_p_edge = pq.top();
            if(random_value > min_p_edge.priority){
                
                PS_adj_l[edge_l[i]].push_back(edge_r[i]);
                sum_degree_l = sum_degree_l + PS_adj_l[edge_l[i]].size()*PS_adj_l[edge_l[i]].size() - (PS_adj_l[edge_l[i]].size() - 1)*(PS_adj_l[edge_l[i]].size() - 1);
                PS_adj_r[edge_r[i]].push_back(edge_l[i]);
                sum_degree_r = sum_degree_r + PS_adj_r[edge_r[i]].size()*PS_adj_r[edge_r[i]].size() - (PS_adj_r[edge_r[i]].size() - 1)*(PS_adj_r[edge_r[i]].size() - 1);
                pq.push(PS_edge(edge_l[i], edge_r[i], random_value));
                
                pq.pop();
                auto it = find(PS_adj_l[min_p_edge.edge_l].begin(), PS_adj_l[min_p_edge.edge_l].end(), min_p_edge.edge_r);
                if (it != PS_adj_l[min_p_edge.edge_l].end()) {
                    PS_adj_l[min_p_edge.edge_l].erase(it);
                }
                sum_degree_l = sum_degree_l + PS_adj_l[min_p_edge.edge_l].size()*PS_adj_l[min_p_edge.edge_l].size() - (PS_adj_l[edge_l[i]].size() + 1)*(PS_adj_l[edge_l[i]].size() + 1);
                
                it = find(PS_adj_r[min_p_edge.edge_r].begin(), PS_adj_r[min_p_edge.edge_r].end(), min_p_edge.edge_l);
                if (it != PS_adj_r[min_p_edge.edge_r].end()) {
                    PS_adj_r[min_p_edge.edge_r].erase(it);
                }
                sum_degree_r = sum_degree_r + PS_adj_r[min_p_edge.edge_r].size()*PS_adj_r[min_p_edge.edge_r].size() - (PS_adj_r[min_p_edge.edge_r].size() + 1)*(PS_adj_r[min_p_edge.edge_r].size() + 1);
                zStar = max(zStar, min_p_edge.priority);
                
                
            }
            
        }
        
                
        
    }
    
    return butterfly_num;
}

uint64_t Graph::Priority_sampling_FURL0(uint32_t BUCKET_BITS){
    this->Reservoir_size = pow(2,BUCKET_BITS)/2;
    uint64_t butterfly_num = 0;
    uint32_t currentEdges = 0;
    random_seed = time(NULL);
    
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
    mt19937_64 generator_FURL(random_seed);
    mt19937 gen(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    
    priority_queue<PS_edge, std::vector<PS_edge>, ComparePriority_min> buffer;
    
    unsigned long long c = dist2(generator_FURL);
    unsigned long long d = dist2(generator_FURL);
    
    vector<vector<uint32_t> > FURL_adj_v, FURL_adj_u;
    uint32_t u, v;
    double *count_u, *count_v;
    
    exactcnt = 1;
    
    for(uint32_t  i = 0; i< m; i++){
        currentEdges++;
        
        if(i%100000 == 0) cout<<i<<endl;
        
        double edge_rank = hash(edge_l[i],edge_r[i],c,d, gen);
//        double edge_rank = EdgeHash_p(edge_l[i],edge_r[i],c,d);
        
        bool sampled = 0;
        
        if(find(FURL_adj_l[edge_l[i]].begin(),FURL_adj_l[edge_l[i]].end(),edge_r[i]) == FURL_adj_l[edge_l[i]].end()){
            if(currentEdges < Reservoir_size){
                buffer.push(PS_edge(edge_l[i],edge_r[i],edge_rank));
                
                FURL_adj_l[edge_l[i]].push_back(edge_r[i]);
                sum_degree_l = sum_degree_l + FURL_adj_l[edge_l[i]].size()*FURL_adj_l[edge_l[i]].size() - (FURL_adj_l[edge_l[i]].size() - 1)*(FURL_adj_l[edge_l[i]].size() - 1);
                FURL_adj_r[edge_r[i]].push_back(edge_l[i]);
                sum_degree_r = sum_degree_r + FURL_adj_r[edge_r[i]].size()*FURL_adj_r[edge_r[i]].size() - (FURL_adj_r[edge_r[i]].size() - 1)*(FURL_adj_r[edge_r[i]].size() - 1);
                sampled = 1;
                
            }else{
                if(exactcnt){
                    exactcnt = 0;
                }
                PS_edge max_p_edge = buffer.top();
                if(edge_rank < max_p_edge.priority){
                    *std::find(FURL_adj_l[max_p_edge.edge_l].begin(), FURL_adj_l[max_p_edge.edge_l].end(),
                               max_p_edge.edge_r) = FURL_adj_l[max_p_edge.edge_l][FURL_adj_l[max_p_edge.edge_l].size() - 1];
                    FURL_adj_l[max_p_edge.edge_l].resize(FURL_adj_l[max_p_edge.edge_l].size() - 1);
                    sum_degree_l = sum_degree_l + FURL_adj_l[max_p_edge.edge_l].size()*FURL_adj_l[max_p_edge.edge_l].size() - (FURL_adj_l[max_p_edge.edge_l].size() + 1)*(FURL_adj_l[max_p_edge.edge_l].size() + 1);
                    
                    *std::find(FURL_adj_r[max_p_edge.edge_r].begin(), FURL_adj_r[max_p_edge.edge_r].end(),
                               max_p_edge.edge_l) = FURL_adj_r[max_p_edge.edge_r][FURL_adj_r[max_p_edge.edge_r].size() - 1];
                    FURL_adj_r[max_p_edge.edge_r].resize(FURL_adj_r[max_p_edge.edge_r].size() - 1);
                    sum_degree_r = sum_degree_r + FURL_adj_r[max_p_edge.edge_r].size()*FURL_adj_r[max_p_edge.edge_r].size() - (FURL_adj_r[max_p_edge.edge_r].size() - 1)*(FURL_adj_r[max_p_edge.edge_r].size() - 1);
                    
                    buffer.pop();
                    
                    buffer.push(PS_edge(edge_l[i],edge_r[i],edge_rank));
                    
                    FURL_adj_l[edge_l[i]].push_back(edge_r[i]);
                    sum_degree_l = sum_degree_l + FURL_adj_l[edge_l[i]].size()*FURL_adj_l[edge_l[i]].size() - (FURL_adj_l[edge_l[i]].size() - 1)*(FURL_adj_l[edge_l[i]].size() - 1);
                    FURL_adj_r[edge_r[i]].push_back(edge_l[i]);
                    sum_degree_r = sum_degree_r + FURL_adj_r[edge_r[i]].size()*FURL_adj_r[edge_r[i]].size() - (FURL_adj_r[edge_r[i]].size() - 1)*(FURL_adj_r[edge_r[i]].size() - 1);
                    sampled = 1;
                    
                    
                }
            }
            
            if(exactcnt){
                if(sum_degree_l < sum_degree_r){
                    FURL_adj_v = FURL_adj_r;
                    FURL_adj_u = FURL_adj_l;
                    v = edge_r[i];
                    u = edge_l[i];
                    count_u = count_l;
                    count_v = count_r;
                }
                    
                else{
                    FURL_adj_v = FURL_adj_l;
                    FURL_adj_u = FURL_adj_r;
                    v = edge_l[i];
                    u = edge_r[i];
                    count_u = count_r;
                    count_v = count_l;
                }
                double weightSum_uv = 0;
                for(uint32_t w : FURL_adj_u[u]){
                    if(w == v) continue;
                    double weightSum_w = 0;
                    for(uint32_t w2 : FURL_adj_v[w]){
                        if(w2 == u) continue;
                        
                        if(find(FURL_adj_v[v].begin(),FURL_adj_v[v].end(),w2) != FURL_adj_v[v].end()){
                            count_u[w2] ++;
                            weightSum_w ++;
                            weightSum_uv ++;
                            butterfly_num++;
                        }
                        
                    }
                    count_v[w] += weightSum_w;
                                    
                }
                count_u[u]+=weightSum_uv;
                count_v[v]+=weightSum_uv;
                    
            }else{
                double qT = 0;
                if(sampled){
                    qT = (((double)Reservoir_size - 4.0) / (double)Reservoir_size) / pow(buffer.top().priority, 4.0);
                    
                    if(sum_degree_l < sum_degree_r){
                        FURL_adj_v = FURL_adj_r;
                        FURL_adj_u = FURL_adj_l;
                        v = edge_r[i];
                        u = edge_l[i];
                        count_u = count_l;
                        count_v = count_r;
                    }
                        
                    else{
                        FURL_adj_v = FURL_adj_l;
                        FURL_adj_u = FURL_adj_r;
                        v = edge_l[i];
                        u = edge_r[i];
                        count_u = count_r;
                        count_v = count_l;
                    }
                    double weightSum_uv = 0;
                    for(uint32_t w : FURL_adj_u[u]){
                        if(w == v) continue;
                        double weightSum_w = 0;
                        for(uint32_t w2 : FURL_adj_v[w]){
                            if(w2 == u) continue;
                            
                            if(find(FURL_adj_v[v].begin(),FURL_adj_v[v].end(),w2) != FURL_adj_v[v].end()){
                                count_u[w2] +=qT;
                                weightSum_w +=qT;
                                weightSum_uv +=qT;
                                
                                butterfly_num+=qT;
                            }
                            
                        }
                        count_v[w] += weightSum_w;
                                        
                    }
                    count_u[u]+=weightSum_uv;
                    count_v[v]+=weightSum_uv;
                }
            }
        }
        
        
    }
    
    return (uint64_t)butterfly_num;
}

uint64_t Graph::Priority_sampling_FURL(uint32_t BUCKET_BITS, double delta, uint32_t J){
    this->Reservoir_size = pow(2,BUCKET_BITS);
    this->delta = delta;
    this->J = J;
    uint64_t butterfly_num = 0;
    uint32_t currentEdges = 0;
    random_seed = time(NULL);
    
    for(uint32_t i = 0;i<l_n;i++){
        FURL_adj_l[i].clear();
    }
    for(uint32_t i = 0;i<r_n;i++){
        FURL_adj_r[i].clear();
        
    }
    
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
    mt19937_64 generator_FURL(random_seed);
    mt19937 gen(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    
    priority_queue<PS_edge, std::vector<PS_edge>, ComparePriority_min> buffer;
    
    unsigned long long c = dist2(generator_FURL);
    unsigned long long d = dist2(generator_FURL);
    
    vector<vector<uint32_t> > FURL_adj_v, FURL_adj_u;
    unordered_map<int, double> &FURL_count_u = FURL_count_l, &FURL_count_v = FURL_count_r;
    uint32_t u, v;
    
    exactcnt = 1;
    TM = 0;
    time_ = 0;
    
    for(uint32_t  i = 0; i< m; i++){
        currentEdges++;
        time_++;
        
        double edge_rank = hash(edge_l[i],edge_r[i], c, d, gen);
//        double edge_rank = EdgeHash_p(edge_l[i],edge_r[i],c,d);
        
        bool sampled = 0;
        
        if( FURL_count_l.find(edge_l[i]) == FURL_count_l.end() ){
            FURL_count_l[edge_l[i]] = 0.0;
        }
        if( FURL_count_r.find(edge_r[i]) == FURL_count_r.end() ){
            FURL_count_r[edge_r[i]] = 0.0;
        }
        
        if((!exactcnt) && ((time_ - TM) % J == 0)){
            weighted_average(delta);

        }
        
        if(find(FURL_adj_l[edge_l[i]].begin(),FURL_adj_l[edge_l[i]].end(),edge_r[i]) == FURL_adj_l[edge_l[i]].end()){
            if(currentEdges < Reservoir_size){
                buffer.push(PS_edge(edge_l[i],edge_r[i],edge_rank));
                
                FURL_adj_l[edge_l[i]].push_back(edge_r[i]);
                sum_degree_l = sum_degree_l + FURL_adj_l[edge_l[i]].size()*FURL_adj_l[edge_l[i]].size() - (FURL_adj_l[edge_l[i]].size() - 1)*(FURL_adj_l[edge_l[i]].size() - 1);
                FURL_adj_r[edge_r[i]].push_back(edge_l[i]);
                sum_degree_r = sum_degree_r + FURL_adj_r[edge_r[i]].size()*FURL_adj_r[edge_r[i]].size() - (FURL_adj_r[edge_r[i]].size() - 1)*(FURL_adj_r[edge_r[i]].size() - 1);
                sampled = 1;
                
            }else{
                if(exactcnt){
                    TM = time_ - 1;
                    weighted_average(0.0);
                    exactcnt = 0;
                    
                }
                PS_edge max_p_edge = buffer.top();
                if(edge_rank < max_p_edge.priority){
                    if(FURL_adj_l[max_p_edge.edge_l].size() == 1) FURL_adj_l[max_p_edge.edge_l].clear();
                    else {
                        *std::find(FURL_adj_l[max_p_edge.edge_l].begin(), FURL_adj_l[max_p_edge.edge_l].end(),
                                   max_p_edge.edge_r) = FURL_adj_l[max_p_edge.edge_l][FURL_adj_l[max_p_edge.edge_l].size() - 1];
                        FURL_adj_l[max_p_edge.edge_l].resize(FURL_adj_l[max_p_edge.edge_l].size() - 1);
                    }
                    sum_degree_l = sum_degree_l + FURL_adj_l[max_p_edge.edge_l].size()*FURL_adj_l[max_p_edge.edge_l].size() - (FURL_adj_l[max_p_edge.edge_l].size() + 1)*(FURL_adj_l[max_p_edge.edge_l].size() + 1);
                    
                    if(FURL_adj_r[max_p_edge.edge_r].size() == 1) FURL_adj_r[max_p_edge.edge_r].clear();
                    else{
                        *std::find(FURL_adj_r[max_p_edge.edge_r].begin(), FURL_adj_r[max_p_edge.edge_r].end(),
                                   max_p_edge.edge_l) = FURL_adj_r[max_p_edge.edge_r][FURL_adj_r[max_p_edge.edge_r].size() - 1];
                        FURL_adj_r[max_p_edge.edge_r].resize(FURL_adj_r[max_p_edge.edge_r].size() - 1);
                    }
                    
                    sum_degree_r = sum_degree_r + FURL_adj_r[max_p_edge.edge_r].size()*FURL_adj_r[max_p_edge.edge_r].size() - (FURL_adj_r[max_p_edge.edge_r].size() - 1)*(FURL_adj_r[max_p_edge.edge_r].size() - 1);
                    
                    buffer.pop();
                    
                    buffer.push(PS_edge(edge_l[i],edge_r[i],edge_rank));
                    
                    FURL_adj_l[edge_l[i]].push_back(edge_r[i]);
                    sum_degree_l = sum_degree_l + FURL_adj_l[edge_l[i]].size()*FURL_adj_l[edge_l[i]].size() - (FURL_adj_l[edge_l[i]].size() - 1)*(FURL_adj_l[edge_l[i]].size() - 1);
                    FURL_adj_r[edge_r[i]].push_back(edge_l[i]);
                    sum_degree_r = sum_degree_r + FURL_adj_r[edge_r[i]].size()*FURL_adj_r[edge_r[i]].size() - (FURL_adj_r[edge_r[i]].size() - 1)*(FURL_adj_r[edge_r[i]].size() - 1);
                    sampled = 1;
                    
                    
                }
            }
            
            if(exactcnt){
                if(sum_degree_l < sum_degree_r){
                    FURL_adj_v = FURL_adj_r;
                    FURL_adj_u = FURL_adj_l;
                    v = edge_r[i];
                    u = edge_l[i];
                    FURL_count_u = FURL_count_l;
                    FURL_count_v = FURL_count_r;
                    
                }
                    
                else{
                    FURL_adj_v = FURL_adj_l;
                    FURL_adj_u = FURL_adj_r;
                    v = edge_l[i];
                    u = edge_r[i];
                    FURL_count_u = FURL_count_r;
                    FURL_count_v = FURL_count_l;
                }
                double weightSum_uv = 0;
                for(uint32_t w : FURL_adj_u[u]){
                    if(w == v) continue;
                    double weightSum_w = 0;
                    for(uint32_t w2 : FURL_adj_v[w]){
                        if(w2 == u) continue;
                        
                        if(find(FURL_adj_v[v].begin(),FURL_adj_v[v].end(),w2) != FURL_adj_v[v].end()){
                            FURL_count_u[w2] ++;
                            weightSum_w ++;
                            weightSum_uv ++;
                            
                        }
                        
                    }
                    FURL_count_v[w] += weightSum_w;
                                    
                }
                FURL_count_u[u]+=weightSum_uv;
                FURL_count_v[v]+=weightSum_uv;
                    
            }else{
                double qT = 0;
                if(sampled){
                    qT = (((double)Reservoir_size - 4.0) / (double)Reservoir_size) / pow(buffer.top().priority, 4.0);
                    
                    if(sum_degree_l < sum_degree_r){
                        FURL_adj_v = FURL_adj_r;
                        FURL_adj_u = FURL_adj_l;
                        v = edge_r[i];
                        u = edge_l[i];
                        FURL_count_u = FURL_count_l;
                        FURL_count_v = FURL_count_r;
                    }
                        
                    else{
                        FURL_adj_v = FURL_adj_l;
                        FURL_adj_u = FURL_adj_r;
                        v = edge_l[i];
                        u = edge_r[i];
                        FURL_count_u = FURL_count_r;
                        FURL_count_v = FURL_count_l;
                    }
                    double weightSum_uv = 0;
                    for(uint32_t w : FURL_adj_u[u]){
                        if(w == v) continue;
                        double weightSum_w = 0;
                        for(uint32_t w2 : FURL_adj_v[w]){
                            if(w2 == u) continue;
                            
                            if(find(FURL_adj_v[v].begin(),FURL_adj_v[v].end(),w2) != FURL_adj_v[v].end()){
                                FURL_count_u[w2] +=qT;
                                weightSum_w +=qT;
                                weightSum_uv +=qT;
                            }
                            
                        }
                        FURL_count_v[w] += weightSum_w;
                                        
                    }
                    FURL_count_u[u]+=weightSum_uv;
                    FURL_count_v[v]+=weightSum_uv;
                }
            }
        }
    }
    
    if((!exactcnt) && (time_ > TM)){
        if((time_ - TM) % J != 0)
            weighted_average(delta);
    }
    
    for(auto& el : FURL_estimations_l){
        butterfly_num += el.second;
    }
    
    return (uint64_t)butterfly_num/2;
}

double Graph::EdgeHash_p(uint32_t u, uint32_t v, unsigned long long a1, unsigned long long b1) {
    uint64_t edge = (hash_func(u^arr[random_seed%8]) << 32) + hash_func(v^arr[(random_seed+3)%8]);
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

double Graph::hash(uint32_t a, uint32_t b, uint32_t c, uint32_t d, mt19937 generator) {
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

void Graph::weighted_average(double delta){
    for(auto &el : FURL_count_l) {
        if(delta == 0.0)
            FURL_estimations_l[el.first] = el.second;
        else{
            double wgt_val = delta * FURL_estimations_l[el.first] + (1.0 - delta) * el.second;
            FURL_estimations_l[el.first] = wgt_val;
        }
    }
    for(auto &el : FURL_count_r) {
        if(delta == 0.0)
            FURL_estimations_r[el.first] = el.second;
        else{
            double wgt_val = delta * FURL_estimations_r[el.first] + (1.0 - delta) * el.second;
            FURL_estimations_r[el.first] = wgt_val;
        }
    }

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

uint64_t Graph::Priority_sampling_PartitionCT(uint32_t BUCKET_BITS){
    this->Reservoir_size = pow(2,BUCKET_BITS);
    double butterfly_num = 0;
    uint32_t currentEdges = 0;
    random_seed = time(NULL);
    
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
    mt19937_64 generator_FURL(random_seed);
    mt19937 gen(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    
    unsigned long long c = dist2(generator_FURL);
    unsigned long long d = dist2(generator_FURL);
    
    uint32_t num_in_S = 0;
    
    vector<vector<uint32_t> > FURL_adj_v, FURL_adj_u;
    uint32_t u, v;
    
    double q = 1.0, q_bf = 1.0;
    S_PartitionCT.resize(Reservoir_size, PS_edge(0, 0, -1));
    
    uint32_t g_max;
    uint32_t y;
    double n_heat = 0, uu;
    double non_empty_bucket_num = 0;
    
    

    
    for(uint32_t  i = 0; i< m; i++){
        currentEdges++;
        g_max = 0;
        
        double edge_rank = hash(edge_l[i],edge_r[i], c, d, gen);
//        double edge_rank = hash_to_double(edge_l[i], edge_r[i]);
        int edge_bucket = hash_bucket(edge_l[i],edge_r[i], c, d, gen);
//        double edge_rank = EdgeHash_p(edge_l[i],edge_r[i],c,d);
        
        if(S_PartitionCT[edge_bucket].priority == -1 || edge_rank < S_PartitionCT[edge_bucket].priority){
            if(S_PartitionCT[edge_bucket].priority != -1){
                UpdateCounts_PartitionCT(0, num_in_S, edge_l[i], edge_r[i], sum_degree_l, sum_degree_r);
                g_max = static_cast<int>(floor(-log2(S_PartitionCT[edge_bucket].priority)));
                
                if(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() == 1) PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].clear();
                else {
                    *std::find(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].begin(), PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].end(),
                               S_PartitionCT[edge_bucket].edge_r) = PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l][PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - 1];
                    PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].resize(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - 1);
                }
                sum_degree_l = sum_degree_l + PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size()*PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - (PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() + 1)*(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() + 1);
                
                if(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() == 1) PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].clear();
                else{
                    *std::find(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].begin(), PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].end(),
                               S_PartitionCT[edge_bucket].edge_l) = PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r][PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1];
                    PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].resize(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1);
                }
                
                sum_degree_r = sum_degree_r + PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size()*PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - (PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1)*(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1);
                
                
            }else non_empty_bucket_num++;
            y = static_cast<int>(floor(-log2(edge_rank)));
            
            double add_bf_num = UpdateCounts_PartitionCT(1, num_in_S, edge_l[i], edge_r[i], sum_degree_l, sum_degree_r);
            
            double p;
            if(n_heat <= 0 || non_empty_bucket_num<=3){
                p = 1;
            }else{
                p = (non_empty_bucket_num/n_heat)*((non_empty_bucket_num-1)/(n_heat-1)) *((non_empty_bucket_num-2)/(n_heat-2))*((non_empty_bucket_num-3)/(n_heat-3));
            }
            butterfly_num += add_bf_num/p;
            
            if(y > g_max){//直接从这里计数三角形数量
//                if(q<0.5){
//                    cout<<222<<endl;
//                }
                
                n_heat += 2.0/q;
                double uu1 = 1.0/pow(2.0, y);
                double uu2 = 1.0/pow(2.0, g_max);
                uu = uu1 - uu2;
                q += (1.0 / Reservoir_size) * uu;
            }
            
            S_PartitionCT[edge_bucket].edge_l = edge_l[i];
            S_PartitionCT[edge_bucket].edge_r = edge_r[i];
            S_PartitionCT[edge_bucket].priority = edge_rank;
            
            
            
            S_PartitionCT[edge_bucket] = PS_edge(edge_l[i], edge_r[i], edge_rank);
            PartitionCT_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + PartitionCT_adj_l[edge_l[i]].size()*PartitionCT_adj_l[edge_l[i]].size() - (PartitionCT_adj_l[edge_l[i]].size() - 1)*(PartitionCT_adj_l[edge_l[i]].size() - 1);
            PartitionCT_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + PartitionCT_adj_r[edge_r[i]].size()*PartitionCT_adj_r[edge_r[i]].size() - (PartitionCT_adj_r[edge_r[i]].size() - 1)*(PartitionCT_adj_r[edge_r[i]].size() - 1);
            
            
        }
        
    }
    
    
    double gamma_4 = ((double)Reservoir_size/n_heat)*((double)(Reservoir_size-1)/(n_heat-1)) *((double)(Reservoir_size-2)/(n_heat-2)) *((double)(Reservoir_size-3)/(n_heat-3));
   
//    double sum = 0;
//    for (int i = 1; i <= 4; ++i) {
//        sum += pow(-1, i - 1) * (choose(3, i) * pow(1 - i / (double)Reservoir_size, n_heat));
//    }
//    
//    
//    double beta_4 = gamma_4 * (1 - sum);
    double estimated_quadrangles = num_in_S / gamma_4;
    cout<<(uint64_t)n_heat<<endl;
    cout<<(uint64_t)num_in_S<<" "<<gamma_4<<" "<<(uint64_t)estimated_quadrangles<<endl;
    return (uint64_t)butterfly_num;
}

double Graph::choose(int n, int k) {
    if (k == 0) return 1;
    return (n * choose(n - 1, k - 1)) / k;
}

double Graph::UpdateCounts_PartitionCT(uint32_t add, uint32_t &num_in_S, uint32_t v_l, uint32_t v_r, uint64_t &sum_degree_l, uint64_t &sum_degree_r){
    vector<vector<uint32_t> > PartitionCT_adj_v, PartitionCT_adj_u;
    uint32_t u, v;
    
    double num = 0.0;
    
    if(sum_degree_l < sum_degree_r){
        PartitionCT_adj_v = PartitionCT_adj_r;
        PartitionCT_adj_u = PartitionCT_adj_l;
        v = v_r;
        u = v_l;
    }
        
    else{
        PartitionCT_adj_v = PartitionCT_adj_l;
        PartitionCT_adj_u = PartitionCT_adj_r;
        v = v_l;
        u = v_r;
    }
    
    for(uint32_t w : PartitionCT_adj_u[u]){
        if(w == v) continue;
        for(uint32_t w2 : PartitionCT_adj_v[w]){
            if(w2 == u) continue;
            
            if(find(PartitionCT_adj_v[v].begin(),PartitionCT_adj_v[v].end(),w2) != PartitionCT_adj_v[v].end()){
                if(add) num_in_S += 1;
                else num_in_S -= 1;
                num +=1.0;
            }
        }
    }
    return num;
}

double Graph::GetCounts_PartitionCT(uint32_t v_l, uint32_t v_r, uint64_t &sum_degree_l, uint64_t &sum_degree_r){
    vector<vector<uint32_t> > PartitionCT_adj_v, PartitionCT_adj_u;
    uint32_t u, v;
    double count_num = 0;
    
    if(sum_degree_l < sum_degree_r){
        PartitionCT_adj_v = PartitionCT_adj_r;
        PartitionCT_adj_u = PartitionCT_adj_l;
        v = v_r;
        u = v_l;
    }
        
    else{
        PartitionCT_adj_v = PartitionCT_adj_l;
        PartitionCT_adj_u = PartitionCT_adj_r;
        v = v_l;
        u = v_r;
    }
    
    for(uint32_t w : PartitionCT_adj_u[u]){
        if(w == v) continue;
        for(uint32_t w2 : PartitionCT_adj_v[w]){
            if(w2 == u) continue;
            
            if(find(PartitionCT_adj_v[v].begin(),PartitionCT_adj_v[v].end(),w2) != PartitionCT_adj_v[v].end()){
                count_num++;
            }
        }
    }
    return count_num;
}

uint64_t Graph::Priority_sampling_PartitionCT_2(uint32_t BUCKET_BITS, uint64_t real_num){
    
    this->Reservoir_size = pow(2,BUCKET_BITS);
    uint64_t butterfly_num = 0;
    uint32_t currentEdges = 0;
    random_seed = time(NULL);
    
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
    for(uint32_t i = 0;i<l_n;i++){
        PartitionCT_adj_l[i].clear();
    }
    for(uint32_t i = 0;i<r_n;i++){
        PartitionCT_adj_r[i].clear();
        
    }

    
    mt19937_64 generator_FURL(random_seed);
    mt19937 gen(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    unsigned long long c = dist2(generator_FURL);
    unsigned long long d = dist2(generator_FURL);
    
    uint32_t num_in_S = 0;
    
    vector<vector<uint32_t> > FURL_adj_v, FURL_adj_u;
    
    double q = 1;
    S_PartitionCT.clear();
    S_PartitionCT.resize(Reservoir_size, PS_edge(0, 0, -1.0));
    
    uint32_t g_max;
    uint32_t y;
    double n_heat = 0, uu;
    
    HyperLogLog *hll = new HyperLogLog(BUCKET_BITS);
    
    double min_p = 1;

    
    for(uint32_t  i = 0; i< m; i++){
        currentEdges++;
        g_max = 0;
        hll->add(edge_l[i], edge_r[i]);
        
//        cout<<hll->count() << " "<< currentEdges<<endl;;
        
//        continue;
        
        double edge_rank = hash(edge_l[i],edge_r[i], c, d, gen);
        int edge_bucket = hash_bucket(edge_l[i],edge_r[i], c, d, gen);
//        double edge_rank = EdgeHash_p(edge_l[i],edge_r[i],c,d);
        
        double p;
        
//        if(hll->count() < Reservoir_size){
//            p = 1;
//        }else{
            p = ((double)Reservoir_size/(double)hll->count())*((double)(Reservoir_size-1)/((double)hll->count()-1)) *((double)(Reservoir_size-2)/((double)hll->count()-2));
//            p = min(1.0, p);
//        }
        
        if(find(PartitionCT_adj_l[edge_l[i]].begin(),PartitionCT_adj_l[edge_l[i]].end(),edge_r[i]) == PartitionCT_adj_l[edge_l[i]].end()){
            
//            cout<<((double)hll->count()/(double)currentEdges)<<" "<<(double)(hll->count() - Reservoir_size)/(double)hll->count()<<endl;
            if(dis(gen)*(double)(hll->count() - Reservoir_size)/(double)hll->count() < ((double)hll->count()/(double)currentEdges))
            butterfly_num += ((uint64_t)(GetCounts_PartitionCT(edge_l[i], edge_r[i], sum_degree_l, sum_degree_r)))/p;
        }
        
        
        
        if(S_PartitionCT[edge_bucket].priority == -1 || edge_rank < S_PartitionCT[edge_bucket].priority){

            if(S_PartitionCT[edge_bucket].priority != -1){
                
                UpdateCounts_PartitionCT(0, num_in_S, edge_l[i], edge_r[i], sum_degree_l, sum_degree_r);
                g_max = static_cast<int>(floor(-log2(S_PartitionCT[edge_bucket].priority)));
                
                if(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() == 1) PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].clear();
                else {
                    *std::find(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].begin(), PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].end(),
                               S_PartitionCT[edge_bucket].edge_r) = PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l][PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - 1];
                    PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].resize(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - 1);
                }
                sum_degree_l = sum_degree_l + PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size()*PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - (PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() + 1)*(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() + 1);
                
                if(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() == 1) PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].clear();
                else{
                    *std::find(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].begin(), PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].end(),
                               S_PartitionCT[edge_bucket].edge_l) = PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r][PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1];
                    PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].resize(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1);
                }
                
                sum_degree_r = sum_degree_r + PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size()*PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - (PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1)*(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1);
                
                
            }
            y = static_cast<int>(floor(-log2(edge_rank)));
            
            if(y > g_max){
                n_heat += 2.0/q;
                double uu1 = 1.0/pow(2.0, y);
                double uu2 = 1.0/pow(2.0, g_max);
                uu = uu1 - uu2;
                q += (1.0 / Reservoir_size) * uu;
            }
            
            UpdateCounts_PartitionCT(1, num_in_S, edge_l[i], edge_r[i], sum_degree_l, sum_degree_r);
            
//            butterfly_num += ((uint64_t)(GetCounts_PartitionCT(edge_l[i], edge_r[i], sum_degree_l, sum_degree_r)))/p;
            
            S_PartitionCT[edge_bucket] = PS_edge(edge_l[i], edge_r[i], edge_rank);
            PartitionCT_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + PartitionCT_adj_l[edge_l[i]].size()*PartitionCT_adj_l[edge_l[i]].size() - (PartitionCT_adj_l[edge_l[i]].size() - 1)*(PartitionCT_adj_l[edge_l[i]].size() - 1);
            PartitionCT_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + PartitionCT_adj_r[edge_r[i]].size()*PartitionCT_adj_r[edge_r[i]].size() - (PartitionCT_adj_r[edge_r[i]].size() - 1)*(PartitionCT_adj_r[edge_r[i]].size() - 1);
            
            
        }
         
        
    }
    
    
    double gamma_4 = ((double)Reservoir_size/(double)hll->count())*((double)(Reservoir_size-1)/((double)hll->count()-1)) *((double)(Reservoir_size-2)/((double)hll->count()-2)) *((double)(Reservoir_size-3)/((double)hll->count()-3));
   
    double sum = 0;
    for (int i = 1; i <= 4; ++i) {
        sum += pow(-1, i - 1) * (choose(3, i) * pow(1 - i / (double)Reservoir_size, (double)hll->count()));
    }
    
    
    double beta_4 = gamma_4 * (1 - sum);
    double estimated_quadrangles = num_in_S / beta_4;
    
//    cout<<beta_4<<" "<<((uint64_t)estimated_quadrangles)/2<<" "<<n_heat<<" "<<hll->count()<<endl;
    
    uint64_t error_num;
    
    if(((uint64_t)estimated_quadrangles) < real_num) error_num = real_num - ((uint64_t)estimated_quadrangles);
    else error_num = ((uint64_t)estimated_quadrangles) - real_num;
    cout<<"Priority_sampling_PartitionCT_num_global: "<<((uint64_t)estimated_quadrangles)<<"  error rate: "<<error_num/(double)real_num<<endl;
    
    
    return butterfly_num/2;
     
//    return hll->count();
}

uint64_t Graph::Priority_sampling_PartitionCT_3(uint32_t BUCKET_BITS){
    this->Reservoir_size = pow(2,BUCKET_BITS);
    double butterfly_num = 0;
    uint32_t currentEdges = 0;
    random_seed = time(NULL);
    
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    
    for(uint32_t i = 0;i<l_n;i++){
        PartitionCT_adj_l[i].clear();
    }
    for(uint32_t i = 0;i<r_n;i++){
        PartitionCT_adj_r[i].clear();
        
    }

    
    mt19937_64 generator_FURL(random_seed);
    mt19937 gen(random_seed);
    uniform_int_distribution<unsigned long long> dist2(Prime2/3,Prime2);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    unsigned long long c = dist2(generator_FURL);
    unsigned long long d = dist2(generator_FURL);
    
    uint32_t num_in_S = 0;
    
    vector<vector<uint32_t> > FURL_adj_v, FURL_adj_u;
    uint32_t u, v;
    
    double q = 1;
    S_PartitionCT.clear();
    S_PartitionCT.resize(Reservoir_size, PS_edge(0, 0, -1.0));
    
    uint32_t g_max;
    uint32_t y;
    double n_heat = 0, uu;

    double min_p = 1;

    uint32_t non_empty_bucket_num = 0;
    
    
    for(uint32_t  i = 0; i< m; i++){
        currentEdges++;
        g_max=0;
        

        if(i%100000 == 0) cout<<i<<endl;


        double edge_rank = hash(edge_l[i],edge_r[i], c, d, gen);
        int edge_bucket = hash_bucket(edge_l[i],edge_r[i], c, d, gen);
        
        double p;
        
        if(n_heat <= 0 || non_empty_bucket_num<=2){
            p = 1;
        }else{
            p = ((double)non_empty_bucket_num/(double)n_heat)*((double)(non_empty_bucket_num-1)/((double)n_heat-1)) *((double)(non_empty_bucket_num-2)/((double)n_heat-2));
        }
        
        if(find(PartitionCT_adj_l[edge_l[i]].begin(),PartitionCT_adj_l[edge_l[i]].end(),edge_r[i]) == PartitionCT_adj_l[edge_l[i]].end()){
            if(n_heat > 0 && dis(gen)*(double)(n_heat - non_empty_bucket_num)/(double)n_heat < ((double)n_heat/(double)currentEdges))
//            if(n_heat > 0 && dis(gen)< ((double)n_heat/(double)currentEdges))
            butterfly_num += (GetCounts_PartitionCT(edge_l[i], edge_r[i], sum_degree_l, sum_degree_r))/p;
        }
        
        
        
        if(S_PartitionCT[edge_bucket].priority == -1 || edge_rank < S_PartitionCT[edge_bucket].priority){

            if(S_PartitionCT[edge_bucket].priority != -1){
                
//                UpdateCounts_PartitionCT(0, num_in_S, edge_l[i], edge_r[i], sum_degree_l, sum_degree_r);
                g_max = static_cast<int>(floor(-log2(S_PartitionCT[edge_bucket].priority)));
                
                if(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() == 1) PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].clear();
                else {
                    *std::find(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].begin(), PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].end(),
                               S_PartitionCT[edge_bucket].edge_r) = PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l][PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - 1];
                    PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].resize(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - 1);
                }
                sum_degree_l = sum_degree_l + PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size()*PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() - (PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() + 1)*(PartitionCT_adj_l[S_PartitionCT[edge_bucket].edge_l].size() + 1);
                
                if(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() == 1) PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].clear();
                else{
                    *std::find(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].begin(), PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].end(),
                               S_PartitionCT[edge_bucket].edge_l) = PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r][PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1];
                    PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].resize(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1);
                }
                
                sum_degree_r = sum_degree_r + PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size()*PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - (PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1)*(PartitionCT_adj_r[S_PartitionCT[edge_bucket].edge_r].size() - 1);
                
                
            }else non_empty_bucket_num++;
            y = static_cast<int>(floor(-log2(edge_rank)));
            
            if(y > g_max){
                n_heat += 2.0/q;
                double uu1 = 1.0/pow(2.0, y);
                double uu2 = 1.0/pow(2.0, g_max);
                uu = uu1 - uu2;
                q += (1.0 / Reservoir_size) * uu;
            }
            
            
//            UpdateCounts_PartitionCT(1, num_in_S, edge_l[i], edge_r[i], sum_degree_l, sum_degree_r);
            
//            butterfly_num += ((uint64_t)(GetCounts_PartitionCT(edge_l[i], edge_r[i], sum_degree_l, sum_degree_r)))/p;
            
            S_PartitionCT[edge_bucket] = PS_edge(edge_l[i], edge_r[i], edge_rank);
            PartitionCT_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + PartitionCT_adj_l[edge_l[i]].size()*PartitionCT_adj_l[edge_l[i]].size() - (PartitionCT_adj_l[edge_l[i]].size() - 1)*(PartitionCT_adj_l[edge_l[i]].size() - 1);
            PartitionCT_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + PartitionCT_adj_r[edge_r[i]].size()*PartitionCT_adj_r[edge_r[i]].size() - (PartitionCT_adj_r[edge_r[i]].size() - 1)*(PartitionCT_adj_r[edge_r[i]].size() - 1);
            
            
        }
        
    }
    
//    double gamma_4 = ((double)Reservoir_size/(double)n_heat)*((double)(Reservoir_size-1)/((double)n_heat-1)) *((double)(Reservoir_size-2)/((double)n_heat-2)) *((double)(Reservoir_size-3)/((double)n_heat-3));
//   
//    double sum = 0;
//    for (int i = 1; i <= 4; ++i) {
//        sum += pow(-1, i - 1) * (choose(3, i) * pow(1 - i / (double)Reservoir_size, (double)n_heat));
//    }
//    
//    
//    double beta_4 = gamma_4 * (1 - sum);
//    double estimated_quadrangles = num_in_S / beta_4;
//    
//    cout<<(uint64_t)n_heat<<endl;
//    cout<<(uint64_t)estimated_quadrangles<<endl;
//    return (uint64_t)butterfly_num*(n_heat)/currentEdges*n_heat/(n_heat-Reservoir_size);
    return (uint64_t)butterfly_num;
}

uint64_t Graph::Reservoir_sampling_de(uint32_t BUCKET_BITS){
    uint64_t butterfly_num = 0;
    this->Reservoir_size = pow(2,BUCKET_BITS);
    double currentEdges = 0;
    double currentEdges2 = 0;
    uint32_t increment = 0;
    double y = 0;
    double Pr = 0;
    uint64_t sum_degree_l = 0, sum_degree_r = 0;
    vector<vector<uint32_t> > RS_adj_v, RS_adj_u;
    uint32_t u, v, del_l, del_r;
    
    srand(static_cast<unsigned int>(time(0)));
    
    std::random_device rd; // 获取随机设备
    std::mt19937 gen(rd()); // 以随机设备生成随机种子
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    RS_sample_edge_l.clear();
    RS_sample_edge_r.clear();
    
    for(uint32_t i = 0;i<l_n;i++){
        RS_adj_l[i].clear();
    }
    for(uint32_t i = 0;i<r_n;i++){
        RS_adj_r[i].clear();
    }
    
    RS_sample_edge_l.reserve(Reservoir_size);
    RS_sample_edge_r.reserve(Reservoir_size);
    
    HyperLogLog *hll = new HyperLogLog(BUCKET_BITS);
    
    for(uint32_t i = 0; i < m; i++){
        currentEdges++;
        
        hll->add(edge_l[i], edge_r[i]);
        
        if(dis(gen) > (double)hll->count()/(double)currentEdges) continue;
        
        
        currentEdges2++;
        
        if(currentEdges2 < Reservoir_size) y = currentEdges2;
        else y = (double)Reservoir_size;
        Pr = (y/currentEdges2) * ((y-1)/(currentEdges2-1)) * ((y-2)/(currentEdges2-2));
        increment = round(1/Pr);
        
        if(sum_degree_l < sum_degree_r){
            RS_adj_v = RS_adj_r;
            RS_adj_u = RS_adj_l;
            v = edge_r[i];
            u = edge_l[i];
        }
            
        else{
            RS_adj_v = RS_adj_l;
            RS_adj_u = RS_adj_r;
            v = edge_l[i];
            u = edge_r[i];
        }
        
        for(uint32_t w : RS_adj_u[u]){
            if(w == v) continue;
            for(uint32_t w2 : RS_adj_v[w]){
                if(w2 == u) continue;
                
                if(find(RS_adj_v[v].begin(),RS_adj_v[v].end(),w2) != RS_adj_v[v].end()){
                    butterfly_num += increment;
                }
                
            }
//            butterfly_num += increment*intersectionSize(RS_adj_v[v], RS_adj_v[w]);
        }
        
        bernoulli_distribution bd((double)Reservoir_size/currentEdges2);
        
        if((uint32_t)currentEdges2 <= Reservoir_size){
            RS_sample_edge_l.push_back(edge_l[i]);
            RS_sample_edge_r.push_back(edge_r[i]);
            
            RS_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + RS_adj_l[edge_l[i]].size()*RS_adj_l[edge_l[i]].size() - (RS_adj_l[edge_l[i]].size() - 1)*(RS_adj_l[edge_l[i]].size() - 1);
            RS_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + RS_adj_r[edge_r[i]].size()*RS_adj_r[edge_r[i]].size() - (RS_adj_r[edge_r[i]].size() - 1)*(RS_adj_r[edge_r[i]].size() - 1);
        }
        else if(bd(gen)){
            
            int del_pos = rand() % Reservoir_size;
            del_l = RS_sample_edge_l[del_pos];
            del_r = RS_sample_edge_r[del_pos];
            RS_sample_edge_l[del_pos] = edge_l[i];
            RS_sample_edge_r[del_pos] = edge_r[i];
            
            *std::find(RS_adj_l[del_l].begin(), RS_adj_l[del_l].end(),
                       del_r) = RS_adj_l[del_l][RS_adj_l[del_l].size() - 1];
            RS_adj_l[del_l].resize(RS_adj_l[del_l].size() - 1);
            
            sum_degree_l = sum_degree_l + RS_adj_l[del_l].size()*RS_adj_l[del_l].size() - (RS_adj_l[del_l].size() + 1)*(RS_adj_l[del_l].size() + 1);
            
            *std::find(RS_adj_r[del_r].begin(), RS_adj_r[del_r].end(),
                       del_l) = RS_adj_r[del_r][RS_adj_r[del_r].size() - 1];
            RS_adj_r[del_r].resize(RS_adj_r[del_r].size() - 1);
            
            sum_degree_r = sum_degree_r + RS_adj_r[del_r].size()*RS_adj_r[del_r].size() - (RS_adj_r[del_r].size() + 1)*(RS_adj_r[del_r].size() + 1);
            
            
            RS_adj_l[edge_l[i]].push_back(edge_r[i]);
            sum_degree_l = sum_degree_l + RS_adj_l[edge_l[i]].size()*RS_adj_l[edge_l[i]].size() - (RS_adj_l[edge_l[i]].size() - 1)*(RS_adj_l[edge_l[i]].size() - 1);
            RS_adj_r[edge_r[i]].push_back(edge_l[i]);
            sum_degree_r = sum_degree_r + RS_adj_r[edge_r[i]].size()*RS_adj_r[edge_r[i]].size() - (RS_adj_r[edge_r[i]].size() - 1)*(RS_adj_r[edge_r[i]].size() - 1);
        }
            
    }
    
    
    return butterfly_num;
}

double Graph::hash_to_double(uint32_t a, uint32_t b){
    double da = static_cast<double>(a);
    double db = static_cast<double>(b);
        
        // 使用 sin 和 cos 函数生成一个近似均匀分布的值
    double result = fabs(sin(da) * cos(db));
        
        // 确保结果在 0 到 1 范围内
    return fmod(result, 1.0);

}
