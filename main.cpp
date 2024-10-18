#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <algorithm>
#include "Graph.hpp"


using namespace std;

int main(int argc, const char * argv[]) {

    string file = "wiki_edit1.txt";

//    file = "Edit-frwiki.txt";
//    file = "livejournal.graph";
    uint32_t Reservoir_size = 12;
    if(argc > 1){
        file = string(argv[1]);
    }
    if(argc > 2) Reservoir_size = stoi(argv[2]);
    double error_num = 0;
    
    
    string file_path = "/Users/milk/test_data/bipartite/all_data/" + file;
    Graph *graph = new Graph(&file_path[0]);
    
    graph->readGraph();
    
//    cout<<graph->Priority_sampling_DEABCPLUS(Reservoir_size)<<endl;
//    cout<<graph->Priority_sampling_DEABCPLUS_2(Reservoir_size, 0)<<endl;
    
    uint64_t real_num = graph->Butterfly_counting();
    if(argc > 3) real_num = stoll(argv[3]);
    cout<<"the number of butterflies: "<<real_num<<endl;
    
    
    uint64_t Reservoir_sampling_num = graph->Reservoir_sampling(Reservoir_size);
    if(Reservoir_sampling_num < real_num) error_num = real_num - Reservoir_sampling_num;
    else error_num = Reservoir_sampling_num - real_num;
    cout<<"Reservoir_sampling_num: "<<Reservoir_sampling_num<<"  error rate: "<<error_num/(double)real_num<<endl;
//    
//    uint64_t Priority_sampling_random_num = graph->Priority_sampling(Reservoir_size, 2);
//    if(Priority_sampling_random_num < real_num) error_num = real_num - Priority_sampling_random_num;
//    else error_num = Priority_sampling_random_num - real_num;
//    cout<<"Priority_sampling_random_num: "<<Priority_sampling_random_num<<"  error rate: "<<error_num/(double)real_num<<endl;
    
//    uint64_t Priority_sampling_hash_FURL_num = graph->Priority_sampling(Reservoir_size, 0);
//    if(Priority_sampling_hash_FURL_num < real_num) error_num = real_num - Priority_sampling_hash_FURL_num;
//    else error_num = Priority_sampling_hash_FURL_num - real_num;
//    cout<<"Priority_sampling_hash_FURL_num: "<<Priority_sampling_hash_FURL_num<<"  error rate: "<<error_num/(double)real_num<<endl;
    
//    uint64_t Priority_sampling_hash_num = graph->Priority_sampling(Reservoir_size, 1);
//    if(Priority_sampling_hash_num < real_num) error_num = real_num - Priority_sampling_hash_num;
//    else error_num = Priority_sampling_hash_num - real_num;
//    cout<<"Priority_sampling_hash_num: "<<Priority_sampling_hash_num<<"  error rate: "<<error_num/(double)real_num<<endl;
    
    uint64_t Priority_sampling_DEABC_num = graph->Priority_sampling_DEABC(Reservoir_size);
    if(Priority_sampling_DEABC_num < real_num) error_num = real_num - Priority_sampling_DEABC_num;
    else error_num = Priority_sampling_DEABC_num - real_num;
    cout<<"Priority_sampling_DEABC_num: "<<Priority_sampling_DEABC_num<<"  error rate: "<<error_num/(double)real_num<<endl;

    uint64_t Priority_sampling_DEABCPLUS_num = graph->Priority_sampling_DEABCPLUS(Reservoir_size);
    if(Priority_sampling_DEABCPLUS_num < real_num) error_num = real_num - Priority_sampling_DEABCPLUS_num;
    else error_num = Priority_sampling_DEABCPLUS_num - real_num;
    cout<<"Priority_sampling_DEABCPLUS_num: "<<Priority_sampling_DEABCPLUS_num<<"  error rate: "<<error_num/(double)real_num<<endl;
    

    return 0;
}

