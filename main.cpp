#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <algorithm>
#include "Graph.hpp"
#include <fstream>

using namespace std;

int main(int argc, const char * argv[]) {

    string file = "mooc_actions";

    uint32_t Reservoir_size = 10;
    if(argc > 1){
        file = string(argv[1]);
    }
    if(argc > 2) Reservoir_size = stoi(argv[2]);
    double error_num = 0;

    struct rusage before, after;
    
    string file_path = file;

    getrusage(RUSAGE_SELF, &before);
    Graph *graph = new Graph(&file_path[0], Reservoir_size);

    vector<double> alg1_error, alg2_error;
    

    graph->readGraph();
    uint64_t real_num;

    getrusage(RUSAGE_SELF, &after);
    // std::cout << "readGraph Memory difference: " << (after.ru_maxrss - before.ru_maxrss)<<endl;
    
//    cout<<graph->Priority_sampling_DEABC_PLUS(Reservoir_size)<<endl;
//    cout<<graph->Priority_sampling_DEABC_PLUS_2(Reservoir_size, 0)<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    
    if(argc > 3) real_num = stoll(argv[3]);
    else{
        cout<<"start Butterfly_counting..."<<endl;
        start = std::chrono::high_resolution_clock::now();
        real_num = graph->Butterfly_counting();
        if(argc > 3) real_num = stoll(argv[3]);
        cout<<"the number of butterflies: "<<real_num<<endl;
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "Function execution time: " << duration.count() << std::endl;
        cout<<endl;
    }

    
    
    
    // cout<<"start Reservoir_sampling..."<<endl;
    // start = std::chrono::high_resolution_clock::now();
    // graph->init();
    // getrusage(RUSAGE_SELF, &before);
    // uint64_t Reservoir_sampling_num = graph->Reservoir_sampling(Reservoir_size);
    // getrusage(RUSAGE_SELF, &after);
    // if(Reservoir_sampling_num < real_num) error_num = real_num - Reservoir_sampling_num;
    // else error_num = Reservoir_sampling_num - real_num;
    // cout<<"Reservoir_sampling_num: "<<Reservoir_sampling_num<<"  error rate: "<<error_num/(double)real_num<<endl;
    // end = std::chrono::high_resolution_clock::now();
    // duration = end - start;
    // std::cout << "Function execution time: " << duration.count() << std::endl;
    // std::cout << "Memory difference: " << after.ru_maxrss <<" "<< before.ru_maxrss <<" "<< (after.ru_maxrss - before.ru_maxrss)<<endl;
    // cout<<endl;

    // start = std::chrono::high_resolution_clock::now();
    // uint64_t Priority_sampling_random_num = graph->Priority_sampling(Reservoir_size, 2);
    // if(Priority_sampling_random_num < real_num) error_num = real_num - Priority_sampling_random_num;
    // else error_num = Priority_sampling_random_num - real_num;
    // cout<<"Priority_sampling_random_num: "<<Priority_sampling_random_num<<"  error rate: "<<error_num/(double)real_num<<endl;
    // end = std::chrono::high_resolution_clock::now();
    // duration = end - start;
    // std::cout << "Function execution time: " << duration.count() << std::endl;
    // cout<<endl;

    

    cout<<"start Priority_sampling_DEABC0..."<<endl;
    start = std::chrono::high_resolution_clock::now();
    graph->init();
    uint64_t Priority_sampling_DEABC0_num = graph->Priority_sampling_DEABC0(Reservoir_size);
    if(Priority_sampling_DEABC0_num < real_num) error_num = real_num - Priority_sampling_DEABC0_num;
    else error_num = Priority_sampling_DEABC0_num - real_num;
    cout<<"Priority_sampling_DEABC0_num: "<<Priority_sampling_DEABC0_num<<"  error rate: "<<error_num/(double)real_num<<endl;
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Function execution time: " << duration.count() << std::endl;


    cout<<"start Priority_sampling_DEABC_PLUS..."<<endl;
    start = std::chrono::high_resolution_clock::now();
    graph->init();
    uint64_t Priority_sampling_DEABC_PLUS_num_3 = graph->Priority_sampling_DEABC_PLUS(Reservoir_size);
    if(Priority_sampling_DEABC_PLUS_num_3 < real_num) error_num = real_num - Priority_sampling_DEABC_PLUS_num_3;
    else error_num = Priority_sampling_DEABC_PLUS_num_3 - real_num;
    cout<<"Priority_sampling_DEABC_PLUS_num_local: "<<Priority_sampling_DEABC_PLUS_num_3<<"  error rate: "<<error_num/(double)real_num<<endl;
    end = std::chrono::high_resolution_clock::now();
    
    duration = end - start;
    std::cout << "Function execution time: " << duration.count() << std::endl;
    cout<<endl;

    
    return 0;
}

