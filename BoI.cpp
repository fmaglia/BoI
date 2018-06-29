#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <time.h>
#include <algorithm>
#include "utils.h"
#include <experimental/filesystem>

using namespace std;
namespace fs = std::experimental::filesystem;

int main(int argc, char* argv[])
{

    //std::cout<<"example: ./BoI SIFT1M 8 100 500 true .."<<endl;
    
    /* Configuration variables*/
    string dataset=argv[1];
    bool debug = false;
    string mode = "multiProbe-lsh-bottleneck";
    bool sublinear = true;
    int hash_dimension = stoi(argv[2]); //LSH -> bit (2^hash_dimension)
    float sigma = 1.0;
    int L = stoi(argv[3]);
    bool fastReRanking = true;
    int vicinato = 3;
    int topN = stoi(argv[4]);

    string lastArgument = argv[5];
    if (lastArgument=="false")
        fastReRanking = false;
    else if (lastArgument=="true")
        fastReRanking = true;
    else {
        cout << "Error on fast re-ranking"<<endl;
        return 0;
        }

    string home = argv[6];

    /* End configuration variables */

    int totQuery = L/4;
    int globalVectorDimension;
    int q_results = 0;

    if (dataset=="SIFT1M" || dataset=="SIFT1B"){
        globalVectorDimension = 128;
        q_results = 10000;
    }
    else if (dataset=="GIST1M") {
        globalVectorDimension = 960;
        q_results = 1000;
    }
    else if (dataset=="DEEP1B"){
        q_results = 10000;
        globalVectorDimension = 96;
    }
    

    cout <<"Variables' state"<<endl;
    cout <<"------------------------------------------------------------------------------------"<<endl;

    if (debug)
        cout <<"Debug ACTIVATED!"<<endl;
    cout <<"Dataset: "<<dataset<<endl;

    cout <<globalVectorDimension<<" d"<<endl <<endl;

    if (mode=="lsh")
        cout <<"BoI LSH indexing.  \n- Created "<<L<<" hash tables of "<<pow(2,hash_dimension)<<" rows"<<endl;
    else if (mode=="multiProbe-lsh")
        cout <<"BoI multi-Probe LSH indexing. \n- Created "<<L<<" hash tables of "<<pow(2,hash_dimension)<<" rows with neighborhood: "<<vicinato<<endl;
    else if (mode=="multiProbe-lsh-bottleneck") {
        cout <<"BoI multi-Probe LSH indexing BOTTLENECK. \n- Created "<<L<<" hash tables of "<<pow(2,hash_dimension)<<" rows with neighborhood: "<<vicinato<<endl;

        if (sublinear)
            cout <<"- Sublinear reduction"<<endl;
        else
            cout <<"- Linear reduction"<<endl;
    }

    cout <<"- TopN: "<<topN<<endl;
    if (fastReRanking)
        cout << "fast re-ranking" << endl;

    cout <<"------------------------------------------------------------------------------------"<<endl;


    vector <string> fileTestSet, fileTrainingSet;
    readTrainingAndTest(home, dataset,fileTrainingSet,fileTestSet);
    vector <vector<int>> lsh_index (pow(2,hash_dimension)*L, vector <int>(0));

    int hashCode = pow(2,hash_dimension);
    vector <vector <int> > query_results(q_results, vector<int>(0));

    auto tEncoding1 = std::chrono::high_resolution_clock::now();

    //precalculate the binary nearest neighbor
    vector <vector<int>> neighbor (0, vector <int>(0));

    for (int i=0; i < pow(2, hash_dimension); ++i) {
        string binary = calculateBinary(i, hash_dimension);
        vector <int> vicini;
        for (int v=1; v <= vicinato; ++v) {
            for (int j=0; j < hash_dimension; ++j) {
                calculateNeighbors(vicini, binary, j, v);
            }
        }
        vicini.insert(vicini.begin(),i);
        neighbor.push_back(vicini);
    }

    int gapBucket = neighbor[0].size();
    int initialGap = gapBucket;

    std::cout << "Gap bucket: "<< gapBucket << endl;

    vector <vector<float>> VLAD_trainingSet;

    //gaussian distribution (mean = 0 and stddev = sigma) for the assignment of value of projection vector
    std::normal_distribution<double> distribution(0.0,sigma);
    std::random_device rd;
    std::mt19937 generator(rd()); //use 0 as a parameter for using VALGRIND (profiler), otherwise use rd()
    //std::mt19937 generator(0); //use 0 as a parameter for using VALGRIND (profiler), otherwise use rd()

    vector <vector<float>> projectionVector (hash_dimension*L, vector<float>(globalVectorDimension));

    for (int i=0; i<projectionVector.size(); ++i) {
        for (int j=0; j<projectionVector[i].size(); ++j) {
            projectionVector[i][j] = distribution(generator);
        }
    }

    int trainingElements = 0;
    /*Reading training/database data (feature detector & descriptor) */
    
    std::ifstream fileStream(fileTrainingSet[0],std::ios::binary);
    int count = 0;
    int counter = 0;
    float f;
    std::vector <float> VLAD_row;
    while (fileStream.read(reinterpret_cast<char*>(&f), sizeof(float))){
        VLAD_row.push_back(f);
        counter++;
        if (counter == globalVectorDimension) {            
            if (fastReRanking)
                VLAD_trainingSet.push_back(VLAD_row);

            for (int hashTables = 0; hashTables < L; ++hashTables) 
                lsh_index[lsh_indexing(hash_dimension, VLAD_row, projectionVector, hashTables)].push_back(count);
            
            count++;
            trainingElements++;
            counter = 0;
            VLAD_row.clear();
        }

    }

        fileStream.close();   


    auto tEncoding2 = std::chrono::high_resolution_clock::now();

    float timeEncoding = (std::chrono::duration_cast<std::chrono::seconds>(tEncoding2 - tEncoding1).count());

    cout << "Encoding of "<<trainingElements<<" database images TERMINATED in "<<timeEncoding<<" s"<< endl;
    cout << "Queries"<<endl;

    double reRankingTime = 0.0;

    std::ifstream DbFileReader(home+"/dataset/"+dataset+"/Db.dat",std::ios::binary);

    auto t1 = std::chrono::high_resolution_clock::now();

    std::ifstream fileStreamQuery(fileTestSet[0],std::ios::binary);
    count = 0;
    float lower_bound = L * 0.2;
    counter = 0;
    std::vector <float> VLAD_testSet;
    while (fileStreamQuery.read(reinterpret_cast<char*>(&f), sizeof(float))){
        VLAD_testSet.push_back(f);
        counter++;

        if (counter == globalVectorDimension) {  
                      
            gapBucket = initialGap;
            vector <int> valueBucket = {68,37,8,1};
            //vector <int> valueBucket = {37,22,8,1};      
            //vector <int> valueBucket = {8, 5, 3, 1};
    
            /*Calcuate the distance between the test set image and the training set images*/
            vector <float> LucaImagePosition(trainingElements,0.0);      
            //auto s1 = std::chrono::high_resolution_clock::now();

            for (int hashTables = 0; hashTables < L; ++hashTables) {

                int offset = hashCode*hashTables;
                int queryRetrieved = lsh_indexing(hash_dimension, VLAD_testSet, projectionVector, hashTables);

                if (sublinear && hashTables >= L/2 && hashTables%totQuery==0 && gapBucket>0){          //sublinear reduction
                    valueBucket.erase(valueBucket.begin());
                    gapBucket = valueBucket[0];
                     //gapBucket -= 2;
                }
                
                searchMultiProbeLSH(neighbor, lsh_index, queryRetrieved, LucaImagePosition, gapBucket, offset, hash_dimension);                
            }

            vector<ranking> imagePosition;
            //auto s3 = std::chrono::high_resolution_clock::now();

            for (int j=0; j<LucaImagePosition.size(); ++j)
                if (LucaImagePosition[j] > lower_bound)
                    imagePosition.push_back(ranking{j,LucaImagePosition[j]});

            std::sort(imagePosition.begin(), imagePosition.end(),
                        [](const ranking& a, const ranking& b) {
                    return a.weight > b.weight;
            });


            //re-ranking according to euclidean distance

            if (imagePosition.size()>topN)
                imagePosition.resize(topN);
                
            auto tReRanking1 = std::chrono::high_resolution_clock::now();

            for (int r=0; r < imagePosition.size(); r++) {

                auto position = imagePosition[r].index;

                if (!fastReRanking) {
                    vector <float> VLAD_row = readIthRow_binary_new(DbFileReader,globalVectorDimension, position);
                    imagePosition[r].weight = l2_norm_2vectors(VLAD_row,VLAD_testSet);
                }
                else
                    imagePosition[r].weight = l2_norm_2vectors(VLAD_trainingSet[position], VLAD_testSet);                
            }

            std::sort(imagePosition.begin(), imagePosition.end(),[](const ranking& a, const ranking& b) {
                    return a.weight < b.weight;
            });

            auto tReRanking2 = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double, std::milli> timeR = tReRanking2 - tReRanking1;

            reRankingTime += (timeR.count());

            if (imagePosition.size() < topN){
                std::cout << "Problem at the query "<<count<<" with only "<<imagePosition.size()<<" elements"<<endl;
                return 0;
            }

            for (int r=0; r < imagePosition.size(); r++) 
                query_results[count].push_back(imagePosition[r].index);
                
            count++;

            if (count%(q_results/3)==0 && count!=0)
                cout << count<<" on "<<q_results<<endl;
        
        counter = 0;
        VLAD_testSet.clear();  
        }
        
    }

    fileStreamQuery.close();     

    cout<<"Matching finished"<<endl;

    auto t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> time = t2 - t1;

    double avgQueryTime = time.count() / q_results;
    reRankingTime /= q_results;

    cout << "AVERAGE QUERY TIME: "<<avgQueryTime<<" ms - avg re-ranking time: "<<reRankingTime<<" ms"<<endl;

    string url;

    if (debug) {
        url = "_DEBUG";
    }
    else {

        url = "BoI"+mode;

        if (mode=="lsh") {
            hash_dimension = pow(2,hash_dimension);
            url+=to_string(hash_dimension)+"_L"+to_string(L);
        }
        else if (mode=="multiProbe-lsh") {
            hash_dimension = pow(2,hash_dimension);
            url+=to_string(hash_dimension)+"_L"+to_string(L);
        }
        else if (mode=="multiProbe-lsh-bottleneck") {
            hash_dimension = pow(2,hash_dimension);
            url+=to_string(hash_dimension)+"_L"+to_string(L);

            if (sublinear)
                url += "_sublinearReduction";
            else
                url += "_linearReduction";

        }

        url+="_top"+to_string(topN);

        if (fastReRanking)
            url += "_fast";
    }

    string result_url = home+"/results/"+dataset+"/"+url;
    writeResults(dataset, result_url, query_results, q_results);
    calcResults(home, dataset, result_url);

    return 0;

}

