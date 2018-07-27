#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
using namespace std;

struct ranking {
    int index;
    float weight;
};

struct ranking_final {
    string name;
    float weight;
};

float l2_norm_2vectors(std::vector<float> &v, std::vector<float> &u) {
    float accum = 0.;
    for (int i = 0; i < v.size(); ++i) {
        accum += pow(v[i] - u[i],2);
    }
    return sqrt(accum);
}

vector <float> readIthRow_binary_new(std::ifstream & re, int dimension, int row){
    //std::ifstream re(url,std::ios::binary);
    // std::cout << "dimension: "<<dimension<<" row: "<<row<<endl;
    re.seekg(int64_t(4) * dimension * row);
    float f;
    int counter = 0;
    vector <float> vettore;
    while (re.read(reinterpret_cast<char*>(&f), sizeof(float))){
        vettore.push_back(f);
        counter++;
        if (counter == dimension)
            break;      
    }

    return vettore;
}

void readTrainingAndTest(const string &home, const string &dataset, vector <string> &fileTrainingSet, vector <string> &fileTestSet) {

    fileTrainingSet.push_back(home+"/dataset/"+dataset+"/Db.dat");
    fileTestSet.push_back(home+"/dataset/"+dataset+"/query.dat");   
}

void writeResults(string dataset, string result_url, const vector <vector <int>> &query_results, const int &q_results){

    cout << "Writing results on file: "<< result_url <<endl;
    if (dataset == "SIFT1M" || dataset=="GIST1M" || dataset=="SIFT1B" || dataset=="DEEP1B"){
        std::ofstream results_file;
        results_file.open(result_url);
        int R = 100;
        if (dataset=="SIFT1B" || dataset=="DEEP1B")
            R = 1000;
        for (int j=0; j < q_results ; j++) {
            for (int element=0; element < R; element ++) {
                /*int value = query_results[j][element];
                if (query_results[j].size() >= element)
                    results_file<< value <<" ";
                else
                    results_file << int(0) <<" ";*/
                results_file<< query_results[j][element] << " ";
                    //std::cout<<query_results[j][element]<<" ";
            }
            //std::cout<<endl;
            results_file<<"\n";
        }
        results_file.close();
    }

}

void calcResults (const string& home, const string &dataset, string result_url) {
    //Launching python or C++ script

    if (dataset=="SIFT1M" || dataset=="GIST1M" || dataset=="SIFT1B" || dataset=="DEEP1B") {
        string command = "python " + home + "/dataset/"+dataset+"/calculate_recall.py "+result_url;
        system(command.c_str());
    }
    else
        cout << "Error on calculating results"<<endl;

}

int lsh_indexing(const int hash_dimension, std::vector<float> &VLAD_row, std::vector<std::vector<float>> &projectionVector, const int iteration) {

    int result=0;
    
    for (int j=0; j < hash_dimension; ++j) {
        float subresult = 0;

        for (int i=0; i < VLAD_row.size(); ++i) {
            subresult += VLAD_row[i]*projectionVector[j+iteration*hash_dimension][i];
        }
        if (subresult > 0) {
            switch (j)
            {
                case 0:  result += 1;
                break;
                case 1: result += 2;
                break;
                case 2: result += 4;
                break;
                case 3: result += 8;
                break;
                case 4: result += 16;
                break;
                case 5: result += 32;
                break;
                case 6: result += 64;
                break;
                case 7: result += 128;
                break;
                case 8: result += 256;
                break;
                case 9: result += 512;
                break;
                case 10: result += 1024;
                break;
                case 11: result += 2048;
                break;
                case 12: result += 4096;
                break;
                case 13: result += 8192;
                break;
                case 14: result += 16384;
                break;
                case 15: result += 32768;
                break;
                case 16: result += 65536;
                break;
                case 17: result += 131072;
                break;
                case 18: result += 262144;
                break;
                case 19: result += 524288;
                break;
                case 20: result += 1048576;
                break;
            }
            //result += pow(2,j);
        }
        //cout << "iteration: "<<iteration<<" hash_dim "<<j<<" subresult: "<<subresult<<" result: "<<result<<endl;
    }

    //update hash for unique vector (using different hash tables)
    //result += iteration*pow(2,hash_dimension);
    switch (hash_dimension) {
        case 4: result += iteration*16;
        break;
        case 5: result += iteration*32;
        break;
        case 6: result += iteration*64;
        break;
        case 7: result += iteration*128;
        break;
        case 8: result += iteration*256;
        break;
        case 9: result += iteration*512;
        break;
        case 10: result += iteration*1024;
        break;
        case 11: result += iteration*2048;
        break;
        case 12: result += iteration*4096;
        break;
        case 13: result += iteration*8192;
        break;
        case 14: result += iteration*16384;
        break;
        case 15: result += iteration*32768;
        break;
        case 16: result += iteration*65536;
        break;
        case 17: result += iteration*131072;
        break;
        case 18: result += iteration*262144;
        break;
        case 19: result += iteration*524288;
        break;
        case 20: result += iteration*1048576;
        break;
    }
    //result += iteration*256;
    
    return result;

}

void searchLSH(const vector <vector<int>> &lsh_index, const int queryRetrieved, vector <float> &imagePosition) {

    for (int elementBucket = 0; elementBucket < lsh_index[queryRetrieved].size(); ++elementBucket)
    {
        ranking img;
        img.index = lsh_index[queryRetrieved][elementBucket];
        img.weight = 1.0;

        imagePosition[img.index] += img.weight;
    }

}

void searchMultiProbeLSH(const vector <vector<int>> &neighbor, const vector <vector<int>> &lsh_index, const int queryRetrieved, vector <float> &imagePosition, const int gapBucket, const int offset, const int hash_dimension, const vector <int> &checkBuckets) {
    int temp = queryRetrieved-offset;
    auto && v_index2 = neighbor[temp];
    
    #pragma omp parallel num_threads(4)
    {
    
    for (int indexBucket=0; indexBucket < gapBucket; ++indexBucket) {          //convert bucket_test[hashTable] to 0...511

        int index2 = v_index2[indexBucket]+offset;

        float weight = 1.0;

        if (indexBucket > 0 && indexBucket <= hash_dimension)
            weight = 0.7;
        else if (indexBucket > hash_dimension && indexBucket <= checkBuckets[1])
            weight = 0.4;
        else if (indexBucket > checkBuckets[1])
            weight = 0.1;

        auto && vec_index = lsh_index[index2];
        auto siz = vec_index.size();
        
        #pragma omp for
        for (int elementBucket = 0; elementBucket < siz; ++elementBucket) {

            auto && index = vec_index[elementBucket];

            imagePosition[index] += weight;

        }
    }
    }
}

string calculateBinary(int number, int hash_dimension) {
    string result = "";
    int counter = 0;

    while (number>0) {
        result.insert(0,to_string(number%2));
        number /= 2;
        counter++;
    }

    while (counter<hash_dimension){
        result.insert(0,"0");
        counter++;
    }

    return result;


}
/*
string calcStringName (const string &start, const int element, const string &dataset){

    if (dataset == "Holidays" || dataset == "Parma" || dataset == "Flickr1M")
        return fs::path(start).stem().string()+".jpg";
    else if (dataset=="ZuBuD")
        return to_string(element+1);
    else if (dataset == "Oxford" || dataset == "Paris") {
        return fs::path(start).stem();
    }
    else if (dataset == "UKB") {
        string photo = fs::path(start).stem();
        return photo.substr(7, photo.size());
    }
    else if (dataset == "SIFT1M" || dataset=="GIST1M" ){
        return fs::path(start).stem();
    }

    cout <<"Error on calculating string"<<endl;

    return "Error";
}*/

int calculateDecimal (string binary) {
    int result = 0;

    for (int i=0; i < binary.size(); ++i) {
        if ('1'==binary[i])
            result += pow(2,binary.size()-i-1);
    }

    return result;
}
string changeBit(string binary, int position) {

    if (binary[position]=='0')
        binary[position]='1';
    else
        binary[position]='0';

    return binary;
}

void calculateNeighbors (vector <int> & vicini, string binary, int position, int vicinato){

    if (vicinato==1) {
        vicini.push_back(calculateDecimal(changeBit(binary, position)));
    }
    else {
        for (int j=0; j<binary.size()-1; ++j){

            string newBinary = changeBit(binary, position);
            for (int v=1; v < vicinato; ++v) {
                if (j+v!=position)
                    newBinary = changeBit(newBinary, j+v);
                else if (j+v-1 >= 0)
                    newBinary = changeBit(newBinary, j+v-1);
                else if (j+v+1 < binary.size())
                    newBinary = changeBit(newBinary, j+v+1);
            }

            int value = calculateDecimal(newBinary);
            if (std::find(vicini.begin(), vicini.end(), value) == vicini.end())
                vicini.push_back(value);
        }

    }
    return;

}

void calculateBuckets(vector<int> &valueBucket, const int hash_dimension){
    if (hash_dimension == 20) {
        valueBucket.insert(valueBucket.end(), {518, 211, 20, 1});
    }
    if (hash_dimension == 19) {
        valueBucket.insert(valueBucket.end(), {464, 191, 19, 1});
    }
    else if (hash_dimension == 18) {
        valueBucket.insert(valueBucket.end(), {413, 172, 18, 1});
    }
    else if (hash_dimension == 17) {
        valueBucket.insert(valueBucket.end(), {365, 154, 17, 1});
    }
    else if (hash_dimension == 16) {
        valueBucket.insert(valueBucket.end(), {320,137,16,1});
    }
    else if (hash_dimension == 15) {
        valueBucket.insert(valueBucket.end(), {278,121,15,1});
    }
    else if (hash_dimension == 14) {
        valueBucket.insert(valueBucket.end(), {239,106,14,1});
    }
    else if (hash_dimension == 13) {
        valueBucket.insert(valueBucket.end(),{203,92,13,1}); 
    }
    else if (hash_dimension == 12) {
        valueBucket.insert(valueBucket.end(),{170,79,12,1}); 
    }
    else if (hash_dimension == 11) {
        valueBucket.insert(valueBucket.end(),{140,67,11,1});
    }
    else if (hash_dimension == 10){
        valueBucket.insert(valueBucket.end(),{113,56,10,1});
    } 
    else if (hash_dimension == 9) {
        valueBucket.insert(valueBucket.end(),{89,46,9,1});
    }
    else if (hash_dimension == 8) {
        valueBucket.insert(valueBucket.end(),{68,37,8,1});
    }
}


#endif // UTILS_H
