# Bag of Indexes (BoI)

[[paper](https://arxiv.org/abs/1806.05946)] [[project](http://implab.ce.unipr.it/?page_id=787)] [[poster](http://implab.ce.unipr.it/wp-content/uploads/2018/07/conference_poster_6.pdf)]

Bag of Indexes (BoI) is a multi-index hashing C++ library, used for solving Approximate Nearest Neighbors (ANN) search problems. In the belove figure, an example of BoI's working is presented.

![BoI_overview](http://implab.ce.unipr.it/wp-content/uploads/2018/06/BoI.png)

LSH and multi-probe LSH are the only, at the moment, projection functions implemented.

## Datasets
The original dataset files are converted in binary through the application of a simple C++ script:
- [SIFT1M](https://drive.google.com/open?id=1w7doJHujhINbR6HXF73_H-udlxL5btiI)
- [GIST1M](https://drive.google.com/open?id=13gnQlNUCdSIdnpZYDr3MUMp4vHQiyo4L) 

After downloading the dat files you need to create a folder called `dataset ` and then put in the uncompressed version.

## Installation
* Requirements:
  * G++ 5.4 or greater
  * openmp
* Build:
` g++ -o BoI BoI.cpp -lstdc++fs -std=c++14 -fopenmp -O2 `

## Test
` ./BoI SIFT1M 16 25 500 true .. `

In the previous command the values of the parameters are:
* δ = 16;
* L = 25;
* topN = 500;
* fast re-ranking = true (loading descriptors from RAM).

## Results

### SIFT1M
The following table represents the results obtained through the application of BoI adaptive multi-probe LSH, trying to change some parameters: 
 * L: number of hash tables;
 * δ: number of buckets per hash table (hash dimension);
 * topN: number of top elements to re-rank based on original distance (usually Euclidean distance).
 
 In every test, the neighborhood is setted to 3.

| Configuration        | 1           | 10  | 100 | avg retrieval time |
| ------------- |:-------------:| :-----:| :---:| ---------:|
| δ = 16, L = 25, top500 | 90.50%   | 91.57% | 91.57% | 6 msec |
| δ = 15, L = 50, top500 | 98.36%   | 99.19% | 99.19% | 18 msec |
| δ = 15, L = 50, top10k | 99.20%   | 100%   | 100%   | 18 msec |

The reduced average retrieval time is obtained through the multi-threading and the cut-off (an unsuperivesd reduction of the BoI structure). For more info see the paper cited in the Reference section.
The reason behind the same retrieval time changing the top elements to re-rank is due to the speed of loading files from RAM memory.

### GIST1M

| Configuration        | 1           | 10  | 100 | avg retrieval time |
| ------------- |:-------------:| :-----:| :---:| ---------:|
| δ = 16, L = 100, top500 | 79.80%   | 80.80% | 80.80% | 60 msec |
| δ = 16, L = 100, top10k | 97.40%   | 98.50%   | 98.50%   | 100 msec |


## Reference

<pre>@article{magliani2018efficient,
  title={Efficient Nearest Neighbors Search for Large-Scale Landmark Recognition},
  author={Magliani, Federico and Fontanini, Tomaso and Prati, Andrea},
  journal={arXiv preprint arXiv:1806.05946},
  year={2018}
}</pre>
