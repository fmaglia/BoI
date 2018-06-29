# Bag of Indexes (BoI)

[[paper](https://arxiv.org/abs/1806.05946)] [[project](http://implab.ce.unipr.it/?page_id=787)]

Bag of Indexes (BoI) is a multi-index hashing C++ library, used for solving Approximate Nearest Neighbors (ANN) search problems. In the belove figure, an example of BoI's working is presented.

![BoI_overview](http://implab.ce.unipr.it/wp-content/uploads/2018/06/BoI.png)

LSH and multi-probe LSH are the only, at the moment, projection functions implemented.

## Datasets
The files are binaries. They are created using a simple C++ script.
- Holidays+Flickr1M:
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
` ./BoI SIFT1M 8 25 500 true .. `

## Reference

<pre>@article{magliani2018efficient,
  title={Efficient Nearest Neighbors Search for Large-Scale Landmark Recognition},
  author={Magliani, Federico and Fontanini, Tomaso and Prati, Andrea},
  journal={arXiv preprint arXiv:1806.05946},
  year={2018}
}</pre>
