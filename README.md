# SWAT

This is the source code of our paper *[SWAT: A System-Wide Approach to Tunable Leakage Mitigation in Encrypted Stores](https://arxiv.org/pdf/2306.16851.pdf)*.

## Introduction
SWAT mitigates system-wide leakage in encrypted data stores via a suite of algorithms with tunable privacy-efficiency trade-offs. We adopt [Intel SGX](https://github.com/intel/linux-sgx) as our trusted execution environment, and [Redis](https://redis.io/) severs as our backend storage engine. The experiemental study shows that SWAT effectively mitigates query-correlation leakage, while sustaining manageable increments in query latency compared to [Pancake](https://www.usenix.org/conference/usenixsecurity20/presentation/grubbs). Meanwhile, it showcases a highly competitive performance compared to [$\mathcal{E}$psolute](https://dl.acm.org/doi/abs/10.1145/3460120.3484786) and [ObliDB](https://dl.acm.org/doi/10.14778/3364324.3364331) with more comprehensive leakage mitigations. In short, SWAT progressively supports the following workloads: 
1. $\theta$-query decorrelation in frequency-smoothed key-value stores;
2. Nearly zero-leakage range queries;
3. Differentially oblivious dynamization.

## Before run
All experiments were conducted on a machine with an Intel(R) Xeon(R) Platinum 8369B CPU @ 2.90GHz of 32 physical cores, with SGXv2 enabled. The machine has 128GB RAM, of which about 64GB is enclaves' protected memory. It operates on Ubuntu 20.04 and uses SGX SDK version 2.19. 

### Prerequisites
We highlight that for SWAT, CPUs have to support SGX and BIOS should be properly configured. One may follow this [doc](https://github.com/intel/linux-sgx?tab=readme-ov-file#introduction) to enable SGX. Meanwhile, please follow the official [documentation](https://github.com/intel/sgx-ra-sample) to enable remote attestation. Besides, the following pacakages are also necessary. 
* cmake-3.15+
* gcc 
* libssl-dev
* redis-server
* ssdb-rocks
* libbzp2-dev
* zlib1g-dev

To install the above requirements on Ubuntu 20.04, run:
```bash
sudo apt-get update
sudo apt-get install cmake gcc libssl-dev redis-server libbz2-dev zlib1g-dev
```


### Compile the code
After installing the requirements, run the following commands to clone this repository and compile the source code: 

```bash
git clone https://github.com/55199789/SWAT.git
cd SWAT
./build.sh
```
After that, the `bin` folder will contain the executables of the client/server proxies as well as the signed enclave as follows: 
```bash
bin
├── clientProxy
├── enclave.signed.so
└── serverProxy
```

## Run

The server proxy has the following options: 
```bash
bin/serverProxy
        --storageType redis # Optional: storage server type (redis, rocksdb), default is redis
        --storageHost 127.0.0.1 # Optional: storage server host name, default is localhost
        --outPort 6379 # Optional: storage server port, default is 6379
        --inPort 12341 # Optional: input port (from the client), default is 12341
        --storageClientNum 32 # Optional: number of clients for the KV store, default is 32
```
To change the size of records (key-value pairs), please modify the macro `RECORD_LENGTH_IN_BYTE` in `build.sh` and re-build the binaries.

The client proxy has the following options: 
```bash
bin/clientProxy
         --serverName localhost # Optional: server proxy host name, default is locahost
         --serverPort 12341 # Optional: Server proxy port, default is 12341
         --seed 11451 # Optional: random seed for reproducing, default is 11451
         --epsilon 1.0 # Optional: differential obliviousness budget, default is 1.0
         --delta 1e-12 # Optional: differential obliviousness failure probability, default is 1e-12
         --k 8 # Optinal: number of layers for dynamization, default is 8
         --type int # Optional: type of keys (int, long long, float, double), default is int
         --storageOverhead 2 # Optional: storage overhead (used in pancake), default is 2
         --batch 3 # Optional: batch size for frequency smoothign, default is 3
         --updatePolicy increment # Optional: weight update policy (constant, increment, multiplication) for query correlation, default is increment
         --domainMin 1 # Optional: minimum value of the domain, default is 1
         --domainMax 10000 # Optional: maximum value of the domain, default is 10000
         --rangeLen 1 # Optional: maximum length of the range query [l, r], i.e., r-l+1, default is 1
         --pointReads 0 # Optional: fraction of point reads, default is 0
         --smallReads 0 # Optional: fraction of small reads, default is 0
         --largeReads 0 # Optional: fraction of large reads, default is 0
         --inserts 100 # Optional: fraction of inserts, default is 100
         --num 0 # Optional: number of initial key-value pairs, default is 0
         --queryNum 100000 # Optional: number of queries, default is 100000
         --threadNum 32 # Optional: number of threads, default is 32
```
The default parameters will run with an initially empty data store and inserting 100,000 key-value pairs, with integer keys uniformly sampled from [1, 10,000]. The seed for generating random numbers is 11451, the insertion will run in a $(1.0, 10^{-12})$-differentially oblivious way, and the dynamization maintains at most $8$ layers.  

In addition, we also provide Python scripts `./exp_scripts.py`, used in our experiments, to run SWAT against a range of parameters. Please note that only a minor portion of the results were presented in our paper. 

We also conduct a separate evaluation of SWAT against Pancake in the `query_decorrelation.py`, focusing on query decorrelation on key-value stores. This is because the full version of SWAT is designed to support range queries, a feature not supported by Pancake.

### About the code
``` bash
clientProxy
├── CMakeLists.txt
└── src
    ├── BinaryIndexedTree.hpp # used by the sampling pool
    ├── client.cpp # the entry (main) file of the client proxy 
    ├── Query.hpp # defines and implements query-related classes
    ├── SamplePool.hpp # defines and implements the sampling pool
    ├── SegmentTree.hpp # used by the sampling pool
    └── UniformQueryFactory.hpp # generates random queries uniformly
common # implements some remote attestation/enclave stuffs
├── ...
include
├── common # defines some remote attestation/enclave stuffs, excepts defs.h
│   ├── defs.h # some macros for default parameters, defines and implements the record (datum)
│   ├── ...
├── enclave
│   ├── Dist.hpp # some random variables used in differential obliviousness
│   ├── DOMerger.hpp # differentially oblivious merge (dynamization)
│   ├── DPInteriorPoint.hpp # differentially private interior point, used by DOMerger
│   ├── DPPrefixSum.hpp # differentially private prefix sum, used by DOMerger
│   ├── Index.hpp # differentially private prefix sum, used by DOMerger
│   ├── ObliviousSort.hpp # oblivious sorter, used by DOMerger
│   └── OHeap.hpp # oblivious heap, intended to replace ObliviousSort but not used in the end
└── libstorage # some definitions and implementations for the backend storage server
    ├── ...
serverProxy
├── CMakeLists.txt
├── enclave
│   ├── enclave.lds
│   ├── private.pem # private key to sign the enclave 
│   ├── server.config.xml # the enclave configuration file
│   ├── server.cpp  # the enclave entry file
│   └── server.edl # defines ecall/ocall functions
└── host
    └── host.cpp # the entry (main) file for the server proxy
```