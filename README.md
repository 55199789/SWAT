# [SWAT: A System-Wide Approach to Tunable Leakage Mitigation in Encrypted Stores](https://arxiv.org/pdf/2306.16851.pdf)

## Prerequisites
* [Intel SGX](https://github.com/intel/linux-sgx)
* [Intel SGX Remote Attestation](https://github.com/intel/sgx-ra-sample)
* cmake-3.15+
* gcc 
* libssl-dev
* redis-server
* ssdb-rocks
* libbzp2-dev
* zlib1g-dev

## Build
After installing the requirements, run: 

```./build.sh```

## Run
```
Client proxy: 
         --serverName: Server host name
         --serverPort: Server port
         --lambda: Security parameter
         --epsilon: Differential obliviousness budget
         --delta: Failure probability
         --seed: Random seed for reproducing
         --k: Number of layers for dynamization
         --type: Type of keys
         --storageOverhead: Storage overhead
         --batch: Batch size for frequency smoothign
         --updatePolicy: Update policy for frequency smoothing
         --domainMin: Minimum value of the domain
         --domainMax: Maximum value of the domain
         --rangeLen: Maximum length of the range
         --pointReads: Number of point reads
         --smallReads: Number of small reads
         --largeReads: Number of large reads
         --inserts: Number of inserts
         --num: Number of initial data
         --threadNum: Number of threads
```

```
Server proxy: 
         --storageHost: Storage server host name
         --storageType: Storage server type (redis, memcached, rocksdb)
         --outPort: Storage server port
         --inPort: Input port
         --storageClientNum: Number of clients for the KV store
```