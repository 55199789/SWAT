import numpy as np
import os
import subprocess
import time
import math
import functools
import redis
import rocksdb

lams = [256, 512, 1024, 2048, 4096]
bucketCapacities = [256, 512, 768, 1024, 1280, 1536, 1792, 2048, 2304, 2560]
nums = [1000, 10000, 100000, 1000000]  # , 10000000]  # , 10000000]
epsilons = [0.001, 0.01, 0.1, 0.25,  0.5, 0.75, 1.0, 2.0, 10.0]
ks = range(4, 13)
# batches = [5, 10, 20, 40, 80, 160, 320, 500, 640, 800, 1000]
threads = [1, 2, 4, 8, 16, 32, 64]
# storageOverheads = [1.1, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
updatePolicies = ["increment", "constant"]  # , "multiplication" is buggy
seeds = [28603, 31849, 17669, 1435837, 532973, 114514]
storages = {"redis": 6379,
            "rocksdb": 3337}
types = ["int", "long long", "float", "double"]
# queryNums = [500, 1000, 5000, 10000, 100000, 1000000]
# domainMins = [1]
domainMaxs = [1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000]
workLoads = {
    "L1": (5.0, 5.0, 0.0, 90.0),   # L1
    "L2": (5.0, 50.0, 5.0, 40.0),  # L2
    "L3": (90.0, 5.0, 5.0, 0.0),   # L3
    "L4": (5.0, 90.0, 5.0, 0.0),   # L4
    "L5": (5.0, 5.0, 90.0, 0.0),    # L5
    "L6": (0.0, 0.0, 0.0, 100.0)    # pure insert
}
selectivities = [0.001, 0.0025, 0.005, 0.01, 0.02]
# rangeLens = [0, 100, 1000, 10000, 100000,
#              1000000, 10000000, 100000000, 1000000000]
recodLengthInBytes = [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

r = redis.Redis()


def Geom(Z, alpha):
    pmf = [(alpha-1)/(alpha+1-2*alpha**(-Z/2))*alpha**(-abs(x-Z/2))
           for x in range(Z+1)]
    return np.array(pmf)


def GeomB(Z, alpha, B):
    base = Geom(Z, alpha)
    cur = Geom(Z, alpha)
    for b in range(B-1):
        cur = np.convolve(cur, base)
    assert (len(cur) == Z*B+1)
    return cur


def getCmp(B, eps, lam):
    alpha = math.exp(eps)
    delta = math.exp(-math.log(lam)**(2))
    lZ = 2
    rZ = (math.ceil(math.log(lam)**5 / eps) | 1) + 1
    while lZ+2 < rZ:
        mZ = (lZ+rZ) // 2
        if mZ % 2 != 0:
            mZ = mZ + 1
        pmf = GeomB(mZ, alpha, B)
        if np.sum(pmf[:mZ]) >= delta:
            lZ = mZ + 2
        else:
            rZ = mZ
    return ((math.floor(math.log(lam)**5 / eps) | 1)+1, lZ)


@functools.lru_cache(maxsize=None)
def getBinCapacity(eps, lam):
    return getCmp(3, eps, lam)[1]


def genFileName(lam, num, epsilon, k, batch, thread, storageOverhead,
                storageType, type, queryNum, domainMin, domainMax, rangeLen,
                updatePolicy, selectivity,
                bucketCapacity, recordLen, omerge,
                seed):
    ret = "lambda" + str(lam) + \
        "_num" + str(num) + \
        "_epsilon" + str(epsilon) + \
        "_k" + str(k) + \
        "_batch" + str(batch) + \
        "_thread" + str(thread) + \
        "_storageOverhead" + str(storageOverhead) + \
        "_storageType" + str(storageType) + \
        "_type" + str(type) + \
        "_queryNum" + str(queryNum) + \
        "_domainMin" + str(domainMin) + \
        "_domainMax" + str(domainMax) + \
        "_rangeLen" + str(rangeLen) + \
        "_updatePolicy" + str(updatePolicy) + \
        "_selectivity" + str(selectivity) + \
        "_bucketCapacity" + str(bucketCapacity) + \
        "_recordLen"+str(recordLen) + \
        "_omerge"+str(omerge) + \
        "_seed" + str(seed) + ".txt"
    return ret


def runExp(lam, num, epsilon, k, batch, thread, storageOverhead,
           storageType, type, queryNum, domainMin, domainMax, rangeLen,
           updatePolicy,
           bucketCapacity, recordLen, omerge,
           outPort,
           selectivity,
           workloadName=None):
    pwd = os.getcwd() + "/"
    print(locals())
    for seed in seeds:
        if workloadName is None:
            fileName = genFileName(lam, num, epsilon, k, batch, thread, storageOverhead,
                                   storageType, type, queryNum, domainMin, domainMax, rangeLen,
                                   updatePolicy, selectivity,
                                   bucketCapacity, recordLen, omerge, seed)
        else:
            _pointRead, _smallRead, _largeRead, _insert = workLoads[workloadName]
            fileName = genFileName(lam, num, epsilon, k, batch, thread, storageOverhead,
                                   storageType, type, queryNum, domainMin, domainMax, workloadName,
                                   updatePolicy, selectivity,
                                   bucketCapacity, recordLen, omerge, seed)
        if os.path.isfile("dataClient/"+fileName):
            continue

        r.flushall()
        os.system("rm /tmp/medsdb")
        os.system(pwd+'closePort.sh')
        g = open("dataServer/" + fileName, "w+")
        content = g.readlines()
        p = subprocess.Popen(
            [pwd+"bin/serverProxy",
                " --storageType=" + str(_storageType) +
                " --outPort=" + str(outPort)
             ],
            stdout=g,
            close_fds=True
        )
        time.sleep(10)

        if workloadName is None:
            os.system("bin/clientProxy" +
                      " --bucketCapacity=" + str(bucketCapacity) +
                      " --lambda=" + str(lam) +
                      " --epsilon="+str(epsilon) +
                      " --k="+str(k) +
                      " --batch="+str(batch) +
                      " --threadNum="+str(thread) +
                      " --storageOverhead=" + str(storageOverhead) +
                      " --seed="+str(seed) +
                      " --num="+str(num) +
                      " --queryNum="+str(queryNum) +
                      " --rangeLen="+str(rangeLen) +
                      " --domainMin="+str(domainMin) +
                      " --domainMax="+str(domainMax) +
                      " --type="+str(type) +
                      " --updatePolicy="+str(updatePolicy) +
                      " >> dataClient/" + fileName
                      )
        else:
            os.system("bin/clientProxy" +
                      " --bucketCapacity=" + str(bucketCapacity) +
                      " --lambda=" + str(lam) +
                      " --epsilon="+str(epsilon) +
                      " --k="+str(k) +
                      " --batch="+str(batch) +
                      " --threadNum="+str(thread) +
                      " --storageOverhead=" + str(storageOverhead) +
                      " --seed="+str(seed) +
                      " --num="+str(num) +
                      " --queryNum="+str(queryNum) +
                      " --domainMin="+str(domainMin) +
                      " --domainMax="+str(domainMax) +
                      " --type="+str(type) +
                      " --updatePolicy="+str(updatePolicy) +
                      " --pointReads="+str(_pointRead) +
                      " --smallReads="+str(_smallRead) +
                      " --largeReads="+str(_largeRead) +
                      " --inserts="+str(_insert) +
                      " >> dataClient/" + fileName
                      )


def computeRangeLen(selectivity, domainMin, domainMax):
    return int(selectivity * (domainMax - domainMin))


def computeBatchSize(num, selectivity, bucketCapacity):
    return math.ceil(num * selectivity / 2 / bucketCapacity * 2.5)


# iterate files in the folder dataClient
# brokenFileCnt = 0
# for dirpath, dnames, fnames in os.walk("dataClient"):
#     print(len(fnames))
#     for fname in fnames:
#         with open("dataClient/"+fname, "r") as f:
#             content = f.readlines()
#             flag = False
#             if len(content) > 26:
#                 for l in content:
#                     if l.find("Client:  success") != -1:
#                         flag = True
#                         break
#                 if flag:
#                     continue
#             brokenFileCnt += 1
#             print("dataClient/"+fname)
#             # print("Input any key to delete the file")
#             # input()
#             # os.system("rm dataClient/"+fname)
#             # os.system("rm dataServer/"+fname)
# print("Wrong file cnt: " + str(brokenFileCnt))
# print("Input any key to continue...")
# input()

cnt = 0
os.system("rm -rf data")
os.system("mkdir data/")
_lam = 512
_k = 8
_epsilon = 1.0
_num = 1000000
_storageType = "redis"
_thread = 32
_outPort = 6379
_storageOverhead = 2.0
_queryNum = 10
_type = "int"
_domainMin = 1
_domainMax = 1000000
_updatePolicy = "constant"
_selectivity = 0.005
_bucketCapacity = 512
_omerge = False
_recordLen = 4096
_rangeLen = computeRangeLen(_selectivity, _domainMin, _domainMax)
_batch = computeBatchSize(_num, _selectivity, _bucketCapacity)


# search over a static set-up data

os.system(
    "cmake -DRECORD_LENGTH_IN_BYTE=%d -DOBLIVIOUSMERGE=OFF -S . -B build " % (_recordLen))
os.system("cmake --build build -j")

runExp(_lam, _num, _epsilon, _k, _batch, _thread, _storageOverhead,
       _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
       _updatePolicy,
       _bucketCapacity, _recordLen, _omerge,
       _outPort,
       _selectivity)

print("Input any key to continue...")
input()
for num in nums:
    batch = computeBatchSize(num, _selectivity, _bucketCapacity)
    runExp(_lam, num, _epsilon, _k, batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity)

for bucketCapacity in bucketCapacities:
    batch = computeBatchSize(_num, _selectivity, bucketCapacity)
    runExp(_lam, _num, _epsilon, _k, batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity)

for selectivity in selectivities:
    rangeLen = computeRangeLen(selectivity, _domainMin, _domainMax)
    batch = computeBatchSize(_num, selectivity, _bucketCapacity)
    runExp(_lam, _num, _epsilon, _k, batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           selectivity)

for thread in threads:
    runExp(_lam, _num, _epsilon, _k, _batch, thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity)

for recordLen in recodLengthInBytes:
    os.system(
        "cmake -DRECORD_LENGTH_IN_BYTE=%d -S . -B build " % (recordLen))
    os.system("cmake --build build -j")

    runExp(_lam, _num, _epsilon, _k, _batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, recordLen, _omerge,
           _outPort,
           _selectivity)

    print("Input any key to continue...")
    input()

os.system(
    "cmake -DRECORD_LENGTH_IN_BYTE=%d -S . -B build " % (_recordLen))
os.system("cmake --build build -j")

for domainMax in domainMaxs:
    rangeLen = computeRangeLen(_selectivity, _domainMin, domainMax)
    runExp(_lam, _num, _epsilon, _k, _batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, domainMax, rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity)


_queryNum = 0
for lam in lams:
    runExp(lam, _num, _epsilon, _k, _batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity)

for epsilon in epsilons:
    runExp(_lam, _num, epsilon, _k, _batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity)

# DO updates
_queryNum = 1000000
for lam in lams:
    runExp(lam, 0, _epsilon, _k, _batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity,
           "L6")

for epsilon in epsilons:
    runExp(_lam, _num, epsilon, _k, _batch, _thread, _storageOverhead,
           _storageType, _type, 0, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity,
           "L6")

for queryNum in nums:
    runExp(_lam, 0, _epsilon, _k, _batch, _thread, _storageOverhead,
           _storageType, _type, queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity,
           "L6")

for k in ks:
    runExp(_lam, 0, _epsilon, k, _batch, _thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity,
           "L6")

for workLoad in workLoads:
    if workLoad == "L6":
        continue
    for storageType, outPort in storages.items():
        runExp(_lam, _num, _epsilon, _k, _batch, _thread, _storageOverhead,
               storageType, _type, 100, _domainMin, _domainMax, _rangeLen,
               _updatePolicy,
               _bucketCapacity, _recordLen, _omerge,
               outPort,
               _selectivity,
               workLoad)

for thread in threads:
    runExp(_lam, 0, _epsilon, _k, _batch, thread, _storageOverhead,
           _storageType, _type, _queryNum, _domainMin, _domainMax, _recordLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           _outPort,
           _selectivity,
           "L6")

for storageType, outPort in storages.items():
    runExp(_lam, 0, _epsilon, _k, _batch, thread, _storageOverhead,
           storageType, _type, _queryNum, _domainMin, _domainMax, _recordLen,
           _updatePolicy,
           _bucketCapacity, _recordLen, _omerge,
           outPort,
           _selectivity,
           "L6")

os.system(
    "cmake -DOBLIVIOUSMERGE=ON -S . -B build ")
os.system("cmake --build build -j")
runExp(_lam, 0, _epsilon, _k, _batch, _thread, _storageOverhead,
       _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
       _updatePolicy,
       _bucketCapacity, _recordLen, True,
       _outPort,
       _selectivity,
       "L6")

_queryNum = 100
runExp(_lam, _num, _epsilon, _k, _batch, _thread, _storageOverhead,
       _storageType, _type, _queryNum, _domainMin, _domainMax, _rangeLen,
       _updatePolicy,
       _bucketCapacity, _recordLen, True,
       _outPort,
       _selectivity)
