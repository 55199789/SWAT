# %%
import numpy as np 
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import random as rm
import queue
import math
from matplotlib.cm import get_cmap
from scipy.linalg import eig
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.ticker as ticker
import tikzplotlib

# %%
A = np.array([[0.3, 0.65, 0.05], [0.9, 0, 0.1], [0.7, 0.3, 0]])
# A = np.array([[0.2, 0.7, 0.1], [0.8, 0, 0.2], [0.7, 0.3, 0]])
N = A.shape[0]  # number of states
# np.random.seed(114514)
np.random.seed(11451)
seeds = np.random.rand(10) * 10000
seeds = [int(i) for i in seeds]
P = A
for i in range(1000):
    P = P @ A
    if np.max(np.abs(P[0] - P[1])) < 1e-6:
        break
p_real = P[0]
start_state = np.random.choice([0, 1, 2], p=p_real)
n = 100000 
queries = [start_state]
initNum = 5

# %%
while len(queries) < n:
    state = np.random.choice([0, 1, 2], p=A[queries[-1]])
    queries.append(state)
cnt = np.zeros(3)
for i in queries: 
    cnt[i] += 1
print(cnt / n)

# %%
name = "tab20c"
cmap = get_cmap(name)  # type: matplotlib.colors.ListedColormap
colors = cmap.colors  # type: list
# colors = cm.rainbow(np.linspace(0, 1, 3))

fig = plt.figure(figsize=(3, 3))
x_labels = [r'$\mathit{k}_1$', r'$\mathit{k}_2$', r'$\mathit{k}_3$']
# creating the bar plot
# plt.axis('off')
plt.bar(x_labels, p_real, color=colors)

# plt.xlabel("Courses offered")
# plt.ylabel("")
plt.yticks([0, 0.25, 0.5])
# plt.xticks(x_labels, {'size': 12, 'weight': 'bold'})
plt.rc('xtick', labelsize=15)
plt.rc('figure', titlesize=15)
ax = plt.gca() 
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.title(r"$f_{real}$", fontsize=16)
# plt.show()
# plt.savefig("independencyFigures/real.pdf", format="pdf", bbox_inches='tight')
tikzplotlib.save("independencyFigures/real.tex")
# plt.close()

# %% [markdown]
# Plotting the 3D Markov Matrix

# %%
class Pancake:
    def __init__(self, p_real, keys_old, customiedQueue=queue.Queue(), initNum = 0):
        self.n = len(keys_old)
        self.p_real = p_real 
        self.keys_old = keys_old
        self.p_fake = []
        assert (len(p_real) == self.n)
        self.keys_new = []
        self.replicas = []
        n_ = 0 
        for i in range(self.n):
            self.replicas.append(int(np.ceil(p_real[i] * self.n)))
            for j in range(self.replicas[i]):
                self.keys_new.append(keys_old[i] + "," + str(j)+"}$")
                self.p_fake.append(
                    (1/self.n-p_real[i]/self.replicas[i]))
            n_ += self.replicas[i]
        for j in range(2*self.n-n_): 
            self.p_fake.append(1/self.n)
            self.keys_new.append(r"$k_{f,"+str(j)+"}$")
        self.delta = 1 / 2
        self.que = customiedQueue
        self.B = 3
        self.cnt = 0
        self.initNum = initNum
        print(self.p_fake)
        
    def batch(self, idx): 
        j = np.random.randint(0, self.replicas[idx])
        self.que.put((self.keys_old[idx] + "," + str(j)+"}$", self.cnt))
        latencies = []
        l = []
        for i in range(self.B):
            type = np.random.choice([0, 1], p=[self.delta, 1-self.delta])
            if type == 0:
                key_ = np.random.choice(self.keys_new, p=self.p_fake)
            else:
                if self.que.qsize() <= self.initNum:
                    idx = np.random.choice(np.arange(self.n), p=p_real)
                    j = np.random.randint(0, self.replicas[idx])
                    key_ = self.keys_old[idx] + "," + str(j)+"}$"
                else:
                    key_, cnt_ = self.que.get()
                    latencies.append(self.cnt - cnt_)
            l.append(key_)
        self.cnt+=1
        return l, latencies
    
    def getDict(self): 
        d = dict()
        for k in self.keys_new: 
            d[k] = dict() 
            for k_ in self.keys_new: 
                d[k][k_] = 0
        return d

# %%
keys_old = [r'$k_{1', r'$k_{2', r'$k_{3']
pancake = Pancake(p_real, keys_old, initNum=initNum)

# %%
def drawFromToHeatMap(A, labels, filaPath, title, yColor='white', min_=0, max_=1):
    fig, axe = plt.subplots(figsize=(8, 8))
    x_labels = labels
    y_labels = labels
    axe.set_xticks(np.arange(len(x_labels)))
    axe.set_yticks(np.arange(len(y_labels)))
    axe.set_xticklabels(x_labels, fontsize=24)
    axe.set_xlabel(r'To $q_{i+1}$', fontsize=30)
    # if yColor != 'white':
    axe.set_yticklabels(y_labels, fontsize=24, color=yColor)
    axe.set_ylabel(r'From $q_i$', fontsize=30, color=yColor)
    # else:
        # axe.yaxis.set_visible(False)
        # axe.spines['left'].set_visible(False)
    im = axe.imshow(A, interpolation='nearest',
                    vmin=min_, vmax=max_)

    # divider = matplotlib.colorbar.make_axes(axe)
    # cax = divider.append_axes("right", size="5%", pad=0.05)

    # if title == 'Independent Queries':
    cbar = fig.colorbar(im, ax=axe)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(24)
    plt.title(title, fontsize=32)
    # plt.show()
    plt.savefig(filaPath, format="pdf", bbox_inches='tight')
    plt.close()

# %%
def sequencesTest(pancake, queries, filename, title, procNum = 0, resultDict=None, yColor='white'):
    freq = dict()
    markovFreq = pancake.getDict()
    prev = -1
    latencies = [] 
    for q in queries:
        batch, latencies_ = pancake.batch(q)
        latencies += latencies_
        for k in batch:
            if prev != -1:
                markovFreq[prev][k] += 1
            prev = k
            if k not in freq:
                freq[k] = 0
            freq[k] += 1
                
    for key in freq:
        freq[key] /= n * 3
    for key in markovFreq:
        for key_ in markovFreq[key]:
            markovFreq[key][key_] /= n * 3

    name = "tab20c"
    cmap = get_cmap(name)  # type: matplotlib.colors.ListedColormap
    colors = cmap.colors  # type: list
    # colors = cm.rainbow(np.linspace(0, 1, 3))

    fig = plt.figure(figsize=(6, 3))
    x_labels = pancake.keys_new
    freqVec = [] 
    for k in x_labels:
        freqVec.append(freq[k])
    # creating the bar plot
    # plt.axis('off')
    plt.yticks([])
    plt.ylabel("Frequency", fontsize=14)
    plt.bar(x_labels, freqVec, color=colors)

    # plt.xlabel("Courses offered")
    # plt.ylabel("")
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # plt.title(r"$$")
    # plt.savefig("independencyFigures/freq_keys" +
    #             filename  + ".pdf", format="pdf", bbox_inches='tight')
    tikzplotlib.save("independencyFigures/freq_keys" +
                filename  + ".tex")
    plt.close() 
    markovFreqMatrix = np.zeros((len(pancake.keys_new), len(pancake.keys_new)))
    m = len(pancake.keys_new)
    for i in range(m):
        for j in range(m):
            markovFreqMatrix[i][j] = markovFreq[pancake.keys_new[i]
                                                ][pancake.keys_new[j]]

    drawFromToHeatMap(markovFreqMatrix, pancake.keys_new,
                      "independencyFigures/heatMap" + filename + ".pdf", title, yColor,
                      1/m * 1/m * 0.8, 1/m * 1/m * 1.16)
    if resultDict is None: 
        return latencies, np.std(markovFreqMatrix) / np.mean(markovFreqMatrix)
    else:
        resultDict[procNum] = (latencies, np.std(markovFreqMatrix)/ np.mean(markovFreqMatrix))

# %%
# keys_old = [r'$k_{1', r'$k_{2', r'$k_{3']
# pancake = Pancake(p_real, keys_old, SamplingPool(), initNum)
# fMin, fMax, fSTD, fRSD = sequencesTest(pancake, queries, "MarkovChainWithSamplingPool" + str(initNum))
# print(fMin, fMax, fSTD, fRSD)
from multiprocess import Process
from multiprocessing import Manager

def pad(x): return x + [np.mean((x))] * (n-len(x))

# %%
keys_old = [r'$k_{1', r'$k_{2', r'$k_{3']
fLatencies = [[] for x in range(10)]
fRSD = [[] for x in range(10)]
for s in seeds:
    np.random.seed(s)
    threads = []
    manager = Manager()
    resultDict = manager.dict()
    for initNum in range(10):
        pancake = Pancake(p_real, keys_old, initNum = initNum)
        threads.append(Process(target=sequencesTest,
                            args=(pancake, queries, 
                                    "MarkovChain", "Correlated Queries", 
                                    initNum, resultDict, 
                                    'black')))
        threads[-1].start()
    for t in threads:
        t.join()
    for procNum in range(10):
        latencies, sRSD_ = resultDict[procNum]
        fRSD[procNum].append(sRSD_)
        fLatencies[procNum].append(pad(latencies))

# %%
# print(np.array(fLatencies).shape)
print(np.array(fRSD).shape)
print(len(seeds))

# %%
keys_old = [r'$k_{1', r'$k_{2', r'$k_{3']
indLatencies = [[] for x in range(10)]
indRSD = [[] for x in range(10)]
for seed in seeds:
    np.random.seed(seed)
    threads = []
    manager = Manager()
    resultDict = manager.dict()
    independentQueries = np.random.choice(np.arange(3), size=n, p=p_real)
    for initNum in range(10):
        pancake = Pancake(p_real, keys_old, initNum=initNum)
        threads.append(Process(target=sequencesTest,
                            args=(pancake, independentQueries, 
                                    "Independent", "Independent Queries",
                                    initNum, resultDict)))
        threads[-1].start()
    for t in threads:
        t.join()
    for procNum in range(10):
        latencies, sRSD_ = resultDict[procNum]
        indRSD[procNum].append(sRSD_)
        indLatencies[procNum].append(pad(latencies))


# %%
class SamplingPool: 
    def __init__(self, initNum = 0, update = None):
        self.pool = [] 
        self.initNum = initNum
        self.update = update
        self.weight = []
    
    def get(self):
        ret = 0
        if self.update is None:
            t = np.random.randint(0, len(self.pool))
            ret = self.pool.pop(t)
        elif self.update == "linear":
            t = np.random.choice(np.arange(len(self.pool)), 
                                 p=np.array(self.weight) / np.sum(self.weight))
            ret = self.pool.pop(t)
            self.weight.pop(t)
            self.weight = [x + 3 for x in self.weight]
        elif self.update == "exp":
            t = np.random.choice(np.arange(len(self.pool)),
                                 p=np.array(self.weight) / np.sum(self.weight))
            ret = self.pool.pop(t)
            self.weight.pop(t)
            self.weight = [x * 5 for x in self.weight]
        return ret 
    def qsize(self):
        return len(self.pool)
    
    def empty(self) -> bool: 
        return len(self.pool) <= self.initNum
    
    def put(self, x): 
        self.pool.append(x) 
        self.weight.append(1)

# %%
keys_old = [r'$k_{1', r'$k_{2', r'$k_{3']
sRSD = [[] for x in range(10)]
lLatencies = [[] for x in range(10)]
for seed in seeds: 
    np.random.seed(seed)
    threads = []
    manager = Manager()
    resultDict = manager.dict()
    for initNum in range(10):
        pancake = Pancake(p_real, keys_old, SamplingPool(initNum), initNum)
        threads.append(Process(target=sequencesTest,
                            args=(pancake, queries, \
                                "MarkovChainWithSamplingPool" + str(initNum), str(initNum)+"-query Decorrelation", 
                                initNum, resultDict)))
        threads[-1].start()
    for t in threads:
        t.join() 
    for procNum in range(10):
        latencies, sRSD_ = resultDict[procNum]
        sRSD[procNum].append(sRSD_)
        lLatencies[procNum].append(pad(latencies))

# %%
sRSDLinear = [[] for x in range(10)]
lLatenciesLinear = [[] for x in range(10)]
for seed in seeds:
    np.random.seed(seed)
    threads = []
    manager = Manager()
    resultDict = manager.dict()
    for initNum in range(10):
        pancake = Pancake(p_real, keys_old, SamplingPool(
            initNum, 'linear'), initNum)
        threads.append(Process(target=sequencesTest,
                            args=(pancake, queries,
                        "MarkovChainWithSamplingPool" +
                        str(initNum)+"Linear", str(initNum) +
                        "-query Decorrelation",
                        initNum, resultDict)))
        threads[-1].start()
        # sRSDLinear.append(sRSD_)
        # lLatenciesLinear.append(lLatencies_)

    for t in threads:
        t.join()
    for procNum in range(10):
        latencies, sRSD_ = resultDict[procNum]
        sRSDLinear[procNum].append(sRSD_)
        lLatenciesLinear[procNum].append(pad(latencies))

# %%
sRSDExp = [[] for x in range(10)]
lLatenciesExp = [[] for x in range(10)]
for seed in seeds:
    np.random.seed(seed)
    threads = []
    manager = Manager()
    resultDict = manager.dict()
    for initNum in range(10):
        pancake = Pancake(p_real, keys_old, SamplingPool(initNum, 'exp'), initNum)
        threads.append(Process(target=sequencesTest,
                            args=(pancake, queries,
                                    "MarkovChainWithSamplingPool" +
                                    str(initNum)+"Exp", str(initNum) +
                                    "-query Decorrelation",
                                    initNum, resultDict)))
        threads[-1].start()
        # sRSDExp.append(sRSD_)
        # lLatenciesLinear.append(lLatencies_)
    for t in threads:
        t.join()
    for procNum in range(10):
        latencies, sRSD_ = resultDict[procNum]
        sRSDExp[procNum].append(sRSD_)
        lLatenciesExp[procNum].append(pad(latencies))

# %%
T = 5
fRSD_ = np.array([100 * np.array(x) for x in fRSD])
sRSD_ = np.array([100 * np.array(x) for x in sRSD])
sRSDLinear_=np.array([100 * np.array(x)for x in sRSDLinear])
sRSDExp_=np.array([100 * np.array(x) for x in sRSDExp])
indRSD_ = np.array([100 * np.array(x) for x in indRSD])
markers = ['o', 's', 'x', 'p', '*']
linestyles = ['--', '-', '-', '-', ':']
# plt.figure(figsize=(9, 6))
# initNums = np.arange(10)[:T]
# plt.xticks(initNums, fontsize=24)
# plt.yticks(fontsize=24)
# plt.plot(initNums, fRSD_[:T], marker='o',
#          linestyle='--',  label="Correlated Queries")
# plt.plot(initNums, sRSDExp_[:T], marker='s',
#          label=r"$\theta$-query Decorrelation (Exp)")
# plt.plot(initNums, sRSDLinear_[:T], marker='x',
#          label=r"$\theta$-query Decorrelation (Linear)")
# plt.plot(initNums, sRSD_[:T], marker='p',
#          label=r"$\theta$-query Decorrelation (Constant)")
# plt.plot(initNums, indRSD_[:T], marker='*', linestyle=':',  label=r"Independent Queries")
# plt.xlabel(r"$\theta$", fontsize=30)
# plt.legend(fancybox=True, loc="best", shadow=True, fontsize=17, framealpha=0.5)
# plt.ylabel(r"RSD%", fontsize=30)
# # plt.title("RSD vs. Different Sampling Pool Sizes", fontsize=16)
# # plt.show()
# plt.savefig('independencyFigures/RSD.pdf', format='pdf', bbox_inches='tight')

# %%
T = 5
initNums = np.arange(10)[:T]
correlatedLatencies = pd.DataFrame({})
indenpendentLatencies = pd.DataFrame({})
constantLatencies = pd.DataFrame({})
linearLatencies = pd.DataFrame({})
expLatencies = pd.DataFrame({})


# def pad(x): return x + [np.mean((x))] * (n-len(x))
# for i in initNums:
#     fLatencies[i] = pad(fLatencies[i])
#     lLatencies[i] = pad(lLatencies[i])
#     lLatenciesLinear[i] = pad(lLatenciesLinear[i])
#     lLatenciesExp[i] = pad(lLatenciesExp[i])
#     indLatencies[i] = pad(indLatencies[i])
#     print(i)
#     print("min: ", min(fLatencies[i]), min(lLatencies[i]), min(lLatenciesLinear[i]), min(lLatenciesExp[i]), min(indLatencies[i]))
#     print("max: ", max(fLatencies[i]), max(lLatencies[i]), max(
#         lLatenciesLinear[i]), max(lLatenciesExp[i]), max(indLatencies[i]))
#     print("mean: ", np.mean(fLatencies[i]), np.mean(lLatencies[i]), np.mean(
#         lLatenciesLinear[i]), np.mean(lLatenciesExp[i]), np.mean(indLatencies[i]))


for i in initNums: 
    correlatedLatencies[i] = np.mean(fLatencies[i], axis=0)
    constantLatencies[i] = np.mean(lLatencies[i], axis=0)
    linearLatencies[i] = np.mean(lLatenciesLinear[i], axis=0)
    expLatencies[i] = np.mean(lLatenciesExp[i], axis=0)
    indenpendentLatencies[i] = np.mean(indLatencies[i], axis=0)

datasets = [correlatedLatencies, expLatencies,
            linearLatencies, constantLatencies,
            indenpendentLatencies]
rsdDatasets = [np.mean(fRSD_, axis=1),
               np.mean(sRSDExp_, axis=1),
               np.mean(sRSDLinear_, axis=1),
               np.mean(sRSD_, axis=1),
               np.mean(indRSD_, axis=1)]
stdDatasets = [np.std(fRSD_, axis=1),
                np.std(sRSDExp_, axis=1),
                np.std(sRSDLinear_, axis=1),
                np.std(sRSD_, axis=1),
                np.std(indRSD_, axis=1)]


# %%
from matplotlib.patches import Patch
x_pos_range = np.arange(len(datasets)) / (len(datasets) - 1)
x_pos = (x_pos_range * 0.5) + 0.75
# Plot

name = "tab20c"
cmap = get_cmap(name)  # type: matplotlib.colors.ListedColormap
colors = cmap.colors  # type: list
groups = [r"Correlated", r"Exp", r"Linear",
                   r"Constant", "Independent"]
fig, ax1 = plt.subplots(figsize=(10, 5))
ax2 = ax1.twinx()

print(x_pos)

for i, data in enumerate(datasets):
    bp = ax1.boxplot(
        np.array(data), notch = True, whis =0.8,  sym='', widths=0.6 / len(datasets),
        labels=list(datasets[0]), patch_artist=True, meanline=True,
        positions=[x_pos[i] + j * 1 for j in range(len(data.T))]
    )
    print(i, [x_pos[i] + j * 1 for j in range(len(data.T))])

    # yData = np.array([fRSD[i], sRSDExp[i], sRSDLinear[i], sRSD[i], indRSD[i]])
    # xData = [i * 1 + x for x in x_pos]
    xData = [np.mean(x_pos) + j * 1 for j in range(len(x_pos))]
    yData = rsdDatasets[i][:T]
    print(" " , i, xData)
    ax2.plot(xData, yData, color=colors[i],  alpha=1,  marker=markers[i], linestyle=linestyles[i],
            label = groups[i],
            # , widths=0.6 / len(datasets),
            # positions=[x_pos[i] + j * 1 for j in range(len(data.T))]
             )
    ax2.fill_between(xData, yData - stdDatasets[i][:T], yData + stdDatasets[i][:T], color=colors[i], alpha=0.1)
    for patch in bp['boxes']:
        patch.set(facecolor=colors[i], alpha=1)

legend_elements = []
legend_elements.append(
    Patch(edgecolor='black', facecolor=None, label="Latency", alpha=1))
# for i in range(len(datasets)):
#     j = i % len(groups)
#     k = i % len(colors)
    # legend_elements.append(
    #     Patch(color=colors[k], linestyle=linestyles[i], marker=markers[i],  label=groups[j], alpha=0.5))
ax2.set_ylabel('Latency', fontsize=20)
ax2.tick_params(axis='y', labelsize=16)
# Titles
# plt.title('')
# ax1.set_yticklabels(fontsize=16)
ax1.set_xlabel(r'$\theta$', fontsize=20)
ax1.set_ylabel('RSD (%)', fontsize=20)
ax1.tick_params(axis='y', labelsize=16)
# plt.ylabel('Latency', fontsize=20)
# Axis ticks and labels
plt.xticks(np.arange(len(list(datasets[0]))) + 1)
plt.gca().xaxis.set_minor_locator(ticker.FixedLocator(
    np.array(range(len(list(datasets[0])) + 1)) + 0.5)
)
# plt.gca().tick_params(axis='x', which='minor', length=4)
# plt.gca().tick_params(axis='x', which='major', length=0)
# Change the limits of the x-axis
plt.xlim([0.5, len(list(datasets[0])) + 0.5])
# ax1.set_yticks(np.linspace(0, 10, num=11, endpoint=True))
legend1 = plt.legend(fancybox=True, loc="upper left", shadow=False, fontsize=15,
           framealpha=0.5,
            labelcolor='linecolor', 
            title="RSD", alignment = 'center', 
                     title_fontsize='15', ncol=2,
        )
        # ,       handles = legend_elements)
plt.legend(handles=legend_elements, loc="upper right",
           shadow=False, fontsize=15,)
plt.gca().add_artist(legend1)

# plt.tick_params(labelsize=14)
plt.rc("xtick", labelsize=22)
# plt.yscale('log', base = 2)
# plt.yticks(np.linspace(0, 8, num=50, endpoint=True), fontsize=24)
# plt.yticks(np.logspace(-10, 3, num=14, endpoint=True, base=2), fontsize=24)
# plt.show()

plt.savefig('independencyFigures/RSDvsLatency.pdf', format='pdf', bbox_inches='tight')

# %%


# %%



