#pragma once 
#include <assert.h>
#include <algorithm>
#include <limits>
#include <list>
#include <random>

#define OHEAP_ROOT_SIZE 8
#define OHEAP_NONROOT_SIZE 3

template <typename Key>
struct Tuple
{
    static const time_t DummyValue = std::numeric_limits<time_t>::max();
    Key key;
    uint32_t pos;
    time_t timeStamp;

    Tuple()
    {
        reset();
    }

    Tuple(const Key& key,
        uint32_t pos,
        uint32_t timeStamp)
        : key(key),
        pos(pos),
        timeStamp(timeStamp)
    {
    }

    Tuple(const Tuple<Key>& o)
        : key(o.key),
        pos(o.pos),
        timeStamp(o.timeStamp)
    {
    }

    Tuple(Tuple<Key>&& o)
        : key(std::move(o.key)),
        pos(std::exchange(o.pos, 0)),
        timeStamp(std::exchange(o.timeStamp, std::numeric_limits<time_t>::max()))
    {
    }

    bool dummy() const
    {
        return timeStamp == DummyValue;
    }

    void reset()
    {
        timeStamp = DummyValue;
    }

    bool operator<(const Tuple& another) const
    {
        if (dummy())
            return false;
        if (another.dummy())
            return true;
        return key < another.key ||
            key == another.key && timeStamp < another.timeStamp;
    }

    bool operator>(const Tuple& another) const
    {
        if (another.dummy())
            return false;
        if (dummy())
            return true;
        return key > another.key ||
            key == another.key && timeStamp > another.timeStamp;
    }

    Tuple<Key>& operator=(const Tuple<Key>& another)
    {
        key = another.key;
        pos = another.pos;
        timeStamp = another.timeStamp;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Tuple<Key>& t)
    {
        if (t.dummy())
            os << "(dummy)";
        else
            os << "(" << t.key << ", " << t.pos << ", " << t.timeStamp << ")";
        return os;
    }
};

template <typename Key>
class Bucket
{
private:
    uint32_t num;
    std::vector<Tuple<Key>> tuples;
    Tuple<Key> subTreeMin;

public:
    Bucket(const uint32_t bucketSize = 0)
        : num(0),
        tuples(bucketSize)
    {
    }

    Bucket(const Bucket& o)
        : num(o.num),
        tuples(o.tuples),
        subTreeMin(o.subTreeMin)
    {
    }

    Bucket(Bucket&& o)
        : num(std::exchange(o.num, 0)),
        tuples(std::move(o.tuples)),
        subTreeMin(std::move(o.subTreeMin))
    {
    }

    Bucket& operator=(const Bucket& o)
    {
        num = o.num;
        tuples = o.tuples;
        subTreeMin = o.subTreeMin;
        return *this;
    }

    const Tuple<Key>& getSubTreeMin() const
    {
        return subTreeMin;
    }

    void add(const Tuple<Key>& tuple)
    {
        tuples[num++] = tuple;
    }

    void del(uint32_t timeStamp)
    {
        for (auto iter = tuples.begin(); iter != tuples.end(); ++iter)
            if (iter->timeStamp == timeStamp)
            {
                iter->reset();
                if (num > (iter - tuples.begin() + 1))
                {
                    std::swap(*iter, tuples[num - 1]);
                }
                --num;
                break;
            }
    }

    void clear()
    {
        for (auto& t : tuples)
            t.reset();
    }

    std::vector<Tuple<Key>>& getTuples()
    {
        return tuples;
    }

    void updateMin()
    {
        subTreeMin.reset();
        for (const auto& t : tuples)
            subTreeMin = std::min(subTreeMin, t);
    }

    void updateMin(const Tuple<Key>& minChild)
    {
        updateMin();
        if (minChild < subTreeMin)
            subTreeMin = minChild;
    }

    void updateNum(uint32_t num)
    {
        this->num = num;
    }

    friend std::ostream& operator<<(std::ostream& os, const Bucket<Key>& bucket)
    {
        os << bucket.subTreeMin << " || ";
        for (const auto& t : bucket.tuples)
            os << t << " ";
        return os;
    }
};

template <typename Key>
class OHeap
{
private:
    uint32_t curTimeStamp;

    const uint32_t levels;
    const uint32_t blockNums;

    std::mt19937 gen;
    std::uniform_int_distribution<uint32_t> unif;

    // Array-like heap, x -> left(x*2), right(x*2 + 1)
    std::vector<Bucket<Key>> buckets;

    // pos is the leaf node index
    void readNRm(uint32_t pos, uint32_t timeStamp)
    {
        while (pos)
        {
            buckets[pos].del(timeStamp);
            pos >>= 1;
        }
    }

    void evictUpdateMin(uint32_t pos)
    {
        std::list<Tuple<Key>> stash;
        uint32_t idx = pos;
        for (uint32_t l = 0; l < levels; l++)
        {
            std::vector<Tuple<Key>>& tuples = buckets[idx].getTuples();
            for (auto& t : tuples)
                if (!t.dummy())
                {
                    stash.push_back(t);
                    t.reset();
                }
                else
                    break;
            idx >>= 1;
        }

        idx = pos;
        for (uint32_t l = 0; l < levels; l++)
        {
            uint32_t lb = idx << l;
            uint32_t ub = lb + (1 << l);
            const uint32_t quota = idx == 1 ? OHEAP_ROOT_SIZE : OHEAP_NONROOT_SIZE;
            uint32_t cnt = 0;
            auto it = stash.begin();
            std::vector<Tuple<Key>>& tuples = buckets[idx].getTuples();
            while (it != stash.end())
            {
                if (it->pos >= lb && it->pos < ub)
                {
                    tuples[cnt++] = *it;
                    it = stash.erase(it);
                    if (cnt >= quota)
                        break;
                }
                else
                    it++;
            }

            buckets[idx].updateNum(cnt);
            idx >>= 1;
        }

        assert(stash.size() == 0);

        buckets[pos].updateMin();
        pos >>= 1;
        while (pos)
        {
            buckets[pos].updateMin(std::min(buckets[pos * 2].getSubTreeMin(), buckets[pos * 2 + 1].getSubTreeMin()));
            pos >>= 1;
        }
    }

public:
    OHeap(const uint32_t blockNums_, const uint32_t seed)
        : curTimeStamp(0),
        levels(std::__lg(blockNums_) + 1 + ((blockNums_ & -blockNums_) != blockNums_)),
        blockNums(1 << (levels - 1)),
        gen(seed),
        unif(0, blockNums - 1)
    {
        std::cout << levels << ", " << blockNums << std::endl;
        buckets.reserve(2u * blockNums);
        // Insert a dummy head
        buckets.emplace_back(0);
        // Insert a root node
        buckets.emplace_back(OHEAP_ROOT_SIZE);
        for (uint32_t i = 2; i < 2u * blockNums; ++i)
            buckets.emplace_back(OHEAP_NONROOT_SIZE);
    }

    const Tuple<Key>& findMin() const
    {
        return buckets[1].getSubTreeMin();
    }

    void insert(const Key& key)
    {
        uint32_t pos = unif(gen) + blockNums;
        buckets[1].add({ key, pos, curTimeStamp++ });

        pos = unif(gen) / 2 + blockNums;
        evictUpdateMin(pos);
        pos = unif(gen) / 2 + blockNums / 2 + blockNums;
        evictUpdateMin(pos);
    }

    Key extractMin()
    {
        Tuple<Key> tmp = findMin();
        readNRm(tmp.pos, tmp.timeStamp);
        evictUpdateMin(tmp.pos);
        return tmp.key;
    }

    uint32_t count()
    {
        uint32_t cnt = 0;
        for (auto& b : buckets)
            for (const auto& t : b.getTuples())
                if (!t.dummy())
                    cnt++;
        return cnt;
    }

    friend std::ostream& operator<<(std::ostream& os, const OHeap<Key>& oheap)
    {
        uint32_t l = 1;
        for (uint32_t i = 1; i < oheap.buckets.size(); i++)
        {
            l = std::__lg(i);
            os << i << ": ";
            for (uint32_t j = 0; j < l; j++)
                os << "   ";
            os << oheap.buckets[i] << std::endl;
        }
        return os;
    }
};