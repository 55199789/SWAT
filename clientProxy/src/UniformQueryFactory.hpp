#pragma once 
#include "Query.hpp"

template <typename Key>
class UniformQueryFactory
{
private:
    std::mt19937 gen;
    Key max;
    using uniform_distribution =
        typename std::conditional<
        std::is_floating_point<Key>::value,
        std::uniform_real_distribution<Key>,
        typename std::conditional<
        std::is_integral<Key>::value,
        std::uniform_int_distribution<Key>,
        void>::type>::type;

    bool rangeQuery;
    Key smallRangeLen, largeRangeLen, rangeLen;
    uniform_distribution pointDist;
    std::discrete_distribution<uint32_t> typeDist;
    UniformQueryFactory() {}    // disable default constructor

    std::unique_ptr<SearchQuery<Key>> createPointRead()
    {
        Key x;
        x = pointDist(gen);
        return std::unique_ptr<SearchQuery<Key>>(new SearchQuery<Key>(x, x));
    }

    std::unique_ptr<SearchQuery<Key>> createSmallRead()
    {
        Key x, y;
        for (;;)
        {
            x = pointDist(gen);
            if (x <= max - smallRangeLen)
                break;
        }
        y = x + smallRangeLen;
        return std::unique_ptr<SearchQuery<Key>>(new SearchQuery<Key>(x, y));
    }

    std::unique_ptr<SearchQuery<Key>> createLargeRead()
    {
        Key x, y;
        for (;;)
        {
            x = pointDist(gen);
            if (x <= max - largeRangeLen)
                break;
        }
        y = x + largeRangeLen;
        return std::unique_ptr<SearchQuery<Key>>(new SearchQuery<Key>(x, y));
    }

    std::unique_ptr<SearchQuery<Key>> createRangeQuery()
    {
        Key x, y;
        for (;;)
        {
            x = pointDist(gen);
            if (x <= max - rangeLen)
                break;
        }
        y = x + rangeLen;
        return std::unique_ptr<SearchQuery<Key>>(new SearchQuery<Key>(x, y));
    }

    std::unique_ptr<InsertQuery<Key>> createInsertQuery()
    {
        Key x = pointDist(gen);
        return std::unique_ptr<InsertQuery<Key>>(new InsertQuery<Key>(x, x));
    }

    std::unique_ptr<DeleteQuery<Key>> createDeleteQuery()
    {
        Key x = pointDist(gen);
        return std::unique_ptr<DeleteQuery<Key>>(new DeleteQuery<Key>(x, x));
    }

public:
    UniformQueryFactory(const UniformQueryFactory<Key>&) = delete;
    UniformQueryFactory<Key>& operator=(const UniformQueryFactory<Key>&) = delete;
    UniformQueryFactory(UniformQueryFactory<Key>&&) = delete;
    UniformQueryFactory<Key>& operator=(UniformQueryFactory<Key>&&) = delete;

    static UniformQueryFactory<Key>& getInstance()
    {
        static UniformQueryFactory<Key> instance;
        return instance;
    }

    void setup(uint32_t seed,
        Key min = std::numeric_limits<Key>::min(), Key max = std::numeric_limits<Key>::max(),
        Key smallRangeLen = 1,
        Key largeRangeLen = std::numeric_limits<Key>::max(),
        const std::initializer_list<double>& weights = { 5.0, 5.0, 90.0, 0.0, 0.0 })
        // pointReads, smallReads, largeReads, insert, delete
    {
        rangeQuery = false;
        this->max = max;
        gen.seed(seed);
        this->smallRangeLen = smallRangeLen;
        this->largeRangeLen = largeRangeLen;
        pointDist = uniform_distribution(min, max);
        typeDist = std::discrete_distribution<uint32_t>(weights);
    }

    void setup(uint32_t seed,
        Key min = std::numeric_limits<Key>::min(), Key max = std::numeric_limits<Key>::max(),
        Key rangeLen = std::numeric_limits<Key>::max(),
        const std::initializer_list<double>& weights = { 100.0, 0.0, 0.0 })
        // range queries, insert, delete
    {
        rangeQuery = true;
        this->max = max;
        gen.seed(seed);
        this->rangeLen = rangeLen;
        pointDist = uniform_distribution(min, max);
        typeDist = std::discrete_distribution<uint32_t>(weights);
    }

    std::unique_ptr<Query<Key>> createQuery()
    {
        uint32_t type = typeDist(gen);
        if (rangeQuery)
        {
            switch (type)
            {
            case 0:
                return createRangeQuery();
                break;
            case 1:
                return createInsertQuery();
                break;
            case 2:
                return createDeleteQuery();
                break;
            default:
                throw std::runtime_error("Invalid query type");
                break;
            }
        }
        else
        {
            switch (type)
            {
            case 0:
                return createPointRead();
                break;
            case 1:
                return createSmallRead();
                break;
            case 2:
                return createLargeRead();
                break;
            case 3:
                return createInsertQuery();
                break;
            case 4:
                return createDeleteQuery();
                break;
            default:
                throw std::runtime_error("Invalid query type");
                break;
            }
        }
    }
    virtual ~UniformQueryFactory() {}

};