#pragma once 
#include <assert.h>
#include <memory>
#include <random>

template <typename Key>
class Query
{
public:
    Query(const Key& left_, const Key& right_): left(left_), right(right_) {}
    virtual OP getOp() = 0;
    std::pair<Key, Key> getQueryData() const
    {
        return std::make_pair(left, right);
    }
    virtual ~Query() {}

protected:
    const Key left;
    const Key right;
};

template <typename Key>
class SearchQuery: public Query<Key>
{
public:
    SearchQuery(const Key& left_, const Key& right_): Query<Key>(left_, right_) {}
    virtual OP getOp() final override
    {
        return OP::SEARCH_OP;
    }
    virtual ~SearchQuery() {}
};

template <typename Key>
class InsertQuery: public Query<Key>
{
public:
    InsertQuery(const Key& left_, const Key& right_): Query<Key>(left_, right_)
    {
        assert(left_ == right_);
    }

    virtual OP getOp() final override
    {
        return OP::INSERT_OP;
    }
    virtual ~InsertQuery() {}
};

template <typename Key>
class DeleteQuery: public Query<Key>
{
public:
    DeleteQuery(const Key& left_, const Key& right_): Query<Key>(left_, right_)
    {
        assert(left_ == right_);
    }

    virtual OP getOp() final override
    {
        return OP::DELETE_OP;
    }
    virtual ~DeleteQuery() {}
};
