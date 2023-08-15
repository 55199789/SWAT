cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(SWAT VERSION 1.0
    DESCRIPTION "Secure DB for range queries with differentially oblivious updates and smoothed search frequency"
    LANGUAGES C CXX)

if(NOT EXISTS ${CMAKE_SOURCE_DIR}/rocksdb/librocksdb.a)
    execute_process(COMMAND make static_lib -j
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/rocksdb
    )
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set(ROCKSDB_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/rocksdb/include)
set(ROCKSDB_LIBARIES
    ${CMAKE_SOURCE_DIR}/rocksdb/librocksdb.a
    snappy
    bz2
    dl
    z
    zstd
    lz4
)
OPTION(OBLIVIOUSMERGE "Option description" ON)

IF(OBLIVIOUSMERGE)
    ADD_DEFINITIONS(-DOBLIVIOUSMERGE)
ENDIF(OBLIVIOUSMERGE)

IF(NOT DEFINED RECORD_LENGTH_IN_BYTE)
    SET(RECORD_LENGTH_IN_BYTE 4096)
ENDIF()

message(RECORD_LENGTH_IN_BYTE=${RECORD_LENGTH_IN_BYTE})

ADD_DEFINITIONS(-DRECORD_LENGTH_IN_BYTE=${RECORD_LENGTH_IN_BYTE})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(SGX REQUIRED)

add_subdirectory(clientProxy)
add_subdirectory(serverProxy)