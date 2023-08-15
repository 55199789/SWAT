set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

include_directories(
    ${CMAKE_SOURCE_DIR}/include/
    ${ROCKSDB_INCLUDE_PATH}
)

set(EDL_SEARCH_PATHS
    server/
    ${CMAKE_SOURCE_DIR}/common/
    ${CMAKE_SOURCE_DIR}/sgx_socket/include/
    ${SGXSSL_INCLUDE_PATH}
)

set(T_SRCS_SGX_SOCKET
    ${CMAKE_SOURCE_DIR}/sgx_socket/htons.c
    ${CMAKE_SOURCE_DIR}/sgx_socket/htonl.c
    ${CMAKE_SOURCE_DIR}/sgx_socket/trusted_sgx_socket.cpp
)

set(T_SRCS_SGX_COMMON
    ${CMAKE_SOURCE_DIR}/common/verify_callback.cpp
    ${CMAKE_SOURCE_DIR}/common/utility.cpp
    ${CMAKE_SOURCE_DIR}/common/openssl_utility.cpp
    ${CMAKE_SOURCE_DIR}/common/tcommon.cpp
    ${CMAKE_SOURCE_DIR}/common/err_msg.cpp
)

set(E_SRCS
    enclave/server.cpp
)

add_trusted_library(trusted_sgx_socket_lib
    SRCS ${T_SRCS_SGX_SOCKET}
    EDL enclave/server.edl
    EDL_SEARCH_PATHS ${EDL_SEARCH_PATHS}
)
add_trusted_library(trusted_common_lib
    SRCS ${T_SRCS_SGX_COMMON}
    EDL enclave/server.edl
    EDL_SEARCH_PATHS ${EDL_SEARCH_PATHS}
)
add_enclave_library(enclave
    SRCS ${E_SRCS}
    TRUSTED_LIBS trusted_sgx_socket_lib trusted_common_lib
    EDL enclave/server.edl
    EDL_SEARCH_PATHS ${EDL_SEARCH_PATHS}
    LDSCRIPT enclave/enclave.lds
)
enclave_sign(enclave
    KEY enclave/private.pem
    CONFIG enclave/server.config.xml
)

set(U_SRCS
    host/host.cpp
    ${CMAKE_SOURCE_DIR}/sgx_socket/untrusted_sgx_socket.cpp
    ${CMAKE_SOURCE_DIR}/common/ucommon.cpp
)

set(DB_LINKS
    cpp_redis
    tacopie
    memcached
    ${ROCKSDB_LIBARIES}
)

add_untrusted_executable(serverProxy
    SRCS ${U_SRCS}
    ENCLAVE enclave
    USER_LINKS ${DB_LINKS}
    EDL enclave/server.edl
    EDL_SEARCH_PATHS ${EDL_SEARCH_PATHS}
)
add_dependencies(serverProxy enclave-sign)