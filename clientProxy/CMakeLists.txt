set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
find_package(OpenSSL REQUIRED)

add_executable(clientProxy
    src/client.cpp

    ${CMAKE_SOURCE_DIR}/common/verify_callback.cpp
    ${CMAKE_SOURCE_DIR}/common/err_msg.cpp
)
target_include_directories(clientProxy
    PUBLIC
    ${SGX_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/include
)
target_compile_definitions(clientProxy
    PRIVATE
    CLIENT_UNTRUSTED
)
target_link_directories(clientProxy
    PUBLIC
    ${SGX_LIBRARY_DIR}
)
target_link_libraries(clientProxy
    Threads::Threads
    OpenMP::OpenMP_CXX
    ${OPENSSL_LIBRARIES}
    ${SGX_URTS_LIB}
    sgx_utls
    sgx_dcap_ql
    sgx_dcap_quoteverify
)