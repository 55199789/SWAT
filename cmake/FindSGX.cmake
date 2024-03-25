# FindPackage cmake file for Intel SGX SDK

cmake_minimum_required(VERSION 3.15) # target_include_directories
include(CMakeParseArguments)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenSSL REQUIRED)

set(SGX_FOUND "NO")
set(SGXSSL_FOUND "NO")

if(EXISTS SGX_DIR)
    set(SGX_PATH ${SGX_DIR})
elseif(EXISTS SGX_ROOT)
    set(SGX_PATH ${SGX_ROOT})
elseif(EXISTS $ENV{SGX_SDK})
    set(SGX_PATH $ENV{SGX_SDK})
elseif(EXISTS $ENV{SGX_DIR})
    set(SGX_PATH $ENV{SGX_DIR})
elseif(EXISTS $ENV{SGX_ROOT})
    set(SGX_PATH $ENV{SGX_ROOT})
else()
    set(SGX_PATH "/opt/intel/sgxsdk")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(SGX_COMMON_CFLAGS -m32)
    set(SGX_LIBRARY_PATH ${SGX_PATH}/lib32)
    set(SGX_ENCLAVE_SIGNER ${SGX_PATH}/bin/x86/sgx_sign)
    set(SGX_EDGER8R ${SGX_PATH}/bin/x86/sgx_edger8r)
else()
    set(SGX_COMMON_CFLAGS -m64)
    set(SGX_LIBRARY_PATH ${SGX_PATH}/lib64)
    set(SGX_ENCLAVE_SIGNER ${SGX_PATH}/bin/x64/sgx_sign)
    set(SGX_EDGER8R ${SGX_PATH}/bin/x64/sgx_edger8r)
endif()

find_path(SGX_INCLUDE_DIR sgx.h "${SGX_PATH}/include" NO_DEFAULT_PATH)
find_path(SGX_LIBRARY_DIR libsgx_urts.so "${SGX_LIBRARY_PATH}" NO_DEFAULT_PATH)

if(SGX_INCLUDE_DIR AND SGX_LIBRARY_DIR)
    set(SGX_FOUND "YES")
    set(SGX_INCLUDE_DIR "${SGX_PATH}/include" CACHE PATH "Intel SGX include directory" FORCE)
    set(SGX_TLIBC_INCLUDE_DIR "${SGX_INCLUDE_DIR}/tlibc" CACHE PATH "Intel SGX tlibc include directory" FORCE)
    set(SGX_LIBCXX_INCLUDE_DIR "${SGX_INCLUDE_DIR}/libcxx" CACHE PATH "Intel SGX libcxx include directory" FORCE)
    set(SGX_INCLUDE_DIRS ${SGX_INCLUDE_DIR} ${SGX_TLIBC_INCLUDE_DIR} ${SGX_LIBCXX_INCLUDE_DIR})
    mark_as_advanced(SGX_INCLUDE_DIR SGX_TLIBC_INCLUDE_DIR SGX_LIBCXX_INCLUDE_DIR SGX_LIBRARY_DIR)
    message(STATUS "Found Intel SGX SDK.")
else()
    message(FATAL_ERROR "NOT found Intel SGX SDK.")
endif()

# if(EXISTS SGXSSL_DIR)
# set(SGXSSL_PATH ${SGXSSL_DIR})
# elseif(EXISTS SGXSSL_ROOT)
# set(SGXSSL_PATH ${SGXSSL_ROOT})
# elseif(EXISTS $ENV{SGXSSL})
# set(SGXSSL_PATH $ENV{SGXSSL})
# elseif(EXISTS $ENV{SGXSSL_DIR})
# set(SGXSSL_PATH $ENV{SGXSSL_DIR})
# elseif(EXISTS $ENV{SGXSSL_ROOT})
# set(SGXSSL_PATH $ENV{SGXSSL_ROOT})
# elseif(EXISTS "/opt/intel/sgxssl")
# set(SGXSSL_PATH "/opt/intel/sgxssl")
# else()
execute_process(COMMAND ${CMAKE_SOURCE_DIR}/prepare_sgxssl.sh)
set(SGXSSL_PATH ${CMAKE_SOURCE_DIR}/sgxssl/Linux/package)
message("NOTE: sgxssl prepared")
# endif()

set(SGXSSL_INCLUDE_PATH ${SGXSSL_PATH}/include)
set(SGXSSL_LIBRARY_PATH ${SGXSSL_PATH}/lib64)
set(SGXSOCKET_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/sgx_socket/include)
find_path(SGXSSL_INCLUDE_DIR tSgxSSL_api.h "${SGXSSL_INCLUDE_PATH}" NO_DEFAULT_PATH)
find_path(SGXSSL_LIBRARY_DIR libsgx_tsgxssl.a "${SGXSSL_LIBRARY_PATH}" NO_DEFAULT_PATH)

if(SGXSSL_INCLUDE_DIR AND SGXSSL_LIBRARY_DIR)
    set(SGXSSL_FOUND "YES")
    message(STATUS "Found Intel SGX SSL.")
else()
    message(STATUS "NOT found Intel SGX SSL.")
endif()

set(SGX_HW ON CACHE BOOL "Run SGX on hardware, OFF for simulation.")
set(SGX_MODE Debug CACHE STRING "SGX build mode: Debug; PreRelease; Release.")

if(SGX_HW)
    set(SGX_URTS_LIB sgx_urts)
    set(SGX_USVC_LIB sgx_uae_service)
    set(SGX_TRTS_LIB sgx_trts)
    set(SGX_TSVC_LIB sgx_tservice)
else()
    set(SGX_URTS_LIB sgx_urts_sim)
    set(SGX_USVC_LIB sgx_uae_service_sim)
    set(SGX_TRTS_LIB sgx_trts_sim)
    set(SGX_TSVC_LIB sgx_tservice_sim)
endif()

if(SGX_MODE STREQUAL "Debug")
    set(SGX_COMMON_CFLAGS "${SGX_COMMON_CFLAGS} -O2 -g -DDEBUG -UNDEBUG -UEDEBUG")
elseif(SGX_MODE STREQUAL "PreRelease")
    set(SGX_COMMON_CFLAGS "${SGX_COMMON_CFLAGS} -O2 -UDEBUG -DNDEBUG -DEDEBUG")
elseif(SGX_MODE STREQUAL "Release")
    set(SGX_COMMON_CFLAGS "${SGX_COMMON_CFLAGS} -O2 -UDEBUG -DNDEBUG -UEDEBUG")
else()
    message(FATAL_ERROR "SGX_MODE ${SGX_MODE} is not Debug, PreRelease or Release.")
endif()

set(ENCLAVE_INC_DIRS "${SGX_INCLUDE_DIR}" "${SGX_TLIBC_INCLUDE_DIR}" "${SGX_LIBCXX_INCLUDE_DIR}")
set(ENCLAVE_C_FLAGS "${SGX_COMMON_CFLAGS} -nostdinc -fvisibility=hidden -fpie -ffunction-sections -fdata-sections -fstack-protector-strong")
set(ENCLAVE_CXX_FLAGS "${ENCLAVE_C_FLAGS} -nostdinc++ -fopenmp")

set(APP_INC_DIRS "${SGX_PATH}/include")
set(APP_C_FLAGS "${SGX_COMMON_CFLAGS} -fPIC -Wno-attributes")
set(APP_CXX_FLAGS "${APP_C_FLAGS}")

function(_build_edl_obj edl edl_search_paths use_prefix)
    get_filename_component(EDL_NAME ${edl} NAME_WE)
    get_filename_component(EDL_ABSPATH ${edl} ABSOLUTE)
    set(EDL_T_C "${CMAKE_CURRENT_BINARY_DIR}/${EDL_NAME}_t.c")
    set(SEARCH_PATHS "")

    foreach(path ${edl_search_paths})
        get_filename_component(ABSPATH ${path} ABSOLUTE)
        list(APPEND SEARCH_PATHS "${ABSPATH}")
    endforeach()

    list(APPEND SEARCH_PATHS "${SGX_PATH}/include")
    string(REPLACE ";" ":" SEARCH_PATHS "${SEARCH_PATHS}")

    if(${use_prefix})
        set(USE_PREFIX "--use-prefix")
    endif()

    add_custom_command(OUTPUT ${EDL_T_C}
        COMMAND ${SGX_EDGER8R} ${USE_PREFIX} --trusted ${EDL_ABSPATH} --search-path ${SEARCH_PATHS}
        MAIN_DEPENDENCY ${EDL_ABSPATH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    add_library(${target}-edlobj OBJECT ${EDL_T_C})
    set_target_properties(${target}-edlobj PROPERTIES COMPILE_FLAGS ${ENCLAVE_C_FLAGS})
    target_include_directories(${target}-edlobj
        PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}
        ${ENCLAVE_INC_DIRS}
        ${SGXSOCKET_INCLUDE_PATH}
    )

    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_CURRENT_BINARY_DIR}/${EDL_NAME}_t.h")
endfunction()

function(add_trusted_library target)
    set(optionArgs USE_PREFIX)
    set(oneValueArgs EDL LDSCRIPT)
    set(multiValueArgs SRCS EDL_SEARCH_PATHS)
    cmake_parse_arguments("SGX" "${optionArgs}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(LDSCRIPT_FLAG "")

    if(NOT "${SGX_LDSCRIPT}" STREQUAL "")
        get_filename_component(LDS_ABSPATH ${SGX_LDSCRIPT} ABSOLUTE)
        set(LDSCRIPT_FLAG "-Wl,--version-script=${LDS_ABSPATH}")
    endif()

    if("${SGX_EDL}" STREQUAL "")
        message("${target}: SGX enclave edl file is not provided; skipping edger8r")
        add_library(${target} STATIC ${SGX_SRCS})
    else()
        if("${SGX_EDL_SEARCH_PATHS}" STREQUAL "")
            message("${target}: SGX enclave edl file search paths are not provided!")
        endif()

        _build_edl_obj(${SGX_EDL} "${SGX_EDL_SEARCH_PATHS}" ${SGX_USE_PREFIX})
        add_library(${target} STATIC ${SGX_SRCS} $<TARGET_OBJECTS:${target}-edlobj>)
    endif()

    set_target_properties(${target} PROPERTIES COMPILE_FLAGS
        "${ENCLAVE_CXX_FLAGS} -include \"tsgxsslio.h\"")
    target_include_directories(${target} PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}
        ${ENCLAVE_INC_DIRS}
        ${SGXSSL_INCLUDE_PATH}
        ${SGXSOCKET_INCLUDE_PATH}
    )

    target_link_libraries(${target}
        "${SGX_COMMON_CFLAGS} \
            -Wl,--no-undefined -nostdlib -nodefaultlibs -nostartfiles \
            -Wl,-z,noexecstack -Wl,-z,relro -Wl,-z,now -pie \
            -L${SGXSSL_LIBRARY_PATH} \
             -Wl,--whole-archive -lsgx_tsgxssl -Wl,--no-whole-archive \
            -lsgx_tsgxssl_ssl -lsgx_tsgxssl_crypto \
            -L${SGX_LIBRARY_PATH} \
            -Wl,--whole-archive -lsgx_tswitchless -l${SGX_TRTS_LIB} -Wl,--no-whole-archive \
            -Wl,--start-group -lsgx_tstdc -lsgx_pthread -lsgx_omp -lsgx_tcxx -lsgx_tcrypto -l${SGX_TSVC_LIB} \
            -lsgx_dcap_tvl -lsgx_ttls -Wl,--end-group \
            -Wl,-Bstatic -Wl,-Bsymbolic -Wl,--no-undefined \
            -Wl,-pie,-eenclave_entry -Wl,--export-dynamic \
             ${LDSCRIPT_FLAG} \
            -Wl,--defsym,__ImageBase=0 -Wl,--gc-sections"
    )
endfunction()

# build enclave shared library
function(add_enclave_library target)
    set(optionArgs USE_PREFIX)
    set(oneValueArgs EDL LDSCRIPT)
    set(multiValueArgs SRCS TRUSTED_LIBS EDL_SEARCH_PATHS)
    cmake_parse_arguments("SGX" "${optionArgs}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if("${SGX_EDL}" STREQUAL "")
        message(FATAL_ERROR "${target}: SGX enclave edl file is not provided!")
    endif()

    if("${SGX_EDL_SEARCH_PATHS}" STREQUAL "")
        message("${target}: SGX enclave edl file search paths are not provided!")
    endif()

    if(NOT "${SGX_LDSCRIPT}" STREQUAL "")
        get_filename_component(LDS_ABSPATH ${SGX_LDSCRIPT} ABSOLUTE)
        set(LDSCRIPT_FLAG "-Wl,--version-script=${LDS_ABSPATH}")
    endif()

    _build_edl_obj(${SGX_EDL} "${SGX_EDL_SEARCH_PATHS}" ${SGX_USE_PREFIX})

    add_library(${target} SHARED ${SGX_SRCS} $<TARGET_OBJECTS:${target}-edlobj>)
    target_compile_definitions(${target}
        PRIVATE
        M_SERVER
    )

    set_target_properties(${target} PROPERTIES COMPILE_FLAGS
        "${ENCLAVE_CXX_FLAGS} -include \"tsgxsslio.h\"")
    target_include_directories(${target}
        PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}
        ${ENCLAVE_INC_DIRS}
        ${SGXSSL_INCLUDE_PATH}
        ${SGXSOCKET_INCLUDE_PATH}
    )

    set(TLIB_LIST "")

    foreach(TLIB ${SGX_TRUSTED_LIBS})
        string(APPEND TLIB_LIST "$<TARGET_FILE:${TLIB}> ")
        add_dependencies(${target} ${TLIB})
    endforeach()

    target_link_libraries(${target}
        "${SGX_COMMON_CFLAGS} \
            -Wl,--no-undefined -nostdlib -nodefaultlibs -nostartfiles \
            -Wl,-z,noexecstack -Wl,-z,relro -Wl,-z,now -pie \
            -L${SGXSSL_LIBRARY_PATH} \
             -Wl,--whole-archive -lsgx_tsgxssl -Wl,--no-whole-archive \
            -lsgx_tsgxssl_ssl -lsgx_tsgxssl_crypto \
            -L${SGX_LIBRARY_PATH} \
            -Wl,--whole-archive -lsgx_tswitchless -l${SGX_TRTS_LIB} -Wl,--no-whole-archive \
            -Wl,--start-group ${TLIB_LIST} -lsgx_tstdc -lsgx_pthread -lsgx_omp -lsgx_tcxx -lsgx_tcrypto -l${SGX_TSVC_LIB} \
            -lsgx_dcap_tvl -lsgx_ttls -Wl,--end-group \
            -Wl,-Bstatic -Wl,-Bsymbolic -Wl,--no-undefined \
            -Wl,-pie,-eenclave_entry -Wl,--export-dynamic \
             ${LDSCRIPT_FLAG} \
            -Wl,--defsym,__ImageBase=0 -Wl,--gc-sections"
    )
endfunction()

# sign the enclave, according to configurations one-step or two-step signing will be performed.
# default one-step signing output enclave name is target.signed.so, change it with OUTPUT option.
function(enclave_sign target)
    set(optionArgs IGNORE_INIT IGNORE_REL)
    set(oneValueArgs KEY CONFIG OUTPUT)
    cmake_parse_arguments("SGX" "${optionArgs}" "${oneValueArgs}" "" ${ARGN})

    if("${SGX_CONFIG}" STREQUAL "")
        message(FATAL_ERROR "${target}: SGX enclave config is not provided! ")
    else()
        get_filename_component(CONFIG_ABSPATH ${SGX_CONFIG} ABSOLUTE)
    endif()

    if("${SGX_KEY}" STREQUAL "")
        if(NOT SGX_HW OR NOT SGX_MODE STREQUAL "Release ")
            message(FATAL_ERROR "${target}: Private key used to sign enclave is not provided! ")
        endif()
    else()
        get_filename_component(KEY_ABSPATH ${SGX_KEY} ABSOLUTE)
    endif()

    if("${SGX_OUTPUT}" STREQUAL "")
        set(OUTPUT_NAME "${target}.signed.so")
    else()
        set(OUTPUT_NAME ${SGX_OUTPUT})
    endif()

    if(${SGX_IGNORE_INIT})
        set(IGN_INIT "-ignore-init-sec-error")
    endif()

    if(${SGX_IGNORE_REL})
        set(IGN_REL "-ignore-rel-error")
    endif()

    if(SGX_HW AND SGX_MODE STREQUAL "Release")
        add_custom_target(${target}-sign ALL
            COMMAND ${SGX_ENCLAVE_SIGNER} gendata
            $<$<NOT:$<STREQUAL:${SGX_CONFIG},>>:-config> $<$<NOT:$<STREQUAL:${SGX_CONFIG},>>:${CONFIG_ABSPATH}>
            -enclave $<TARGET_FILE:${target}> -out $<TARGET_FILE_DIR:${target}>/${target}_hash.hex ${IGN_INIT} ${IGN_REL}
            COMMAND ${CMAKE_COMMAND} -E cmake_echo_color
            --cyan "SGX production enclave first step signing finished, \ use ${CMAKE_CURRENT_BINARY_DIR}/${target}_hash.hex for second step "
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    else()
        add_custom_target(${target}-sign ALL ${SGX_ENCLAVE_SIGNER} sign -key ${KEY_ABSPATH}
            $<$<NOT:$<STREQUAL:${SGX_CONFIG},>>:-config> $<$<NOT:$<STREQUAL:${SGX_CONFIG},>>:${CONFIG_ABSPATH}>
            -enclave $<TARGET_FILE:${target}>
            -out $<TARGET_FILE_DIR:${target}>/${OUTPUT_NAME}
            ${IGN_INIT} ${IGN_REL}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    add_custom_command(TARGET ${target}-sign POST_BUILD
        COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT_NAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Copy enclave.sign.so..."
    )
    set(CLEAN_FILES "$<TARGET_FILE_DIR:${target}>/${OUTPUT_NAME};$<TARGET_FILE_DIR:${target}>/${target}_hash.hex;${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${OUTPUT_NAME}")
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES} ")
endfunction()

function(add_untrusted_executable target)
    set(optionArgs USE_PREFIX)
    set(oneValueArgs ENCLAVE)
    set(multiValueArgs SRCS USER_LINKS EDL EDL_SEARCH_PATHS)
    cmake_parse_arguments("SGX" "${optionArgs}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if("${SGX_EDL}" STREQUAL "")
        message(FATAL_ERROR "${target}: SGX enclave edl file is not provided! ")
    endif()

    if("${SGX_EDL_SEARCH_PATHS}" STREQUAL "")
        message("${target}: SGX enclave edl file search paths are not provided! ")
    endif()

    set(EDL_U_SRCS "")

    foreach(EDL ${SGX_EDL})
        get_filename_component(EDL_NAME ${EDL} NAME_WE)
        get_filename_component(EDL_ABSPATH ${EDL} ABSOLUTE)
        set(EDL_U_C "${CMAKE_CURRENT_BINARY_DIR}/${EDL_NAME}_u.c")
        set(EDL_U_H "${CMAKE_CURRENT_BINARY_DIR}/${EDL_NAME}_u.h")
        set(SEARCH_PATHS "")

        foreach(path ${SGX_EDL_SEARCH_PATHS})
            get_filename_component(ABSPATH ${path} ABSOLUTE)
            list(APPEND SEARCH_PATHS "${ABSPATH}")
        endforeach()

        list(APPEND SEARCH_PATHS "${SGX_PATH}/include")
        string(REPLACE ";" ":" SEARCH_PATHS "${SEARCH_PATHS}")

        if(${SGX_USE_PREFIX})
            set(USE_PREFIX "--use-prefix")
        endif()

        add_custom_command(OUTPUT ${EDL_U_C}
            COMMAND ${SGX_EDGER8R} ${USE_PREFIX} --untrusted ${EDL_ABSPATH} --search-path "${SEARCH_PATHS}"
            MAIN_DEPENDENCY ${EDL_ABSPATH}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

        list(APPEND EDL_U_SRCS ${EDL_U_C})
        list(APPEND EDL_U_HDRS ${EDL_U_H})
    endforeach()

    add_executable(${target} ${SGX_SRCS} ${EDL_U_SRCS})
    set_target_properties(${target} PROPERTIES COMPILE_FLAGS
        ${APP_CXX_FLAGS}
    )
    target_compile_definitions(${target}
        PRIVATE
        ENCLAVE_PATH="${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${SGX_ENCLAVE}.signed.so")

    target_include_directories(${target}
        PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}
        ${APP_INC_DIRS}
        ${OPENSSL_INCLUDE_DIR}
    )

    target_link_libraries(${target}
        Threads::Threads
        OpenMP::OpenMP_CXX
        ${SGX_COMMON_CFLAGS}
        "-L${SGXSSL_LIBRARY_PATH} \
            -Wl,--whole-archive -lsgx_usgxssl -lsgx_uswitchless -Wl,--no-whole-archive \
        -L${SGX_LIBRARY_PATH} -l${SGX_URTS_LIB}"
        sgx_utls
        sgx_dcap_ql
        sgx_dcap_quoteverify
        ${OPENSSL_LIBRARIES}
        ${SGX_USER_LINKS}
    )
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${EDL_U_HDRS})
endfunction()