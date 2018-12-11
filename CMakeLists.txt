cmake_minimum_required(VERSION 3.12)

project(primer3  VERSION      2.05.00
        DESCRIPTION "primer 3"
        HOMEPAGE_URL https://github.com/qPCR4vir/primer3-ori
        LANGUAGES CXX         )



function(install_exe dir files)
    foreach ( test ${${files}})
        add_executable(${test} "${${dir}}${test}.cpp")
        set_property( TARGET ${test} PROPERTY CXX_STANDARD 14 )
        target_link_libraries(${test} ${primer3_lib} ${primer3_LINKS})
        install(TARGETS ${test} RUNTIME DESTINATION "${primer3_INSTALL_PATH}bin")
        message("... to build:  ${${dir}}${test}.cpp" )
    endforeach( test ${${files}})
endfunction(install_exe)


set(PRIMER_EXE     primer3_core)
set(NTDPAL_EXE     ntdpal)
set(NTTHAL_EXE     ntthal)
set(OLIGOTM_EXE    oligotm)
set(LONG_SEQ_EXE   long_seq_tm_test)
set(MASKER_EXE     primer3_masker)


set(LIBOLIGOTM      liboligotm.a)
set(LIBOLIGOTM_DYN  liboligotm.so.1.2.0)
set(LIBDPAL         libdpal.a)
set(LIBDPAL_DYN     libdpal.a.so.1.0.0)
set(LIBTHAL         libthal.a)
set(LIBTHAL_DYN     libthal.a.so.1.0.0)
set(LIBPRIMER3      libprimer3.a)
set(LIBPRIMER3_DYN  libprimer3.so.1.0.0)
set(LIBMASKER       libmasker.a)
set(LIBMASKER_DYN   libmasker.a.so.1.0.0)


#add_library(dpal src/dpal.c src/dpal.h)
add_library(thal)
set_property( TARGET thal PROPERTY CXX_STANDARD 14 )
target_link_libraries(thal PUBLIC "-std=c++14 -Wall -static-libgcc -static-libstdc++  -pthread -lm" )
target_sources(thal PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src/thal.cc ${CMAKE_CURRENT_LIST_DIR}/src/thal.hpp)

include(${CMAKE_CURRENT_LIST_DIR}/select_filesystem.cmake)

set(EXES     ${PRIMER_EXE}
             ${NTDPAL_EXE} 
             ${NTTHAL_EXE}
             ${OLIGOTM_EXE}
             ${LONG_SEQ_EXE}
             ${MASKER_EXE}    )

set (all ${EXES} ${LIBRARIES})

add_executable(${PRIMER_EXE}

        src/primer3_boulder_main.c
        src/format_output.c
        src/read_boulder.c
        src/print_boulder.c

        src/libprimer3.cc
        src/p3_seq_lib.c
        src/dpal.cc
        src/thal.cc
        src/oligotm.c
        src/masker.cc
        src/mman.c
        )

set_property( TARGET ${PRIMER_EXE} PROPERTY CXX_STANDARD 14 )
target_include_directories(${PRIMER_EXE} PUBLIC src)
target_link_libraries(${PRIMER_EXE} "-std=c++14 -Wall -static-libgcc -static-libstdc++  -pthread -lm" thal)
install(TARGETS ${PRIMER_EXE} RUNTIME DESTINATION "C:/Prog/ExtLib/primer3-ori/bin")


add_executable(${NTDPAL_EXE}
        src/ntdpal_main.cc
        src/dpal.cc
        )

set_property( TARGET ${NTDPAL_EXE} PROPERTY CXX_STANDARD 14 )
target_include_directories(${NTDPAL_EXE} PUBLIC src)
target_link_libraries(${NTDPAL_EXE} "-std=c++14 -Wall -static-libgcc -static-libstdc++  -pthread -lm" )
install(TARGETS ${NTDPAL_EXE} RUNTIME DESTINATION "C:/Prog/ExtLib/primer3-ori/bin")


add_executable(${NTTHAL_EXE}
        src/thal_main.c
        src/thal.c
        )

set_property( TARGET ${NTTHAL_EXE} PROPERTY CXX_STANDARD 14 )
target_include_directories(${NTTHAL_EXE} PUBLIC src)
target_link_libraries(${NTTHAL_EXE} "-std=c++14 -Wall -static-libgcc -static-libstdc++  -pthread -lm" )
install(TARGETS ${NTTHAL_EXE} RUNTIME DESTINATION "C:/Prog/ExtLib/primer3-ori/bin")


add_executable(${OLIGOTM_EXE}
        src/oligotm_main.c
        src/oligotm.c
        )

set_property( TARGET ${OLIGOTM_EXE} PROPERTY CXX_STANDARD 14 )
target_include_directories(${OLIGOTM_EXE} PUBLIC src)
target_link_libraries(${OLIGOTM_EXE} "-std=c++14 -Wall -static-libgcc -static-libstdc++  -pthread -lm" )
install(TARGETS ${OLIGOTM_EXE} RUNTIME DESTINATION "C:/Prog/ExtLib/primer3-ori/bin")

add_executable(${LONG_SEQ_EXE}
        src/long_seq_tm_test_main.c
        src/oligotm.c
        )

set_property( TARGET ${LONG_SEQ_EXE} PROPERTY CXX_STANDARD 14 )
target_include_directories(${LONG_SEQ_EXE} PUBLIC src)
target_link_libraries(${LONG_SEQ_EXE} "-std=c++14 -Wall -static-libgcc -static-libstdc++  -pthread -lm" )
install(TARGETS ${LONG_SEQ_EXE} RUNTIME DESTINATION "C:/Prog/ExtLib/primer3-ori/bin")

add_executable(${MASKER_EXE}
        src/masker_main.c

        src/libprimer3.cc
        src/p3_seq_lib.c
        src/dpal.c
        src/thal.c
        src/oligotm.c
        src/masker.cc
        src/mman.c
        )

set_property( TARGET ${MASKER_EXE} PROPERTY CXX_STANDARD 14 )
target_include_directories(${MASKER_EXE} PUBLIC src)
target_link_libraries(${MASKER_EXE} "-std=c++14 -Wall -static-libgcc -static-libstdc++  -pthread -lm -ffloat-store" )
install(TARGETS ${MASKER_EXE} RUNTIME DESTINATION "C:/Prog/ExtLib/primer3-ori/bin")
