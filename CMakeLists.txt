cmake_minimum_required(VERSION 3.12)

project(primer3  VERSION      2.05.00
        DESCRIPTION "primer 3"
        HOMEPAGE_URL https://github.com/qPCR4vir/primer3-ori
        LANGUAGES CXX         )

function(Set_Properties targ)
    set_property( TARGET ${targ} PROPERTY CXX_STANDARD 17 )
	target_include_directories(${targ} PUBLIC src)
    target_link_libraries(${targ} PUBLIC "-std=c++17 -Wall -static-libgcc -static-libstdc++  -pthread -lm" )
    install(TARGETS ${test} RUNTIME DESTINATION "${primer3_INSTALL_PATH}bin")
    message("... to build:  ${targ}" )
endfunction(Set_Properties)

add_library   ( dpal   src/dpal.cc src/dpal.hpp )
Set_Properties( dpal )

add_library   ( thal   src/thal.cc src/thal.hpp)
Set_Properties( thal )

add_library   ( oligo_tm src/oligotm.cc)
Set_Properties( oligo_tm )

add_library   ( primer3  src/libprimer3.cc)
Set_Properties( primer3 )

add_library   ( masker   src/masker.cc)
Set_Properties( masker )

add_library   ( p3_seq src/masker.cc)
Set_Properties( p3_seq )

include(${CMAKE_CURRENT_LIST_DIR}/select_filesystem.cmake)

add_executable      ( primer3_core src/primer3_boulder_main.cc
		src/format_output.cc
									src/read_boulder.c
		src/print_boulder.cc
									src/mman.c       )
Set_Properties       ( primer3_core )
target_link_libraries( primer3_core PUBLIC thal dpal  oligo_tm  primer3 masker p3_seq)

add_executable       ( ntdpal  src/ntdpal_main.cc      )
Set_Properties       ( ntdpal )
target_link_libraries( ntdpal PUBLIC dpal)

add_executable       ( ntthal src/thal_main.c    )
Set_Properties       ( ntthal )
target_link_libraries( ntthal PUBLIC thal)

add_executable       ( oligotm src/oligotm_main.cc)
Set_Properties       ( oligotm )
target_link_libraries( oligotm PUBLIC oligo_tm)

add_executable       ( long_seq_tm_test src/long_seq_tm_test_main.c  )
Set_Properties       ( long_seq_tm_test )
target_link_libraries( long_seq_tm_test PUBLIC oligo_tm)

add_executable       ( primer3_masker src/masker_main.cc src/mman.c )
Set_Properties       ( primer3_masker )
target_link_libraries( primer3_masker PUBLIC  thal dpal oligo_tm primer3 masker p3_seq)