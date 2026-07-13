add_executable(pdsch_bler_simulation pdsch_bler_simulation.cpp)
target_link_libraries(pdsch_bler_simulation PRIVATE nr_phy)
target_include_directories(pdsch_bler_simulation PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(test_bler test_bler.cpp)
target_link_libraries(test_bler PRIVATE nr_phy)
target_include_directories(test_bler PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(pdsch_bler_qpsk pdsch_bler_qpsk.cpp)
target_link_libraries(pdsch_bler_qpsk PRIVATE nr_phy)
target_include_directories(pdsch_bler_qpsk PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(pdsch_bler_qam16 pdsch_bler_qam16.cpp)
target_link_libraries(pdsch_bler_qam16 PRIVATE nr_phy)
target_include_directories(pdsch_bler_qam16 PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(quick_test quick_test.cpp)
target_link_libraries(quick_test PRIVATE nr_phy)
target_include_directories(quick_test PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(debug_params debug_params.cpp)
target_link_libraries(debug_params PRIVATE nr_phy)
target_include_directories(debug_params PRIVATE ${PROJECT_SOURCE_DIR}/include)

if(OpenMP_CXX_FOUND)
    target_link_libraries(pdsch_bler_simulation PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(test_bler PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_qpsk PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_qam16 PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(quick_test PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(debug_params PRIVATE OpenMP::OpenMP_CXX)
endif()
