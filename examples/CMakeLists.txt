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

add_executable(pdsch_bler_mcs27 pdsch_bler_mcs27.cpp)
target_link_libraries(pdsch_bler_mcs27 PRIVATE nr_phy)
target_include_directories(pdsch_bler_mcs27 PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(test_chain test_chain.cpp)
target_link_libraries(test_chain PRIVATE nr_phy)
target_include_directories(test_chain PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(pdsch_bler_mcs27_fine pdsch_bler_mcs27_fine.cpp)
target_link_libraries(pdsch_bler_mcs27_fine PRIVATE nr_phy)
target_include_directories(pdsch_bler_mcs27_fine PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(test_tdl test_tdl.cpp)
target_link_libraries(test_tdl PRIVATE nr_phy)
target_include_directories(test_tdl PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(pdsch_bler_mcs27_tdla pdsch_bler_mcs27_tdla.cpp)
target_link_libraries(pdsch_bler_mcs27_tdla PRIVATE nr_phy)
target_include_directories(pdsch_bler_mcs27_tdla PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(pdsch_bler_all_mcs pdsch_bler_all_mcs.cpp)
target_link_libraries(pdsch_bler_all_mcs PRIVATE nr_phy)
target_include_directories(pdsch_bler_all_mcs PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(print_mcs_params print_mcs_params.cpp)
target_link_libraries(print_mcs_params PRIVATE nr_phy)
target_include_directories(print_mcs_params PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(pdsch_bler_cbs_test pdsch_bler_cbs_test.cpp)
target_link_libraries(pdsch_bler_cbs_test PRIVATE nr_phy)
target_include_directories(pdsch_bler_cbs_test PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(test_features test_features.cpp)
target_link_libraries(test_features PRIVATE nr_phy)
target_include_directories(test_features PRIVATE ${PROJECT_SOURCE_DIR}/include)

if(OpenMP_CXX_FOUND)
    target_link_libraries(pdsch_bler_simulation PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(test_bler PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_qpsk PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_qam16 PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(quick_test PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(debug_params PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_mcs27 PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(test_chain PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_mcs27_fine PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(test_tdl PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_mcs27_tdla PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_all_mcs PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(print_mcs_params PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(pdsch_bler_cbs_test PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(test_features PRIVATE OpenMP::OpenMP_CXX)
endif()
