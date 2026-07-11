add_executable(pdsch_bler_simulation pdsch_bler_simulation.cpp)
target_link_libraries(pdsch_bler_simulation PRIVATE nr_phy)
target_include_directories(pdsch_bler_simulation PRIVATE ${PROJECT_SOURCE_DIR}/include)

if(OpenMP_CXX_FOUND)
    target_link_libraries(pdsch_bler_simulation PRIVATE OpenMP::OpenMP_CXX)
endif()
