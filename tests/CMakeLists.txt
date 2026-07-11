add_executable(test_crc test_crc.cpp)
target_link_libraries(test_crc PRIVATE nr_phy)
add_test(NAME CRC COMMAND test_crc)

add_executable(test_modulator test_modulator.cpp)
target_link_libraries(test_modulator PRIVATE nr_phy)
add_test(NAME Modulator COMMAND test_modulator)
