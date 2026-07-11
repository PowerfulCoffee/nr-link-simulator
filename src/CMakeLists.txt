set(NR_SOURCES
    common/NrTables.cpp
    common/Utils.cpp
    phy/CrcEncoder.cpp
    phy/LdpcCodec.cpp
    phy/RateMatcher.cpp
    phy/Scrambler.cpp
    phy/Modulator.cpp
    phy/LayerMapper.cpp
    phy/Precoder.cpp
    phy/DmrsGenerator.cpp
    phy/ResourceMapper.cpp
    phy/OfdmModulator.cpp
    phy/ChannelEstimator.cpp
    phy/Equalizer.cpp
    channel/AwgnChannel.cpp
    channel/TdlChannel.cpp
    channel/CdlChannel.cpp
    phy/PdschProcessor.cpp
)

add_library(nr_phy STATIC ${NR_SOURCES})
target_link_libraries(nr_phy PUBLIC armadillo)
target_include_directories(nr_phy PUBLIC ${PROJECT_SOURCE_DIR}/include)

if(NR_USE_MKL)
    find_package(MKL REQUIRED)
    target_link_libraries(nr_phy PUBLIC MKL::MKL)
    target_compile_definitions(nr_phy PUBLIC NR_HAS_MKL)
endif()
