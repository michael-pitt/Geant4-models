#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(CaloR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find HepMC (required package)
#
#include_directories(${HEPMC_INCLUDE_DIR})
#find_package(HepMC REQUIRED)
#----------------------------------------------------------------------------
# Find Pythia6 (optional package)
#
#find_package(Pythia6 QUIET)
#if(PYTHIA6_FOUND)
#  message(STATUS "G4 Examples: Pythia6 found. --> HepMCEx02 example with Pythia6 enabled.") 
#  add_definitions(-DG4LIB_USE_PYTHIA) 
#else()
#  set(PYTHIA6_LIBRARIES "")  
#endif()

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include("${ROOT_USE_FILE}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include_directories(${ROOT_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Find Pythia8 (required package)
#
#find_package(Pythia8 REQUIRED)
#if(Pythia8_FOUND)
#  message(STATUS "G4 Examples: Pythia8 found. --> CaloR example with Pythia8 enabled.") 
#  add_definitions(-DG4LIB_USE_PYTHIA8) 
#else()
#  set(PYTHIA8_LIBRARIES "")
  #set(PYTHIA8_HEPMC_LIBRARY "")  
#endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
#                    ${HEPMC_INCLUDE_DIR}
#                    ${PYTHIA8_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR}
                    ${ROOTSYS}
                   )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Generate dictionary for input/output via HepMC ROOT files
#

set_target_properties(ROOT::Core PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}")


#ROOT_GENERATE_DICTIONARY(HepMC LINKDEF include/LinkDef.h )
#add_library(HepMCdict ${sources} HepMC.cxx)
#target_link_libraries(HepMCdict ${HEPMC_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(CaloR CaloR.cc ${sources} ${headers})
target_link_libraries(CaloR ${Geant4_LIBRARIES})

#target_link_libraries(CaloR ${Geant4_LIBRARIES}
#                      ${HEPMC_LIBRARIES} ${HEPMC_FIO_LIBRARIES}
#                      ${PYTHIA6_LIBRARIES}
#                      ${PYTHIA8_LIBRARIES}
#                      HepMCdict 
#                      ${ROOT_LIBRARIES}
#                      #HepMCdict
#                     )

# if pythia is compiled with g77, link with -lg2c instead.
#target_link_libraries(CaloR ${Geant4_LIBRARIES}
#                      ${HEPMC_LIBRARIES} ${HEPMC_FIO_LIBRARIES}
#                      ${PYTHIA6_LIBRARIES} g2c)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build CaloR. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(CaloR_SCRIPTS
    run_gun.in vis.mac init_vis.mac gui.mac
  )

foreach(_script ${CaloR_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS CaloR DESTINATION bin)

