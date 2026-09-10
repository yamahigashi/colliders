cmake_minimum_required(VERSION 3.10)

project (yddColliders)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package( Maya REQUIRED )

set(sources
  sources/pluginIdentity.h
  sources/main.cpp 
  sources/bellCollider.cpp
  sources/bellCollider.h
  sources/colliderDrawGeometry.h
  sources/bellColliderSolver.cpp
  sources/bellColliderSolver.h
  sources/planeCollider.cpp
  sources/planeCollider.h
  sources/skirtBellCollider.cpp
  sources/skirtBellCollider.h
  sources/skirtRingFrames.h
  sources/skirtCollideDeformer.cpp
  sources/skirtCollideDeformer.h
  sources/skirtWaveDeformer.cpp
  sources/skirtWaveDeformer.h
  sources/utils.hpp)

add_library(yddColliders SHARED ${sources})

MAYA_PLUGIN( yddColliders )
install(TARGETS yddColliders
        RUNTIME DESTINATION modules
        LIBRARY DESTINATION modules
        ARCHIVE DESTINATION modules)
