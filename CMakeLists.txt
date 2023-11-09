cmake_minimum_required (VERSION 3.12.0)

include_guard()

set (VERSION "1.11.0")
project(base-accessors VERSION ${VERSION} LANGUAGES C CXX)

find_package(Git REQUIRED)

include(CheckCXXCompilerFlag)
include (GNUInstallDirs)
# Load the version IDs of all the dependency sub repos
include(dependencies.cmake)

message (STATUS "Forming build label")
execute_process(
   COMMAND git describe --tags --always --dirty
   OUTPUT_VARIABLE BUILD_STR
   OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Build label is ${BUILD_STR}")

message (STATUS "DEV_OVERRIDES = $ENV{DEV_OVERRIDES}")
if(DEFINED ENV{DEV_OVERRIDES})
  message (STATUS "cmake config overrides file specified, $ENV{DEV_OVERRIDES}")
  include($ENV{DEV_OVERRIDES} OPTIONAL)
else()
  message (STATUS "No cmake config overrides file specified, using defaults")
endif()

message(STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}")

# Declare all of sub repos and tools
include(FetchContent)
# Put our askap cmake modules in the cmake path
# We need this before any other Fetch calls so we have the find_package routines available
# -------- askap cmake modules --------
FetchContent_Declare(askap-cmake
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/askap-cmake.git
  GIT_TAG           ${ASKAP_CMAKE_TAG}
  SOURCE_DIR        askap-cmake
)

message (STATUS "Fetching askap-cmake files")
FetchContent_GetProperties(askap-cmake)
if(NOT askap-cmake_POPULATED)
  FetchContent_Populate(askap-cmake)
endif()

list(APPEND CMAKE_MODULE_PATH "${askap-cmake_SOURCE_DIR}")
#message(STATUS "askap-cmake path is ${askap-cmake_SOURCE_DIR}")
message(STATUS "askap-cmake path is ${CMAKE_MODULE_PATH}")

# -------- lofar common --------
FetchContent_Declare(lofar-common
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/lofar-common.git
  GIT_TAG           ${LOFAR_COMMON_TAG}
  SOURCE_DIR        lofar-common
)

# -------- lofar blob --------
FetchContent_Declare(lofar-blob
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/lofar-blob.git
  GIT_TAG           ${LOFAR_BLOB_TAG}
  SOURCE_DIR        lofar-blob
)

# -------- base-askap --------
FetchContent_Declare(base-askap
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/base-askap.git
  GIT_TAG           ${BASE_ASKAP_TAG}
  SOURCE_DIR        base-askap
)

# -------- base-askapparallel --------
FetchContent_Declare(base-askapparallel
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/base-askapparallel.git
  GIT_TAG           ${BASE_ASKAPPARALLEL_TAG}
  SOURCE_DIR        base-askapparallel
)

# -------- base-imagemath --------
FetchContent_Declare(base-imagemath
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/base-imagemath.git
  GIT_TAG           ${BASE_IMAGEMATH_TAG}
  SOURCE_DIR        base-imagemath
)

# -------- base-scimath --------
FetchContent_Declare(base-scimath
  GIT_REPOSITORY    ${ASKAPSDP_GIT_URL}/base-scimath.git
  GIT_TAG           ${BASE_SCIMATH_TAG}
  SOURCE_DIR        base-scimath
)

message (STATUS "Fetching sub repos")
FetchContent_MakeAvailable( lofar-common lofar-blob
                           base-askap base-imagemath base-askapparallel base-scimath
)
message (STATUS "Done - Fetching sub repos")

configure_file(askap_accessors.cc.in ${CMAKE_CURRENT_BINARY_DIR}/askap_accessors.cc)

# Required Version of ASKAP dependencies is the MAJOR and MINOR version
# of this package. This allows the PATCH versions to change on everything
# but the version number cannot be less or more than the current version
#
set(REQUIRED_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}")


option (CXX11 "Compile as C++11 if possible" YES)
option (ENABLE_SHARED "Build shared libraries" YES)
option (ENABLE_RPATH "Include rpath in executables and shared libraries" YES)
option (ENABLE_OPENMP "Build with OPENMP Support" NO)

# find packages
#find_package(lofar-common REQUIRED)
#find_package(lofar-blob REQUIRED)
#find_package(askap-askap 1.2 REQUIRED)
#find_package(askap-scimath 1.2 REQUIRED)
#find_package(askap-parallel 1.1 REQUIRED)
#find_package(askap-askap REQUIRED)
#find_package(askap-scimath REQUIRED)
#find_package(askap-parallel REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem system thread chrono)
find_package(log4cxx REQUIRED)
find_package(Casacore REQUIRED COMPONENTS  ms images mirlib coordinates fits lattices measures scimath scimath_f tables casa)
find_package(XercesC REQUIRED)
find_package(CPPUnit)
find_package(MPI)

# uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

if (ENABLE_OPENMP)
	find_package(OpenMP)
    if (OPENMP_FOUND)
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        else( )
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif( )
    else (OPENMP_FOUND)
        message(WARNING "Cannot fullfill USE_OPENMP, compiler does not support it")
    endif (OPENMP_FOUND)
endif()


if (ENABLE_SHARED)
option (BUILD_SHARED_LIBS "" YES)
    if (ENABLE_RPATH)
        # Set RPATH to use for installed targets; append linker search path
        set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
        set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif (ENABLE_RPATH)
endif(ENABLE_SHARED)

if (CASACORE3 OR CXX11)
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

add_library(accessors SHARED
	${CMAKE_CURRENT_BINARY_DIR}/askap_accessors.cc
)

macro(add_sources_to_accessors)
	foreach(arg IN ITEMS ${ARGN})
		target_sources(accessors PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${arg})
	endforeach()
endmacro()

add_subdirectory(askap/imageaccess)
add_subdirectory(askap/dataaccess)
add_subdirectory(askap/calibaccess)
add_subdirectory(askap/votable)
add_subdirectory(apps)

target_compile_definitions(accessors PUBLIC
	casa=casacore
	HAVE_AIPSPP
	HAVE_BOOST
	HAVE_LOG4CXX
  ADIOS2_USE_MPI
)
# add some more tests and sub-directories

set_target_properties(accessors PROPERTIES
	OUTPUT_NAME askap_accessors
)

target_link_libraries(accessors
	lofar::Blob
	lofar::Common
	askap::parallel
	askap::scimath
	askap::askap
	${CASACORE_LIBRARIES}
	${log4cxx_LIBRARY}
	${XercesC_LIBRARY}
    ${Boost_LIBRARIES}
)

if (OPENMP_FOUND)
    target_link_libraries(accessors OpenMP::OpenMP_CXX)
endif (OPENMP_FOUND)

if (CPPUNIT_FOUND)
    target_link_libraries(accessors
        ${CPPUNIT_LIBRARY})
    target_include_directories(accessors PUBLIC
        ${CPPUNIT_INCLUDE_DIR}
    )
endif ()

target_include_directories(accessors PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/askap>
	$<INSTALL_INTERFACE:include>
	$<INSTALL_INTERFACE:include/askap>
	${Boost_INCLUDE_DIRS}
	${log4cxx_INCLUDE_DIRS}
	${CASACORE_INCLUDE_DIRS}
	${XercesC_INCLUDE_DIR}
)

add_library(askap::accessors ALIAS accessors)

install (
	TARGETS accessors
	EXPORT askap-accessors-targets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)
install (FILES
	askap_accessors.h
	DESTINATION include/askap/
	)

include(CTest)
if (BUILD_TESTING)
  if (NOT CPPUNIT_FOUND)
    message(FATAL_ERROR "CPPUnit is needed for building and running unit tests")
  endif()
  enable_testing()
  add_subdirectory(tests/calibaccess)
  add_subdirectory(tests/dataaccess)
  add_subdirectory(tests/votable)
  add_subdirectory(tests/imageaccess)
endif ()

# look for dataOnly member in Table.h
# (tried and failed to get CheckCXXSymbolExists and CheckCXXSourceCompiles working)
file (READ "${CASACORE_INCLUDE_DIR}/casacore/tables/Tables/Table.h" TMPTXT)
string (FIND "${TMPTXT}" "dataOnly" FOUND)
message(STATUS "dataonly flag: ${FOUND}")
unset (TMPTXT)
if (FOUND GREATER 0)
    message(STATUS "Found casacore::Table::dataOnly: adding HAVE_CASACORE_TABLE_DATAONLY" )
    add_compile_definitions(HAVE_CASACORE_TABLE_DATAONLY)
endif()


#include(yanda_export)
#yanda_export(askap)
