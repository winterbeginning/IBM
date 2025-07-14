cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(IBMFoam LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG")

if(DEFINED ENV{WM_PROJECT_DIR})
	MESSAGE(STATUS "OpenFOAM: " $ENV{WM_PROJECT_DIR})
else()
    message(FATAL_ERROR "Cannot find OpenFOAM")
endif(DEFINED ENV{WM_PROJECT_DIR})

set(FOAMVERSION $ENV{WM_PROJECT_VERSION}) 
set(FOAMDIR $ENV{WM_PROJECT_DIR})
set(FOAMLIB $ENV{FOAM_LIBBIN})
set(FOAMSRC $ENV{FOAM_SRC})
set(FOAMMPI $ENV{FOAM_MPI})
set(FOAMAPP $ENV{FOAM_APP})
set(FOAM_USER_APPBIN $ENV{FOAM_USER_APPBIN})
set(FOAM_USER_LIBBIN $ENV{FOAM_USER_LIBBIN})
set(WM_PROJECT_USER_DIR $ENV{WM_PROJECT_USER_DIR})

set(FOAMFLAG "-std=c++17 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wno-unused-parameter -Wno-overloaded-virtual -Wold-style-cast
-Wnon-virtual-dtor -Wno-unused-variable -Wno-invalid-offsetof -DNoRepository -ftemplate-depth-100 -fPIC")

# 验证必需的环境变量
foreach(var IN ITEMS FOAM_USER_APPBIN FOAM_USER_LIBBIN FOAMSRC WM_PROJECT_USER_DIR)
    if(NOT ${var})
        message(FATAL_ERROR "Environment variable not set: ${var}")
    endif()
endforeach()

set(LIB_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/Include)

set(SRC_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src/Include)

set(FOAMINC
    ${FOAMSRC}/transportModels/lnInclude
    ${FOAMSRC}/MomentumTransportModels/momentumTransportModels/lnInclude
    ${FOAMSRC}/MomentumTransportModels/incompressible/lnInclude
    ${FOAMSRC}/finiteVolume/lnInclude
    ${FOAMSRC}/dynamicFvMesh/lnInclude
    ${FOAMSRC}/meshTools/lnInclude
    ${FOAMSRC}/sampling/lnInclude
    ${FOAMSRC}/OpenFOAM/lnInclude
    ${FOAMSRC}/OSspecific/POSIX/lnInclude
    ${LIB_INCLUDE_PATH}
    ${SRC_INCLUDE_PATH}
)

# 添加链接目录
link_directories(
    ${FOAM_USER_LIBBIN}
    ${FOAMLIB}
    ${FOAMLIB}/${FOAMMPI}
)

add_subdirectory(lib)

add_subdirectory(src)