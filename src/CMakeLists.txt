# 添加头文件搜索路径
include_directories(${FOAMINC} ${CMAKE_CURRENT_SOURCE_DIR})

# 创建可执行文件
add_executable(${PROJECT_NAME} IBMFoam.C)

add_definitions(${FOAMFLAG})
# 设置输出目录
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${FOAM_USER_APPBIN}
)

set (INCLUDE_SRC     
    ${FOAMSRC}/finiteVolume/lnInclude 
    ${FOAMSRC}/meshTools/lnInclude 
    ${FOAMSRC}/dynamicMesh/lnInclude 
    ${FOAMSRC}/dynamicFvMesh/lnInclude 
    ${LIB_INCLUDE_PATH}
    ${SRC_INCLUDE_PATH}
)

target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_SRC})

# 添加链接库
target_link_libraries(${PROJECT_NAME} PUBLIC
    finiteVolume 
    meshTools
    dynamicMesh
    dynamicFvMesh
    ImmersedBoundaryMethod
    OpenFOAM
    Pstream
)
