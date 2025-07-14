# 添加头文件搜索路径
include_directories(${FOAMINC} ${CMAKE_CURRENT_SOURCE_DIR})

set(LIB_NAME ImmersedBoundaryMethod)

set (IBM_LIB_SRC
    eulerMesh/eulerMesh.C
    immersedBoundaryMethod/immersedBoundaryMethod.C
    IBObject/IBObject.C
    IBObject/IBParticle/IBParticle.C
    IBOStream/IBOStream.C
    IBModel/IBModel.C
    IBModel/directForcing/directForcing.C
    IBMotions/IBMotion/IBMotion.C
    IBMotions/IBMotion/sixDoFMotion/sixDoFMotion.C
    IBMotions/IBMotion/fixedLinearMotion/fixedLinearMotion.C
    IBMotions/IBMotionSolver/IBMotionSolver.C
    IBMotions/IBMotionSolver/newton/newton.C
    IBMotions/IBMotionSolver/modifiedNewton/modifiedNewton.C
)

# 创建可执行文件
add_library(${LIB_NAME} ${IBM_LIB_SRC})

add_definitions(${FOAMFLAG})
# 设置输出目录
set_target_properties(${LIB_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${FOAM_USER_LIBBIN}
)

set (INCLUDE_LIB     
    ${FOAMSRC}/finiteVolume/lnInclude 
    ${FOAMSRC}/meshTools/lnInclude 
    ${FOAMSRC}/dynamicMesh/lnInclude 
    ${FOAMSRC}/dynamicFvMesh/lnInclude 
    ${FOAMAPP}/utilities/mesh/generation/foamyMesh/conformalVoronoiMesh/vectorTools
    ${LIB_INCLUDE_PATH}
)

target_include_directories(${LIB_NAME} PUBLIC ${INCLUDE_LIB})

# 添加链接库
target_link_libraries(${LIB_NAME} PUBLIC
    finiteVolume 
    meshTools
    dynamicMesh
    dynamicFvMesh
    topoChangerFvMesh
)
