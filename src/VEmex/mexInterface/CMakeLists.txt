PROJECT(MEXINTERFACE)

INCLUDE_DIRECTORIES(../VESimulator/include include)

SET(COMMON_MEX src/extractStructs.cpp
                src/matlab.def
    )

SET(MEX_SRC VETransportCPU)

SET(MEX_OUTPUT_PATH ${CMAKE_BINARY_DIR}/mex)

IF(WIN32)
    SET(PREFIX "")
    SET(CMAKE_SHLIB_SUFFIX ".dll")
ELSE(WIN32)
    SET(PREFIX "lib")
    SET(CMAKE_SHLIB_SUFFIX ".so")
ENDIF(WIN32)

IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    IF(WIN32)
        SET(MEX_ENDING "mexw32")
    ELSE(WIN32)
        SET(MEX_ENDING "mexglx")
    ENDIF(WIN32)
ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
    IF(${FORCE_32_BIT})
        IF(WIN32)
            SET(MEX_ENDING "mexw32")
        ELSE(WIN32)
            SET(MEX_ENDING "mexglx")
        ENDIF(WIN32)
    ELSE(${FORCE_32_BIT})
        SET(MEX_ENDING "mexa64")
    ENDIF(${FORCE_32_BIT})
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)

FOREACH(EXE ${MEX_SRC})
    ADD_LIBRARY(${EXE} SHARED ${COMMON_MEX} src/${EXE}.cpp)
    TARGET_LINK_LIBRARIES(${EXE} ${VESimulatorCPU} ${MATLAB_LIBRARIES})

    SET_TARGET_PROPERTIES(${EXE} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
        LINK_FLAGS "${OpenMP_CXX_FLAGS}")

    IF(${FORCE_32_BIT})
        SET_TARGET_PROPERTIES(${EXE} PROPERTIES COMPILE_FLAGS "-m32"
            LINK_FLAGS "-m32")
    ENDIF(${FORCE_32_BIT})

    #Now move the file
    SET(INNAME "${PREFIX}${EXE}${CMAKE_SHLIB_SUFFIX}")
    SET(OUTNAME "${EXE}.${MEX_ENDING}")

    ADD_CUSTOM_COMMAND( TARGET ${EXE}
                        SOURCE ${EXE}
                        COMMAND ${CMAKE_COMMAND}
                        ARGS -E
                        copy
                ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/${INNAME}
                ${MEX_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/${OUTNAME}
                )

    ADD_CUSTOM_COMMAND( TARGET ${EXE}
                        COMMAND ${CMAKE_COMMAND}
                        ARGS -E
                        remove
                ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INDIR}/${INNAME}
            )

ENDFOREACH(EXE ${MEX_SRC})
