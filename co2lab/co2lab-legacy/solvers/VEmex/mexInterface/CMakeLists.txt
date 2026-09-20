PROJECT(MEXINTERFACE)

INCLUDE_DIRECTORIES(../VESimulator/include include)

SET(COMMON_MEX src/extractStructs.cpp
                src/matlab.def
    )

SET(MEX_SRC VETransportCPU)

SET(MEX_OUTPUT_PATH ${CMAKE_BINARY_DIR}/mex)

IF(WIN32)
    SET(PREFIX "")
ELSE(WIN32)
    SET(PREFIX "lib")
ENDIF(WIN32)

IF(CMAKE_SHARED_LIBRARY_SUFFIX)
    SET(SHLIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
ELSE(CMAKE_SHARED_LIBRARY_SUFFIX)
    SET(SHLIB_SUFFIX ${CMAKE_SHLIB_SUFFIX})
ENDIF(CMAKE_SHARED_LIBRARY_SUFFIX)

IF(MATLAB_MEX_EXTENSION)
    SET(MEX_ENDING "${MATLAB_MEX_EXTENSION}")
ELSE(MATLAB_MEX_EXTENSION)
    IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
        IF(WIN32)
            SET(MEX_ENDING "mexw32")
        ELSEIF(APPLE)
            SET(MEX_ENDING "mexmaci")
        ELSE(WIN32)
            SET(MEX_ENDING "mexglx")
        ENDIF(WIN32)
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
        IF(${FORCE_32_BIT})
            IF(WIN32)
                SET(MEX_ENDING "mexw32")
            ELSEIF(APPLE)
                SET(MEX_ENDING "mexmaci")
            ELSE(WIN32)
                SET(MEX_ENDING "mexglx")
            ENDIF(WIN32)
        ELSE(${FORCE_32_BIT})
            IF(WIN32)
                SET(MEX_ENDING "mexw64")
            ELSEIF(APPLE)
                IF("$ENV{MATLAB_ARCH}" MATCHES "maca64" OR
                   CMAKE_OSX_ARCHITECTURES MATCHES "(^|;)arm64($|;)")
                    SET(MEX_ENDING "mexmaca64")
                ELSE()
                    SET(MEX_ENDING "mexmaci64")
                ENDIF()
            ELSE(WIN32)
                SET(MEX_ENDING "mexa64")
            ENDIF(WIN32)
        ENDIF(${FORCE_32_BIT})
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
ENDIF(MATLAB_MEX_EXTENSION)

FOREACH(EXE ${MEX_SRC})
    ADD_LIBRARY(${EXE} SHARED ${COMMON_MEX} src/${EXE}.cpp)
    TARGET_LINK_LIBRARIES(${EXE} VESimulatorCPU ${MATLAB_LIBRARIES})

    IF(TARGET OpenMP::OpenMP_CXX)
        TARGET_LINK_LIBRARIES(${EXE} OpenMP::OpenMP_CXX)
    ELSE(TARGET OpenMP::OpenMP_CXX)
        SET_TARGET_PROPERTIES(${EXE} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
            LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    ENDIF(TARGET OpenMP::OpenMP_CXX)

    IF(${FORCE_32_BIT})
        SET_TARGET_PROPERTIES(${EXE} PROPERTIES COMPILE_FLAGS "-m32"
            LINK_FLAGS "-m32")
    ENDIF(${FORCE_32_BIT})

    #Now move the file
    SET(INNAME "${PREFIX}${EXE}${SHLIB_SUFFIX}")
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
                ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR}/${INNAME}
            )

ENDFOREACH(EXE ${MEX_SRC})
