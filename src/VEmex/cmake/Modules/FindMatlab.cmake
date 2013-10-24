# - this module looks for Matlab
# Defines:
#  MATLAB_INCLUDE_DIR: include path for mex.h, engine.h
#  MATLAB_LIBRARIES:   required libraries: libmex, etc
#  MATLAB_MEX_LIBRARY: path to libmex.lib
#  MATLAB_MX_LIBRARY:  path to libmx.lib
#  MATLAB_ENG_LIBRARY: path to libeng.lib

#=============================================================================
# Copyright 2005-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

SET(MATLAB_FOUND 0)
IF(WIN32)
  SET(MATLAB_PATHS
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.10;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.9;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.8;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.7;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.6;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.5;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.4;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.3;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.2;MATLABROOT]"
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.1;MATLABROOT]"
  )
  # Will crash badly if none of the MATLAB versions is present.
  FIND_PATH(MATLAB_DIR "license.txt" ${MATLAB_PATHS} NO_DEFAULT_PATH)

  IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
    SET(MATLAB_ROOT "${MATLAB_DIR}/extern/lib/win32/microsoft/msvc60")
  ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
    IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
      # Assume people are generally using 7.1,
      # if using 7.0 need to link to: ../extern/lib/win32/microsoft/msvc70
      SET(MATLAB_ROOT "${MATLAB_DIR}/extern/lib/win32/microsoft/msvc71")
    ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
      IF(${CMAKE_GENERATOR} MATCHES "Borland")
        # Same here, there are also: bcc50 and bcc51 directories
        SET(MATLAB_ROOT "${MATLAB_DIR}/extern/lib/win32/microsoft/bcc54")
      ELSE(${CMAKE_GENERATOR} MATCHES "Borland")
          IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 10")
              SET(MATLAB_ROOT "${MATLAB_DIR}/extern/lib/win32/microsoft")
          ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 10")
            IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 9")
                SET(MATLAB_ROOT "${MATLAB_DIR}/extern/lib/win32/microsoft")
            ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 9")
                IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 8")
                    SET(MATLAB_ROOT "${MATLAB_DIR}/extern/lib/win32/microsoft")
                ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 8")
                    IF(MATLAB_FIND_REQUIRED)
                        MESSAGE(FATAL_ERROR "Generator not compatible: ${CMAKE_GENERATOR}")
                    ENDIF(MATLAB_FIND_REQUIRED)
                ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 8")
            ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 9")
          ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 10")
      ENDIF(${CMAKE_GENERATOR} MATCHES "Borland")
    ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
  ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
  FIND_LIBRARY(MATLAB_MEX_LIBRARY
    libmex
    ${MATLAB_ROOT}
    )
  FIND_LIBRARY(MATLAB_MX_LIBRARY
    libmx
    ${MATLAB_ROOT}
    )
  FIND_LIBRARY(MATLAB_ENG_LIBRARY
    libeng
    ${MATLAB_ROOT}
    )

  FIND_PATH(MATLAB_INCLUDE_DIR
    "mex.h"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.10;MATLABROOT]/extern/include"
    )
ELSE( WIN32 )
  IF(NOT $ENV{MATLAB_ARCH} EQUAL "")
    # Maybe x86 on a AMD64 machine
    SET(MATLAB_ROOT
      /usr/local/matlab-7sp1/bin/$ENV{MATLAB_ARCH}/
      /opt/matlab-7sp1/bin/$ENV{MATLAB_ARCH}/
      $ENV{HOME}/matlab-7sp1/bin/$ENV{MATLAB_ARCH}/
      $ENV{HOME}/redhat-matlab/bin/$ENV{MATLAB_ARCH}/
      $ENV{MATLAB_PATH}/bin/$ENV{MATLAB_ARCH}/
      )
  ELSEIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # Regular x86
    SET(MATLAB_ROOT
      /usr/local/matlab-7sp1/bin/glnx86/
      /opt/matlab-7sp1/bin/glnx86/
      $ENV{HOME}/matlab-7sp1/bin/glnx86/
      $ENV{HOME}/redhat-matlab/bin/glnx86/
      $ENV{MATLAB_PATH}/bin/glnx86/
      )
  ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # AMD64:
    SET(MATLAB_ROOT
      /usr/local/matlab-7sp1/bin/glnxa64/
      /opt/matlab-7sp1/bin/glnxa64/
      $ENV{HOME}/matlab7_64/bin/glnxa64/
      $ENV{HOME}/matlab-7sp1/bin/glnxa64/
      $ENV{HOME}/redhat-matlab/bin/glnxa64/
      $ENV{MATLAB_PATH}/bin/glnxa64/
      )
  ENDIF(NOT $ENV{MATLAB_ARCH} EQUAL "")
  FIND_LIBRARY(MATLAB_MEX_LIBRARY
    mex
    ${MATLAB_ROOT}
    )
  FIND_LIBRARY(MATLAB_MX_LIBRARY
    mx
    ${MATLAB_ROOT}
    )
  FIND_LIBRARY(MATLAB_ENG_LIBRARY
    eng
    ${MATLAB_ROOT}
    )
  FIND_PATH(MATLAB_INCLUDE_DIR
    "mex.h"
    "/usr/local/matlab-7sp1/extern/include/"
    "/opt/matlab-7sp1/extern/include/"
    "$ENV{HOME}/matlab-7sp1/extern/include/"
    "$ENV{HOME}/redhat-matlab/extern/include/"
    "$ENV{MATLAB_PATH}/extern/include"
    )
ENDIF(WIN32)

# This is common to UNIX and Win32:
SET(MATLAB_LIBRARIES
  ${MATLAB_MEX_LIBRARY}
  ${MATLAB_MX_LIBRARY}
  ${MATLAB_ENG_LIBRARY}
)

IF(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARIES)
  SET(MATLAB_FOUND 1)
ENDIF(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARIES)

MARK_AS_ADVANCED(
  MATLAB_LIBRARIES
  MATLAB_MEX_LIBRARY
  MATLAB_MX_LIBRARY
  MATLAB_ENG_LIBRARY
  MATLAB_INCLUDE_DIR
  MATLAB_FOUND
  MATLAB_ROOT
)

