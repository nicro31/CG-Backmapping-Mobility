# Install script for directory: /home/x_nicro/LOE-CTP-FRAG/QC_Tools

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calc_J" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calc_J")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calc_J"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/bin/calc_J")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calc_J" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calc_J")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calc_J"
         OLD_RPATH "/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/ARGUMENTS/PROPERTIES:/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/ARGUMENTS:/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO/FILE_READERS:/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/IO:/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/MATRIX:/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/STRING_SUPPORT:/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/PARAMETERS:/home/x_nicro/LOE-CTP-FRAG/QC_Tools/src/QC_FUNCTIONS:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/software/sse/easybuild/prefix/software/binutils/2.30-GCCcore-7.3.0/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calc_J")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/MATRIX/cmake_install.cmake")
  include("/home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/STRING_SUPPORT/cmake_install.cmake")
  include("/home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/IO/cmake_install.cmake")
  include("/home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/PARAMETERS/cmake_install.cmake")
  include("/home/x_nicro/LOE-CTP-FRAG/QC_Tools/QC_Tools_Build/src/QC_FUNCTIONS/cmake_install.cmake")

endif()

