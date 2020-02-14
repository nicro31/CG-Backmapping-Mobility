# Install script for directory: /home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/lib

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so.2.1.0" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so.2.1.0")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so.2.1.0"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/libApproxMVBB.so.2.1.0")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so.2.1.0" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so.2.1.0")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/software/sse/easybuild/prefix/software/binutils/2.30-GCCcore-7.3.0/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so.2.1.0")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/libApproxMVBB.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/software/sse/easybuild/prefix/software/binutils/2.30-GCCcore-7.3.0/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libApproxMVBB.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake/approxmvbb-config.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake/approxmvbb-config.cmake"
         "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/CMakeFiles/Export/share/ApproxMVBB/cmake/approxmvbb-config.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake/approxmvbb-config-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake/approxmvbb-config.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake" TYPE FILE FILES "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/CMakeFiles/Export/share/ApproxMVBB/cmake/approxmvbb-config.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake" TYPE FILE FILES "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/CMakeFiles/Export/share/ApproxMVBB/cmake/approxmvbb-config-noconfig.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake" TYPE FILE FILES "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/cmake/approxmvbb-config-dependencies.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake/modules" TYPE FILE FILES "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/cmake/FindMeta.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB/cmake" TYPE FILE FILES "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/ApproxMVBB_Build/lib/cmake/approxmvbb-config-version.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ApproxMVBB" TYPE FILE PERMISSIONS OWNER_READ FILES
    "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/COPYING"
    "/home/x_nicro/LOE-CTP-FRAG/Approx_MVBB/ApproxMVBB/README.md"
    )
endif()

