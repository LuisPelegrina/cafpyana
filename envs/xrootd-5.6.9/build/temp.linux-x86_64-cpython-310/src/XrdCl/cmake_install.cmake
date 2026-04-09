# Install script for directory: /home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
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

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310/pyxrootd/libXrdCl.so.3.0.0"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310/pyxrootd/libXrdCl.so.3"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libXrdCl.so.3.0.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libXrdCl.so.3"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310/pyxrootd/libXrdCl.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libXrdCl.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libXrdCl.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libXrdCl.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/xrootd/XrdCl" TYPE FILE FILES
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClAnyObject.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClBuffer.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClConstants.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClCopyProcess.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClDefaultEnv.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClEnv.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFile.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFileSystem.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFileSystemUtils.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClMonitor.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClStatus.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClURL.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClXRootDResponses.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClOptional.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPlugInInterface.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPropertyList.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClLog.hh"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/xrootd/private/XrdCl" TYPE FILE FILES
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClJobManager.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClMessage.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPlugInManager.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPostMaster.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPostMasterInterfaces.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClTransportManager.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClResponseJob.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClSyncQueue.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClZipArchive.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClZipCache.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClOperations.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClOperationHandlers.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClOperationTimeout.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClArg.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClCtx.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFwd.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClParallelOperation.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFileOperations.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFileSystemOperations.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFinalOperation.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClUtils.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClXRootDTransport.hh"
    "/home/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClZipOperations.hh"
    )
endif()

