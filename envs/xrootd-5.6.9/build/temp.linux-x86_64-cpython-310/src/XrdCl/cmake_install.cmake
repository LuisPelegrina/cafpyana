# Install script for directory: /exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl

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
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310/pyxrootd/libXrdCl.so.3.0.0"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310/pyxrootd/libXrdCl.so.3"
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310/pyxrootd/libXrdCl.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libXrdCl.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libXrdCl.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libXrdCl.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/xrootd/XrdCl" TYPE FILE FILES
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClAnyObject.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClBuffer.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClConstants.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClCopyProcess.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClDefaultEnv.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClEnv.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFile.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFileSystem.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFileSystemUtils.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClMonitor.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClStatus.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClURL.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClXRootDResponses.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClOptional.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPlugInInterface.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPropertyList.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClLog.hh"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/xrootd/private/XrdCl" TYPE FILE FILES
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClJobManager.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClMessage.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPlugInManager.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPostMaster.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClPostMasterInterfaces.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClTransportManager.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClResponseJob.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClSyncQueue.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClZipArchive.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClZipCache.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClOperations.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClOperationHandlers.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClOperationTimeout.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClArg.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClCtx.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFwd.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClParallelOperation.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFileOperations.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFileSystemOperations.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClFinalOperation.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClUtils.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClXRootDTransport.hh"
    "/exp/sbnd/data/users/lpelegri/cafpyana/envs/xrootd-5.6.9/src/XrdCl/XrdClZipOperations.hh"
    )
endif()

