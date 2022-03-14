cmake_minimum_required(VERSION 3.5)

# project name and the language is c++
project(ECS_Client LANGUAGES CXX)

# c++ standard
set(CMAKE_CXX_STANDARD 17)
# c++ enforce c++ standard
set(CMAKE_CXX_STANDARD_REQUIRED True)

# header files root path
set(Header_path include parser/include)
# executable name
set(EXECUTEABLE_NAME feature_detect)

# source files
set(SRC_list parser/src/arguparser.cpp src/main.cpp)

# determine the number of processes to perform the link time optimization
include(ProcessorCount)
ProcessorCount(N)
message(ParallelCount: ${N})
if(NOT N EQUAL 0)
  set(LINK_OPT_PARALLEL -flto=${N})
else()
  set(LINK_OPT_PARALLEL -flto)
endif()

# common flags used shared by both CXXFLAG and LDFLAG
set(FLAGSCOMMON
    -Wall # all warnings about constructions that some consider questionable
    -Wextra # turn on warnings that is not turned on by -Wall
    -pedantic # Issue all the warnings demanded by strict ISO C and ISO C++
    -Wno-deprecated-declarations
    # -Wundef
    -Wformat=2
    -Wformat-truncation
    # -Wconversion -Wold-style-cast
    -O3 # optimization level 3
    ${LINK_OPT_PARALLEL} # turn on link time optimization in parallel
    -fdevirtualize-at-ltrans # aggressive link-time devirtualization
    -fipa-pta # interprocedural pointer, modification and reference analysis
    # Use pipes rather than temporary files for communication between the
    # various stages of compilation to speed up compiling
    -pipe
    # -fstack-protector-strong # emit extra code to check for buffer overflows
    # https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html
)

# define a variable to hold pre-compile definitions
set(COMPILE_DEF_COMMON PCL_NO_PRECOMPILE # For custom PCL point type to work
                                         # properly
)

# Add Verbose debug macros
if(DEFINED VERBOSE_DEBUG)

endif()

# define a variable to hold compile options
set(CXXFLAGS ${FLAGSCOMMON})

# define a variable to hold link options
set(LDFLAGS ${FLAGSCOMMON})

# define a variable to hold the pre-compile definitions for debugging
set(DBGPRINTDEF __ECS_DEBUG__)

# define a variable to hold the flag for compiling the binary with GDB symbols
set(GDBFLAG -g3)

# define a variable to hold the debug flag for compiling with sanitizer
set(MEMCHECKFLAG
    # maximum number of hash table slots to use during variable tracking
    # dataflow analysis of any function
    --param=max-vartrack-size=60000000
    -fsanitize=address,undefined,bounds-strict,pointer-compare,pointer-subtract,leak,alignment
)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
  # don't ignore CMAKE_REQUIRED_LIBRARIES see cmake --help-policy CMP0075 for
  # more details
  cmake_policy(SET CMP0075 NEW)
endif()

# find Qt5 Widgets
find_package(Qt5 REQUIRED COMPONENTS Widgets Core Gui)

# find PCL
find_package(PCL REQUIRED)

# find libboost
find_package(
  Boost REQUIRED
  COMPONENTS system # basic boost functionality
             filesystem # required by PCL
             date_time # required by PCL
             iostreams # required by PCL
             serialization # required by PCL
             regex # required by PCL
             program_options # for program command line options parser
)

# find VTK find_package(VTK REQUIRED)

# find libeigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# fine libopencv find_package(OpenCV REQUIRED)

# find pthread
set(CMAKE_THREAD_PREFER_PTHREAD True) # prefer pthread
set(THREADS_PREFER_PTHREAD_FLAG True)
find_package(Threads REQUIRED)
if(NOT ${CMAKE_USE_PTHREADS_INIT}) # terminate if pthread not found
  message(FATAL_ERROR "pthread is not found")
endif()

# find tcmalloc (optional)
unset(FOUND_TCMALLOC CACHE)
find_library(FOUND_TCMALLOC tcmalloc)
if(FOUND_TCMALLOC)
  message("libtcmalloc found (lib: ${FOUND_TCMALLOC})")
  # use tcmalloc instead of built-in malloc
  set(ADDITIONAL_LIBS ${ADDITIONAL_LIBS} tcmalloc)
  # disable built-in malloc
  set(ADDITIONAL_FLAGS -fno-builtin-malloc -fno-builtin-calloc
                       -fno-builtin-realloc -fno-builtin-free)
else()
  message("libtcmalloc not found, use default malloc")
endif()

# flag to automatically generate Qt moc files
set(CMAKE_AUTOMOC ON)
# flag to automatically generate header files from Qt ui file
set(CMAKE_AUTOUIC ON)

include(${VTK_USE_FILE})

include_directories(. ${Header_path} ${Boost_INCLUDE_DIRS} ${PCL_INCLUDE_DIRS}
                    ${VTK_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS})

link_directories(${Boost_LIBRARY_DIRS} ${PCL_LIBRARY_DIRS})

# define the targes for building executables
add_executable(${EXECUTEABLE_NAME} # normal executable
               ${SRC_list} parser/include/line_parser.h parser/include/point_type.h)

if(${MULTI_TARGETS})
  add_executable(# executable with more prints
                 ${EXECUTEABLE_NAME}_print ${SRC_list})
  add_executable(# executable with GDB symbols
                 ${EXECUTEABLE_NAME}_gdb ${SRC_list})
  add_executable(# executable with debug sanitizer
                 ${EXECUTEABLE_NAME}_memcheck ${SRC_list})

  # add shortcut targets for executables
  add_custom_target(print)
  add_custom_target(gdb)
  add_custom_target(memcheck)
  # connect shortcut to the origin executable targets
  add_dependencies(print ${EXECUTEABLE_NAME}_print)
  add_dependencies(gdb ${EXECUTEABLE_NAME}_gdb)
  add_dependencies(memcheck ${EXECUTEABLE_NAME}_memcheck)
  # do not build those executables by default
  set_target_properties(${EXECUTEABLE_NAME}_print PROPERTIES EXCLUDE_FROM_ALL
                                                             TRUE)
  set_target_properties(${EXECUTEABLE_NAME}_gdb PROPERTIES EXCLUDE_FROM_ALL
                                                           TRUE)
  set_target_properties(${EXECUTEABLE_NAME}_memcheck PROPERTIES EXCLUDE_FROM_ALL
                                                                TRUE)
endif()

# target-specific pre-compile definitions
set(COMPILE_DEF_SHARE ${Qt5Widgets_COMPILE_DEFINITIONS} ${VTK_DEFINITIONS}
                      ${PCL_DEFINITIONS} ${COMPILE_DEF_COMMON})

target_compile_definitions(${EXECUTEABLE_NAME} PRIVATE ${COMPILE_DEF_SHARE})
if(${MULTI_TARGETS})
  target_compile_definitions(${EXECUTEABLE_NAME}_print
                             PRIVATE ${COMPILE_DEF_SHARE} ${DBGPRINTDEF})
  target_compile_definitions(${EXECUTEABLE_NAME}_gdb
                             PRIVATE ${COMPILE_DEF_SHARE})
  target_compile_definitions(${EXECUTEABLE_NAME}_memcheck
                             PRIVATE ${COMPILE_DEF_SHARE})
endif()

# target-specific link-time libraries
set(TARGET_LINK_LIBS
    modbus
    Threads::Threads
    ${Boost_LIBRARIES}
    ${PCL_LIBRARIES}
    Eigen3::Eigen
    ${VTK_LIBRARIES}
    ${Qt5Widgets_LIBRARIES}
    ${OpenCV_LIBS}
    ${ADDITIONAL_LIBS})

target_link_libraries(${EXECUTEABLE_NAME} PRIVATE ${TARGET_LINK_LIBS})
if(${MULTI_TARGETS})
  target_link_libraries(${EXECUTEABLE_NAME}_print PRIVATE ${TARGET_LINK_LIBS})
  target_link_libraries(${EXECUTEABLE_NAME}_gdb PRIVATE ${TARGET_LINK_LIBS})
  target_link_libraries(${EXECUTEABLE_NAME}_memcheck
                        PRIVATE ${TARGET_LINK_LIBS})
endif()

# target-specific compile options
target_compile_options(${EXECUTEABLE_NAME} PRIVATE ${CXXFLAGS}
                                                   ${ADDITIONAL_FLAGS})
if(${MULTI_TARGETS})
  target_compile_options(${EXECUTEABLE_NAME}_print PRIVATE ${CXXFLAGS}
                                                           ${ADDITIONAL_FLAGS})
  target_compile_options(${EXECUTEABLE_NAME}_gdb PRIVATE ${CXXFLAGS} ${GDBFLAG}
                                                         ${ADDITIONAL_FLAGS})
  target_compile_options(
    ${EXECUTEABLE_NAME}_memcheck PRIVATE ${CXXFLAGS} ${MEMCHECKFLAG}
                                         ${ADDITIONAL_FLAGS})
endif()

# target-specific link-time options
set(LDFLAGS ${LDFLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}
            ${ADDITIONAL_FLAGS})
if(CMAKE_VERSION VERSION_LESS 3.13)
  # cmake < 3.13 doesn't have command target_link_options(), so we use an
  # alternative way to assign linker options cmake list --> string replace ";"
  # with space
  string(REGEX REPLACE ";" " " LDFLAGS "${LDFLAGS}")
  # set link time flag
  set_target_properties(${EXECUTEABLE_NAME} PROPERTIES LINK_FLAGS "${LDFLAGS}")
  if(${MULTI_TARGETS})
    set(GDBLDFLAGS ${LDFLAGS} ${GDBFLAG})
    set(MEMCHECKLDFLAGS ${LDFLAGS} ${MEMCHECKFLAG})
    # replace ";" with space
    string(REGEX REPLACE ";" " " GDBLDFLAGS "${GDBLDFLAGS}")
    # replace ";" with space
    string(REGEX REPLACE ";" " " MEMCHECKLDFLAGS "${MEMCHECKLDFLAGS}")
    set_target_properties(${EXECUTEABLE_NAME}_print PROPERTIES LINK_FLAGS
                                                               "${LDFLAGS}")
    set_target_properties(${EXECUTEABLE_NAME}_gdb PROPERTIES LINK_FLAGS
                                                             "${GDBLDFLAGS}")
    set_target_properties(${EXECUTEABLE_NAME}_memcheck
                          PROPERTIES LINK_FLAGS "${MEMCHECKLDFLAGS}")
  endif()
  if(CMAKE_VERSION VERSION_LESS 3.7.0) # For Qt uic to work properly with cmake
                                       # < 3.7.0
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
  endif()
else()
  target_link_options(${EXECUTEABLE_NAME} PRIVATE ${LDFLAGS})
  if(${MULTI_TARGETS})
    target_link_options(${EXECUTEABLE_NAME}_print PRIVATE ${LDFLAGS})
    target_link_options(${EXECUTEABLE_NAME}_gdb PRIVATE ${LDFLAGS} ${GDBFLAG})
    target_link_options(${EXECUTEABLE_NAME}_memcheck PRIVATE ${LDFLAGS}
                        ${MEMCHECKFLAG})
  endif()
endif()

# [OPTIONAL] The following area is only to display debug info
if(${CMAKELIST_DEBUG_MODE})
  execute_process(
    COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST
  )# Get all propreties that cmake supports

  # Convert command output into a CMake list
  string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
  string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

  # function to print all support property of cmake
  function(print_properties)
    message("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
  endfunction(print_properties)

  # function to print all properties of target ${tgt}
  function(print_target_properties tgt)
    if(NOT TARGET ${tgt})
      message("There is no target named '${tgt}'")
      return()
    endif()

    foreach(prop ${CMAKE_PROPERTY_LIST})
      string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
      # Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-
      # the- location-property-may-not-be-read-from-target-error-i
      if(prop STREQUAL "LOCATION"
         OR prop MATCHES "^LOCATION_"
         OR prop MATCHES "_LOCATION$")
        continue()
      endif()
      # message ("Checking ${prop}")
      get_property(
        propval
        TARGET ${tgt}
        PROPERTY ${prop}
        SET)
      if(propval)
        get_target_property(propval ${tgt} ${prop})
        message("${tgt} ${prop} = ${propval}")
      endif()
    endforeach(prop)
  endfunction(print_target_properties)

  # display found lib variables
  message(Qt5Widgets_INCLUDE_DIRS = ${Qt5Widgets_INCLUDE_DIRS})
  message(Qt5Widgets_DEFINITIONS = ${Qt5Widgets_DEFINITIONS})
  message(Qt5Widgets_COMPILE_DEFINITIONS = ${Qt5Widgets_COMPILE_DEFINITIONS})
  message(Qt5Widgets_LIBRARIES = ${Qt5Widgets_LIBRARIES})
  message(Qt5Widgets_EXECUTABLE_COMPILE_FLAGS =
          ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})
  message(VTK_USE_FILE = ${VTK_USE_FILE})
  message(VTK_DEFINITIONS = ${VTK_DEFINITIONS})
  message(VTK_LIBRARIES = ${VTK_LIBRARIES})
  message(VTK_INCLUDE_DIRS = ${VTK_INCLUDE_DIRS})
  message(VTK_INSTALL_PREFIX = ${VTK_INSTALL_PREFIX})
  message(VTK_LIBRARY_DIRS = ${VTK_LIBRARY_DIRS})
  message(Boost_LIBRARIES = ${Boost_LIBRARIES})
  message(PCL_DEFINITIONS = ${PCL_DEFINITIONS})
  message(PCL_INCLUDE_DIRS = ${PCL_INCLUDE_DIRS})
  message(PCL_LIBRARIES = ${PCL_LIBRARIES})

  print_target_properties(Qt5::Widgets)
  print_target_properties(Qt5::Core)
  print_target_properties(Qt5::Gui)
endif()

# set(CMAKE_VERBOSE_MAKEFILE ON)
