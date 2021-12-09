cmake_minimum_required(VERSION 3.13)

# Aurora Examples Project
project(Aurora_Examples)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++14")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

find_library(LIBSNDFILE_LIBRARY NAMES sndfile libsndfile-1 libsndfile)
find_path(LIBSNDFILE_INCLUDE_DIRECTORY sndfile.h)

message(STATUS "LIBSNDFILE FOUND: ${LIBSNDFILE_LIBRARY}
${LIBSNDFILE_INCLUDE_DIRECTORY}")

# function to build Aurora examples
function(add_program NAME)
  add_executable(${NAME} examples/${NAME}.cpp)
  target_include_directories(${NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
  if(${ARGC} GREATER 1)
    target_link_libraries(${NAME} ${ARGV1})
  if(${ARGC} GREATER 2)
    target_include_directories(${NAME} PUBLIC ${ARGV2})
  endif()
  endif()
endfunction()

# function to build C utilities
function(add_cutil NAME LIBS INCLUDES)
  add_executable(${NAME} utilities/${NAME}.c)
  target_include_directories(${NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(${NAME} ${LIBS})
  target_include_directories(${NAME} PUBLIC ${INCLUDES})
endfunction()

# these have no dependencies
add_program(custom)
add_program(drive)
add_program(stackedfm)
add_program(stackedpm)
add_program(operatorpm)
add_program(lopass)
add_program(lpwave)
add_program(wave)
add_program(svfdrive)
add_program(oscil)

# these have libsndfile dependency
if(LIBSNDFILE_LIBRARY)
 add_program(delay ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
 add_program(filter ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
 add_program(reverb ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
 add_program(flanger ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
 add_program(chorus ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
 add_program(follow ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
 add_program(equaliser ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
 add_cutil(towav ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
 add_cutil(wtgen ${LIBSNDFILE_LIBRARY} ${LIBSNDFILE_INCLUDE_DIRECTORY})
endif() 
