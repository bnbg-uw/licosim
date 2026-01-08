cmake_minimum_required(VERSION 3.24)
project(licosim VERSION 1.0)

if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})

include(licosim.cmake)

add_executable(licosim ${LICOSIM_SOURCES})
target_include_directories(licosim PRIVATE ${LICOSIM_EXTERNAL_INCLUDES})
target_link_libraries(licosim PRIVATE ${LICOSIM_EXTERNAL_LINKS})
copy_proj_db_after_build(licosim)
copy_licosim_resources_after_build(licosim)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT licosim)

set_target_properties(licosim PROPERTIES
    VS_DEBUGGER_COMMAND_ARGUMENTS "-p F:/inyo_test/test.shp -u F:/inyo_test/test.shp -t 1 -l //172.25.182.82/pfc-frl/01_LiDAR_data/Processed_Lapis/Sierra_B22/del13/ -o G:/inyo_licosim/ --dbhmin 10.16 --dbhmax 76.2 --writeunits -r E:/Dropbox/Inyo_WBI/inyowbiref.csv --allom_slope 0.758902 --allom_intercept 1.69211 --allom_rsq 0.705867 --allom_transform power --allom_inunits feet --allom_outunits inches"
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/licosim.cmake
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/.gitignore
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/README.md
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (MSVC)
	target_compile_options(licosim PRIVATE /W3 /WX)
else()
	target_compile_options(licosim PRIVATE -Wall -Wextra -Werror)
endif()