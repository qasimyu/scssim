cmake_minimum_required(VERSION 2.8)
project(SCSsim)

SET(CMAKE_BUILD_TYPE Release)

add_definitions(-std=c++11)

add_subdirectory(lib)

include_directories (
	${SCSsim_SOURCE_DIR}/lib/amplicon
	${SCSsim_SOURCE_DIR}/lib/config
	${SCSsim_SOURCE_DIR}/lib/fastahack
	${SCSsim_SOURCE_DIR}/lib/fragment
	${SCSsim_SOURCE_DIR}/lib/genome
	${SCSsim_SOURCE_DIR}/lib/malbac
	${SCSsim_SOURCE_DIR}/lib/matrix
	${SCSsim_SOURCE_DIR}/lib/mydefine
	${SCSsim_SOURCE_DIR}/lib/profile
	${SCSsim_SOURCE_DIR}/lib/seqwriter
	${SCSsim_SOURCE_DIR}/lib/snp
	${SCSsim_SOURCE_DIR}/lib/split
	${SCSsim_SOURCE_DIR}/lib/threadpool
	${SCSsim_SOURCE_DIR}/lib/vcfparser
)

add_executable(scssim src/scssim.cpp)
#add_executable(simuVars src/simuVars.cpp)
#add_executable(learnProfile src/learnProfile.cpp)
#add_executable(genReads src/genReads.cpp)

target_link_libraries(scssim config genome mydefine profile malbac)
#target_link_libraries(simuVars config genome mydefine)
#target_link_libraries(learnProfile config genome profile mydefine)
#target_link_libraries(genReads config genome malbac profile mydefine)

SET_TARGET_PROPERTIES(scssim PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SCSsim_BINARY_DIR}/bin)
#SET_TARGET_PROPERTIES(simuVars PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SCSsim_BINARY_DIR}/bin)
#SET_TARGET_PROPERTIES(learnProfile PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SCSsim_BINARY_DIR}/bin)
#SET_TARGET_PROPERTIES(genReads PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SCSsim_BINARY_DIR}/bin)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
	set_property(TARGET scssim PROPERTY COMPILE_OPTIONS "-pthread")
	set_property(TARGET scssim PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
	#set_property(TARGET simuVars PROPERTY COMPILE_OPTIONS "-pthread")
	#set_property(TARGET simuVars PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
	#set_property(TARGET learnProfile PROPERTY COMPILE_OPTIONS "-pthread")
	#set_property(TARGET learnProfile PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
	#set_property(TARGET genReads PROPERTY COMPILE_OPTIONS "-pthread")
	#set_property(TARGET genReads PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
	target_link_libraries(scssim "${CMAKE_THREAD_LIBS_INIT}")
	#target_link_libraries(simuVars "${CMAKE_THREAD_LIBS_INIT}")
	#target_link_libraries(learnProfile "${CMAKE_THREAD_LIBS_INIT}")
	#target_link_libraries(genReads "${CMAKE_THREAD_LIBS_INIT}")
endif()

INSTALL(TARGETS scssim
	RUNTIME DESTINATION bin
)
