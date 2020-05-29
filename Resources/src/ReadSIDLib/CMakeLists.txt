cmake_minimum_required(VERSION 2.8.9)
project(readSIDlib)

#all source and header files
file(GLOB SOURCES "*.cpp")
file(GLOB HEADER "*.h")

add_definitions(-DEXPORTLIB)

add_library(readSIDlib SHARED
    ${SOURCES}
	${HEADER}
)
