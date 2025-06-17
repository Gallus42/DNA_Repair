cmake_minimum_required(VERSION 3.19)
project(Dna_Reparation LANGUAGES CXX)

cmake_policy(SET CMP0167 OLD)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt5 REQUIRED COMPONENTS Core Widgets)
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(ROOT REQUIRED CONFIG)
find_package(Boost)

#TODO исправить статичный путь
# include_directories("C:/boost_1_88_0")
include(${ROOT_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)

add_executable(
    Dna_Reparation #WIN32
    main.cpp
    mainwindow.cpp
    mainwindow.h
    ${sources}
    ${headers}
)

# Инклуд-пути
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}
        ${Geant4_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${Boost_INCLUDE_DIR}
)

# Линкуем с Qt и Geant4
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        Qt::Core
        Qt::Widgets
        ${Geant4_LIBRARIES}
        ${ROOT_LIBRARIES}
        # DNARepair
)

message(STATUS "Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}")
message(STATUS "Boost_LIBRARY_DIR = ${Boost_LIBRARY_DIR}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")

include(GNUInstallDirs)
