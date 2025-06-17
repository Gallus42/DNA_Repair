cmake_minimum_required(VERSION 3.19)
project(Dna_Reparation LANGUAGES CXX)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt5 REQUIRED COMPONENTS Core Widgets)
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(ROOT REQUIRED CONFIG)

#TODO исправить статичный путь
include_directories("C:/boost_1_88_0")
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

message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
message(STATUS "ROOT include dirs: ${ROOT_INCLUDE_DIRS}")

include(GNUInstallDirs)
