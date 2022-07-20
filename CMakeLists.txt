cmake_minimum_required(VERSION 3.22)
project(untitled)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_EXE main.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(SOURCE_LIB extra_tools.cpp)			# Тоже самое, но для библиотеки

add_library(extra STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем foo

add_executable(untitled ${SOURCE_EXE})	# Создает исполняемый файл с именем main

target_link_libraries(untitled extra)		# Линковка программы с библиотекой


