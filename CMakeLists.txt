cmake_minimum_required(VERSION 3.22)
project(untitled)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS "-lm -lfftw3")

#set(CMAKE_C_FLAGS   "-lm -lfftw3")

set(SOURCE_LIB extra_tools.cpp base_algorithm.cpp base_algorithm.h)			# Тоже самое, но для библиотеки

set(SOURCE_EXE main.cpp)		# Установка переменной со списком исходников для исполняемого файла

#INCLUDE_DIRECTORIES(/usr/include)
#LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu)
#add_library(fftw3 STATIC IMPORTED)
#set(CMAKE_C_OUTPUT_EXTENSION_REPLACE 1)
#set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

add_library(extra STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем extra

add_executable(untitled ${SOURCE_EXE})	# Создает исполняемый файл с именем untitled

#target_link_libraries(untitled fftw3 extra)		# Линковка программы с библиотекой

find_library(FFTW_LIBRARY
        NAMES fftw3 fftw)
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
if(UNIX AND NOT WIN32)
    find_library(FFTW_libm_LIBRARY
            NAMES m)
    list(APPEND FFTW_LIBRARIES "${FFTW_libm_LIBRARY}")
endif()

#target_link_libraries(untitled ${FFTW_LIBRARIES} extra)
target_link_libraries(extra ${FFTW_LIBRARIES})
target_link_libraries(untitled extra)