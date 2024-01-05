project(mmio LANGUAGES C)

add_library(mmio STATIC mmio.c)

target_include_directories(mmio PUBLIC include)