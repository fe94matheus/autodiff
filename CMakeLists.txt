cmake_minimum_required(VERSION 3.10)

# Nome do projeto
project(AutoDiff)

# Configuração do padrão C++
set(CMAKE_CXX_STANDARD 23)

 

# Incluir diretórios de cabeçalhos
include_directories(include)

# Adicionar a flag de otimização
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpfr")

# Adicionar o executável
add_executable(AutoDiff src/autodiff.cpp src/main.cpp)

target_link_libraries(AutoDiff mpfr gmp)

# Tarefa customizada para exibir o compilador durante o build
add_custom_target(show_compiler
    COMMAND ${CMAKE_COMMAND} -E echo "Using C++ Compiler: ${CMAKE_CXX_COMPILER}"
)

# Faz com que o alvo "show_compiler" seja executado antes do build completo
add_dependencies(AutoDiff show_compiler)
