PETSC_DIR=/usr/local/Cellar/petsc/3.5.3-debug
SLEPC_DIR=/usr/local/Cellar/slepc/3.5.3-debug

include ${SLEPC_DIR}/conf/slepc_common

CPP_FLAGS=-Iinclude/ -I/Users/spott/Code/c++/petsc_cpp_wrapper/include/ -I. -std=c++1y -g -O0
LDFLAGS= -L/Users/spott/Code/c++/petsc_cpp_wrapper/lib/ -lpetsc_cpp -lgsl -lboost_program_options-mt -lboost_iostreams-mt

basis_src=src/test/basis_test.cpp
hamiltonian_src=src/test/dipole_test.cpp
parameters_src=src/parameters/basis.cpp src/parameters/hamiltonian.cpp
utilities_src=src/utilities/types.cpp src/utilities/math.cpp
basis_objects=$(basis_src:.cpp=.o)
hamiltonian_objects=$(hamiltonian_src:.cpp=.o)
parameters_objects=$(parameters_src:.cpp=.o)
utilities_objects=$(utilities_src:.cpp=.o)
SOURCES=${basis_src} ${hamiltonian_src} ${parameters_src} ${utilities_src}
HEADERS=include/time_independent/* include/utilities/* include/parameters/*
OBJECTS=$(SOURCES:.cpp=.o)
executables=testing/test_basis testing/test_hamiltonian
DEFAULT=all

all: format $(SOURCES) $(executables)

testing/test_basis: ${basis_objects} ${parameters_objects} ${utilities_objects} chkopts
	-${CLINKER} -o $@ ${basis_objects} ${parameters_objects} ${utilities_objects} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

testing/test_hamiltonian: ${hamiltonian_objects} ${parameters_objects} ${utilities_objects} chkopts
	-${CLINKER} -o $@ ${hamiltonian_objects} ${parameters_objects} ${utilities_objects} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}


syntax_check: chkopts
	clang++ -fsyntax-only ${SOURCES} -I${SLEPC_DIR}/include/ -I./include/ -std=c++1y -I${PETSC_DIR}/include/

format:
	clang-format -style=file -i ${SOURCES}
	clang-format -style=file -i ${HEADERS}

cleanup:
	${RM} ${OBJECTS}
#clean:
#	${RM} ${OBJECTS}
#
#clean:
#	rm *.o;
#	rm src/petsc_cpp/*.o
