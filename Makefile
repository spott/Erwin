PETSC_CPP_HOME=${PWD}/../petsc_cpp/

include ${SLEPC_DIR}/conf/slepc_common


RELEASE_FLAGS=-Wall -Wpedantic -Wextra
DEBUG_FLAGS=-fdiagnostics-show-template-tree -Wall -Wpedantic -Wextra -Werror -Wbind-to-temporary-copy -Weverything -Wno-c++98-compat-pedantic -Wno-old-style-cast -Wno-error=padded
CPP_FLAGS=-Iinclude/ -I${PETSC_CPP_HOME}/include/ -I. -std=c++1y -g -O0
LDFLAGS= -L${PETSC_CPP_HOME}/lib/ -lpetsc_cpp -lgsl -lboost_program_options-mt -lboost_iostreams-mt

basis_src=src/test/basis_test.cpp
hamiltonian_src=src/test/dipole_test.cpp
propagate_src=src/test/propagate_test.cpp
parameters_src=src/parameters/basis.cpp src/parameters/hamiltonian.cpp src/parameters/laser.cpp
utilities_src=src/utilities/types.cpp src/utilities/math.cpp
basis_objects=$(basis_src:.cpp=.o)
hamiltonian_objects=$(hamiltonian_src:.cpp=.o)
propagate_objects=$(propagate_src:.cpp=.o)
parameters_objects=$(parameters_src:.cpp=.o)
utilities_objects=$(utilities_src:.cpp=.o)
SOURCES=${basis_src} ${hamiltonian_src} ${propagate_src} ${parameters_src} ${utilities_src}
HEADERS=include/time_independent/* include/utilities/* include/parameters/*
OBJECTS=$(SOURCES:.cpp=.o)
executables=testing/test_basis testing/test_hamiltonian testing/test_propagate
ifeq ($(UNAME), Linux)
	DEFAULT=release
endif
ifeq ($(UNAME), Darwin)
	DEFAULT=debug
endif

debug: CPP_FLAGS += ${DEBUG_FLAGS}
debug: format $(SOURCES) $(executables)

release: CPP_FLAGS += ${RELEASE_FLAGS}
release: $(SOURCES) $(executables)

testing/test_basis: ${basis_objects} ${parameters_objects} ${utilities_objects} chkopts
	-${CLINKER} -o $@ ${basis_objects} ${parameters_objects} ${utilities_objects} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

testing/test_hamiltonian: ${hamiltonian_objects} ${parameters_objects} ${utilities_objects} chkopts
	-${CLINKER} -o $@ ${hamiltonian_objects} ${parameters_objects} ${utilities_objects} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

testing/test_propagate: ${propagate_objects} ${parameters_objects} ${utilities_objects} chkopts
	-${CLINKER} -o $@ ${propagate_objects} ${parameters_objects} ${utilities_objects} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

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
