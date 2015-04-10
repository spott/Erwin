PETSC_CPP_HOME=../petsc-cpp/

include ${SLEPC_DIR}/lib/slepc-conf/slepc_variables

release_cpp_flags = -Wall -Wpedantic -Wextra -fdiagnostics-color=auto
clang_cpp_flags   = -fdiagnostics-show-template-tree -Wbind-to-temporary-copy       \
										-Weverything -Wno-c++98-compat-pedantic -Wno-error=weak-vtables \
										-Wno-error=exit-time-destructors -DDEBUG
debug_cpp_flags   = ${release_cpp_flags} -Werror -Wno-old-style-cast -Wno-padded \
									  -Wno-deprecated-declarations

CPP_FLAGS         = -I./include/ -I${PETSC_CPP_HOME}include -std=c++1y ${SLEPC_CC_INCLUDES} ${PETSC_CC_INCLUDES}

boost							= -lboost_program_options-mt -lboost_iostreams-mt
gsl								= -lgsl
LD_FLAGS          = -L${PETSC_CPP_HOME}lib/ -lpetsc_cpp ${boost} ${gsl} ${SLEPC_LIB} ${PETSC_LIB} 

#Directories
source     = src
includes   = include
parameters = parameters
utilities  = utilities
build      = build
testing    = testing
test       = test

basis_src          = basis_test.cpp
hamiltonian_src    = dipole_test.cpp
propagate_src      = propagate_test.cpp
output_src         = check_prototype.cpp
dipole_test_src    = check_dipole.cpp

parameters_src     = basis.cpp hamiltonian.cpp laser.cpp propagate.cpp absorber.cpp dipole.cpp eigenstates.cpp
parameters_objects = ${patsubst %.cpp, ${build}/${parameters}/%.o, ${parameters_src}}

utilities_src      = types.cpp math.cpp
utilities_objects  = ${patsubst %.cpp, ${build}/${utilities}/%.o, ${utilities_src}}

executables        = ${patsubst %.cpp, ${testing}/%, ${basis_src} ${hamiltonian_src} ${propagate_src} ${output_src} ${dipole_test_src}}

clang: CXX=clang++
clang: CPP_FLAGS += ${clang_cpp_flags} ${debug_cpp_flags}
clang: $(executables)

gpp: OMPI_CXX=g++-4.9
gpp: CPP_FLAGS += ${debug_cpp_flags}
gpp: $(executables)

release: CPP_FLAGS += ${release_cpp_flags}
release: $(executables)

${testing}/%: ${build}/${test}/%.o ${parameters_objects} ${utilities_objects}
	mpicxx -o $@ $^ ${LD_FLAGS}

${build}/${test}/%.o: ${source}/${test}/%.cpp
	@mkdir -p ${dir $@}
	mpicxx -o $@ -c $< ${CPP_FLAGS}

${build}/${utilities}/%.o: ${source}/${utilities}/%.cpp
	@mkdir -p ${dir $@}
	mpicxx -o $@ -c $< ${CPP_FLAGS}

${build}/${parameters}/%.o: ${source}/${parameters}/%.cpp
	@mkdir -p ${dir $@}
	mpicxx -o $@ -c $< ${CPP_FLAGS}

${source}/${test}/%.cpp: ${includes}/time_independent/*.hpp
	-clang-format -style=file -i $@

${source}/${utilities}/%.cpp: ${includes}/${utilities}/%.hpp
	-clang-format -style=file -i $@

${source}/${parameters}/%.cpp: ${includes}/${parameters}/%.hpp
	-clang-format -style=file -i $@

%.hpp:
	-clang-format -style=file -i $@

syntax_check: chkopts
	mpicxx -fsyntax-only ${SOURCES} ${CPP_FLAGS_} ${CLANG_ONLY_FLAGS} ${DEBUG_FLAGS} -I${SLEPC_DIR}/include/ -I${PETSC_DIR}/include/

clean:
	rm -rf ${build}
	rm -f ${executables}

print-%: ; @echo $*=$($*)

.PHONEY: clean print-%
