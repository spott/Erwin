PETSC_DIR=/usr/local/Cellar/petsc/3.5.3-debug
SLEPC_DIR=/usr/local/Cellar/slepc/3.5.3-debug

include ${SLEPC_DIR}/conf/slepc_common

CPP_FLAGS=-Iinclude/ -I/Users/spott/Code/c++/petsc_cpp_wrapper/include/ -I. -std=c++1y
LDFLAGS= -L/Users/spott/Code/c++/petsc_cpp_wrapper/lib/ -lpetsc_cpp -lgsl -lboost_program_options-mt

SOURCES=src/test/ti_test.cpp
HEADERS=include/time_independent/* include/utilities/*
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=testing/test
DEFAULT=all

all: format $(SOURCES) $(EXECUTABLE)

${EXECUTABLE}: ${OBJECTS}  chkopts
	-${CLINKER} -o ${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

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
