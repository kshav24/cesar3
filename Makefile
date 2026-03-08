CC=gcc
LDLIBS=-lm -lc -lpthread
LDFLAGS=${LDLIBS}
FLAGS_DEBUG=-g -Wall -DDEBUG -DLOGLEVEL=0 -Wextra -Wno-unused-function

# CESAR3: OpenMP for parallel Viterbi state loop (--threads N at runtime)
# CESAR3: CHECKPOINT_INTERVAL controls traceback memory footprint (default 1000)
# To change checkpoint interval: make CHECKPOINT_INTERVAL=500
CHECKPOINT_INTERVAL ?= 1000

CFLAGS=-O3 -flto -march=native -std=c11 -Wall -Wno-unused-result -DCHECKPOINT_INTERVAL=$(CHECKPOINT_INTERVAL)
CFLAGS_OMP=$(CFLAGS)

srcfiles=$(wildcard src/*.c)
objects=$(srcfiles:%.c=%.o)
%.o: %.c
	${CC} ${CFLAGS_OMP} -c -o $@ $<


# compile (fast) production version with OpenMP
# CESAR3: compile directly to 'cesar' — avoids 'mv CESAR cesar' failing on
# case-insensitive filesystems (macOS HFS+/APFS) where CESAR == cesar.
default: cesar
.PHONY: cesar
cesar: ${objects}
	${CC} ${CFLAGS} -o cesar $^ ${LDLIBS}

.PHONY: cesar_omp
cesar_omp:
	${CC} ${CFLAGS_OMP} -o cesar_omp ${srcfiles} ${LDLIBS}


# compile without OpenMP (single-threaded, for systems without OpenMP support)
.PHONY: cesar_noomp
cesar_noomp:
	${CC} ${CFLAGS_NOOMP} -o cesar_noomp ${srcfiles} ${LDLIBS}


# compile debugging version
.PHONY: debug
debug: ${objects}
	${CC} ${FLAGS_DEBUG} -fopenmp -o $@ $^ ${LDLIBS}


# compile with AddressSanitizer (for bug detection)
.PHONY: asan
asan:
	${CC} -O1 -std=c11 -fsanitize=address -fno-omit-frame-pointer -fopenmp \
	      -DCHECKPOINT_INTERVAL=$(CHECKPOINT_INTERVAL) \
	      -o cesar_asan ${srcfiles} ${LDLIBS}


# dev setup
workspace:
	virtualenv venv
	pip install -r requirements.txt


# dev doc
.PHONY: doc
doc:
	mkdir -p doc/
	doxygen doxygen.conf


.PHONY: test
test: ${objects} cesar
	tests.sh ./cesar


.PHONY: test/valgrind
test/valgrind: ${objects} cesar
	test/run.sh cesar
	make -C test valgrind

# cleaners
.PHONY: clean
clean: cleanCode cleanTest
	rm -rf doc
	mkdir doc

cleanCode:
	rm -fv *.o
	rm -fv src/*.o
	rm -fv src/*.gch
	rm -f cesar
	rm -f cesar_noomp
	rm -f cesar_asan
	rm -f debug

cleanTest:
	make -i -C test clean

.PHONY: cesar3
cesar3:
	g++ -O3 -std=c++17 -o cesar3 cesar3.cpp -lpthread
