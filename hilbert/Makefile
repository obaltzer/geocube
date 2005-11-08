CFLAGS := -pg -funroll-all-loops -funroll-loops -O3 -Wall -ansi -I. -DWITH_FP -DWITHOUT_MP -I/usr/local/include
LDFLAGS := -lm -L/usr/local/lib

all:	test_fp	sort_test tree_test_fp 
# test

tree_test_fp:	hilbert.o quicksort.o sort.o tree.o tree_test.c util.o
	gcc -o tree_test_fp tree_test.c hilbert.o \
		sort.o quicksort.o tree.o util.o $(CFLAGS) $(LDFLAGS)

sort_test:	hilbert.o quicksort.o sort.o sort_main.c util.o norm.o
	gcc -o sort_test sort_main.c hilbert.o \
		sort.o quicksort.o util.o norm.o $(CFLAGS) $(LDFLAGS)

test_fp:	hilbert.o hilbert_test.c
	gcc -o test_fp hilbert_test.c hilbert.o $(CFLAGS) $(LDFLAGS)

test_mp:	hilbert.o hilbert_test.c
	gcc -o test_mp hilbert_test.c hilbert.o \
		$(CFLAGS) $(LDFLAGS)

hilbert.o:	hilbert.c
	gcc -c hilbert.c $(CFLAGS)

quicksort.o:	quicksort.c quicksort.h
	gcc -c quicksort.c $(CFLAGS)

sort.o:	sort.c quicksort.h hilbert.h
	gcc -c sort.c $(CFLAGS)

tree.o:	tree.c sort.h hilbert.h quicksort.h
	gcc -c tree.c $(CFLAGS)

util.o:	util.c sort.h
	gcc -c util.c $(CFLAGS)

norm.o:	hilbert.h sort.h quicksort.h norm.c
	gcc -c norm.c $(CFLAGS)

test:
	@cd tests ; \
	$(MAKE)

test_clean:
	cd tests ; \
	$(MAKE) clean

clean-dir:	
	rm -f test_fp test_mp hilbert.o sort.o quicksort.o \
		tree.o tree_test_fp util.o norm.o

clean:	clean-dir test_clean
