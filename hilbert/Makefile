CFLAGS := -O3 -g -pg -Wall -ansi -I. -DWITH_FP
LDFLAGS := -lgmp -lm

all:	test_fp	test_mp sort_test_fp tree_test_fp

tree_test_fp:	hilbert.o quicksort.o sort.o tree.o tree_test.c
	gcc -o tree_test_fp tree_test.c hilbert.o \
		sort.o quicksort.o tree.o $(CFLAGS) $(LDFLAGS)

sort_test_fp:	hilbert.o quicksort.o sort.o sort_test.c
	gcc -o sort_test_fp sort_test.c hilbert.o \
		sort.o quicksort.o $(CFLAGS) $(LDFLAGS)

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

clean:
	rm -f test_fp test_mp sort_test_fp hilbert.o sort.o quicksort.o \
		tree.o tree_test_fp
