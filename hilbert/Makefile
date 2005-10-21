all:	test_fp	test_mp sort_test_fp tree_test_fp

tree_test_fp:	hilbert.o quicksort.o sort.o tree.o tree_test.c
	gcc -pg -g -Wall -ansi -o tree_test_fp tree_test.c hilbert.o \
		sort.o quicksort.o tree.o -I. -DWITH_FP -lgmp

sort_test_fp:	hilbert.o quicksort.o sort.o sort_test.c
	gcc -pg -g -Wall -ansi -o sort_test_fp sort_test.c hilbert.o \
		sort.o quicksort.o -I. -DWITH_FP -lgmp

test_fp:	hilbert.o hilbert_test.c
	gcc -pg -Wall -ansi -o test_fp hilbert_test.c hilbert.o -I. -DWITH_FP -lgmp

test_mp:	hilbert.o hilbert_test.c
	gcc -pg -Wall -ansi -o test_mp hilbert_test.c hilbert.o -I. -lgmp

hilbert.o:	hilbert.c
	gcc -pg -g -Wall -ansi -c hilbert.c -I.

quicksort.o:	quicksort.c quicksort.h
	gcc -pg -g -Wall -ansi -c quicksort.c

sort.o:	sort.c quicksort.h hilbert.h
	gcc -pg -g -Wall -ansi -c sort.c -I.

tree.o:	tree.c sort.h hilbert.h quicksort.h
	gcc -pg -g -Wall -ansi -c tree.c -I.

clean:
	rm -f test_fp test_mp sort_test_fp hilbert.o sort.o quicksort.o \
		tree.o tree_test_fp
