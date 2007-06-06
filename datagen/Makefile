CFLAGS = -Wall -O3 -ansi
CC = gcc
all: datagen

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c

datagen: datagen.o
	$(CC) -o datagen datagen.o -lc
	
clean:	datagen-clean

datagen-clean:
	rm -f datagen
