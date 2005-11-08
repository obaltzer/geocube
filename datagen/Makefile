CFLAGS = -Wall -O3 -ansi

all: datagen

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c

datagen: datagen.o
	$(CC) -o datagen datagen.o
	
clean:	datagen-clean

datagen-clean:
	rm -f datagen
