CC = gcc
CFLAGS = -c -O3 -g -Wall -Wextra -march=native
LDFLAGS =
LIBFLAGS = -lm

SOURCES = main.c \
		  input.c \
		  input.h \
		  global.h \
		  macros.h \
		  print.h \
		  print.c \
		  wdrop.h \
		  wdrop.c \

OBJECTS = $(SOURCES:.c=.o)

all: wdrop

wdrop: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< $(LIBFLAGS) -o $@

clean:
	rm *.o wdrop

install:
	cp wdrop $(HOME)/bin


