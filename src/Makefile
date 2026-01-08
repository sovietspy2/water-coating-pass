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
		  pass.h \
		  pass.c \

OBJECTS = $(SOURCES:.c=.o)

all: pass

pass: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< $(LIBFLAGS) -o $@

clean:
	rm *.o pass

install:
	cp pass $(HOME)/bin


