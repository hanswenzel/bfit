CC = gcc
CFLAGS = -c -O      

BASE_DIR = .

LIBRARY = $(BASE_DIR)/libconf.a

OBJS = libconf.o

all: $(LIBRARY) conf_example

$(LIBRARY): $(OBJS)
	@echo Updating library libconf.a ...
	$(AR) -r $(LIBRARY) $(OBJS)

libconf.o: libconf.c libconf.h
	@echo Compiling libconf.c ...
	$(CC) $(CFLAGS) libconf.c

conf_example: conf_example.c libconf.a
	@echo Building example program
	$(CC) conf_example.c -L. -lconf -o conf_example

clean: 
	@rm -f *.o *~ libconf.a conf_example
