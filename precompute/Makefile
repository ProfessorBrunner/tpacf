include ../make.config

CC = $(CCSERIAL)
CFLAGS = $(CFLAGSSERIAL)
INCLDIRS = -I /usr/include -I ../
LFLAGS = -lm # -lsvml
SOURCES = driver.c pcsource.c tree.c
OBJECTS = $(SOURCES:.c=.o)
HEADERS = precompute.h ../params.h
EXEC = precompute

.c.o:
	$(CC) $(CFLAGS) $(INCLDIRS) -c $<

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) $(LFLAGS)

$(OBJECTS): $(HEADERS) Makefile ../make.config

clean:
	rm -f $(EXEC) *~ ../include/*~ $(OBJECTS)

