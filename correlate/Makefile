include ../make.config

CC = $(CC_COR)
CFLAGS = $(CFLAGS_COR)
INCLDIRS = -I /usr/include -I ../
LFLAGS = -lm ${LOTHER}
SOURCES = driver.c pcsource.c tree.c bins.c correlate.c
OBJECTS = $(SOURCES:.c=.o)
HEADERS = globals.h correlation.h timing.h ../params.h
EXEC = correlate

.c.o:
	$(CC) $(CFLAGS) $(INCLDIRS) -c $<

all: $(EXEC) 

$(EXEC): $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) $(LFLAGS)

$(OBJECTS): $(HEADERS) Makefile ../make.config

clean:
	rm -f $(EXEC) $(OBJECTS)

