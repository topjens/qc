CC=gcc
WEAVE=cweave
TANGLE=ctangle
TEX=pdftex

CFLAGS=
LIBS=-lm

TARGET=0m 0m.pdf

all: $(TARGET)

0m: 0m.o
	$(CC) $(CFLAGS) $(LIBS) -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.c: %.w
	$(TANGLE) $<

%.pdf: %.tex
	$(TEX) -interaction nonstopmode $<

%.tex: %.w
	$(WEAVE) $<

clean:
	rm -f *.tex *.log *.idx *.scn *.toc *.c *.o $(TARGET)
