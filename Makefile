TARGET=intersec.exe

CFLAGS=-Wall -Werror
CFLAGS+=-g -O0

MV=mv

LDFLAGS=-lm

all:$(TARGET)

check: $(TARGET)
	./$(TARGET)

intersec.exe: intersec
	$(MV) -f $< $@

clean:
	$(RM) $(TARGET)
clobber: clean
	$(RM) *~ *.tga
