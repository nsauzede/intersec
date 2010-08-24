TARGET=intersec.exe

CFLAGS=-Wall -Werror

OPTIM=1
ifdef OPTIM
CFLAGS+=-O2
else
CFLAGS+=-g -O0
endif

#USE32=1
ifdef USE32
CFLAGS+=-m32
LDFLAGS+=-m32
endif

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
