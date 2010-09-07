TARGET=intersec.exe

CFLAGS=-Wall -Werror

OPTIM=1
ifdef OPTIM
CFLAGS+=-O2
else
CFLAGS+=-g -O0
endif

USE32=1
ifdef USE32
CFLAGS+=-m32
LDFLAGS+=-m32
endif

USE_SKY=1
ifdef USE_SKY
CFLAGS+=-DUSE_SKY=$(USE_SKY)
endif

#DEBUG=1
ifdef DEBUG
CFLAGS+=-DDEBUG=$(DEBUG)
endif

USETIMEOFDAY=1
ifdef USETIMEOFDAY
CFLAGS+=-DUSETIMEOFDAY=$(USETIMEOFDAY)
endif

MV=mv

LDFLAGS=-lm

all:$(TARGET)

check: $(TARGET)
	./$(TARGET)

intersec.exe: intersec
	$(MV) -f $< $@

sdlr3:CFLAGS=`sdl-config --cflags` -g -O2
sdlr3:LDFLAGS+=`sdl-config --libs`
sdlr3:sdlr3.c
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

clean:
	$(RM) $(TARGET)
clobber: clean
	$(RM) *~ *.tga
