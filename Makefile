TARGET=intersec.exe
TARGET+=plane.exe
TARGET+=fake3d.exe
TARGET+=bench.exe
TARGET+=fbench.exe

CFLAGS=-Wall -Werror

USE_SDL=1
#USE_CROSS=1

ifdef USE_CROSS
CCCROSS=arm-none-linux-gnueabi-
SDLCROSS=~/nico/install/webos/
endif
CC=$(CCCROSS)gcc
SDLCONFIG=$(SDLCROSS)sdl-config

CFLAGS+=-g
OPTIM=1
ifdef OPTIM
#CFLAGS+=-O2
CFLAGS+=-O3
else
CFLAGS+=-O0
endif

#USE32=1
ifdef USE32
ifndef USE_CROSS
CFLAGS+=-m32
LDFLAGS+=-m32
endif
endif

#USE_SKY=1
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

ifdef USE_DOUBLE
CFLAGS+=-DUSE_DOUBLE
endif

MV=mv

LDFLAGS=-lm

all:$(TARGET)

check: $(TARGET)
	./$(TARGET)

ifdef USE_SDL
intersec.exe: CFLAGS+=-DUSE_SDL
intersec.exe:CFLAGS+=`$(SDLCONFIG) --cflags`
#intersec.exe:LDFLAGS+=`$(SDLCONFIG) --libs` -mno-windows
intersec.exe:LDFLAGS+=`$(SDLCONFIG) --libs`
fake3d.exe:CFLAGS+=`$(SDLCONFIG) --cflags`
fake3d.exe:LDFLAGS+=`$(SDLCONFIG) --libs`

#%.exe:CFLAGS+=`$(SDLCONFIG) --cflags`
#%.exe:LDFLAGS+=`$(SDLCONFIG) --libs`

#plane.exe: CFLAGS+=-DUSE_SDL
#plane.exe:CFLAGS+=`$(SDLCONFIG) --cflags`
#plane.exe:LDFLAGS+=`$(SDLCONFIG) --libs`
endif

intersec.exe: intersec.o vec.o
fake3d.exe: fake3d.o vec.o
bench.exe: bench.o vec.o
fbench.exe: fbench.o fvec.o

%.exe: %.o
	$(CC) -o $@ $^ $(LDFLAGS)

sdlr3:CFLAGS=`$(SDLCONFIG) --cflags`
sdlr3:LDFLAGS+=`$(SDLCONFIG) --libs`
sdlr3:sdlr3.c
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

clean:
	$(RM) $(TARGET) *.o
clobber: clean
	$(RM) *~ *.tga *.exe
