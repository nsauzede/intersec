TARGET=intersec.exe

LDFLAGS=-lm

all:$(TARGET)

check: $(TARGET)
	./$(TARGET)

intersec.exe: intersec

clean:
	$(RM) $(TARGET)
