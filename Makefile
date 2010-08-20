TARGET=intersec.exe

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
	$(RM) *~
