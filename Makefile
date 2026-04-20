ifdef SANITIZE
     CFLAGS += -fsanitize=address -fsanitize=undefined
     LDFLAGS += -fsanitize=address -fsanitize=undefined
endif


CFLAGS += -Wall -O3 -g
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  CFLAGS += -I/opt/local/include
  LDFLAGS  += -L/opt/local/lib
  SOFLAGS = -bundle -Wl,-undefined,dynamic_lookup
else
  LDLIBS += -lbsd
  SOFLAGS = -shared
endif


all: fmt

# install in user's home directory
install:
	cp -i fmt ~/bin/

fmt: fmt.o load_wav.o
	cc $(LDFLAGS) -o fmt fmt.o load_wav.o $(LDLIBS) -lfftw3f -lm

clean:
	rm -f fmt *.o

fmt.o: fmt.c

load_wav.o: load_wav.c
