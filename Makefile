CC = gcc

PROJECT = alvar

FILE = alvar.c

ifndef RELEASE
 # this is triggered by 'make RELEASE=x', with x being whatever.
 CFLAGS = -O0 -g -Wall -Wextra -pedantic -lm -std=c99 -fopenmp
 BFLAGS =
else
 CFLAGS = -O3 -Wall -lm -std=c99 -fopenmp -ffast-math
 BFLAGS = -Wl,--strip-all
endif

ifdef windir
 delcall = del
 binfix = .exe
 pathfix = $(subst /,\,$1)
else
 ifeq ($(shell uname), Linux)
  delcall = rm -f
  binfix =
  pathfix = $1
 else
  ifneq ($(shell uname), Linux)
   $(error OS not identified)
  endif
 endif
endif

$(PROJECT): $(FILE)
	$(CC) $(CFLAGS) $(BFLAGS) -o $@ $^ -lm

clean:
	$(delcall) $(PROJECT)$(binfix)
