CC = gcc
CFLAGS = -Os -std=gnu99 -Wall -D_XOPEN_SOURCE=600 -I/usr/local/trmm/include -I/usr/local/include
LDFLAGS = -L/usr/local/trmm/lib -L/usr/local/lib -lrsl -lnetcdf

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
LDFLAGS += -lhdf5_hl -lhdf5
endif
LDFLAGS += -lz -lm

PROGS = radar2nc

ALL: $(PROGS)

radar2nc: radar2nc.c
	$(CC) $(CFLAGS) -o $@ $@.c $(LDFLAGS)

clear:
	rm -rf $(PROGS)
