## -*- Makefile -*-

BINARY          := chain

PREFIX=/usr/local
LOCALPREFIX=~
EULERBIN=~/eulerbin

SRCDIR=./src

#### Compiler and tool definitions shared by all build targets #####
CCC = g++
CXX = g++

# Flags
CFLAGS          := 
#CXXFLAGS        :=  -march=native -O3
#CXXFLAGS        := -pedantic -Wall -Wno-sign-compare -O3  -std=c++0x -fopenmp
CXXFLAGS        := -pedantic -Wall -Wno-sign-compare -O3 -g  -std=c++0x -fopenmp -O0
#CXXFLAGS        :=   -g -std=c++0x
LDFLAGS         :=  -fopenmp
GGOFLAGS        :=  -N -e --unamed-opts --output-dir=$(SRCDIR) --src-output-dir=$(SRCDIR)


# Commands
GGO                     := gengetopt
RM                      := rm -f

# Rules
all:	$(BINARY) 

$(BINARY):	$(SRCDIR)/$(BINARY).cpp $(SRCDIR)/$(BINARY)cmdline.c $(SRCDIR)/$(BINARY)cmdline.h $(SRCDIR)/$(tilingsSRC) $(SRCDIR)/$(linfactorsSRC) $(SRCDIR)/$(imageandfileSRC)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(SRCDIR)/$(BINARY).cpp $(SRCDIR)/$(BINARY)cmdline.c -o $(BINARY)

$(SRCDIR)/$(BINARY)cmdline.c: $(SRCDIR)/$(BINARY).ggo	
	-$(GGO) $(GGOFLAGS) -i $(SRCDIR)/$(BINARY).ggo --file-name=$(BINARY)cmdline
#NOTE: THE ABOVE LINE WILL BE SKIPPED OVER IF IT FAILS.  this should be fine if you haven't modified the ggo file.


clean:
	$(RM) $(BINARY) $(SRCDIR)/$(BINARY)cmdline.c $(SRCDIR)/$(BINARY)cmdline.h 
install: $(BINARY)
	install $(BINARY) $(PREFIX)/bin

uninstall:
	-rm $(PREFIX)/bin/$(BINARY)

linstall: $(BINARY) 
	install $(BINARY) $(LOCALPREFIX)/bin

luninstall:
	-rm $(LOCALPREFIX)/bin/$(BINARY)


.PHONY: install
