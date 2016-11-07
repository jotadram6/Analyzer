ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXX = g++
CXXFLAGS += -Wall -O2 $(ROOTCFLAGS) -I./ 

LD = g++
LDFLAGS += -Wall -O2 $(ROOTLIBS) -lGenVector

LIBS=

SRCDIR = src
OBJDIR = obj
EXE = Analyzer

#------------------------------------------------------------------------------
SOURCES = $(wildcard src/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

#------------------------------------------------------------------------------

all: $(OBJECTS) 
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJECTS) $(LIBS) 

Analyzer: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) 

obj/main.o: src/main.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(SRCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

%: $(OBJDIR)/%.o
	$(LD) -o $@ $(LDFLAGS) $<  $(LIBS) 


clean :
	rm $(OBJDIR)/*

