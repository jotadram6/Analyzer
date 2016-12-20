ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXX = g++
CXXFLAGS += -Wall $(ROOTCFLAGS) -I./ -O3

LD = g++
LDFLAGS += -Wall $(ROOTLIBS) -lGenVector -O3

LIBS=

SRCDIR = src
OBJDIR = obj
EXE = Analyzer

#------------------------------------------------------------------------------
SOURCES = $(wildcard src/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

#------------------------------------------------------------------------------

all: $(EXE)

$(EXE): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) 

obj/main.o: src/main.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(SRCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

%: $(OBJDIR)/%.o
	$(LD) -o $@ $(LDFLAGS) $<  $(LIBS) 


clean :
	rm $(OBJDIR)/*

