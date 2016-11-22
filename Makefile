ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXX = g++
CXXFLAGS += -Wall $(ROOTCFLAGS) -I./ -g

LD = g++
LDFLAGS += -Wall $(ROOTLIBS) -lGenVector -g

LIBS=

SRCDIR = src
OBJDIR = obj
EXE = AnalyzerNew

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

