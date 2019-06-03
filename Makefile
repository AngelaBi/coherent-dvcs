ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I$(CLAS12TOOL)/Hipo3 -I$(CLAS12TOOL)/Clas12Banks3
HIPOLIBS    := -L$(CLAS12TOOL)/lib -lHipo3 -lClas12Banks3

LZ4LIBS     := -L$(CLAS12TOOL)/Lz4/lib -llz4
LZ4INCLUDES := -I$(CLAS12TOOL)/Lz4/lib

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)


all: readBanks readParticles analysis clas12event_example mesonexevent_example

readBanks: readBanks.o
	$(CXX) -o readBanks $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

readParticles: readParticles.o
	$(CXX) -o readParticles $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

analysis: analysis.o
	$(CXX) -o analysis $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clas12event_example: clas12event_example.o
	$(CXX) -o clas12event_example $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

mesonexevent_example: mesonexevent_example.o
	$(CXX) -o mesonexevent_example $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o readParticles readBanks analysis clas12event_example mesonexevent_example *~

%.o: %.cc
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
