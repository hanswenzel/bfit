ObjSuf		= o
SrcSuf		= cc
ExeSuf		= 
DllSuf		= so
OutPutOpt	= -o
HeadSuf		= h


ROOTCFLAGS 	= $(shell root-config --cflags)
ROOTLIBS      	= $(shell root-config --libs) -lMinuit
ROOTGLIBS     	= $(shell root-config --glibs) -lMinuit


CXX           = g++
CXXFLAGS      = -O -Wall -fPIC -I. -I../config_file
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) $(LIBCONF)
LIBCONF       = -L../config_file -lconf
GLIBS         = $(ROOTGLIBS) $(SYSLIBS) $(LIBCONF)
CONFIGFILE    = ../config_file/libconf.a
#------------------------------------------------------------------------------

BDECAYO       = BDecay.$(ObjSuf) BDecayDict.$(ObjSuf) 
BDECAYS       = BDecay.$(SrcSuf) BDecayDict.$(SrcSuf) 
BDECAYSO      = BDecay.$(DllSuf) 

MEASUREMENTO       = Measurement.$(ObjSuf) MeasurementDict.$(ObjSuf) 
MEASUREMENTS       = Measurement.$(SrcSuf) MeasurementDict.$(SrcSuf) 
MEASUREMENTSO      = Measurement.$(DllSuf)


BEXCLFITO       = bexclfit.$(ObjSuf) bexclfit_fcts.$(ObjSuf)
BEXCLFITS       = bexclfit.$(SrcSuf) bexclfit_fcts.$(SrcSuf)
BEXCLFIT        = bexclfit$(ExeSuf)

TESTO       = test.$(ObjSuf) bexclfit_fcts.$(ObjSuf)
TESTS       = test.$(SrcSuf) bexclfit_fcts.$(SrcSuf)
TEST        = test$(ExeSuf)

MKBDECAYO	= mkbdecay.$(ObjSuf) bexclfit_fcts.$(ObjSuf)
MKBDECAYS	= mkbdecay.$(SrcSuf) bexclfit_fcts.$(SrcSuf)
MKBDECAY	= mkbdecay$(ExeSuf)

H1000O		= h1000.$(ObjSuf)
H1000S		= h1000.$(SrcSuf) 
H1000SO		= libh1000.$(DllSuf)


OBJS	=   $(MEASUREMENTO)  $(BDECAYO) $(BEXCLFITO) $(TESTO) $(H1000O)

INCLUDES = BDecay.h Measurement.h  bexclfit_fcts.h global.h

PROGRAMS=   $(MEASUREMENTSO) $(BDECAYSO) $(BEXCLFIT) $(TEST) $(MKBDECAY)

all : $(CONFIGFILE) $(PROGRAMS) 

$(CONFIGFILE):
		(cd ../config_file; make)

$(MKBDECAY):    $(MKBDECAYO)  $(H1000SO) $(MEASUREMENTO) $(BDECAYSO) $(BDECAYO)
		$(LD) $(LDFLAGS) $(MKBDECAYO) $(MEASUREMENTO) $(BDECAYO) $(H1000O) $(LIBS) $(OutPutOpt) $(MKBDECAY)
		@cp mkbdecay ..
		@echo "$(MKBDECAY) done"

$(H1000SO):	$(H1000O)
		$(LD) $(SOFLAGS) $(LDFLAGS) $(H1000O) $(OutPutOpt) $(H1000SO)
		@echo "$(H1000SO) done"


$(BEXCLFIT):  	$(BEXCLFITO) $(BDECAYSO) $(MEASUREMENTSO) $(MEASUREMENTO) $(BDECAYO)
	      	$(LD) $(LDFLAGS) $(BEXCLFITO) $(MEASUREMENTO) $(BDECAYO) $(GLIBS) $(OutPutOpt) $(BEXCLFIT)
		@cp bexclfit ..
	      	@echo "$(BEXCLFIT) done"

$(TEST):  	$(TESTO) $(BDECAYSO) $(MEASUREMENTSO) $(MEASUREMENTO) $(BDECAYO)
	      	$(LD) $(LDFLAGS) $(TESTO) $(MEASUREMENTO) $(BDECAYO) $(GLIBS) $(OutPutOpt) $(TEST)
		@cp test ..
	      	@echo "$(TEST) done"

$(BDECAYSO):      $(BDECAYO) $(BDECAYS)  $(MEASUREMENTO) $(MEASUREMENTS) $(INCLUDES)
	        $(LD) $(SOFLAGS) $(LDFLAGS) $(BDECAYO) $(MEASUREMENTO) $(OutPutOpt) $(BDECAYSO)
	        @echo "$(BDECAYSO) done"

$(MEASUREMENTSO):      $(MEASUREMENTO) $(MEASUREMENTS) $(INCLUDES)
	       $(LD) $(SOFLAGS) $(LDFLAGS) $(MEASUREMENTO) $(OutPutOpt) $(MEASUREMENTSO)
	       @echo "$(MEASUREMENTSO) done"


clean:
	@rm -f $(OBJS) $(PROGRAMS) *Dict.* core *.o *.so

tidy:
	@rm -f $(OBJS)  *Dict.* core *.o *.so *~ $(PROGRAMS) mkbdecay
.SUFFIXES: .$(SrcSuf) .$(ObjSuf)



BDecay.$(ObjSuf): BDecay.h 
BDecayDict.$(SrcSuf): BDecay.h
	@echo "Generating dictionary BDecayDict..."
	@$(ROOTSYS)/bin/rootcint -f  BDecayDict.$(SrcSuf) -c  BDecay.h

BDecay.$(ObjSuf): Measurement.h
MeasurementDict.$(SrcSuf): Measurement.h
	@echo "Generating dictionary MeasurementDict..."
	@$(ROOTSYS)/bin/rootcint -f  MeasurementDict.$(SrcSuf) -c Measurement.h

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<





