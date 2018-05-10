HOME = .

include $(HOME)/make.inc

SRC = $(HOME)/src
UTILS = $(HOME)/utils
WRITER = $(UTILS)/VisitWriter
CHTIMER = $(UTILS)/chtimer
EIGEN = $(UTILS)/eigen

# import settings from make.inc 
CFLAGS += -DDIM=$(DIM) $(OPT)

# openmp 
ifdef OMP 
CFLAGS += -fopenmp -DOMP 
endif

# use FFTW_MEASURE to optimize 1D FFT calls 
ifdef MEASURE 
CFLAGS += -DMEASURE
endif

# write VTK files in text 
ifdef VISIT_ASCII
CFLAGS += -DVISIT_ASCII=0 
else 
CFLAGS += -DVISIT_ASCII=1
endif 

ifdef NWRITE 
CFLAGS += -DNWRITE 
endif

# zero highest mode out 
ifdef ZERO 
CFLAGS += -DZERO
endif

# use paper cutter method v regular LU 
ifdef PaperCutter
CFLAGS += -DPC 
endif

# turn on tau's
ifdef UseTau
CFLAGS += -DTAU 
endif

# enable CH_Timer 
ifdef CH_NTIMER 
CFLAGS += -DCH_NTIMER
endif

# enable checks 
ifdef ErrorCheck 
CFLAGS += -DErrorCheck=$(ErrorCheck) 
endif

# enable warnings 
ifdef WARNINGS
CFLAGS += -DWARNINGS 
endif

# FFTW setup 
FFTW_INC = -I$(FFTW_HOME)/include 
FFTW_LIB = -L$(FFTW_HOME)/lib -lfftw3 
ifdef FFTW_HOME
CFLAGS += -DFFTW $(FFTW_INC)
LIBS += $(FFTW_LIB) 
endif 

# look for source files in
VPATH = $(SRC) $(WRITER) $(CHTIMER) $(EIGEN)/src 
# look for includes in 
CFLAGS += -I$(SRC) -I$(WRITER) -I$(CHTIMER) -I$(EIGEN)

OBJDIR = $(HOME)/obj
DEPDIR = $(HOME)/dep

SRCFILES = $(notdir $(wildcard $(SRC)/*.cpp $(WRITER)/*.cpp $(CHTIMER)/*.cpp)) 
OBJS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(SRCFILES))
DEPS = $(patsubst $(OBJDIR)/%.o, $(DEPDIR)/%.d, $(OBJS))

# don't delete intermediate files (.o's) 
.SECONDARY :

$(OBJDIR)/%.o : %.cpp $(HOME)/Makefile $(HOME)/make.inc
	mkdir -p $(OBJDIR); $(CXX) -c $(CFLAGS) $(LIBS) $(UPC) $< -o $@
	mkdir -p $(DEPDIR) 
	$(CXX) -MM $(CFLAGS) $(LIBS) $(UPC) $< | sed -e '1s@^@$(OBJDIR)\/@' > $*.d
	mv $*.d $(DEPDIR)

clean :
	rm -f *.exe *.vtk *.visit time.table *.3D
cleantree :
	rm -rf $(DEPDIR) $(OBJDIR) $(HOME)/test/*.vtk $(HOME)/test/*.visit \
		$(HOME)/test/*.exe $(HOME)/exec/*.vtk $(HOME)/exec/*.exe \
		$(HOME)/exec/*.visit $(HOME)/exec/time.table $(HOME)/test/time.table

-include $(DEPS)

.PHONY : docs
docs :
	cd $(HOME)/docs; doxygen Doxyfile 
objs :
	@echo $(OBJS)
flags : 
	@echo $(CFLAGS)
