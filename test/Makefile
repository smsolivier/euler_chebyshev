override HOME := ./..
include $(HOME)/Makefile 

TARGETS = $(wildcard *.cpp) 
TARGETS := $(filter-out 2d.cpp steady_heat.cpp, $(TARGETS))
TESTS = $(patsubst %.cpp, %.exe, $(TARGETS))

test : $(TESTS)
	@for exe in $(TESTS); do \
		./$$exe ; done 

%.exe : %.cpp $(OBJS) 
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) $< -o $@ 

targets :
	@echo $(TARGETS)
