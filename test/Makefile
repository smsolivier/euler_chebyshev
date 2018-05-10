override HOME := ./..

TARGETS = $(wildcard *.cpp) 
TESTS = $(patsubst %.cpp, %.exe, $(TARGETS))

test : $(TESTS)
	@for exe in $(TESTS); do \
		./$$exe ; done 

%.exe : %.cpp $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) $< -o $@ 

include $(HOME)/Makefile 
