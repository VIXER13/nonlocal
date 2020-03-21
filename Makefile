TARGET = a.out
HDRS = include \
	   include/containers \
	   include/finite_elements \
	   Eigen

SRC = src/main.cpp
	  

.PHONY: all clean

all: $(SRCS)
		g++ -std=c++17 -fopenmp -O2 -Wall -Wextra -m64 -I${MKLROOT}/include $(addprefix -I,$(HDRS)) -o $(TARGET) $(CXXFLAGS) $(SRC) \
		-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

clean:
