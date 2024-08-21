CXX = g++ -std=c++2a	
CXXFLAGS = -fopenmp -pthread -march=native -mtune=native -O3 -lgurobi_c++ -lgurobi110 -lm -lz -lpthread -ldl
GUROBI_HOME = /home/ghanshyam/opt/gurobi1101/linux64
INLCLUDES = -I$(GUROBI_HOME)/include
LIBS = -L$(GUROBI_HOME)/lib

all: AlphaASM

src_dir := src

OBJS = $(src_dir)/main.o $(src_dir)/gfa-io.o $(src_dir)/gfa-base.o \
		$(src_dir)/options.o $(src_dir)/kalloc.o \
		$(src_dir)/misc.o $(src_dir)/sys.o $(src_dir)/ILP_index.o \
		$(src_dir)/MurmurHash3.o

AlphaASM: $(OBJS)
	$(CXX) $^ -o $@ $(INLCLUDES) $(LIBS) $(CXXFLAGS)

$(src_dir)/%.o: $(src_dir)/%.cpp
	$(CXX) -c $< -o $@ $(INLCLUDES) $(LIBS) $(CXXFLAGS)

clean:
	rm -f $(src_dir)/*.o AlphaASM