CXX = g++ -std=c++2a	
CXXFLAGS = -ggdb3 -fopenmp -pthread -march=native -O3 -ljemalloc -lgurobi_c++ -lgurobi110 -lm -lz -lpthread -ldl -lstdc++fs
INLCLUDES = -I/home/ghanshyam/opt/gurobi1101/linux64/include
LIBS = -L/home/ghanshyam/opt/gurobi1101/linux64/lib

all: AlphaASM

src_dir := src

OBJS = $(src_dir)/main.o $(src_dir)/sketch.o $(src_dir)/gfa-ed.o \
		$(src_dir)/gfa-io.o $(src_dir)/sketch.o $(src_dir)/gfa-base.o \
		$(src_dir)/kthread.o $(src_dir)/options.o $(src_dir)/kalloc.o \
		$(src_dir)/misc.o $(src_dir)/gfa-bbl.o \
		$(src_dir)/sys.o $(src_dir)/ILP_index.o

AlphaASM: $(OBJS)
	$(CXX) $^ -o $@ $(INLCLUDES) $(LIBS) $(CXXFLAGS)

$(src_dir)/%.o: $(src_dir)/%.cpp
	$(CXX) -c $< -o $@ $(INLCLUDES) $(LIBS) $(CXXFLAGS)

clean:
	rm -f $(src_dir)/*.o AlphaASM