#include <iostream>
#include <jemalloc/jemalloc.h>
#include <omp.h>
#include "AApriv.h"
#include <assert.h>
#include <math.h>
#include "khashl.h"
#include "kalloc.h"
#include <zlib.h>
#include "kseq.h"
#include "kvec-km.h"
#include "kthread.h"
#include "sys.h"

// CPP
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <future>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <queue>
#include <stack>
#include <functional>
#include <tuple>
#include <map>

// Gurobi
#include "gurobi_c++.h"


class ILP_index {
    public:

        // Data Structures
        gfa_t *g; // This is graph
        std::vector<std::vector<uint32_t>> adj_list;
        std::vector<uint32_t> node_len;
        std::vector<std::string> node_seq;
        std::vector<std::vector<uint32_t>> paths;
        std::vector<std::vector<uint32_t>> haps;

        // Support Variables
        int32_t lin_ref = 0;
        uint32_t n_vtx = 0;
        bool is_gfa_v12 = false;
        uint32_t num_walks = 0;
        int32_t num_bubbles = 0;
        int32_t num_threads = 4;
        std::string hap_file = "";
        std::string hap_name = "";
        bool debug = false;
        int32_t k_mer;
        int32_t window;
        int32_t bucket_bits;
        int32_t max_occ;

        // Constructor
        ILP_index(gfa_t *g);	// This is constructor

        // Functions
        void read_gfa();
        void read_ip_reads(std::vector<std::pair<std::string, std::string>> &ip_reads, std::string ip_reads_file);
        void ILP_function(std::vector<std::pair<std::string, std::string>> &ip_reads);

        /* Please add your ILP functions here */
};