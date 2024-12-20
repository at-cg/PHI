#include <iostream>
#include <stdlib.h>
#include "gfa-priv.h"
#include "PHIpriv.h"
#include "ketopt.h"
#include <zlib.h>
#include "ILP_index.h"
#include <map>
#include <set>
#include "sys.h"
#include <omp.h>

// CPP
#include <vector>
#include <unordered_map>


static ko_longopt_t long_options[] = {
    { const_cast<char*>("version"), ko_no_argument, 300 },
    { 0, 0, 0 }
};

int main(int argc, char *argv[]) {
    // Increase stack size
    struct rlimit rl;
    // Set the new stack size limit
    rl.rlim_cur = std::pow(10, 10);  // 10 GB
    if (rl.rlim_cur > rl.rlim_max) {
        rl.rlim_max = rl.rlim_cur;  // Ensure the hard limit is at least as large as the soft limit
    }

    // check the current status of the stack size
    if (setrlimit(RLIMIT_STACK, &rl) == -1) {
        perror("setrlimit");
        exit(EXIT_FAILURE);
    }

    const char *opt_str = "x:d:c:l:s:m:R:q:T:N:m:h:k:w:t:g:r:o:DSc";
    ketopt_t o = KETOPT_INIT;
	mg_mapopt_t opt;
	mg_idxopt_t ipt;
    bool debug = false;
    int32_t recombination = 100;
    int32_t is_qclp = 1;
    int32_t is_naive_exp = 0;
    float threshold = 1.0f;
    bool is_mixed = true;

    int i, c, ret;
	FILE *fp_help = stderr;
    int32_t help = 0;
    mg_verbose = 3;
    int32_t max_occ = 5000;

	mg_opt_set(0, &ipt, &opt);
	o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'w') ipt.w = atoi(o.arg);
		else if (c == 'k') ipt.k = atoi(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
        else if (c == 'm') is_mixed = atoi(o.arg);
        else if (c == 'g') opt.gfa_file = o.arg;
        else if (c == 'R') recombination = atoi(o.arg);
        else if (c == 'q') is_qclp = atoi(o.arg);
        else if (c == 'N') is_naive_exp = atoi(o.arg);
        else if (c == 'T') threshold = atof(o.arg);
        else if (c == 'r') opt.reads_file = o.arg;
        else if (c == 'o') opt.hap_file = o.arg;
        else if (c == 'c') max_occ = atoi(o.arg);
        else if (c == 'd') debug = atoi(o.arg);
        else if (c == 'h') help = 1;
        else if (c == 300) {
            fprintf(fp_help, "PHI version: %s\n", PHI_VERSION);
            return 0;
        }
	}

	if (argv[1] == NULL || opt.gfa_file == "" || opt.reads_file == "" || opt.hap_file == "" || help == 1 || fp_help == stdout) {
		fprintf(fp_help, "Usage: PHI -g <target.gfa> -r <reads.fa> -o <haplotype.fasta> \n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "    -k INT       K-mer size [%d]\n", ipt.k);
		fprintf(fp_help, "    -w INT       Minimizer window size [%d]\n", ipt.w);
        fprintf(fp_help, "    -R INT       Recombination penalty [%d]\n", recombination);
        fprintf(fp_help, "    -q INT       Mode QP/ILP (default IQP i.e q1, use q0 for ILP) [%d]\n", is_qclp);
        // fprintf(fp_help, "    -N INT       Mode OPT/Naive expanded graph (default Optimized i.e N0, use N1 for Naive) [%d]\n", is_naive_exp);
        fprintf(fp_help, "    -m INT       Mixed/Interger programming (default Mixed i.e -m1, use -m0 for Integer) [%d]\n", is_mixed);
        fprintf(fp_help, "    -T FLOAT     Threshold for minimizer filtering [%.3f]\n", threshold);
        fprintf(fp_help, "    -t INT       Threads [%d]\n", opt.n_threads);
        fprintf(fp_help, "    -g INT       GFA file [%s]\n", opt.gfa_file.c_str());
        fprintf(fp_help, "    -r INT       Read [%s]\n", opt.reads_file.c_str());
        fprintf(fp_help, "    -o INT       Output haplotype [%s]\n", opt.hap_file.c_str());
        fprintf(fp_help, "    -d bool      Debug mode [%d]\n", debug);
		return fp_help == stdout? 0 : 1;
	};

    // start time 
    mg_realtime0 = realtime();

    std::string reads = opt.reads_file;
    gfa_t *g = gfa_read(opt.gfa_file.c_str());
    if (g == 0) {
        fprintf(stderr, "[E::%s] failed to load the GFA file\n", __func__);
        return 1;
    }
    fprintf(stderr, "[M::%s::%.3f*%.2f] Loaded graph from: %s\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), opt.gfa_file.c_str());

    // Get haplotype name (used as an id for the haplotype)
    std::string hap_name = "";
    get_hap_name((char *)opt.gfa_file.c_str(), (char *)opt.reads_file.c_str(), hap_name);


    // Index the graph
    ILP_index *ILP_handle = new ILP_index(g);
    ILP_handle->read_gfa();

    // Read params
    ILP_handle->num_threads = opt.n_threads; // number of threads
    ILP_handle->hap_file = opt.hap_file; // haplotype file to be written
    ILP_handle->debug = debug; // debug mode
    ILP_handle->hap_name = hap_name; // haplotype name to be written as id of the haplotype
    ILP_handle->k_mer = ipt.k; // k-mer size
    ILP_handle->window = ipt.w; // window size
    ILP_handle->bucket_bits = 14; // bucket bits
    ILP_handle->max_occ = max_occ; // maximum k-mer occurence
    ILP_handle->recombination = recombination; // recombination penalty
    ILP_handle->max_occ = max_occ; // maximum k-mer occurence
    ILP_handle->is_qclp = is_qclp; // mode QCLP/ILP
    ILP_handle->is_naive_exp = is_naive_exp; // naive mode
    ILP_handle->threshold = threshold; // threshold for k-mer filtering
    ILP_handle->is_mixed = is_mixed; // mixed mode



    // Read the reads from "-r" file
    std::vector<std::pair<std::string, std::string>> ip_reads; //ip_reads[idx] = (read_id, sequence)
    ILP_handle->read_ip_reads(ip_reads, reads); // read the reads from the file

    // execute the ILP function
    ILP_handle->ILP_function(ip_reads);


    // Print runtime statistics
    fprintf(stderr, "[M::%s] PHI Version: %s\n", __func__, PHI_VERSION);
    fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - mg_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);

    return 0;
}