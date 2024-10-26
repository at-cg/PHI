#ifndef PHI_H
#define PHI_H

#include <stdint.h>
#include "gfa.h"
#include <string>
#include <vector>

#define PHI_VERSION "1.0"

typedef struct { uint64_t x, y; } mg128_t;
typedef struct { size_t n, m; mg128_t *a; } mg128_v;
typedef struct { int32_t n, m; uint32_t *a; } mg32_v;
typedef struct { int32_t n, m; uint64_t *a; } mg64_v;

typedef struct {
	int w, k;
	int bucket_bits;
} mg_idxopt_t;

typedef struct {
	uint64_t flag;
	int64_t mini_batch_size;
	std::string gfa_file;
	std::string reads_file;
	std::string hap_file;
	int64_t cap_kalloc;
	int32_t n_threads;
} mg_mapopt_t;

typedef struct {
	const gfa_t *g;
	gfa_edseq_t *es;
	int32_t b, w, k, flag, n_seg;
	struct mg_idx_bucket_s *B; // index (hidden)
} mg_idx_t;

typedef struct {
	int32_t len, n_off, *off;
	char *ds;
} mg_ds_t;

typedef struct mg_tbuf_s mg_tbuf_t;
extern int mg_verbose, mg_dbg_flag;
extern double mg_realtime0;

#ifdef __cplusplus
extern "C" {
#endif

// options
int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo);

// index operations
mg_idx_t *mg_index(gfa_t *g, const mg_idxopt_t *io, int n_threads, mg_mapopt_t *mo, std::vector<int32_t> &bubble_vertices); // combine mg_index_core() and mg_opt_update()
void mg_idx_destroy(mg_idx_t *gi);

#ifdef __cplusplus
}
#endif

#endif