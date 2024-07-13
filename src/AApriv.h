#ifndef MGPRIV_H
#define MGPRIV_H

#include <stdlib.h>
#include "AlphaASM.h"

// by  GS
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cassert>
#include <sys/resource.h>
#include <errno.h>
#include <math.h>

bool compare_T(const mg128_t &a, const mg128_t &b);

#define MG_SEED_IGNORE     (1ULL<<41)
#define MG_SEED_TANDEM     (1ULL<<42)
#define MG_SEED_FIXED      (1ULL<<43)

#define MG_MAX_SEG        255
#define MG_SEED_SEG_SHIFT  48
#define MG_SEED_SEG_MASK   (0xffULL<<(MG_SEED_SEG_SHIFT))
#define mg_seg_id(a) ((int32_t)(((a).y&MG_SEED_SEG_MASK) >> MG_SEED_SEG_SHIFT))

#define MG_SEED_OCC_SHIFT  56

#define MG_MAX_SHORT_K  15

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	unsigned l, m;
	char *s;
} kstring_t;
#endif

// shortest path
typedef struct {
	// input
	uint32_t v;
	int32_t target_dist;
	uint32_t target_hash;
	uint32_t meta:30, check_hash:1, inner:1;
	int32_t qlen;
	// output
	uint32_t n_path:31, is_0:1;
	int32_t path_end;
	int32_t dist;
	uint32_t hash;
} mg_path_dst_t;

typedef struct {
	uint32_t v, d;
	int32_t pre;
} mg_pathv_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline float mg_log2(float x) // NB: this doesn't work when x<2
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}

extern unsigned char seq_nt4_table[256];

void mg_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, mg128_v *p);

void *mg_idx_a2h(void *km, int32_t n_a, mg128_t *a, int suflen, uint64_t **q_, int32_t *n_);
const uint64_t *mg_idx_hget(const void *h_, const uint64_t *q, int suflen, uint64_t minier, int *n);
void mg_idx_hfree(void *h_);

const uint64_t *mg_idx_get(const mg_idx_t *gi, uint64_t minier, int *n);
void mg_idx_cal_quantile(const mg_idx_t *gi, int32_t m, float f[], int32_t q[]);
void mg_update_anchors(int32_t n_a, mg128_t *a, int32_t n, const int32_t *mini_pos);

void mg_sprintf_lite(kstring_t *s, const char *fmt, ...);
void mg_sprintf_km(void *km, kstring_t *s, const char *fmt, ...);
void mg_str_write(void *km, kstring_t *s, int32_t len, char *str);
void mg_str_reserve(void *km, kstring_t *s, int32_t len);

void radix_sort_128x(mg128_t *beg, mg128_t *end);
void radix_sort_gfa64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);


// Custom function by GS
void get_hap_name(char *gfa_name, char *reads_name, std::string &hap_name);
std::vector<std::pair<int, int>> find_bridges(std::vector<std::vector<int>> &adj_list);

#ifdef __cplusplus
}
#endif

#endif