#include "ILP_index.h"

#define MG_M_NO_DIAG      0x400000

// Constructor
ILP_index::ILP_index(gfa_t *g) {
    this->g = g;
}

KSEQ_INIT(gzFile, gzread)

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASHL_MAP_INIT(KH_LOCAL, idxhash_t, mg_hidx, uint64_t, uint64_t, idx_hash, idx_eq)

typedef struct mg_idx_bucket_s {
	mg128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mg_idx_bucket_t;

mg_idx_t *mg_idx_init(int k, int w, int b)
{
	mg_idx_t *gi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	KCALLOC(0, gi, 1);
	gi->w = w, gi->k = k, gi->b = b;
	KCALLOC(0, gi->B, 1<<b);
	return gi;
}

void mg_idx_destroy(mg_idx_t *gi)
{
	uint32_t i;
	if (gi == 0) return;
	if (gi->B) {
		for (i = 0; i < 1U<<gi->b; ++i) {
			free(gi->B[i].p);
			free(gi->B[i].a.a);
			mg_hidx_destroy((idxhash_t*)gi->B[i].h);
		}
		free(gi->B);
	}
	gfa_edseq_destroy(gi->n_seg, gi->es);
	free(gi);
}

/****************
 * Index access *
 ****************/

const uint64_t *mg_idx_hget(const void *h_, const uint64_t *q, int suflen, uint64_t minier, int *n)
{
	khint_t k;
	const idxhash_t *h = (const idxhash_t*)h_;
	*n = 0;
	if (h == 0) return 0;
	k = mg_hidx_get(h, minier>>suflen<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &q[kh_val(h, k)>>32];
	}
}

const uint64_t *mg_idx_get(const mg_idx_t *gi, uint64_t minier, int *n)
{
	int mask = (1<<gi->b) - 1;
	mg_idx_bucket_t *b = &gi->B[minier&mask];
	return mg_idx_hget(b->h, b->p, gi->b, minier, n);
}

void mg_idx_cal_quantile(const mg_idx_t *gi, int32_t m, float f[], int32_t q[])
{
	int32_t i;
	uint64_t n = 0;
	khint_t *a, k;
	for (i = 0; i < 1<<gi->b; ++i)
		if (gi->B[i].h) n += kh_size((idxhash_t*)gi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = 0, n = 0; i < 1<<gi->b; ++i) {
		idxhash_t *h = (idxhash_t*)gi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	for (i = 0; i < m; ++i)
		q[i] = ks_ksmall_uint32_t(n, a, (size_t)((1.0 - (double)f[i]) * n));
	free(a);
}

/***************
 * Index build *
 ***************/

static void mg_idx_add(mg_idx_t *gi, int n, const mg128_t *a)
{
	int i, mask = (1<<gi->b) - 1;
	for (i = 0; i < n; ++i) {
		mg128_v *p = &gi->B[a[i].x>>8&mask].a;
		kv_push(mg128_t, 0, *p, a[i]);
	}
}

void mg_idx_hfree(void *h_)
{
	idxhash_t *h = (idxhash_t*)h_;
	if (h == 0) return;
	mg_hidx_destroy(h);
}

void *mg_idx_a2h(void *km, int32_t n_a, mg128_t *a, int suflen, uint64_t **q_, int32_t *n_)
{
	int32_t N, n, n_keys;
	int32_t j, start_a, start_q;
	idxhash_t *h;
	uint64_t *q;

	*q_ = 0, *n_ = 0;
	if (n_a == 0) return 0;

	// sort by minimizer
	radix_sort_128x(a, a + n_a);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, N = 0; j <= n_a; ++j) {
		if (j == n_a || a[j].x>>8 != a[j-1].x>>8) {
			++n_keys;
			if (n > 1) N += n;
			n = 1;
		} else ++n;
	}
	h = mg_hidx_init2(km);
	mg_hidx_resize(h, n_keys);
	KCALLOC(km, q, N);
	*q_ = q, *n_ = N;

	// create the hash table
	for (j = 1, n = 1, start_a = start_q = 0; j <= n_a; ++j) {
		if (j == n_a || a[j].x>>8 != a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mg128_t *p = &a[j-1];
			itr = mg_hidx_put(h, p->x>>8>>suflen<<1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					q[start_q + k] = a[start_a + k].y;
				radix_sort_gfa64(&q[start_q], &q[start_q + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_q<<32 | n;
				start_q += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	assert(N == start_q);
	return h;
}

static void worker_post(void *g, long i, int tid)
{
	mg_idx_t *gi = (mg_idx_t*)g;
	mg_idx_bucket_t *b = &gi->B[i];
	if (b->a.n == 0) return;
	b->h = (idxhash_t*)mg_idx_a2h(0, b->a.n, b->a.a, gi->b, &b->p, &b->n);
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}

int mg_gfa_overlap(const gfa_t *g)
{
	int64_t i;
	for (i = 0; i < g->n_arc; ++i) // non-zero overlap
		if (g->arc[i].ov != 0 || g->arc[i].ow != 0)
			return 1;
	return 0;
}

// For finding anchors
struct mg_tbuf_s {
	void *km;
	int frag_gap;
};

mg_tbuf_t *mg_tbuf_init(void)
{
	mg_tbuf_t *b;
	b = (mg_tbuf_t*)calloc(1, sizeof(mg_tbuf_t));
	if (true) b->km = km_init();
	return b;
}

void mg_tbuf_destroy(mg_tbuf_t *b)
{
	if (b == 0) return;
	if (b->km) km_destroy(b->km);
	free(b);
}

void *mg_tbuf_get_km(mg_tbuf_t *b)
{
	return b->km;
}

static void collect_minimizers(void *km, const mg_mapopt_t *opt, const mg_idx_t *gi, int n_segs, const int *qlens, const char **seqs, mg128_v *mv)
{
	int i, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		mg_sketch(km, seqs[i], qlens[i], gi->w, gi->k, i, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mg128_t, heap_lt)

typedef struct {
	uint32_t n;
	uint32_t q_pos, q_span;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} mg_match_t;

static mg_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mg_idx_t *gi, const mg128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, int32_t **mini_pos)
{
	int rep_st = 0, rep_en = 0, n_m;
	size_t i;
	mg_match_t *m;
	*n_mini_pos = 0;
	KMALLOC(km, *mini_pos, mv->n);
	m = (mg_match_t*)kmalloc(km, mv->n * sizeof(mg_match_t));
	for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
		const uint64_t *cr;
		mg128_t *p = &mv->a[i];
		uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
		int t;
		cr = mg_idx_get(gi, p->x>>8, &t);
		if (t >= max_occ) {
			int en = (q_pos >> 1) + 1, st = en - q_span;
			if (st > rep_en) {
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			mg_match_t *q = &m[n_m++];
			q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
			q->is_tandem = 0;
			if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
			if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
			*n_a += q->n;
			(*mini_pos)[(*n_mini_pos)++] = q_pos>>1;
		}
	}
	*rep_len += rep_en - rep_st;
	*_n_m = n_m;
	return m;
}

static mg128_t *collect_seed_hits(void *km, const mg_mapopt_t *opt, int max_occ, const mg_idx_t *gi, const char *qname, const mg128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, int32_t **mini_pos)
{
	int i, n_m;
	mg_match_t *m;
	mg128_t *a;
	m = collect_matches(km, &n_m, max_occ, gi, mv, n_a, rep_len, n_mini_pos, mini_pos);
    a = (mg128_t*)malloc(*n_a * sizeof(mg128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mg_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		uint32_t k;
		for (k = 0; k < q->n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mg128_t *p;
			if (qname && (opt->flag & MG_M_NO_DIAG)) {
				const gfa_seg_t *s = &gi->g->seg[r[k]>>32];
				const char *gname = s->snid >= 0 && gi->g->sseq? gi->g->sseq[s->snid].name : s->name;
				int32_t g_pos;
				if (s->snid >= 0 && gi->g->sseq)
					gname = gi->g->sseq[s->snid].name, g_pos = s->soff + (uint32_t)r[k];
				else
					gname = s->name, g_pos = (uint32_t)r[k];
				if (g_pos == q->q_pos && strcmp(qname, gname) == 0)
					continue;
			}
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) // forward strand
				p->x = r[k]>>32<<33 | rpos;
			else // reverse strand
				p->x = r[k]>>32<<33 | 1ULL<<32 | (gi->g->seg[r[k]>>32].len - (rpos + 1 - q->q_span) - 1);
			p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			p->y |= (uint64_t)q->seg_id << MG_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MG_SEED_TANDEM;
			p->y |= (uint64_t)(q->n < 255? q->n : 255) << MG_SEED_OCC_SHIFT;
		}
	}
	kfree(km, m);
	radix_sort_128x(a, a + (*n_a));
	return a;
}

std::vector<mg128_t> find_anchors(mg_idx_t *gfa_idx, std::string read_seq, int32_t max_occ)
{
    mg_mapopt_t opt;
    mg_idxopt_t ipt;
    mg_opt_set(0, &ipt, &opt);
    mg_tbuf_t *b = mg_tbuf_init();
    void *km = mg_tbuf_get_km(b);
    int i, n, rep_len, n_mini_pos;
    int32_t *mini_pos;
    mg128_t *a;
    mg128_v mv;
    int64_t n_a;
    mv.n = 0, mv.m = 0, mv.a = 0;
    int32_t *read_len = (int32_t*)malloc(sizeof(int32_t));
    read_len[0] = read_seq.size();
    const char **seqs = (const char**)malloc(sizeof(char*));
    seqs[0] = read_seq.c_str();
    collect_minimizers(km, &opt, gfa_idx, 1, read_len, seqs, &mv);
    a = collect_seed_hits(km, &opt, max_occ, gfa_idx, 0, &mv, read_seq.size(), &n_a, &rep_len, &n_mini_pos, &mini_pos);
    
    // Free Memory
    mg_tbuf_destroy(b);
    free(read_len);
    free(seqs);

    std::vector<mg128_t> a_;
    for (int i = 0; i < n_a; i++)
    {
        a_.push_back(a[i]);
    }
    free(a);

    // return the anchors
    return a_;
}

void::ILP_index::read_gfa() 
{
    uint v;
    n_vtx = gfa_n_vtx(g);
    // Resize node_len
    node_len.resize(n_vtx/2, 0);
    node_seq.resize(n_vtx/2, "");
    std::vector<std::vector<uint32_t>> adj_;
    adj_.resize(n_vtx); 
    /* Node_len */
    for (int v = 0; v < n_vtx/2; v++)
    {
        gfa_seg_t segment = (g)->seg[v];
        int len =  segment.len;
        node_len[v] = len;
        node_seq[v] = segment.seq;
    }
    // look for all the edges , if the sum of all the 
    // edges are zero then, that's a linear reference
    u_int num_edges = 0;
    for (v = 0; v < n_vtx; v++)
    {
        gfa_arc_t *av = gfa_arc_a(g, v);
        int n_edges = gfa_arc_n(g, v);
        int v_ = av->v_lv >> 32;
        // std::cerr << " node : " << v << " node_len : " << node_len[v] << std::std::endl; 
        for (int i = 0; i < n_edges; i++)
        {
            num_edges++;
        }
    }


    for (v = 0; v < n_vtx; v++)
    {
        gfa_arc_t *av = gfa_arc_a(g, v);
        int n_edges = gfa_arc_n(g, v);
        if (num_edges == 0) // Linear Reference
        {
            lin_ref = 1; // Mark as a linear reference
        }else
        {
            int v_ = av->v_lv >> 32;
            // std::cerr << " node : " << v << " node_len : " << node_len[v] << std::std::endl; 
            for (int i = 0; i < n_edges; i++)
            {
                uint w = av[i].w;
                adj_[v_].push_back(w);
            }
        }
    }

    // Copy the adjacency list only for forward strand
    n_vtx /= 2;
    adj_list.resize(n_vtx);
    for (int i = 0; i < 2 * n_vtx; i++)
    {
        if (i % 2 == 0) // Only forward strand
        {
            for (auto v : adj_[i])
            {
                adj_list[i/2].push_back(v/2);
            }
        }
    }

    adj_.clear();

    // Read walks
    if (g->n_walk > 0) is_gfa_v12 = true;
    num_walks = g->n_walk; // count number of walks
    haps.resize(n_vtx);
    paths.resize(num_walks);

    // Fill the paths
    for (size_t w = 0; w < g->n_walk; w++)
    {
        int32_t idx = 0;
        for (size_t n = 0; n < g->walk[w].n_v; n++)
        {
            int v = g->walk[w].v[n];
            if (v%2 != 0) {
                fprintf(stderr, "Error: Walk %d has reverse strand vertices %d\n", w, v);
                exit(1);
            }
			v /= 2; // for forward strand
            haps[v].push_back(w); // Map forward walk to haplotype
            paths[w].push_back(v); // Map forward walk to path
        }
    }
}

void printConstraints(GRBModel& model) {
    GRBConstr* constraints = model.getConstrs();
    int numConstraints = model.get(GRB_IntAttr_NumConstrs);

    for (int i = 0; i < numConstraints; ++i) {
        std::string constrName = constraints[i].get(GRB_StringAttr_ConstrName);
        GRBLinExpr constrExpr = model.getRow(constraints[i]);
        double rhs = constraints[i].get(GRB_DoubleAttr_RHS);
        char sense = constraints[i].get(GRB_CharAttr_Sense);

        std::cout << "Constraint " << constrName << ": ";

        for (int j = 0; j < constrExpr.size(); ++j) {
            GRBVar var = constrExpr.getVar(j);
            double coeff = constrExpr.getCoeff(j);
            std::string varName = var.get(GRB_StringAttr_VarName);

            std::cout << coeff << " * " << varName;
            if (j < constrExpr.size() - 1) {
                std::cout << " + ";
            }
        }

        if (sense == GRB_EQUAL) {
            std::cout << " == ";
        } else if (sense == GRB_LESS_EQUAL) {
            std::cout << " <= ";
        } else if (sense == GRB_GREATER_EQUAL) {
            std::cout << " >= ";
        }

        std::cout << rhs << std::endl;
    }

    delete[] constraints;
}

void printNonZeroVariables(GRBModel& model) {
    int numVars = model.get(GRB_IntAttr_NumVars);
    GRBVar* vars = model.getVars();

    std::cout << "Non-zero variables:\n";
    for (int i = 0; i < numVars; ++i) {
        double val = vars[i].get(GRB_DoubleAttr_X);
        if (val != 0.0) {
            std::string varName = vars[i].get(GRB_StringAttr_VarName);
            std::cout << varName << " = " << val << std::endl;
        }
    }

    delete[] vars;
}

// Read the reads
void ILP_index::read_ip_reads(std::vector<std::pair<std::string, std::string>> &ip_reads, std::string ip_reads_file)
{
    // Read the reads with gzip
    gzFile fp;
    kseq_t *seq;
    int32_t l;

    fp = gzopen(ip_reads_file.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        ip_reads.push_back(std::make_pair(seq->name.s, seq->seq.s));
    }

    kseq_destroy(seq);
    gzclose(fp);
}

void ILP_index::ILP_function(std::vector<std::pair<std::string, std::string>> &ip_reads)
{
    /* 
    * 1. Graph -> adj_list, node_len, node_seq, paths, haps
    * 2. Reads -> ip_reads
    * 3. Haplotype -> hap_file, hap_name [This is id for the haplotype]
    */

   // Print the graph stats
   fprintf(stderr, "[M::%s::%.3f*%.2f] Graph has %d vertices, %d walks and read has %d reads\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), n_vtx, num_walks, ip_reads.size());

   /*
        1) Get the haplotypes as a sequence.
        2) Get the reads as a sequence. 
   */
    std::vector<int32_t> hap_sizes(num_walks, 0);
    mg128_v a = {0,0,0};
    mg_idx_t *gi = mg_idx_init(k_mer, window, bucket_bits);
    for (size_t h = 0; h < num_walks; h++)
    {
        std::string haplotype = "";
        for (size_t i = 0; i < paths[h].size(); i++)
        {
            haplotype += node_seq[paths[h][i]];
        }
        hap_sizes[h] = haplotype.size();
        a.n = 0;
        mg_sketch(0, haplotype.c_str(), haplotype.size(), k_mer, window, h, &a);
        mg_idx_add(gi, a.n, a.a);
    }
    free(a.a);
	kt_for(num_threads, worker_post, gi, 1<<gi->b);

    // sketching haplotypes
    fprintf(stderr, "[M::%s::%.3f*%.2f] Haplotypes sketched\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

    int32_t num_reads = ip_reads.size();
    std::vector<std::vector<mg128_t>> read_hits(num_reads);

    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_reads; i++)
    {
        read_hits[i] = find_anchors(gi, ip_reads[i].second, max_occ);
        // fprintf(stderr, "Read : %d, Anchors : %d\n", i, read_hits[i].size());
    }

    std::vector<std::set<std::vector<int32_t>>> unique_kmers(num_walks);
    for (int i = 0; i < num_reads; i++)
    {
        for (int j = 0; j < read_hits[i].size(); j++)
        {
            int32_t hap = read_hits[i][j].x >> 32;
            std::vector<int32_t> k_mers;
            if (hap % 2 == 0) // FWD strand
            {
                hap = hap/2;
                int32_t end_x = (int32_t)read_hits[i][j].x;
                int32_t start_x = end_x - k_mer + 1;
                // std::cout << "Hap: " << hap <<" Start : " << start_x << " End : " << end_x << std::endl;
                for (int32_t k = start_x; k <= end_x; k++)
                {
                    k_mers.push_back(k);
                }
            }else // REV strand
            {
                hap = (hap - 1)/2; // REV strand -> FWD strand
                int32_t end_x = hap_sizes[hap] - 2 - (int32_t)read_hits[i][j].x + k_mer;
                int32_t start_x = end_x - k_mer + 1;
                // std::cout << "Hap: " << hap <<" Start : " << start_x << " End : " << end_x << std::endl;
            }
            unique_kmers[hap].insert(k_mers);
        }
        read_hits[i].clear();
    }

    // For every path print the unique kmers
    int32_t num_kmers = 0;
    for (int i = 0; i < num_walks; i++)
    {
        num_kmers += unique_kmers[i].size();
        // fprintf(stderr, "Hap : %d, Kmers : %d\n", i, unique_kmers[i].size());
    }
    std::vector<std::vector<std::vector<int32_t>>> Kmers(num_walks);

    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_walks; i++)
    {
        for (auto it = unique_kmers[i].begin(); it != unique_kmers[i].end(); it++)
        {
            for (int j = 0; j < it->size(); j++)
            {
                Kmers[i].push_back((*it));
            }
        }
    }
    unique_kmers.clear();
    kfree(0, gi->B);
    free(gi);

    fprintf(stderr, "[M::%s::%.3f*%.2f] %d K-mers found\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), num_kmers);

    std::vector<std::vector<int32_t>> in_nodes(n_vtx);
    for (int32_t i = 0; i < n_vtx; i++)
    {
        for (auto v : adj_list[i])
        {
            in_nodes[v].push_back(i);
        }
    }

    // Write an ILP with Gurobi
    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_Threads, num_threads);
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        std::map<std::string, GRBVar> vars;

        // Kmer constraints
        for (int32_t i = 0; i < num_walks; i++) {
            for (int32_t j = 0; j < Kmers[i].size(); j++) {
                GRBLinExpr kmer_expr;
                for (int32_t k = 1; k < Kmers[i][j].size(); k++) {
                    int32_t u = Kmers[i][j][k-1];
                    int32_t v = Kmers[i][j][k];
                    std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(i);
                    vars[var_name] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                    kmer_expr += vars[var_name];
                }
                std::string exra_var = "z_" + std::to_string(i) + "_" + std::to_string(j);
                GRBVar kmer_expr_var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, exra_var);
                vars[exra_var] = kmer_expr_var;
                std::string constraint_name = "Kmer_constraints_" + std::to_string(i) + "_" + std::to_string(j);
                int32_t kmer_weight = Kmers[i][j].size() - 1;
                model.addConstr(kmer_expr == kmer_weight * kmer_expr_var, constraint_name);
            }
        }

        fprintf(stderr, "[M::%s::%.3f*%.2f] Kmer constraints added to the model\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        // Create the objective function
        GRBLinExpr obj;

        GRBLinExpr start_expr;
        GRBLinExpr end_expr;
        for (int32_t i = 0; i < num_walks; i++) {
            int32_t u_start = paths[i][0];
            std::string var_name_start = "s_" + std::to_string(u_start) + "_" + std::to_string(i);
            GRBVar var_start = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name_start);
            vars[var_name_start] = var_start;
            start_expr += var_start;

            int32_t u_end = paths[i][paths[i].size() - 1];
            std::string var_name_end = std::to_string(u_end) + "_" + std::to_string(i) + "_e";
            GRBVar var_end = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name_end);
            vars[var_name_end] = var_end;
            end_expr += var_end;
        }

        // set start_expr <= 1 and end_expr <= 1
        model.addConstr(start_expr == 1, "Start_expr");
        model.addConstr(end_expr == 1, "End_expr");

        // Add vertex constraints from adj_list
        GRBLinExpr vtx_expr;
        for (auto u = 0; u < n_vtx; u++) {
            for (auto v : adj_list[u]) {
                for (auto i: haps[u]) {
                    for (auto j: haps[v]) {
                        std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(j);
                        if (i != j) {
                            GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                            vars[var_name] = var;
                            vtx_expr += 2 * var;
                        } else {
                            if (vars.find(var_name) == vars.end()) { // Variable does not exist
                                GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                                vars[var_name] = var;
                                vtx_expr += var;
                            }else {
                                vtx_expr += 0 * vars[var_name]; // Already variable exists in Kmers hence no need to add
                            }
                        }
                    }
                }
            }
        }

        obj = start_expr + vtx_expr + end_expr;

        fprintf(stderr, "[M::%s::%.3f*%.2f] Objective function added to the model\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        // Flow constraints for each vertex
        for (auto u = 0; u < n_vtx; u++) {
            if (in_nodes[u].size() == 0 || adj_list[u].size() == 0) continue;

            GRBLinExpr in_expr;
            for (auto v: in_nodes[u]) {
                for (auto i: haps[v]) {
                    for (auto j: haps[u]) {
                        std::string var_name = std::to_string(v) + "_" + std::to_string(i) + "_" + std::to_string(u) + "_" + std::to_string(j);
                        in_expr += vars[var_name];
                    }
                }
            }

            GRBLinExpr out_expr;
            for (auto v: adj_list[u]) {
                for (auto i: haps[u]) {
                    for (auto j: haps[v]) {
                        std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(j);
                        out_expr += vars[var_name];
                    }
                }
            }

            std::string constraint_name = "Flow_conservation_" + std::to_string(u);
            model.addConstr(in_expr == out_expr, constraint_name);
        }

        // Flow constraints for source nodes
        for (int32_t i = 0; i < num_walks; i++) {
            int32_t u = paths[i][0];
            GRBLinExpr s_expr;
            s_expr += vars["s_" + std::to_string(u) + "_" + std::to_string(i)];
            for (auto v: adj_list[u]) {
                for (auto j: haps[v]) {
                    s_expr -= vars[std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(j)];
                }
            }
            std::string constraint_name = "Source_conservation_" + std::to_string(u) + "_" + std::to_string(i);
            model.addConstr(s_expr == 0, constraint_name);
        }

        // Flow constraints for sink nodes
        for (int32_t i = 0; i < num_walks; i++) {
            int32_t u = paths[i][paths[i].size() - 1];
            GRBLinExpr e_expr;
            for (auto v: in_nodes[u]) {
                for (auto j: haps[v]) {
                    e_expr += vars[std::to_string(v) + "_" + std::to_string(j) + "_" + std::to_string(u) + "_" + std::to_string(i)];
                }
            }
            e_expr += -1 * vars[std::to_string(u) + "_" + std::to_string(i) + "_e"];
            std::string constraint_name = "Sink_conservation_" + std::to_string(u) + "_" + std::to_string(i);
            model.addConstr(e_expr == 0, constraint_name);
        }

        fprintf(stderr, "[M::%s::%.3f*%.2f] Flow conservation constraints added to the model\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        model.setObjective(obj, GRB_MINIMIZE);
        // Optimize model
        model.optimize();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Model optimized\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        // Print constraints
        bool debug = false;
        if (debug)
        {
            printConstraints(model);
            printNonZeroVariables(model);
        }

    } catch (GRBException e) {
        std::cerr << "Error code = " << e.getErrorCode() << std::endl;
        std::cerr << e.getMessage() << std::endl;
    } catch (...) {
        std::cerr << "Exception during optimization" << std::endl;
    }

    
    // write haplotype as to a file as fasta from the path
    std::string path_str = "ATCG"; // replace this with the actual string spelled by the path
    std::ofstream hap_file_stream(hap_file, std::ios::out);
    hap_file_stream << ">" << hap_name << " LN:" << path_str.size() << std::endl;
    // write the path_str to the file 80 characters per line
    for (size_t i = 0; i < path_str.size(); i += 80) {
        hap_file_stream << path_str.substr(i, 80) << std::endl;
    }
    hap_file_stream.close();

    fprintf(stderr, "[M::%s::%.3f*%.2f] Haplotype written to: %s\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), hap_file.c_str());
}




// try {
//     // Create an environment
//     GRBEnv env = GRBEnv(true);
//     env.set("LogFile", "gurobi.log");
//     env.start();

//     // Create an empty model
//     GRBModel model = GRBModel(env);

//     // Create binary variables
//     GRBVar su1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "su1");
//     GRBVar u1v1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "u1v1");
//     GRBVar su2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "su2");
//     GRBVar v1x1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "v1x1");
//     GRBVar w2x2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "w2x2");
//     GRBVar u2w2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "u2w2");
//     GRBVar u1w2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "u1w2");
//     GRBVar v1x2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "v1x2");
//     GRBVar w2x1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "w2x1");
//     GRBVar x1z2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x1z2");
//     GRBVar x2y1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x2y1");
//     GRBVar x1y1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x1y1");
//     GRBVar x2z2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x2z2");
//     GRBVar y1t = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y1t");
//     GRBVar z2t = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z2t");
//     GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");


//     // Set objective: maximize -2x + 3y + 4z
//     GRBLinExpr objective = su1 + 0 * u1v1 + su2 + u2w2 + 2* u1w2 + 2*v1x2 + 2*w2x1 + 2*x1z2 + 2*x2y1 + x1y1 + x2z2 + w2x2 + 0 * v1x1 + y1t + z2t;
//     model.setObjective(objective, GRB_MINIMIZE);

//     // Add constraints to enforce x + y ∈ {0, 2}
//     model.addConstr(u1v1 + v1x1 <= 2 * z, "c0");  // If a == 0, then x + y == 0
//     model.addConstr(u1v1 + v1x1 >= 2 * z, "c1");  // If a == 1, then x + y == 2

//     model.addConstr(su1 + su2 == 1, "c2");  // flow conservation at u1
//     model.addConstr(y1t + z2t == 1, "c3");  // flow conservation at t
//     model.addConstr(su1 - u1v1 - u2w2 == 0, "c4");  // flow conservation at u2
//     model.addConstr(su2 - u2w2 == 0, "c5");  // flow conservation at w2
//     model.addConstr(u1v1 - v1x1 - v1x2 == 0, "c6");  // flow conservation at v1
//     model.addConstr(u1w2 + u2w2 - w2x1 - w2x2 == 0, "c7");  // flow conservation at w2
//     model.addConstr(v1x1 + w2x1 - x1y1 - x1z2 == 0, "c8");  // flow conservation at x1
//     model.addConstr(v1x2 + w2x2 - x2y1 - x2z2 == 0, "c9");  // flow conservation at x2
//     model.addConstr(x1y1 + x2y1 - y1t == 0, "c10");  // flow conservation at y1
//     model.addConstr(x2z2 + x1z2 - z2t == 0, "c11");  // flow conservation at z2

//     // Optimize model
//     model.optimize();
        
//     // print the objective value and path i.e. non zero variables
//     std::cout << "Objective value: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
//     for (int i = 0; i < model.get(GRB_IntAttr_NumVars); i++) {
//         GRBVar var = model.getVar(i);
//         if (var.get(GRB_DoubleAttr_X) != 0) {
//             std::cout << var.get(GRB_StringAttr_VarName) << " = " << var.get(GRB_DoubleAttr_X) << std::endl;
//         }
//     }

// } catch (GRBException e) {
//     std::cout << "Error code = " << e.getErrorCode() << std::endl;
//     std::cout << e.getMessage() << std::endl;
// } catch (...) {
//     std::cout << "Exception during optimization" << std::endl;
// }