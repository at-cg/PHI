#include "ILP_index.h"

// Constructor
ILP_index::ILP_index(gfa_t *g) {
    this->g = g;
}

KSEQ_INIT(gzFile, gzread)

uint64_t fnv1a_hash_64(const std::string& str) {
    // FNV-1a 64-bit parameters
    const uint64_t FNV_prime = 0x100000001b3;
    const uint64_t offset_basis = 0xcbf29ce484222325;

    // Initialize hash with offset basis
    uint64_t hash = offset_basis;

    // Process each byte of the string
    for (char c : str) {
        hash ^= static_cast<uint64_t>(c);
        hash *= FNV_prime;
    }

    return hash;
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
    in_paths.resize(num_walks, std::vector<int32_t>(n_vtx, 0));

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
            in_paths[w][v] = 1;
        }
    }

    // compute tpological order with kahns algorithm
    std::vector<int32_t> in_degree(n_vtx, 0);
    for (int32_t i = 0; i < n_vtx; i++)
    {
        for (auto v: adj_list[i])
        {
            in_degree[v]++;
        }
    }

    std::queue<int32_t> q;
    for (int32_t i = 0; i < n_vtx; i++)
    {
        if (in_degree[i] == 0)
        {
            q.push(i);
        }
    }

    while (!q.empty())
    {
        int32_t u = q.front();
        q.pop();
        top_order.push_back(u);
        for (auto v: adj_list[u])
        {
            in_degree[v]--;
            if (in_degree[v] == 0)
            {
                q.push(v);
            }
        }
    }

    // create map for topological order
    top_order_map.resize(n_vtx);
    for (int32_t i = 0; i < top_order.size(); i++)
    {
        top_order_map[top_order[i]] = i;
    }
}

void printQuadraticConstraints(GRBModel& model) {
    GRBQConstr* qconstraints = model.getQConstrs();
    int numQConstraints = model.get(GRB_IntAttr_NumQConstrs);

    for (int i = 0; i < numQConstraints; ++i) {
        std::string qconstrName = qconstraints[i].get(GRB_StringAttr_QCName);
        GRBQuadExpr qconstrExpr = model.getQCRow(qconstraints[i]);
        double qrhs = qconstraints[i].get(GRB_DoubleAttr_QCRHS);
        char qsense = qconstraints[i].get(GRB_CharAttr_QCSense);

        std::cout << "Quadratic Constraint " << qconstrName << ": ";

        // Print linear terms
        for (int j = 0; j < qconstrExpr.getLinExpr().size(); ++j) {
            GRBVar var = qconstrExpr.getLinExpr().getVar(j);
            double coeff = qconstrExpr.getLinExpr().getCoeff(j);
            std::string varName = var.get(GRB_StringAttr_VarName);

            std::cout << coeff << " * " << varName;
            if (j < qconstrExpr.getLinExpr().size() - 1 || qconstrExpr.size() > 0) {
                std::cout << " + ";
            }
        }

        // Print quadratic terms
        for (int j = 0; j < qconstrExpr.size(); ++j) {
            GRBVar var1 = qconstrExpr.getVar1(j);
            GRBVar var2 = qconstrExpr.getVar2(j);
            double qcoeff = qconstrExpr.getCoeff(j);
            std::string varName1 = var1.get(GRB_StringAttr_VarName);
            std::string varName2 = var2.get(GRB_StringAttr_VarName);

            std::cout << qcoeff << " * " << varName1 << " * " << varName2;
            if (j < qconstrExpr.size() - 1) {
                std::cout << " + ";
            }
        }

        if (qsense == GRB_EQUAL) {
            std::cout << " == ";
        } else if (qsense == GRB_LESS_EQUAL) {
            std::cout << " <= ";
        } else if (qsense == GRB_GREATER_EQUAL) {
            std::cout << " >= ";
        }

        std::cout << qrhs << std::endl;
    }

    delete[] qconstraints;
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

void printObjectiveFunction(GRBModel& model) {
    GRBQuadExpr objExpr = model.getObjective();
    int sense = model.get(GRB_IntAttr_ModelSense);

    std::cout << "Objective Function: ";

    // Print linear terms
    bool firstTerm = true; // To handle the '+' sign placement
    for (int i = 0; i < objExpr.getLinExpr().size(); ++i) {
        GRBVar var = objExpr.getLinExpr().getVar(i);
        double coeff = objExpr.getLinExpr().getCoeff(i);
        std::string varName = var.get(GRB_StringAttr_VarName);

        if (!firstTerm) {
            std::cout << " + ";
        }
        std::cout << coeff << " * " << varName;
        firstTerm = false;
    }

    // Print quadratic terms
    for (int i = 0; i < objExpr.size(); ++i) {
        GRBVar var1 = objExpr.getVar1(i);
        GRBVar var2 = objExpr.getVar2(i);
        double qcoeff = objExpr.getCoeff(i);
        std::string varName1 = var1.get(GRB_StringAttr_VarName);
        std::string varName2 = var2.get(GRB_StringAttr_VarName);

        if (!firstTerm) {
            std::cout << " + ";
        }
        std::cout << qcoeff << " * " << varName1 << " * " << varName2;
        firstTerm = false;
    }

    // Print constant term if it exists
    double constant = model.get(GRB_DoubleAttr_ObjCon);
    if (constant != 0) {
        if (!firstTerm) {
            std::cout << " + ";
        }
        std::cout << constant;
    }

    if (sense == GRB_MINIMIZE) {
        std::cout << " (Minimize)" << std::endl;
    } else {
        std::cout << " (Maximize)" << std::endl;
    }
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

std::string reverse_strand(std::string seq)
{
    std::string rev_seq = "";
    for (int i = seq.size() - 1; i >= 0; i--)
    {
        if (seq[i] == 'A' || seq[i] == 'a')
        {
            rev_seq += 'T';
        }
        else if (seq[i] == 'T' || seq[i] == 't')
        {
            rev_seq += 'A';
        }
        else if (seq[i] == 'C' || seq[i] == 'c')
        {
            rev_seq += 'G';
        }
        else if (seq[i] == 'G' || seq[i] == 'g')
        {
            rev_seq += 'C';
        }else
        {
            rev_seq += seq[i];
        }
        
    }
    return rev_seq;
}

std::vector<std::pair<uint64_t, Anchor>> ILP_index::index_kmers(int32_t hap)
{
    std::vector<std::pair<uint64_t, Anchor>> kmer_index;
    std::string haplotype;
    for (size_t i = 0; i < paths[hap].size(); i++) {
        haplotype += node_seq[paths[hap][i]];
    }
    // transform lower case letters to upper case
    std::transform(haplotype.begin(), haplotype.end(), haplotype.begin(), ::toupper);

    if (haplotype.size() < window + k_mer - 1) return kmer_index;

    int32_t hap_size = haplotype.size();
    std::vector<int32_t> idx_vtx_map(hap_size, -1);
    int32_t count_idx = 0;
    for (size_t i = 0; i < paths[hap].size(); i++) {
        for (size_t j = 0; j < node_seq[paths[hap][i]].size(); j++) {
            idx_vtx_map[count_idx] = paths[hap][i]; // base_idx -> vertex_idx
            count_idx++;
        }
    }

    int32_t count_kmers = window + k_mer - 1;
    
    uint64_t prev_hash = 0;
    for (int32_t i = 0; i <= haplotype.size() - count_kmers; i++) {
        std::string fwd_kmer = std::string(k_mer, 'Z');
        std::string rev_kmer = std::string(k_mer, 'Z');
        int32_t start_idx_fwd = i;
        int32_t start_idx_rev = i;

        for (size_t j = i; j < i + window; j++) {
            std::string kmer = haplotype.substr(j, k_mer);
            std::string local_fwd_kmer = kmer;
            std::string local_rev_kmer = reverse_strand(kmer);

            if (local_fwd_kmer < fwd_kmer) {
                fwd_kmer = local_fwd_kmer;
                start_idx_fwd = j;
            }

            if (local_rev_kmer < rev_kmer) {
                rev_kmer = local_rev_kmer;
                start_idx_rev = j;
            }
        }

        if (fwd_kmer == rev_kmer) continue; // Cannonical k-mer

        uint64_t hash = fnv1a_hash_64(fwd_kmer);
        int32_t start_idx = start_idx_fwd;

        if (fwd_kmer > rev_kmer) {
            hash = fnv1a_hash_64(rev_kmer);
            start_idx = start_idx_rev;
        }

        if (prev_hash == hash) continue; // Skip the same k-mer
        prev_hash = hash; // handles first k-mer as well

        Anchor anchor;
        anchor.h = hap;
        std::unordered_set<int32_t> set;
        std::vector<int32_t> unique_vtxs_vec;
        for (size_t j = start_idx; j < start_idx + k_mer; j++) {
            int32_t value = idx_vtx_map[j];
            if (set.find(value) == set.end()) {
                set.insert(value);
                unique_vtxs_vec.push_back(value);
            }
        }
        set.clear();
        // push the anchor
        for (auto v: unique_vtxs_vec) {
            anchor.k_mers.push_back(v);
        }
        kmer_index.push_back(std::make_pair(hash, anchor));
    }

    bool only_unique = false;
    if (only_unique)
    {
        std::map<uint64_t, Anchor> kmer_index_set;
        for (auto kmer: kmer_index) {
            kmer_index_set[kmer.first] = kmer.second;
        }
        kmer_index.clear();
        for (auto kmer: kmer_index_set) {
            kmer_index.push_back(kmer);
        }
    }

    return kmer_index;
}

std::set<uint64_t> ILP_index::compute_hashes(std::string &read_seq)
{
    // create upper case
    std::transform(read_seq.begin(), read_seq.end(), read_seq.begin(), ::toupper);
    // Find the minimizers in the read and match with the haplotype and return the anchors
    std::set<uint64_t> read_hashes;
    int32_t count_kmers = window + k_mer - 1;
    if (read_seq.size() < count_kmers) return read_hashes;

    uint64_t prev_hash = 0;
    for (int32_t i = 0; i <= read_seq.size() - count_kmers; i++) {
        std::string fwd_kmer = std::string(k_mer, 'Z'); // ZZ...Z
        std::string rev_kmer = std::string(k_mer, 'Z');

        for (size_t j = i; j < i + window; j++) {
            std::string kmer = read_seq.substr(j, k_mer);
            std::string local_fwd_kmer = kmer;
            std::string local_rev_kmer = reverse_strand(kmer);

            if (local_fwd_kmer < fwd_kmer) {
                fwd_kmer = local_fwd_kmer;
            }

            if (local_rev_kmer < rev_kmer) {
                rev_kmer = local_rev_kmer;
            }
        }

        if (fwd_kmer == rev_kmer) continue; // Cannonical k-mer

        uint64_t hash = fnv1a_hash_64(fwd_kmer);

        if (fwd_kmer > rev_kmer) {
            hash = fnv1a_hash_64(rev_kmer);
        }

        if (prev_hash == hash) continue; // Skip the same k-mer
        prev_hash = hash; // handles first k-mer as well

        read_hashes.insert(hash);
    }

    return read_hashes;
}

std::vector<std::vector<std::vector<int32_t>>> ILP_index::compute_anchors(std::vector<std::pair<uint64_t, Anchor>> &minimizers, std::map<uint64_t, int32_t> &read_hashes)
{
    std::vector<std::vector<std::vector<int32_t>>> anchors;
    std::vector<std::vector<std::pair<int32_t, std::vector<int32_t>>>> local_anchors(num_threads);
    #pragma omp parallel for num_threads(num_threads)
    for (int64_t i = 0; i < minimizers.size(); i++)
    {
        int32_t tid = omp_get_thread_num();
        auto minimizer  = minimizers[i];
        auto hash = minimizer.first;
        if (read_hashes.find(hash) != read_hashes.end()) // Found a match
        {
            std::vector<int32_t> anchor;
            for (size_t j = 0; j < minimizer.second.k_mers.size(); j++)
            {
                anchor.push_back(minimizer.second.k_mers[j]);
            }
            local_anchors[tid].push_back(std::make_pair(read_hashes[hash], anchor)); // id, anchor
        }
    }

    anchors.resize(read_hashes.size());
    for (int32_t i = 0; i < num_threads; i++)
    {
        for (auto anchor: local_anchors[i])
        {
            anchors[anchor.first].push_back(anchor.second);  // (id -> anchor_path)
        }
    }
    local_anchors.clear();
    return anchors;
}

// void create_upper_case(std::string &str)
// {
//     // if any index has a lower case then convert the whole string to upper case
//     // only a -> A, c -> C, g -> G, t -> T and n -> N
//     // std::transform(str_.begin(), str_.end(), str_.begin(), ::toupper);
// }

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
    std::vector<int32_t> hap_sizes(num_walks);

    #pragma omp parallel for num_threads(num_threads)
    for (size_t h = 0; h < num_walks; h++)
    {
        std::string haplotype = "";
        for (size_t i = 0; i < paths[h].size(); i++) haplotype += node_seq[paths[h][i]];
        hap_sizes[h] = haplotype.size();
    }

    int32_t num_reads = ip_reads.size();

    // Index the kmers
    std::vector<std::vector<std::pair<uint64_t, Anchor>>> kmer_index(num_walks);
    #pragma omp parallel for num_threads(num_threads)
    for (int32_t h = 0; h < num_walks; h++)
    {
        kmer_index[h] = index_kmers(h);
        fprintf(stderr, "Hap : %d, Kmers : %d\n", h, kmer_index[h].size());
    }
    fprintf(stderr, "[M::%s::%.3f*%.2f] Haplotypes sketched\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

    // Compute the anchors
    int64_t num_kmers = 0;
    std::vector<std::set<uint64_t>> Read_hashes(num_reads);
    std::map<uint64_t, int32_t> Sp_R;
    #pragma omp parallel for num_threads(num_threads)
    for (int32_t r = 0; r < num_reads; r++)
    {
        Read_hashes[r] = compute_hashes(ip_reads[r].second);
    }
    // push all the unique read hashes to a set
    for (int32_t r = 0; r < num_reads; r++)
    {
        for (auto hash: Read_hashes[r])
        {
            Sp_R[hash]++; // duplicate hashes are not allowed
        }
    }
    // reset Sp_R values from 0 to max_Sp_R
    int32_t count_sp_r = 0;
    for (auto &hash: Sp_R)
    {
        hash.second = count_sp_r++; // unique hash -> map_id
    }
    assert(Sp_R.size() == count_sp_r);
    // clear the read hashes
    Read_hashes.clear();

    // print Indexed reads with spectrum size: Sp_R.size()
    fprintf(stderr, "[M::%s::%.3f*%.2f] Indexed reads with spectrum size: %d\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), Sp_R.size());

    std::vector<std::vector<std::vector<std::vector<int32_t>>>> Anchor_hits(Sp_R.size(), std::vector<std::vector<std::vector<int32_t>>>(num_walks));
    // compute the anchors
    for (int32_t h = 0; h < num_walks; h++)
    {
        std::vector<std::vector<std::vector<int32_t>>> loc_match = compute_anchors(kmer_index[h], Sp_R); // parallel execution
        for (int32_t r = 0; r < Sp_R.size(); r++)
        {
            for (auto anchor: loc_match[r])
            {
                Anchor_hits[r][h].push_back(anchor);
            }
        } 
    }

    // find number of kmers
    for (int32_t h = 0; h < num_walks; h++)
    {
        int32_t loc_count = 0;
        for (int32_t r = 0; r < Sp_R.size(); r++)
        {
            loc_count += Anchor_hits[r][h].size();
        }
        num_kmers += loc_count;
        printf("Hap : %d, Anchors : %d\n", h, loc_count);
    }
    Sp_R.clear();

    fprintf(stderr, "[M::%s::%.3f*%.2f] %d unique k-mers matches found\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), num_kmers);

    std::vector<std::vector<int32_t>> in_nodes(n_vtx);
    for (int32_t i = 0; i < n_vtx; i++)
    {
        for (auto v : adj_list[i])
        {
            in_nodes[v].push_back(i);
        }
    }


    // For backtracking the haplotype
    std::string haplotype;
    // Write an ILP with Gurobi
    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_Threads, num_threads);
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Set parameters to speed up the model
        model.set("PreSparsify", "1"); // Sparsify the model
        model.set("Heuristics", "0.50"); // Spent 50% time on heuristics
        model.set("NodefileStart", "0.5"); // 0.5 GB nodefile start
        model.set("Presolve", "2"); // Aggressive presolve to reduce the model size

        // create map to store variables
        std::map<std::string, GRBVar> vars;
        std::vector<GRBVar> Zvars;

        // print started ILP model
        fprintf(stderr, "[M::%s::%.3f*%.2f] ILP model started\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        int32_t c_1 = recombination; // INF no recombination


        bool is_ilp = true; // ILP
        if(is_qclp) is_ilp = false; // QCLP
        int32_t count_kmer_matches = 0;

        if (is_ilp)
        {
            // Kmer constraints
            for (int32_t i = 0; i < count_sp_r; i++) {
                // GRBQuadExpr kmer_expr;
                GRBLinExpr z_expr;
                int32_t temp = 0;
                for (int32_t j = 0; j < num_walks; j++) {
                    for (int32_t k = 0; k < Anchor_hits[i][j].size(); k++) {
                        GRBLinExpr kmer_expr; // kmer-expression
                        std::string extra_var = "z_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                        GRBVar kmer_expr_var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, extra_var);
                        if (Anchor_hits[i][j][k].size() - 1 == 0) continue; // ignore matches with only one vertex
                        // int32_t weight = (k_mer - 1) - (Anchor_hits[i][j][k].size() - 1);
                        // kmer_expr += weight * kmer_expr_var; // weight * z_{i,j,k}
                        // kmer_expr += weight; // z_{i,j,k}
                        for (int32_t l = 1; l < Anchor_hits[i][j][k].size(); l++) {
                            int32_t u = Anchor_hits[i][j][k][l - 1];
                            int32_t v = Anchor_hits[i][j][k][l];
                            std::string var_name = std::to_string(u) + "_" + std::to_string(j) + "_" + std::to_string(v) + "_" + std::to_string(j);
                            if (vars.find(var_name) == vars.end()) // Variable does not exist
                            {
                                vars[var_name] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                            }
                            // kmer_expr += vars[var_name] * kmer_expr_var;
                            kmer_expr += vars[var_name];
                        }
                        int32_t weight = Anchor_hits[i][j][k].size() - 1;
                        model.addConstr(kmer_expr >= weight * kmer_expr_var, "Kmer_constraints_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
                        z_expr += kmer_expr_var;
                        temp += 1;
                        count_kmer_matches++;
                    }
                }
                if (temp != 0)
                {
                    std::string constraint_name = "Kmer_constraints_" + std::to_string(i);
                    int32_t kmer_weight = k_mer - 1;
                    std::string z_var = "z_" + std::to_string(i);
                    GRBVar z_var_r = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, z_var);
                    Zvars.push_back(z_var_r);
                    // model.addQConstr(kmer_expr == kmer_weight * z_var_r, constraint_name);
                    model.addConstr(z_expr == z_var_r, "Z_constraint_" + std::to_string(i));
                }
            }
        } else
        {
            // Kmer constraints
            for (int32_t i = 0; i < count_sp_r; i++) {
                GRBQuadExpr kmer_expr;
                GRBLinExpr z_expr;
                int32_t temp = 0;
                for (int32_t j = 0; j < num_walks; j++) {
                    for (int32_t k = 0; k < Anchor_hits[i][j].size(); k++) {
                        std::string extra_var = "z_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                        GRBVar kmer_expr_var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, extra_var);
                        if (Anchor_hits[i][j][k].size() - 1 == 0) continue; // ignore matches with only one vertex
                        int32_t weight = (k_mer - 1) - (Anchor_hits[i][j][k].size() - 1);
                        kmer_expr += weight * kmer_expr_var; // weight * z_{i,j,k}
                        for (int32_t l = 1; l < Anchor_hits[i][j][k].size(); l++) {
                            int32_t u = Anchor_hits[i][j][k][l - 1];
                            int32_t v = Anchor_hits[i][j][k][l];
                            std::string var_name = std::to_string(u) + "_" + std::to_string(j) + "_" + std::to_string(v) + "_" + std::to_string(j);
                            if (vars.find(var_name) == vars.end()) // Variable does not exist
                            {
                                vars[var_name] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                            }
                            kmer_expr += vars[var_name] * kmer_expr_var;
                        }
                        z_expr += kmer_expr_var;
                        temp += 1;
                    }
                }
                if (temp != 0)
                {
                    std::string constraint_name = "Kmer_constraints_" + std::to_string(i);
                    int32_t kmer_weight = k_mer - 1;
                    std::string z_var = "z_" + std::to_string(i);
                    GRBVar z_var_r = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, z_var);
                    Zvars.push_back(z_var_r);
                    model.addQConstr(kmer_expr == kmer_weight * z_var_r, constraint_name);
                    model.addConstr(z_expr == z_var_r, "Z_constraint_" + std::to_string(i));
                }
            }
        }

        // print count_sp_r_ilp/count_sp_r * 100% kmer matches are in ilp 
        fprintf(stderr, "[M::%s::%.3f*%.2f] %.2f%% k-mer matches are in ILP\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), (count_kmer_matches * 100.0) / num_kmers);

        // clear memory
        for (int32_t i = 0; i < num_walks; i++)
        {
            kmer_index[i].clear();
        }
        kmer_index.clear();
        Anchor_hits.clear();

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

            int32_t u_end = paths[i].back();
            std::string var_name_end = std::to_string(u_end) + "_" + std::to_string(i) + "_e";
            GRBVar var_end = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name_end);
            vars[var_name_end] = var_end;
            end_expr += var_end;
        }

        // set start_expr <= 1 and end_expr <= 1
        model.addConstr(start_expr == 1, "Start_expr");
        model.addConstr(end_expr == 1, "End_expr");

        // Add vertex constraints from paths
        GRBLinExpr vtx_expr;
        // GRBLinExpr recomb_expr;

        std::map<std::string, std::vector<std::string>> new_adj;

        // w/o recombination
        for (int32_t i = 0; i < num_walks; i++)
        {
            for (int32_t idx = 0; idx < paths[i].size() - 1; idx++)
            {
                int32_t u = paths[i][idx];
                int32_t v = paths[i][idx + 1];
                std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(i);
                
                // New adjacency list
                new_adj[std::to_string(u) + "_" + std::to_string(i)].push_back(std::to_string(v) + "_" + std::to_string(i));
                
                if (vars.find(var_name) == vars.end()) { // Variable does not exist
                    GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                    vars[var_name] = var;
                    vtx_expr += 0 * var; // no need without recombination
                }
            }
        }

        std::map<std::string, bool> visited_vtx;
        // with recombination
        for (int32_t u = 0; u < adj_list.size(); u++)
        {
            for (auto v : adj_list[u])
            {
                std::string new_vtx = "w_" + std::to_string(u) + "_" + std::to_string(v);

                std::set<int32_t> haps_u(haps[u].begin(), haps[u].end());
                std::set<int32_t> haps_v(haps[v].begin(), haps[v].end());
                std::set<int32_t> not_common_haps;
                std::set_difference(haps_u.begin(), haps_u.end(), haps_v.begin(), haps_v.end(), std::inserter(not_common_haps, not_common_haps.begin()));
                if (not_common_haps.size() != 0) // new vertex does not exist
                {
                    for (auto j : not_common_haps)
                    {
                        std::string var_name_1 = std::to_string(u) + "_" + std::to_string(j) + "_" + new_vtx;
                        new_adj[std::to_string(u) + "_" + std::to_string(j)].push_back(new_vtx);
                        if (vars.find(var_name_1) == vars.end()) // Variable does not exist
                        {
                            GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name_1);
                            vars[var_name_1] = var;
                        }
                        vtx_expr += c_1 * vars[var_name_1];
                    }

                    for (auto i : haps[v])
                    {
                        std::string var_name_2 = new_vtx + "_" + std::to_string(v) + "_" + std::to_string(i);
                        new_adj[new_vtx].push_back(std::to_string(v) + "_" + std::to_string(i));

                        if (vars.find(var_name_2) == vars.end()) // Variable does not exist
                        {
                            GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name_2);
                            vars[var_name_2] = var;
                        }
                        vtx_expr += c_1 * vars[var_name_2];
                    }
                }
            }
        }
        // int32_t max_recomb = 100;
        // model.addConstr(recomb_expr <= max_recomb, "Recombination_constraint_less");
        // model.addConstr(recomb_expr >= max_recomb, "Recombination_constraint_greater");

        // add (1-z_{i}) constraints
        GRBLinExpr z_expr;
        for (int32_t i = 0; i < Zvars.size(); i++) {
            z_expr += (1 - Zvars[i]);
        }

        obj =  vtx_expr + z_expr;

        // Create the reverse adjacency list
        std::map<std::string, std::vector<std::string>> in_nodes_new;
        for (auto v: new_adj) {
            for (auto u: v.second) {
                in_nodes_new[u].push_back(v.first);
            }
        }

        fprintf(stderr, "[M::%s::%.3f*%.2f] Objective function added to the model\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        // paths based flow constraints
        for (int32_t i = 0; i < num_walks; i++)
        {
            for (int32_t idx = 0; idx < paths[i].size(); idx++)
            {
                if (idx == 0 || idx == paths[i].size() - 1) continue; // skip source and sink nodes
                GRBLinExpr in_expr;
                GRBLinExpr out_expr;

                int32_t v = paths[i][idx];
                std::string vtx = std::to_string(v) + "_" + std::to_string(i);
                for (auto u: in_nodes_new[vtx]) {
                    std::string edge = u + "_" + vtx;
                    in_expr += vars[edge];
                }
                for (auto u: new_adj[vtx]) {
                    std::string edge = vtx + "_" + u;
                    out_expr += vars[edge];
                }

                std::string constraint_name = "Flow_conservation_" + std::to_string(v) + "_" + std::to_string(i);
                model.addConstr(in_expr == out_expr, constraint_name);
            }
        }

        for (int32_t u = 0; u < n_vtx; u++)
        {
            for (auto v : adj_list[u])
            {
                GRBLinExpr in_expr;
                GRBLinExpr out_expr;
                // for w_u_v vertices
                std::string w_vtx = "w_" + std::to_string(u) + "_" + std::to_string(v);
                if (new_adj.find(w_vtx) != new_adj.end()) // w_vtx exists
                {
                    for (auto u: new_adj[w_vtx]) {
                        std::string edge = w_vtx + "_" + u;
                        out_expr += vars[edge];
                    }
                    for (auto u: in_nodes_new[w_vtx]) {
                        std::string edge = u + "_" + w_vtx;
                        in_expr += vars[edge];
                    }
                    std::string constraint_name = "Flow_conservation_" + w_vtx;
                    model.addConstr(in_expr == out_expr, constraint_name);
                }
            }
        }

        // // print new adj list
        // for (auto v: new_adj) {
        //     std::cout << v.first << " -> ";
        //     for (auto u: v.second) {
        //         std::cout << u << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // exit(0);

        // Flow constraints for source nodes
        for (int32_t i = 0; i < num_walks; i++) {
            int32_t u = paths[i][0];
            GRBLinExpr s_expr;
            s_expr += vars["s_" + std::to_string(u) + "_" + std::to_string(i)];
            std::string vtx = std::to_string(u) + "_" + std::to_string(i);
            for (auto v: new_adj[vtx]) {
                std::string edge = vtx + "_" + v;
                s_expr -= vars[edge];
            }
            std::string constraint_name = "Source_conservation_" + std::to_string(u) + "_" + std::to_string(i);
            model.addConstr(s_expr == 0, constraint_name);
        }

        // Flow constraints for sink nodes
        for (int32_t i = 0; i < num_walks; i++) {
            int32_t u = paths[i].back();
            GRBLinExpr e_expr;
            std::string vtx = std::to_string(u) + "_" + std::to_string(i);
            for (auto v: in_nodes_new[vtx]) {
                std::string edge = v + "_" + vtx;
                e_expr += vars[edge];
            }
            std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_e";
            e_expr += -1 * vars[var_name];
            std::string constraint_name = "Sink_conservation_" + std::to_string(u) + "_" + std::to_string(i);
            model.addConstr(e_expr == 0, constraint_name);
        }

        // clear vars
        vars.clear();
        new_adj.clear();
        in_nodes_new.clear();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Flow conservation constraints added to the model\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        model.setObjective(obj, GRB_MINIMIZE);

        // Check the default optimality tolerance
        double defaultTol = model.get(GRB_DoubleParam_OptimalityTol);
        std::cout << "Default Optimality Tolerance: " << defaultTol << std::endl;
        // Set optimality tolerance
        model.set(GRB_DoubleParam_OptimalityTol, 1.00e-08);
        // Optimize model
        model.optimize();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Model optimized\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        // Print constraints
        bool debug = false;
        if (debug)
        {
            printObjectiveFunction(model);
            printConstraints(model);
            printQuadraticConstraints(model);
            printNonZeroVariables(model);
        }

        // Vector to store the names of non-zero binary variables
        std::vector<std::string> path_strs;

        // Get the list of variables in the model
        GRBVar* variables = model.getVars();
        int num_vars = model.get(GRB_IntAttr_NumVars);

        for (int i = 0; i < num_vars; ++i) {
            GRBVar var = variables[i];
            // Check if the variable is binary and non-zero
            if (var.get(GRB_CharAttr_VType) == GRB_BINARY && var.get(GRB_DoubleAttr_X) > 0.0) {
                // first letter is or last letter is e then skip
                std::string var_name = var.get(GRB_StringAttr_VarName);
                if (var_name[0] == 's' || var_name[0] == 'z'  || var_name[var_name.size() - 1] == 'e') {
                    continue;
                }
                path_strs.push_back(var_name);
                if (debug) std::cerr << "var name : " << var_name << std::endl;
            }
        }

        // print paths strs
        std::vector<std::string> path_strs_vec;
        std::set<int32_t> path_vertices;
        int32_t recombination_count = 0;
        for (int i = 0; i < path_strs.size(); i++)
        {
            // std::cout << path_strs[i] << std::endl;
            std::stringstream ss (path_strs[i]);
            std::vector<std::string> tokens;
            std::string item;
            while (std::getline(ss, item, '_')) {
                // Convert the string item to an integer and add to the vector
                tokens.push_back(item);
            }

            std::string hap_1;
            std::string hap_2;
            int32_t u;
            int32_t v;

            if (tokens.size() == 4)
            {
                u = std::stoi(tokens[0]);
                hap_1 = tokens[1];
                v = std::stoi(tokens[2]);
                hap_2 = tokens[3];
                path_strs_vec.push_back(hap_1);
                path_strs_vec.push_back(hap_2);
                path_vertices.insert(u);
                path_vertices.insert(v);
                if (debug) std::cerr << "(vtx, hap) => " << "(" << u << "," << hap_1 << ")" << std::endl;
                if (debug) std::cerr << "(vtx, hap) => " << "(" << v << "," << hap_2 << ")" << std::endl;
            }else
            {
                if (tokens[2] == "w")
                {
                    u = std::stoi(tokens[0]);
                    hap_1 = tokens[1];
                    path_strs_vec.push_back(hap_1);
                    path_vertices.insert(u);
                    if (debug) std::cerr << "(vtx, hap) => " << "(" << u << "," << hap_1 << ")" << std::endl;
                    recombination_count++; // only one count for each w vertex
                }else {
                    v = std::stoi(tokens[3]);
                    hap_2 = tokens[4];
                    path_strs_vec.push_back(hap_2);
                    path_vertices.insert(v);
                    if (debug) std::cerr << "(vtx, hap) => " << "(" << v << "," << hap_2 << ")" << std::endl;
                }
            }
            // pritn token[0] -> token[3]
        }
        path_strs.clear();
        std::cout << "Recombination count: " << recombination_count << std::endl;

        // generate a set of vertices from the path edges
        std::vector<int32_t> hap_path(path_vertices.begin(), path_vertices.end());
        std::sort(hap_path.begin(), hap_path.end(), [&](int32_t a, int32_t b) {
            return top_order_map[a] < top_order_map[b];
        }); // To ensure the path is in topological order
        
        // verify the path vertices by checking if there exist and edge between the vertices
        if (debug) std::cout << "(" << "s" << "," << hap_path[0] << ")" << "->";
        for (int i = 1; i < hap_path.size(); i++)
        {
            int32_t u = hap_path[i - 1];
            int32_t v = hap_path[i];
            bool exist_edge = false;
            for (auto w: adj_list[u])
            {
                if (w == v)
                {
                    exist_edge = true;
                    break;
                }
            }
            if (!exist_edge)
            {
                fprintf(stderr, "Error: No edge between %d and %d\n", u, v);
                exit(1);
            }
            if (debug) std::cout << "(" << u << "," << v << ")" << "->";
        }
        if (debug) std::cout << "(" << hap_path.back() << "," << "e" << ")" << std::endl;

        // Get the path string and store in haplotype
        for (int i = 0; i < hap_path.size(); i++)
        {
            haplotype += node_seq[hap_path[i]];
        }

    } catch (GRBException e) {
        std::cerr << "Error code = " << e.getErrorCode() << std::endl;
        std::cerr << e.getMessage() << std::endl;
    } catch (...) {
        std::cerr << "Exception during optimization" << std::endl;
    }
    
    // write haplotype as to a file as fasta from the path
    std::string path_str = haplotype;
    std::ofstream hap_file_stream(hap_file, std::ios::out);
    hap_file_stream << ">" << hap_name << " LN:" << path_str.size() << std::endl;
    // write the path_str to the file 80 characters per line
    for (size_t i = 0; i < path_str.size(); i += 80) {
        hap_file_stream << path_str.substr(i, 80) << std::endl;
    }
    hap_file_stream.close();

    fprintf(stderr, "[M::%s::%.3f*%.2f] Haplotype of size: %d written to: %s\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), path_str.size(), hap_file.c_str());
}