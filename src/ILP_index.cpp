#include "ILP_index.h"

// Constructor
ILP_index::ILP_index(gfa_t *g) {
    this->g = g;
}

KSEQ_INIT(gzFile, gzread)

uint64_t hash128_to_64(const std::string &str) {
    uint64_t hash_output[2]; // To store the 128-bit output

    // Hash the string using MurmurHash3_x64_128
    MurmurHash3_x64_128(str.data(), str.length(), 0, hash_output);

    // Combine the two 64-bit parts of the 128-bit hash into a single 64-bit hash
    return hash_output[0] ^ hash_output[1]; // Simple XOR combination
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
        std::string walk_name = std::string(g->walk[w].sample) + "." + std::to_string(g->walk[w].hap);
        hap_id2name.push_back(walk_name);
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

std::vector<std::pair<uint64_t, Anchor>> ILP_index::index_kmers(int32_t hap) {
    std::vector<std::pair<uint64_t, Anchor>> kmer_index;
    std::string haplotype;

    // Concatenate sequences in the path to form the haplotype
    for (size_t i = 0; i < paths[hap].size(); i++) {
        haplotype += node_seq[paths[hap][i]];
    }

    // Convert to uppercase
    std::transform(haplotype.begin(), haplotype.end(), haplotype.begin(), ::toupper);
    size_t hap_size = haplotype.size();

    if (hap_size < window + k_mer - 1) return kmer_index;

    // Precompute index to vertex map
    std::vector<int32_t> idx_vtx_map(hap_size, -1);
    int32_t count_idx = 0;
    for (size_t i = 0; i < paths[hap].size(); i++) {
        for (size_t j = 0; j < node_seq[paths[hap][i]].size(); j++) {
            idx_vtx_map[count_idx++] = paths[hap][i]; // base_idx -> vertex_idx
        }
    }

    uint64_t prev_hash = UINT64_MAX; // Initialize to the maximum possible value

    // Deque to store k-mers and their positions
    std::deque<std::pair<std::string, int32_t>> window_deque;

    for (int32_t i = 0; i <= hap_size - k_mer; ++i) {
        // Extract the forward and reverse k-mers
        std::string fwd_kmer = haplotype.substr(i, k_mer);
        std::string rev_kmer = reverse_strand(fwd_kmer);

        // Choose the lexicographically smallest k-mer
        std::string min_kmer = std::min(fwd_kmer, rev_kmer);

        // Remove elements from the back of the deque if they are larger than the current k-mer
        while (!window_deque.empty() && window_deque.back().first >= min_kmer) {
            window_deque.pop_back();
        }

        // Add current element at the back of the deque
        window_deque.emplace_back(min_kmer, i);

        // Remove elements from the front if they are out of the window
        if (!window_deque.empty() && window_deque.front().second <= i - window) {
            window_deque.pop_front();
        }

        // Process the minimizer once the window is fully populated
        if (i >= window - 1) {  
            std::string best_kmer = window_deque.front().first;
            uint64_t best_hash = hash128_to_64(best_kmer.c_str());
            if (best_hash != prev_hash) {
                prev_hash = best_hash;

                // Create anchor
                Anchor anchor;
                anchor.h = hap;
                std::unordered_set<int32_t> unique_vtx_set;
                std::vector<int32_t> unique_vtxs_vec;

                // Map k-mer positions to vertex indices
                int32_t best_start_idx = window_deque.front().second;
                for (int32_t j = best_start_idx; j < best_start_idx + k_mer; j++) {
                    int32_t vtx_idx = idx_vtx_map[j];
                    if (unique_vtx_set.insert(vtx_idx).second) {
                        unique_vtxs_vec.push_back(vtx_idx);
                    } else if (vtx_idx != unique_vtxs_vec.back()) {
                        unique_vtxs_vec.push_back(vtx_idx);
                    }
                }

                // Assign unique vertices to the anchor
                anchor.k_mers = std::move(unique_vtxs_vec);
                kmer_index.emplace_back(best_hash, anchor);
            }
        }
    }

    return kmer_index;
}

std::set<uint64_t> ILP_index::compute_hashes(std::string &read_seq) {
    // Convert to uppercase
    std::transform(read_seq.begin(), read_seq.end(), read_seq.begin(), ::toupper);

    std::set<uint64_t> read_hashes;
    int32_t seq_size = read_seq.size();
    if (seq_size < window + k_mer - 1) return read_hashes;

    uint64_t prev_hash = UINT64_MAX; // Initialize to maximum possible value

    // Precompute the reverse complement of the entire sequence
    std::string rev_seq = reverse_strand(read_seq);

    // Deque to store k-mers and their positions
    std::deque<std::pair<std::string, int32_t>> window_deque;

    for (int32_t i = 0; i <= seq_size - k_mer; ++i) {
        // Extract the forward and reverse k-mers
        std::string fwd_kmer = read_seq.substr(i, k_mer);
        std::string rev_kmer = rev_seq.substr(seq_size - i - k_mer, k_mer);

        // Choose the lexicographically smallest k-mer
        std::string min_kmer = std::min(fwd_kmer, rev_kmer);

        // Remove elements from the back of the deque if they are larger than the current k-mer
        while (!window_deque.empty() && window_deque.back().first >= min_kmer) {
            window_deque.pop_back();
        }

        // Add current element at the back of the deque
        window_deque.emplace_back(min_kmer, i);

        // Remove elements from the front if they are out of the window
        if (!window_deque.empty() && window_deque.front().second <= i - window) {
            window_deque.pop_front();
        }

        // Only insert the hash if it is different from the previous hash
        if (i >= window - 1) {  // Once we reach the window size
            std::string best_kmer = window_deque.front().first;
            uint64_t best_hash = hash128_to_64(best_kmer.c_str());
            if (best_hash != prev_hash) {
                read_hashes.insert(best_hash);
                prev_hash = best_hash;
            }
        }
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
    std::cerr << "Number of Minimizers" << std::endl;
    std::vector<std::vector<std::pair<uint64_t, Anchor>>> kmer_index(num_walks);
    #pragma omp parallel for num_threads(num_threads)
    for (int32_t h = 0; h < num_walks; h++)
    {
        kmer_index[h] = index_kmers(h);
        std::string hap_name = hap_id2name[h];
        fprintf(stderr, "%s : %d\n", hap_name.c_str(), kmer_index[h].size());
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
    Sp_R.clear();

    // find number of kmers
    int32_t num_kmers_tot = 0;
    for (int32_t h = 0; h < num_walks; h++)
    {
        int32_t loc_count = 0;
        for (int32_t r = 0; r < count_sp_r; r++)
        {
            loc_count += Anchor_hits[r][h].size();
        }
        num_kmers_tot += loc_count;
    }

    std::vector<std::vector<std::vector<std::vector<int32_t>>>> Anchor_hits_1(
    count_sp_r, std::vector<std::vector<std::vector<int32_t>>>(num_walks));

    #pragma omp parallel for num_threads(num_threads)
    for (int32_t r = 0; r < count_sp_r; r++) {
        std::map<std::string, std::pair<int32_t, std::vector<std::pair<int32_t, std::vector<int32_t>>>>> Anchor_hits_map;

        for (int32_t h = 0; h < num_walks; h++) {
            for (int32_t k = 0; k < Anchor_hits[r][h].size(); k++) {
                std::string anchor_str;
                for (auto v : Anchor_hits[r][h][k]) {
                    anchor_str += std::to_string(v) + "_";
                }

                if (Anchor_hits_map.find(anchor_str) == Anchor_hits_map.end()) { // does not exist
                    Anchor_hits_map[anchor_str].first = 1;
                    Anchor_hits_map[anchor_str].second.push_back(std::make_pair(h, Anchor_hits[r][h][k]));
                } else {
                    Anchor_hits_map[anchor_str].first++;
                    Anchor_hits_map[anchor_str].second.push_back(std::make_pair(h, Anchor_hits[r][h][k]));
                }
            }
        }

        // Check for all horizontal matches
        for (const auto &anchor : Anchor_hits_map) {
            if (anchor.second.first < threshold * num_walks) {
                for (const auto &hap_anchor : anchor.second.second) {
                    Anchor_hits_1[r][hap_anchor.first].push_back(hap_anchor.second);
                }
            }
        }
    }

    // Replace Anchor_hits with Anchor_hits_1
    Anchor_hits = std::move(Anchor_hits_1);
    Anchor_hits_1.clear();

    // find number of kmers
    std::cerr << "Number of Anchors" << std::endl;
    for (int32_t h = 0; h < num_walks; h++)
    {
        int32_t loc_count = 0;
        for (int32_t r = 0; r < count_sp_r; r++)
        {
            loc_count += Anchor_hits[r][h].size();
        }
        num_kmers += loc_count;
        std::string hap_name = hap_id2name[h];
        fprintf(stderr, "%s : %d\n", hap_name.c_str(), loc_count);
    }

    // pritn num_kmers/num_kmers_tot * 100 % are part of the haplotype
    fprintf(stderr, "[M::%s::%.3f*%.2f] Total/Filtered anchors: %d/%d, Fraction of anchors retained: %.4f\n",
        __func__, 
        realtime() - mg_realtime0, 
        cputime() / (realtime() - mg_realtime0), 
        num_kmers_tot, 
        num_kmers_tot - num_kmers, 
        (float)num_kmers / (float)num_kmers_tot);

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
        model.set("Method", "3"); // Concurrent method

        // create map to store variables
        std::map<std::string, GRBVar> vars;
        std::vector<GRBVar> Zvars;
        int32_t c_1 = recombination; // INF no recombination

        bool is_ilp = true; // ILP
        if(is_qclp) is_ilp = false; // QCLP
        int32_t count_kmer_matches = 0;

        if (is_ilp)
        {
            fprintf(stderr, "[M::%s::%.3f*%.2f] ILP model started\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
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
                    count_kmer_matches++;
                }
            }
        } else
        {
            fprintf(stderr, "[M::%s::%.3f*%.2f] QCLP model started\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
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
                    count_kmer_matches++;
                }
            }
        }

        // print count_sp_r_ilp/count_sp_r * 100% kmer matches are in ilp 
        fprintf(stderr, "[M::%s::%.3f*%.2f] %.2f%% Kmers are in ILP\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), (count_kmer_matches * 100.0) / count_sp_r);

        // clear memory
        for (int32_t i = 0; i < num_walks; i++)
        {
            kmer_index[i].clear();
        }
        kmer_index.clear();
        Anchor_hits.clear();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Kmer constraints added to the model\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        if (is_naive_exp)
        {
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

            // Add vertex constraints from paths
            GRBLinExpr vtx_expr;
            // GRBLinExpr recomb_expr;

            // w/o recombination
            for (int32_t i = 0; i < num_walks; i++)
            {
                for (int32_t idx = 0; idx < paths[i].size() - 1; idx++)
                {
                    int32_t u = paths[i][idx];
                    int32_t v = paths[i][idx + 1];
                    std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(i);
                    if (vars.find(var_name) == vars.end()) { // Variable does not exist
                        GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                        vars[var_name] = var;
                        vtx_expr += 0 * var; // no need without recombination
                    }
                }
            }

            // with recombination
            for (int32_t i = 0; i < num_walks; i++)
            {
                for (int32_t idx = 0; idx < paths[i].size() - 1; idx++)
                {
                    int32_t u = paths[i][idx];
                    for (auto v : adj_list[u])
                    {
                        for (auto j : haps[v])
                        {
                            if (i != j)
                            {
                                std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(j);
                                if (vars.find(var_name) == vars.end()) // Variable does not exist
                                {
                                    GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                                    vars[var_name] = var;
                                    vtx_expr += c_1 * var;
                                    // recomb_expr += var;
                                }
                            }
                        }
                    }
                }
            }

            // add (1-z_{i}) constraints
            GRBLinExpr z_expr;
            for (int32_t i = 0; i < Zvars.size(); i++) {
                z_expr += (1 - Zvars[i]);
            }

            obj =  vtx_expr + z_expr;
            model.setObjective(obj, GRB_MINIMIZE);

            // paths based flow constraints
            for (int32_t i = 0; i < num_walks; i++)
            {
                for (int32_t idx = 0; idx < paths[i].size(); idx++)
                {
                    if (idx == 0 || idx == paths[i].size() - 1) continue; // skip source and sink nodes
                    GRBLinExpr in_expr;
                    GRBLinExpr out_expr;
                    int32_t v = paths[i][idx];
                    int32_t v_in = paths[i][idx - 1];
                    int32_t v_out = paths[i][idx + 1];
                    std::string var_name = std::to_string(v_in) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(i);
                    if (vars.find(var_name) == vars.end()) { // Variable does not exist
                        GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                        vars[var_name] = var;
                    }
                    in_expr += vars[var_name];

                    var_name = std::to_string(v) + "_" + std::to_string(i) + "_" + std::to_string(v_out) + "_" + std::to_string(i);
                    if (vars.find(var_name) == vars.end()) { // Variable does not exist
                        GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                        vars[var_name] = var;
                    }
                    out_expr += vars[var_name];
                    
                    // In expression
                    for (auto u : in_nodes[v])
                    {
                        for (auto j : haps[u])
                        {
                            if (i != j)
                            {
                                var_name = std::to_string(u) + "_" + std::to_string(j) + "_" + std::to_string(v) + "_" + std::to_string(i);
                                if (vars.find(var_name) == vars.end()) { // Variable does not exist
                                    GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                                    vars[var_name] = var;
                                }
                                in_expr += vars[var_name];
                            }
                        }
                    }

                    // Out expression
                    for (auto u : adj_list[v])
                    {
                        for (auto j : haps[u])
                        {
                            if (i != j)
                            {
                                var_name = std::to_string(v) + "_" + std::to_string(i) + "_" + std::to_string(u) + "_" + std::to_string(j);
                                if (vars.find(var_name) == vars.end()) { // Variable does not exist
                                    GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                                    vars[var_name] = var;
                                }
                                out_expr += vars[var_name];
                            }
                        }
                    }
                    std::string constraint_name = "Flow_conservation_" + std::to_string(v) + "_" + std::to_string(i);
                    model.addConstr(in_expr == out_expr, constraint_name);
                }
            }
            

            // Flow constraints for source nodes
            for (int32_t i = 0; i < num_walks; i++) {
                int32_t u = paths[i][0];
                GRBLinExpr s_expr;
                s_expr += vars["s_" + std::to_string(u) + "_" + std::to_string(i)];
                for (auto v: adj_list[u]) {
                    for (auto j: haps[v]) {
                        std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(j);
                        if (vars.find(var_name) == vars.end()) { // Variable does not exist
                            GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                            vars[var_name] = var;
                        }
                        s_expr -= vars[var_name];
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
                        std::string var_name = std::to_string(v) + "_" + std::to_string(j) + "_" + std::to_string(u) + "_" + std::to_string(i);
                        if (vars.find(var_name) == vars.end()) { // Variable does not exist
                            GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                            vars[var_name] = var;
                        }
                        e_expr += vars[var_name];
                    }
                }
                std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_e";
                if (vars.find(var_name) == vars.end()) { // Variable does not exist
                    GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                    vars[var_name] = var;
                }
                e_expr += -1 * vars[var_name];
                std::string constraint_name = "Sink_conservation_" + std::to_string(u) + "_" + std::to_string(i);
                model.addConstr(e_expr == 0, constraint_name);
            }

            // clear vars
            vars.clear();
            Zvars.clear();
            fprintf(stderr, "[M::%s::%.3f*%.2f] Naive expanded graph constructed\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
        }else
        {
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

            // Populate hash table to quickly get vertex index in haplotype
            std::vector<std::unordered_map<int, int>> elementIndexMaps(paths.size());
            for(int h = 0; h < paths.size(); h++)
            {
                std::unordered_map<int, int> elementIndexMap;
                for (int i = 0; i < paths[h].size(); ++i) 
                {
                    elementIndexMap[paths[h][i]] = i;
                }
                elementIndexMaps[h] = elementIndexMap;
            }

            // Add recombination vertices and edges
            for (int32_t u = 0; u < adj_list.size(); u++)
            {
                for (auto v : adj_list[u])
                {
    
                    std::string new_vtx = "w_" + std::to_string(u) + "_" + std::to_string(v);

                    bool new_vertex_used = false;
                    for(auto h : haps[u])
                    {
                        // check if the next entry paths[h] after u is v
                        int index = elementIndexMaps[h][u];

                        if(index == paths[h].size()-1 || paths[h][index+1] != v)
                        {

                            if(!new_vertex_used)
                            {
                                new_vertex_used = true;
                            }
                            std::string var_name_1 = std::to_string(u) + "_" + std::to_string(h) + "_" + new_vtx;
                            new_adj[std::to_string(u) + "_" + std::to_string(h)].push_back(new_vtx);
                            
                            if (vars.find(var_name_1) == vars.end()) // Variable does not exist
                            {
                                GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name_1);
                                vars[var_name_1] = var;
                            }
                            vtx_expr += c_1 * vars[var_name_1];

                        }
                    }

                    if(new_vertex_used)
                    {
                        for(auto h : haps[v])
                        {
                            std::string var_name_2 = new_vtx + "_" + std::to_string(v) + "_" + std::to_string(h);
                            new_adj[new_vtx].push_back(std::to_string(v) + "_" + std::to_string(h));

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
            elementIndexMaps.clear();

            // add (1-z_{i}) constraints
            GRBLinExpr z_expr;
            for (int32_t i = 0; i < Zvars.size(); i++) {
                z_expr += (1 - Zvars[i]);
            }

            obj =  vtx_expr + z_expr;

            model.setObjective(obj, GRB_MINIMIZE);

            // Create the reverse adjacency list
            std::map<std::string, std::vector<std::string>> in_nodes_new;
            for (auto v: new_adj) {
                for (auto u: v.second) {
                    in_nodes_new[u].push_back(v.first);
                }
            }

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

            fprintf(stderr, "[M::%s::%.3f*%.2f] Optimized expanded graph constructed\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
        }


        // Check the default optimality tolerance
        double defaultTol = model.get(GRB_DoubleParam_OptimalityTol);
        std::cout << "Default Optimality Tolerance: " << defaultTol << std::endl;
        // Set optimality tolerance
        model.set(GRB_DoubleParam_OptimalityTol, 1.00e-08);
        // Optimize model
        model.optimize();

        fprintf(stderr, "[M::%s::%.3f*%.2f] Model optimized\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

        // Print constraints
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
        std::set<std::pair<int32_t, int32_t>> path_vertices_hap;
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
                path_vertices_hap.insert({u, std::stoi(hap_1)});
                path_vertices_hap.insert({v, std::stoi(hap_2)});
                if (debug) std::cerr << "(vtx, hap) => " << "(" << u << "," << hap_1 << ")" << std::endl;
                if (debug) std::cerr << "(vtx, hap) => " << "(" << v << "," << hap_2 << ")" << std::endl;
            }else
            {
                if (tokens[2] == "w")
                {
                    u = std::stoi(tokens[0]);
                    hap_1 = tokens[1];
                    path_vertices_hap.insert({u, std::stoi(hap_1)});
                    if (debug) std::cerr << "(vtx, hap) => " << "(" << u << "," << hap_1 << ")" << std::endl;
                    int32_t u_int = std::stoi(tokens[0]);
                    int32_t hap_int = std::stoi(tokens[1]);
                }else {
                    v = std::stoi(tokens[3]);
                    hap_2 = tokens[4];
                    path_vertices_hap.insert({v, std::stoi(hap_2)});
                    if (debug) std::cerr << "(vtx, hap) => " << "(" << v << "," << hap_2 << ")" << std::endl;
                    int32_t v_int = std::stoi(tokens[3]);
                    int32_t hap_int = std::stoi(tokens[4]);
                }
            }
        }
        std::vector<std::pair<int32_t, int32_t>> path_vertices_hap_vec(path_vertices_hap.begin(), path_vertices_hap.end());
        // sort path_vertices_hap_vec based on first element
        std::sort(path_vertices_hap_vec.begin(), path_vertices_hap_vec.end(), [&](std::pair<int32_t, int32_t> a, std::pair<int32_t, int32_t> b) {
            return top_order_map[a.first] < top_order_map[b.first];
        });

        int32_t prev_hap = path_vertices_hap_vec[0].second;
        int32_t prev_str_id = 0;
        int32_t str_id = 0;
        std::vector<std::string> hap_st_en_vec;
        str_id += node_seq[path_vertices_hap_vec[0].first].size();
        for (int32_t i = 1; i < path_vertices_hap_vec.size(); ++i)
        {
            str_id += node_seq[path_vertices_hap_vec[i].first].size();

            if (prev_hap != path_vertices_hap_vec[i].second) // only prints prev_hap not current hap hence additionally we need to print the last segment
            {
                recombination_count++;
                std::string str = ">(" + hap_id2name[prev_hap] + ",[" + std::to_string(prev_str_id) + "," + std::to_string(str_id - 1) + "])";
                hap_st_en_vec.push_back(str);
                prev_hap = path_vertices_hap_vec[i].second;
                prev_str_id = str_id;
            }
        }

        // Capture the last segment after the loop
        std::string str = ">(" + hap_id2name[path_vertices_hap_vec.back().second] + ",[" + std::to_string(prev_str_id) + "," + std::to_string(str_id - 1) + "])";
        hap_st_en_vec.push_back(str);


        std::cerr << "Recombination count: " << recombination_count << std::endl;
        if (recombination_count > 0)
        {
            std::cerr << "Recombined haplotypes: ";
            for (int i = 0; i < hap_st_en_vec.size(); i++)
            {
                std::cerr << hap_st_en_vec[i];
            }
            std::cerr << std::endl;
        }else
        {
            std::cerr << "Recombined haplotypes: ";
            int32_t sum_str = 0;
            for (int i = 0; i < path_vertices_hap_vec.size(); i++)
            {
                sum_str += node_seq[path_vertices_hap_vec[i].first].size();
            }
            std::cerr << ">(" << hap_id2name[prev_hap] << ",[" << 0 << "," << sum_str - 1 << "])" << std::endl;
        }
        
        
        // verify the path vertices by checking if there exist and edge between the vertices
        if (debug) std::cout << "(" << "s" << "," << path_vertices_hap_vec[0].first << ")" << "->";
        for (int i = 1; i < path_vertices_hap_vec.size(); i++)
        {
            int32_t u = path_vertices_hap_vec[i - 1].first;
            int32_t v = path_vertices_hap_vec[i].first;
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
                fprintf(stderr, "Error: No edge between %d and %d\n", u , v);
                exit(1);
            }
            if (debug) std::cout << "(" << u << "," << v << ")" << "->";
        }
        if (debug) std::cout << "(" << path_vertices_hap_vec.back().first << "," << "e" << ")" << std::endl;

        // Get the path string and store in haplotype
        for (int i = 0; i < path_vertices_hap_vec.size(); i++)
        {
            haplotype += node_seq[path_vertices_hap_vec[i].first];
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