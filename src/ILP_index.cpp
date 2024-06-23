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

std::unordered_map<uint64_t, Anchor> ILP_index::index_kmers(int32_t hap)
{
    std::unordered_map<uint64_t, Anchor> kmer_index;
    std::string haplotype;
    for (size_t i = 0; i < paths[hap].size(); i++) {
        haplotype += node_seq[paths[hap][i]];
    }

    int32_t count_kmers = window + k_mer - 1;

    for (int32_t i = 0; i <= haplotype.size() - count_kmers; i++) {
        uint64_t fwd_hash = std::numeric_limits<uint64_t>::max();
        uint64_t rev_hash = std::numeric_limits<uint64_t>::max();
        int32_t start_idx_fwd = i;
        int32_t start_idx_rev = i;

        for (size_t j = i; j < i + window; j++) {
            std::string kmer = haplotype.substr(j, k_mer);
            uint64_t local_fwd_hash = fnv1a_hash_64(kmer);
            uint64_t local_rev_hash = fnv1a_hash_64(reverse_strand(kmer));

            if (local_fwd_hash < fwd_hash) {
                fwd_hash = local_fwd_hash;
                start_idx_fwd = j;
            }

            if (local_rev_hash < rev_hash) {
                rev_hash = local_rev_hash;
                start_idx_rev = j;
            }
        }

        if (fwd_hash == rev_hash) continue; // Cannonical k-mer

        uint64_t hash = fwd_hash;
        int32_t start_idx = start_idx_fwd;

        if (fwd_hash > rev_hash) {
            hash = rev_hash;
            start_idx = start_idx_rev;
        }

        Anchor anchor;
        anchor.h = hap;
        for (size_t j = start_idx; j < start_idx + k_mer; j++) {
            anchor.k_mers.push_back(paths[hap][j]);
        }
        kmer_index[hash] = anchor; // unique k-mer
    }

    return kmer_index;
}

std::set<uint64_t> ILP_index::compute_hashes(std::string &read_seq)
{
    // Find the minimizers in the read and match with the haplotype and return the anchors
    std::set<uint64_t> read_hashes;
    int32_t count_kmers = window + k_mer - 1;
    for (int32_t i = 0; i <= read_seq.size() - count_kmers; i++) {
        uint64_t fwd_hash = std::numeric_limits<uint64_t>::max();
        uint64_t rev_hash = std::numeric_limits<uint64_t>::max();

        for (size_t j = i; j < i + window; j++) {
            std::string kmer = read_seq.substr(j, k_mer);
            uint64_t local_fwd_hash = fnv1a_hash_64(kmer);
            uint64_t local_rev_hash = fnv1a_hash_64(reverse_strand(kmer));

            if (local_fwd_hash < fwd_hash) {
                fwd_hash = local_fwd_hash;
            }

            if (local_rev_hash < rev_hash) {
                rev_hash = local_rev_hash;
            }
        }

        if (fwd_hash == rev_hash) continue; // Cannonical k-mer

        uint64_t hash = fwd_hash;

        if (fwd_hash > rev_hash) {
            hash = rev_hash;
        }

        read_hashes.insert(hash);
    }

    return read_hashes;
}

std::vector<Anchor> ILP_index::compute_anchors(std::unordered_map<uint64_t, Anchor> &minimizers, std::set<uint64_t> &read_hashes)
{
    std::vector<Anchor> anchors;
    for (auto hash: read_hashes)
    {
        if (minimizers.find(hash) != minimizers.end())
        {
            anchors.push_back(minimizers[hash]);
        }
    }
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
    std::vector<std::unordered_map<uint64_t, Anchor>> kmer_index(num_walks);
    #pragma omp parallel for num_threads(num_threads)
    for (int32_t h = 0; h < num_walks; h++)
    {
        kmer_index[h] = index_kmers(h);
        // fprintf(stderr, "Hap : %d, Kmers : %d\n", h, kmer_index[h].size());
    }
    fprintf(stderr, "[M::%s::%.3f*%.2f] Haplotypes sketched\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

    // Compute the anchors
    int64_t num_kmers = 0;
    std::vector<std::set<uint64_t>> Read_hashes(num_reads);
    std::set<uint64_t> Unique_read_hashes;
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
            Unique_read_hashes.insert(hash);
        }
    }
    // clear the read hashes
    Read_hashes.clear();

    std::vector<std::vector<Anchor>> Anchor_hits(num_walks);
    std::vector<std::vector<std::vector<int32_t>>> Kmers(num_walks);
    // compute the anchors
    #pragma omp parallel for num_threads(num_threads)
    for (int32_t h = 0; h < num_walks; h++)
    {
        Anchor_hits[h] = compute_anchors(kmer_index[h], Unique_read_hashes);
        for (auto anchor: Anchor_hits[h])
        {
            Kmers[h].push_back(anchor.k_mers);
        }
    }
    // find number of kmers
    for (int32_t h = 0; h < num_walks; h++)
    {
        num_kmers += Kmers[h].size();
        printf("Hap : %d, Anchors : %d\n", h, Kmers[h].size());
    }

    fprintf(stderr, "[M::%s::%.3f*%.2f] %d unique k-mers matches found\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), num_kmers);
    exit(0);

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
        model.set("PreSparsify", "1");
        model.set(GRB_DoubleParam_NodefileStart, 128.0); // Beyond this memory, it will write to disk

        // create map to store variables
        std::map<std::string, GRBVar> vars;

        // // Kmer constraints
        // for (int32_t i = 0; i < num_walks; i++) {
        //     for (int32_t j = 0; j < Kmers[i].size(); j++) {
        //         GRBLinExpr kmer_expr;
        //         for (int32_t k = 1; k < Kmers[i][j].size(); k++) {
        //             int32_t u = Kmers[i][j][k-1];
        //             int32_t v = Kmers[i][j][k];
        //             std::string var_name = std::to_string(u) + "_" + std::to_string(i) + "_" + std::to_string(v) + "_" + std::to_string(i);
        //             vars[var_name] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        //             kmer_expr += vars[var_name];
        //         }
        //         std::string exra_var = "z_" + std::to_string(i) + "_" + std::to_string(j);
        //         GRBVar kmer_expr_var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, exra_var);
        //         vars[exra_var] = kmer_expr_var;
        //         std::string constraint_name = "Kmer_constraints_" + std::to_string(i) + "_" + std::to_string(j);
        //         int32_t kmer_weight = Kmers[i][j].size() - 1;
        //         model.addConstr(kmer_expr == kmer_weight * kmer_expr_var, constraint_name);
        //     }
        // }

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
                if (var_name[0] == 's' || var_name[var_name.size() - 1] == 'e') {
                    continue;
                }
                path_strs.push_back(var_name);
            }
        }

        // print paths strs
        std::vector<std::pair<int32_t, int32_t>> path_edges;
        for (int i = 0; i < path_strs.size(); i++)
        {
            // std::cout << path_strs[i] << std::endl;
            std::stringstream ss (path_strs[i]);
            std::vector<int32_t> tokens;
            std::string item;
            while (std::getline(ss, item, '_')) {
                // Convert the string item to an integer and add to the vector
                tokens.push_back(std::stoi(item));
            }
            int32_t u = tokens[0];
            int32_t v = tokens[2];
            path_edges.push_back(std::make_pair(u, v));
        }
        path_strs.clear();

        // pritn the path edges
        for (int i = 0; i < path_edges.size(); i++)
        {
            std::cout << path_edges[i].first << " " << path_edges[i].second << std::endl;
        }

        // generate a set of vertices from the path edges
        std::set<int32_t> path_vertices;
        for (int i = 0; i < path_edges.size(); i++)
        {
            path_vertices.insert(path_edges[i].first);
            path_vertices.insert(path_edges[i].second);
        }
        std::vector<int32_t> hap_path(path_vertices.begin(), path_vertices.end());
        // verify the path vertices by checking if there exist and edge between the vertices
        for (int i = 1; i < path_edges.size(); i++)
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
        }

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

    fprintf(stderr, "[M::%s::%.3f*%.2f] Haplotype written to: %s\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), hap_file.c_str());
}