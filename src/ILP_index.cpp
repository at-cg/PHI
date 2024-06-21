#include "ILP_index.h"

// Constructor
ILP_index::ILP_index(gfa_t *g) {
    this->g = g;
}

KSEQ_INIT(gzFile, gzread)

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