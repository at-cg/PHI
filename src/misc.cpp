#include <stdlib.h>
#include "PHIpriv.h"
#include "ksort.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cassert>
#include <string>

void dfs(int u, std::vector<int>& visited, std::vector<int>& parent, std::vector<int>& low,
    std::vector<int>& disc, int& time, std::vector<std::vector<int>>& adj_list,
    std::vector<std::pair<int, int>>& bridges) {
    visited[u] = 1;
    low[u] = disc[u] = ++time;

    for (int v : adj_list[u]) {
        if (!visited[v]) {
            parent[v] = u;
            dfs(v, visited, parent, low, disc, time, adj_list, bridges);
            low[u] = std::min(low[u], low[v]);
            if (low[v] > disc[u]) {
                bridges.push_back({u, v});
            }
        } else if (v != parent[u]) {
            low[u] = std::min(low[u], disc[v]);
        }
    }
}

std::vector<std::pair<int, int>> find_bridges(std::vector<std::vector<int>> &adj_list) {
    std::vector<std::pair<int, int>> bridges;
    int n = adj_list.size();
    std::vector<int> visited(n, 0);
    std::vector<int> parent(n, -1);
    std::vector<int> low(n, 0);
    std::vector<int> disc(n, 0);
    int time = 0;

    std::cerr << "n = " << n << std::endl;

    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            dfs(i, visited, parent, low, disc, time, adj_list, bridges);
        }
    }

    // reverse the order of the bridges
    std::reverse(bridges.begin(), bridges.end());

    return bridges;
}

bool compare_T(const mg128_t &a, const mg128_t &b) {
	return a.x < b.x;
};

void get_hap_name(char *gfa_name, char *reads_name, std::string &hap_name)
{
	// split opt.gfa_file.c_str() to get the hap_name by last '/' and .gfa
    std::string gfa_file = gfa_name;
    size_t found = gfa_file.find_last_of("/\\");
    if (found != std::string::npos) {
        hap_name = gfa_file.substr(found+1);
    } else {
        hap_name = gfa_file;
    }
    found = hap_name.find_last_of(".");
    if (found != std::string::npos) {
        hap_name = hap_name.substr(0, found);
    }

    // split opt.reads_file to get the hap_name by last '/' and .fa and add to hap_name
    std::string reads_file = reads_name;
    found = reads_file.find_last_of("/\\");
    if (found != std::string::npos) {
        hap_name += "_";
        hap_name += reads_file.substr(found+1);
    } else {
        hap_name += "_";
        hap_name += reads_file;
    }
    found = hap_name.find_last_of(".");
    if (found != std::string::npos) {
        hap_name = hap_name.substr(0, found);
    }
}

int mg_verbose = 1;
int mg_dbg_flag = 0;
double mg_realtime0;


#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mg128_t, sort_key_128x, 8) 

KSORT_INIT_GENERIC(uint32_t)