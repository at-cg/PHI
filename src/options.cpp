#include <string.h>
#include "PHIpriv.h"

void mg_idxopt_init(mg_idxopt_t *io)
{
	memset(io, 0, sizeof(mg_idxopt_t));
	io->k = 31;
	io->w = 25;
	io->bucket_bits = 14;
}

void mg_mapopt_init(mg_mapopt_t *mo)
{
	memset(mo, 0, sizeof(mg_mapopt_t));
	mo->cap_kalloc = 1000000000;
    mo->n_threads = 4;
}

int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo)
{
	if (preset == 0) {
		mg_idxopt_init(io);
		mg_mapopt_init(mo);
	}
	return 0;
}