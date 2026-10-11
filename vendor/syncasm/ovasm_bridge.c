/* OrganelleVerse MIT-licensed binding. Upstream files retain their original licenses. */
#include <stddef.h>
#include "misc.h"
#include "syncasm.h"

int VERBOSE = 0;
int syncasm(char **file_in, int n_file, size_t m_data, int k, int s,
            int bubble_size, int tip_size, int min_k_cov, double min_a_cov_f,
            double weak_cross, int do_ec, int do_unzip, int n_threads,
            char *out, scg_meta_t *meta, int verbose);

int ov_syncasm(const char *const *reads, int n, int k, int s, int cutoff,
               int threads, const char *prefix, uint64_t *read_count, uint64_t *base_count)
{
    scg_meta_t meta = {0};
    *read_count = *base_count = 0;
    sys_init();
    int result = syncasm((char **) reads, n, 0, k, s, 100000, 10000, cutoff,
                        0.35, 0.3, 1, 3, threads, (char *) prefix, &meta, 0);
    if (result == 0) {
        *read_count = meta.sr_db->n;
        uint64_t bases = 0;
        /* The reader records raw length before HPC or syncmer correction.
         * Include reads without a syncmer and every ambiguous base. */
        for (size_t i = 0; i < meta.sr_db->n; ++i) {
            bases += meta.sr_db->a[i].original_length;
        }
        *base_count = bases;
    }
    scg_meta_clean(&meta);
    return result;
}
