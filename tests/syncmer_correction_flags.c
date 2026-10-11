/* Regression for strand-independent original sequence support after correction. */
#include <assert.h>
#include "../vendor/syncasm/syncerr.c"

int main(void)
{
    uint64_t kmers[] = {0, 3, 4, 5}; /* original0, corrected1, original2, corrected2 */
    uint32_t positions[] = {1, 0, 2, 3}; /* reverse, forward, forward, reverse */
    sr_t read = {.sid = 0, .n = 4, .k_mer = kmers, .m_pos = positions};
    sr_db_t reads = {.n = 1, .a = &read};
    syncmer_t nodes[3] = {{0}};
    syncmer_db_t database = {.n = 3, .a = nodes};
    for (int strand = 0; strand < 2; ++strand) {
        update_syncmer_db(&reads, &database);
        assert(nodes[0].cov == 1 && !nodes[0].del);
        assert(nodes[1].cov == 1 && nodes[1].del);
        assert(nodes[2].cov == 2 && !nodes[2].del);
        for (int i = 0; i < 4; ++i) positions[i] ^= 1;
    }
    for (int i = 0; i < 3; ++i) free(nodes[i].m_pos);
    return 0;
}
