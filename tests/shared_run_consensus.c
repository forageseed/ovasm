/* Reproduce disagreeing run estimates and duplicate read-position weighting. */
#include <assert.h>
#include <string.h>
#include "../vendor/syncasm/syncasm.c"

static void packed_read(sr_t *read, const char *hp, uint32_t extra_position, uint32_t extra)
{
    read->hoco_l = strlen(hp);
    read->hoco_s = calloc((read->hoco_l + 3) / 4, 1);
    read->ho_rl = calloc(read->hoco_l, 1);
    for (uint32_t i = 0; i < read->hoco_l; ++i)
        read->hoco_s[i/4] |= seq_nt4_table[(unsigned char) hp[i]] << (((i & 3) ^ 3) << 1);
    read->ho_rl[extra_position] = extra < 255 ? extra : 255;
    if (extra >= 255) {
        read->ho_l_rl = malloc(sizeof(uint32_t));
        read->ho_l_rl[0] = extra;
    }
}

static void check(uint32_t a, uint32_t b, uint32_t c, uint32_t expected, int reverse_read, int reverse_unitigs, int branch)
{
    sr_t read[4] = {{0}};
    packed_read(&read[0], reverse_read ? "CGTACGT" : "ACGTACG", reverse_read ? 4 : 2, a);
    packed_read(&read[1], reverse_read ? "TACGT" : "ACGTA", 2, b);
    packed_read(&read[2], reverse_read ? "CGTAC" : "GTACG", reverse_read ? 4 : 0, c);
    packed_read(&read[3], reverse_read ? "CGTACGT" : "ACGTACG", reverse_read ? 4 : 2, 255);
    uint64_t kmers0[] = {reverse_read ? 2 : 0, reverse_read ? 0 : 2};
    uint64_t kmers1[] = {0}, kmers2[] = {2};
    uint64_t corrected[] = {reverse_read ? 3 : 1, reverse_read ? 1 : 3};
    uint32_t positions0[] = {reverse_read ? 1 : 0, reverse_read ? 5 : 4};
    uint32_t positions1[] = {reverse_read}, positions2[] = {reverse_read};
    read[0].k_mer = kmers0; read[0].m_pos = positions0; read[0].n = 2;
    read[1].k_mer = kmers1; read[1].m_pos = positions1; read[1].n = 1;
    read[2].k_mer = kmers2; read[2].m_pos = positions2; read[2].n = 1;
    read[3].k_mer = corrected; read[3].m_pos = positions0; read[3].n = 2;
    sr_db_t reads = {.a = read, .n = 4, .k = 5};
    uint64_t observations0[] = {reverse_read ? 3 : 0, (1ULL << 32) | reverse_read, (3ULL << 32) | (reverse_read ? 3 : 0)};
    uint64_t observations1[] = {reverse_read ? 1 : 2, (2ULL << 32) | reverse_read, (3ULL << 32) | (reverse_read ? 1 : 2)};
    syncmer_t nodes[] = {{.cov = 3, .m_pos = observations0}, {.cov = 3, .m_pos = observations1}};
    syncmer_db_t database = {.a = nodes, .n = 2};
    uint64_t mers0[] = {reverse_unitigs}, mers1[] = {2 | reverse_unitigs};
    asmg_vtx_t vertices[] = {{.a = mers0, .n = 1, .cov = 2}, {.a = mers1, .n = 1, .cov = 2},
                             {.a = mers1, .n = 1, .cov = 2}};
    asmg_arc_t arcs[] = {{.v = reverse_unitigs, .w = 2 | reverse_unitigs, .cov = 1},
                        {.v = reverse_unitigs, .w = 4 | reverse_unitigs, .cov = 1, .del = !branch},
                        {.v = 3 ^ reverse_unitigs, .w = 1 ^ reverse_unitigs, .cov = 1, .comp = 1},
                        {.v = 5 ^ reverse_unitigs, .w = 1 ^ reverse_unitigs, .cov = 1, .comp = 1, .del = !branch}};
    uint64_t index[6] = {0}, number[6] = {0};
    index[arcs[0].v] = 0; index[arcs[2].v] = 2; index[arcs[3].v] = 3;
    number[arcs[0].v] = 2; number[arcs[2].v] = number[arcs[3].v] = 1;
    asmg_t graph = {.vtx = vertices, .n_vtx = branch ? 3 : 2, .arc = arcs, .n_arc = 4, .idx_p = index, .idx_n = number};
    scg_t scg = {.scm_db = &database, .utg_asmg = &graph};
    scg_consensus(&reads, &scg, 0, 1, NULL);
    assert(arcs[0].ls == expected + 3 && arcs[2].ls == arcs[0].ls);
    if (branch) {
        assert(arcs[1].ls == arcs[0].ls && arcs[3].ls == arcs[0].ls);
        assert(!strcmp(vertices[1].seq, vertices[2].seq));
    }
    assert(vertices[0].len == expected + 5 && vertices[1].len == expected + 5);
    if (!reverse_unitigs) {
        assert(vertices[0].seq[0] == 'A' && vertices[0].seq[1] == 'C');
        assert(vertices[1].seq[expected + 3] == 'C' && vertices[1].seq[expected + 4] == 'G');
        assert(!memcmp(vertices[0].seq + 2, vertices[1].seq, arcs[0].ls));
    } else {
        assert(!memcmp(vertices[1].seq + 2, vertices[0].seq, arcs[0].ls));
    }
    for (int i = 0; i < (branch ? 3 : 2); ++i) free(vertices[i].seq);
    for (int i = 0; i < 4; ++i) {
        free(read[i].hoco_s); free(read[i].ho_rl); free(read[i].ho_l_rl);
    }
}

int main(void)
{
    for (int read_strand = 0; read_strand < 2; ++read_strand)
        for (int graph_strand = 0; graph_strand < 2; ++graph_strand)
          for (int branch = 0; branch < 2; ++branch) {
            check(0, 0, 3, 1, read_strand, graph_strand, branch); /* independent means 1 vs 3; pooled run 2 */
            check(3, 0, 0, 1, read_strand, graph_strand, branch); /* counting the shared read twice would give run 3 */
            check(255, 0, 0, 85, read_strand, graph_strand, branch); /* long runs use the side array */
        }
    return 0;
}
