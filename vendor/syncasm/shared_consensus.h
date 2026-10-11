/* OVASM extension to the MIT-licensed Syncasm kernel.
 * Restore shared HPC graph positions with one consensus of distinct original
 * read positions. Independent syncmer consensuses otherwise disagree in run
 * length and cannot legitimately be exported as exact base-coordinate links.
 * Included after calc_syncmer_overlap; not a separate compilation unit. */

typedef struct {
    uint64_t next, unitig, graph_start;
    uint32_t read_start, reverse;
} shared_run_placement_t;

static uint64_t shared_run_root(uint64_t *parent, uint64_t x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static uint64_t shared_run_coordinate(uint64_t offset, uint64_t length, uint64_t oriented, uint64_t p)
{
    return offset + ((oriented & 1) ? length - 1 - p : p);
}

static char shared_run_base(char base, uint64_t oriented)
{
    if (!(oriented & 1)) return base;
    switch (base) {
        case 'A': return 'T'; case 'C': return 'G';
        case 'G': return 'C'; case 'T': return 'A';
        default: return base;
    }
}

/* On entry the arc lengths are in HPC coordinates, and vertex sequences are
 * the existing per-unitig base consensus. Only shared coordinates are pooled;
 * interior positions retain that consensus. No edge or observation is removed. */
static void scg_shared_consensus(sr_db_t *reads, scg_t *scg)
{
    asmg_t *g = scg->utg_asmg;
    syncmer_t *nodes = scg_a_scm(scg);
    kstring_t *hp;
    uint64_t *offset, *parent, *sum, *count, *runs;
    uint8_t *shared;
    MYCALLOC(hp, g->n_vtx);
    MYCALLOC(offset, g->n_vtx + 1);
    for (uint64_t u = 0; u < g->n_vtx; ++u) {
        asmg_vtx_t *v = &g->vtx[u];
        if (!v->del)
            scg_unitig_consensus(reads, v->a, v->n, nodes, &hp[u], 1);
        offset[u + 1] = offset[u] + hp[u].l;
    }
    uint64_t n = offset[g->n_vtx];
    MYMALLOC(parent, n);
    MYCALLOC(shared, n);
    MYCALLOC(sum, n);
    MYCALLOC(count, n);
    MYMALLOC(runs, n);
    for (uint64_t p = 0; p < n; ++p) parent[p] = p;

    for (uint64_t i = 0; i < g->n_arc; ++i) {
        asmg_arc_t *a = &g->arc[i];
        if (a->del || a->comp) continue;
        uint64_t u = a->v >> 1, v = a->w >> 1, l = a->ls;
        assert(l <= hp[u].l && l <= hp[v].l);
        for (uint64_t j = 0; j < l; ++j) {
            uint64_t x = shared_run_coordinate(offset[u], hp[u].l, a->v, hp[u].l - l + j);
            uint64_t y = shared_run_coordinate(offset[v], hp[v].l, a->w, j);
            assert(shared_run_base(hp[u].s[x-offset[u]], a->v) ==
                   shared_run_base(hp[v].s[y-offset[v]], a->w));
            shared[x] = shared[y] = 1;
            uint64_t rx = shared_run_root(parent, x), ry = shared_run_root(parent, y);
            if (rx < ry) parent[ry] = rx;
            else parent[rx] = ry;
        }
    }

    /* Bucket syncmer placements by original read. A duplicate observation can
     * only belong to that same read, so distinct-position bookkeeping need not
     * retain every base observation from the entire input at once. */
    uint64_t *heads;
    MYMALLOC(heads, reads->n);
    for (uint64_t sid = 0; sid < reads->n; ++sid) heads[sid] = UINT64_MAX;
    shared_run_placement_t *placements = NULL;
    size_t number = 0, capacity = 0;
    for (uint64_t u = 0; u < g->n_vtx; ++u) {
        asmg_vtx_t *v = &g->vtx[u];
        if (v->del) continue;
        /* Recover the existing run lengths without conflating padding N bases
         * with actual homopolymer positions. */
        uint64_t b = 0;
        uint64_t *prefix;
        MYCALLOC(prefix, hp[u].l + 1);
        for (uint64_t j = 0; j < hp[u].l; ++j) {
            char base = hp[u].s[j];
            assert(b < v->len && v->seq[b] == base);
            uint64_t start = b++;
            if (base != 'N') while (b < v->len && v->seq[b] == base) ++b;
            runs[offset[u] + j] = b - start;
            prefix[j + 1] = prefix[j] + shared[offset[u] + j];
        }
        assert(b == v->len);
        if (!prefix[hp[u].l]) { free(prefix); continue; }

        int64_t position = 0;
        for (uint64_t i = 0; i < v->n; ++i) {
            uint64_t node = v->a[i];
            if (i) position += calc_syncmer_overlap(reads, &nodes[v->a[i-1]>>1], v->a[i-1]&1,
                                                    &nodes[node>>1], node&1, NULL);
            assert(position >= 0 && position + reads->k <= hp[u].l);
            if (prefix[position + reads->k] == prefix[position]) continue;
            syncmer_t *mer = &nodes[node >> 1];
            for (uint32_t m = 0; m < mer->cov; ++m) {
                uint64_t sid = mer->m_pos[m] >> 32;
                sr_t *read = &reads->a[sid];
                uint32_t index = mer->m_pos[m] >> 1 & MAX_RD_SCM;
                if (read->k_mer[index] & 1) continue;
                uint32_t start = read->m_pos[index] >> 1;
                uint32_t reverse = (read->m_pos[index] & 1) ^ (node & 1);
                if (number == capacity) {
                    capacity = capacity ? capacity * 2 : 1024;
                    placements = realloc(placements, capacity * sizeof(*placements));
                    assert(placements);
                }
                placements[number] = (shared_run_placement_t) {
                    heads[sid], u, offset[u] + position, start, reverse
                };
                heads[sid] = number++;
            }
        }
        free(prefix);
    }
    /* Reuse Syncasm's integer-key map for the exact (graph position, read
     * position) pair. Clearing it between reads preserves the original
     * distinct-read-position estimator and bounds its lifetime to one read. */
    kh_u128_t *observations = kh_u128_init();
    for (uint64_t sid = 0; sid < reads->n; ++sid) {
        if (heads[sid] == UINT64_MAX) continue;
        kh_u128_m_clear(observations);
        sr_t *read = &reads->a[sid];
        for (uint64_t i = heads[sid]; i != UINT64_MAX; i = placements[i].next) {
            shared_run_placement_t *a = &placements[i];
            uint32_t large = 0;
            if (read->ho_l_rl)
                for (uint32_t j = 0; j < a->read_start; ++j) large += read->ho_rl[j] == 255;
            for (uint32_t j = 0; j < (uint32_t) reads->k; ++j) {
                uint32_t rp = a->read_start + j, extra = read->ho_rl[rp];
                if (extra == 255) extra = read->ho_l_rl[large++];
                uint64_t p = a->graph_start + (a->reverse ? reads->k - 1 - j : j);
                if (!shared[p]) continue;
                char base = char_nt4_table[(read->hoco_s[rp/4] >> (((rp & 3) ^ 3) << 1)) & 3];
                assert(shared_run_base(base, a->reverse) == hp[a->unitig].s[p-offset[a->unitig]]);
                uint64_t root = shared_run_root(parent, p);
                int absent;
                khint_t slot = kh_u128_put(observations, (uint128_t) root << 64 | rp, &absent);
                if (absent) {
                    kh_val(observations, slot) = extra;
                    sum[root] += extra;
                    ++count[root];
                } else assert(kh_val(observations, slot) == extra);
            }
        }
    }
    kh_u128_destroy(observations);
    free(placements); free(heads);
    for (uint64_t p = 0; p < n; ++p) {
        if (!shared[p]) continue;
        uint64_t root = shared_run_root(parent, p);
        assert(count[root] > 0);
        runs[p] = 1 + (uint64_t) lround((double) sum[root] / count[root]);
    }
    for (uint64_t u = 0; u < g->n_vtx; ++u) {
        asmg_vtx_t *v = &g->vtx[u];
        if (v->del) continue;
        kstring_t seq = {0, 0, 0};
        for (uint64_t j = 0; j < hp[u].l; ++j)
            for (uint64_t r = 0; r < runs[offset[u] + j]; ++r) kputc_(hp[u].s[j], &seq);
        kputc_('\0', &seq);
        free(v->seq);
        v->seq = seq.s;
        v->len = seq.l - 1;
    }
    for (uint64_t i = 0; i < g->n_arc; ++i) {
        asmg_arc_t *a = &g->arc[i];
        if (a->del || a->comp) continue;
        uint64_t u = a->v >> 1, l = a->ls, restored = 0;
        for (uint64_t j = 0; j < l; ++j)
            restored += runs[shared_run_coordinate(offset[u], hp[u].l, a->v, hp[u].l-l+j)];
        a->ls = restored;
        asmg_arc(g, a->w^1, a->v^1)->ls = restored;
    }
    for (uint64_t u = 0; u < g->n_vtx; ++u) free(hp[u].s);
    free(hp); free(offset); free(parent); free(shared); free(sum); free(count); free(runs);
}
