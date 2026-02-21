#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define HASH_SIZE (1u << 22)

typedef struct {
    int idx;
    uint64_t mask;
    mpz_t val;
    unsigned int gen;
    char initialized;
} Cache;

static Cache *cache;
static unsigned int current_gen = 1;
static int adj[64][64];
static int adj_len[64];
static int row_order[64];

static inline uint64_t get_hash(int idx, uint64_t mask) {
    uint64_t x = mask ^ ((uint64_t)(unsigned)idx * 0x9e3779b97f4a7c15ULL);
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x & (HASH_SIZE - 1);
}

static void next_generation(void) {
    current_gen++;
    if (current_gen == 0) {
        for (uint64_t i = 0; i < HASH_SIZE; i++) cache[i].gen = 0;
        current_gen = 1;
    }
}

static int cache_get(int idx, uint64_t mask, mpz_t out) {
    uint64_t h = get_hash(idx, mask);
    for (uint64_t step = 0; step < HASH_SIZE; step++) {
        Cache *slot = &cache[(h + step) & (HASH_SIZE - 1)];
        if (slot->gen != current_gen) return 0;
        if (slot->idx == idx && slot->mask == mask) {
            mpz_set(out, slot->val);
            return 1;
        }
    }
    return 0;
}

static void cache_put(int idx, uint64_t mask, const mpz_t val) {
    uint64_t h = get_hash(idx, mask);
    for (uint64_t step = 0; step < HASH_SIZE; step++) {
        Cache *slot = &cache[(h + step) & (HASH_SIZE - 1)];
        if (slot->gen != current_gen || (slot->idx == idx && slot->mask == mask)) {
            if (!slot->initialized) {
                mpz_init(slot->val);
                slot->initialized = 1;
            }
            slot->idx = idx;
            slot->mask = mask;
            slot->gen = current_gen;
            mpz_set(slot->val, val);
            return;
        }
    }
}

static void count_matching(int idx, uint64_t mask, int n, mpz_t res) {
    if (idx == n) {
        mpz_set_ui(res, 1);
        return;
    }

    if (cache_get(idx, mask, res)) return;

    mpz_set_ui(res, 0);
    mpz_t temp;
    mpz_init(temp);

    int u = row_order[idx];
    for (int i = 0; i < adj_len[u]; i++) {
        int v = adj[u][i];
        if (!(mask & (1ULL << v))) {
            count_matching(idx + 1, mask | (1ULL << v), n, temp);
            mpz_add(res, res, temp);
        }
    }

    cache_put(idx, mask, res);
    mpz_clear(temp);
}

static void prepare(int n) {
    int deg[64] = {0};
    for (int i = 0; i < n; i++) {
        adj_len[i] = 0;
        row_order[i] = i;
        for (int j = 0; j < n; j++) {
            int u = i + 1;
            int v = j + 1;
            if (u % v == 0 || v % u == 0) {
                adj[i][adj_len[i]++] = j;
                deg[i]++;
            }
        }
    }
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (deg[row_order[i]] > deg[row_order[j]]) {
                int t = row_order[i];
                row_order[i] = row_order[j];
                row_order[j] = t;
            }
        }
    }
}

int main(int argc, char **argv) {
    int nmax = 32;
    if (argc >= 2) {
        nmax = atoi(argv[1]);
        if (nmax < 0 || nmax > 32) {
            fprintf(stderr, "usage: %s [nmax<=32]\n", argv[0]);
            return 1;
        }
    }

    cache = (Cache *)calloc(HASH_SIZE, sizeof(Cache));
    if (!cache) {
        fprintf(stderr, "cache allocation failed\n");
        return 1;
    }

    FILE *out = fopen("b320843.txt", "w");
    if (!out) {
        fprintf(stderr, "failed to open output file: b320843.txt\n");
        free(cache);
        return 1;
    }

    fprintf(out, "0 1\n");
    fflush(out);

    for (int n = 1; n <= nmax; n++) {
        next_generation();
        prepare(n);
        mpz_t res;
        mpz_init(res);
        count_matching(0, 0ULL, n, res);
        char *s = mpz_get_str(NULL, 10, res);
        if (!s) {
            fprintf(stderr, "mpz_get_str failed\n");
            mpz_clear(res);
            fclose(out);
            free(cache);
            return 1;
        }
        fprintf(out, "%d %s\n", n, s);
        free(s);
        mpz_clear(res);
        fflush(out);
    }

    for (uint64_t i = 0; i < HASH_SIZE; i++) {
        if (cache[i].initialized) mpz_clear(cache[i].val);
    }
    fclose(out);
    free(cache);
    return 0;
}
