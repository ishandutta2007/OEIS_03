#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_PRIMES 6

static const uint64_t PRIMES[MAX_PRIMES] = {
    2305843009213693951ULL, /* 2^61 - 1 */
    2305843009213693921ULL,
    2305843009213693907ULL,
    2305843009213693723ULL,
    2305843009213693693ULL,
    2305843009213693669ULL
};

static inline uint64_t gcd_u64(uint64_t a, uint64_t b) {
    while (b) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static inline uint64_t lcm_u64(uint64_t a, uint64_t b) {
    if (a == 0 || b == 0) return 0;
    return (a / gcd_u64(a, b)) * b;
}

static inline uint64_t add_mod(uint64_t a, uint64_t b, uint64_t p) {
    uint64_t s = a + b;
    return (s >= p) ? (s - p) : s;
}

static inline uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t p) {
    return (a >= b) ? (a - b) : (a + (p - b));
}

static inline uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t p) {
    return (uint64_t)((__uint128_t)a * (__uint128_t)b % p);
}

static uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t p) {
    uint64_t r = 1;
    while (e) {
        if (e & 1ULL) r = mul_mod(r, a, p);
        a = mul_mod(a, a, p);
        e >>= 1ULL;
    }
    return r;
}

static inline uint64_t mod_inv_prime(uint64_t a, uint64_t p) {
    return pow_mod(a, p - 2, p);
}

static inline uint64_t i64_to_mod(int64_t x, uint64_t p) {
    if (x >= 0) return (uint64_t)x % p;
    uint64_t t = (uint64_t)(-x) % p;
    return (t == 0) ? 0 : (p - t);
}

static void build_row_masks_and_coldeg(int n, uint64_t *row_masks, int16_t *coldeg) {
    for (int j = 0; j < n; ++j) coldeg[j] = 0;
    for (int i = 0; i < n; ++i) {
        int u = i + 1;
        uint64_t mask = 0;
        for (int j = 0; j < n; ++j) {
            int v = j + 1;
            if (lcm_u64((uint64_t)u, (uint64_t)v) <= (uint64_t)n) {
                mask |= (1ULL << j);
                ++coldeg[j];
            }
        }
        row_masks[i] = mask;
    }
}

/* Glynn formula with delta_0 fixed to +1 and Gray code over remaining n-1 signs.
   per(A) = 2^(1-n) * sum_{delta_0=1, delta_i in {-1,+1}}
            (prod_i delta_i) * prod_j (sum_i delta_i a_{ij})
*/
static uint64_t permanent_mod_glynn_parallel(const uint64_t *row_masks,
                                             const int16_t *coldeg,
                                             int n, uint64_t p) {
    if (n == 0) return 1;
    if (n == 1) return (row_masks[0] & 1ULL) ? 1ULL : 0ULL;
    if (n > 63) {
        fprintf(stderr, "n=%d exceeds supported limit for bit masks\n", n);
        exit(1);
    }

    const int vars = n - 1;
    const uint64_t total = 1ULL << vars;
    const int colsum_offset = n;
    const int colsum_span = 2 * n + 1;

    uint64_t *colsum_mod = (uint64_t *)malloc((size_t)colsum_span * sizeof(uint64_t));
    if (!colsum_mod) {
        fprintf(stderr, "allocation failed\n");
        exit(1);
    }
    for (int v = -n; v <= n; ++v) colsum_mod[v + colsum_offset] = i64_to_mod((int64_t)v, p);

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    int chunks = nthreads * 8;
    if ((uint64_t)chunks > total) chunks = (int)total;
    if (chunks < 1) chunks = 1;

    uint64_t *chunk_acc = (uint64_t *)calloc((size_t)chunks, sizeof(uint64_t));
    if (!chunk_acc) {
        fprintf(stderr, "allocation failed\n");
        exit(1);
    }

#pragma omp parallel for schedule(static) if(chunks > 1)
    for (int c = 0; c < chunks; ++c) {
        uint64_t start = ((uint64_t)c * total) / (uint64_t)chunks;
        uint64_t end = ((uint64_t)(c + 1) * total) / (uint64_t)chunks;
        if (start >= end) {
            chunk_acc[c] = 0;
            continue;
        }

        int8_t *delta = (int8_t *)malloc((size_t)n);
        int16_t *colsum = (int16_t *)malloc((size_t)n * sizeof(int16_t));
        if (!delta || !colsum) {
            fprintf(stderr, "allocation failed\n");
            exit(1);
        }

        for (int i = 0; i < n; ++i) delta[i] = 1;
        for (int j = 0; j < n; ++j) colsum[j] = coldeg[j];

        uint64_t gray = start ^ (start >> 1);
        uint64_t tmp = gray;
        while (tmp) {
            int b = __builtin_ctzll(tmp);
            int row = b + 1;
            delta[row] = -1;
            uint64_t m = row_masks[row];
            while (m) {
                int j = __builtin_ctzll(m);
                colsum[j] -= 2;
                m &= (m - 1);
            }
            tmp &= (tmp - 1);
        }

        int odd = (__builtin_popcountll(gray) & 1U);
        uint64_t acc = 0;

        for (uint64_t t = start; t < end; ++t) {
            uint64_t prod = 1;
            for (int j = 0; j < n; ++j) {
                prod = mul_mod(prod, colsum_mod[colsum[j] + colsum_offset], p);
                if (prod == 0) break;
            }

            if (!odd) {
                acc = add_mod(acc, prod, p);
            } else {
                acc = sub_mod(acc, prod, p);
            }

            if (t + 1 < end) {
                int b = __builtin_ctzll(t + 1);
                int row = b + 1;
                delta[row] = (int8_t)-delta[row];
                int16_t change = (int16_t)(2 * (int)delta[row]);
                uint64_t m = row_masks[row];
                while (m) {
                    int j = __builtin_ctzll(m);
                    colsum[j] = (int16_t)(colsum[j] + change);
                    m &= (m - 1);
                }
                odd ^= 1;
            }
        }

        free(delta);
        free(colsum);
        chunk_acc[c] = acc;
    }

    uint64_t sum = 0;
    for (int c = 0; c < chunks; ++c) sum = add_mod(sum, chunk_acc[c], p);
    free(chunk_acc);
    free(colsum_mod);

    uint64_t scale = mod_inv_prime(pow_mod(2, (uint64_t)(n - 1), p), p);
    return mul_mod(sum, scale, p);
}

static double log2_factorial(int n) {
    static int inited = 0;
    static double pref[64];
    if (!inited) {
        pref[0] = 0.0;
        for (int i = 1; i < 64; ++i) pref[i] = pref[i - 1] + log2((double)i);
        inited = 1;
    }
    if (n < 0) return 0.0;
    if (n < 64) return pref[n];
    double s = pref[63];
    for (int i = 64; i <= n; ++i) s += log2((double)i);
    return s;
}

static double log2_bregman_upper_bound(const uint64_t *row_masks, int n) {
    double bits = 0.0;
    for (int i = 0; i < n; ++i) {
        int r = __builtin_popcountll(row_masks[i]);
        if (r == 0) return 0.0;
        bits += log2_factorial(r) / (double)r;
    }
    return bits;
}

static int choose_prime_count_for_bits(double bits_upper) {
    int need = (int)(bits_upper + 8.0);
    int got = 0;
    for (int i = 0; i < MAX_PRIMES; ++i) {
        got += 61;
        if (got >= need) return i + 1;
    }
    return -1;
}

static void crt_reconstruct(const uint64_t *residues, int pcnt, mpz_t out) {
    mpz_t x, mod;
    mpz_init_set_ui(x, 0);
    mpz_init_set_ui(mod, 1);

    for (int i = 0; i < pcnt; ++i) {
        uint64_t p = PRIMES[i];
        uint64_t x_mod_p = mpz_fdiv_ui(x, p);
        uint64_t mod_mod_p = mpz_fdiv_ui(mod, p);
        uint64_t rhs = (residues[i] >= x_mod_p) ? (residues[i] - x_mod_p)
                                                : (residues[i] + (p - x_mod_p));
        uint64_t inv = mod_inv_prime(mod_mod_p, p);
        uint64_t t = mul_mod(rhs, inv, p);
        mpz_addmul_ui(x, mod, t);
        mpz_mul_ui(mod, mod, p);
    }

    mpz_set(out, x);
    mpz_clear(x);
    mpz_clear(mod);
}

int main(void) {
    FILE *out = fopen("b354756.txt", "w");
    if (!out) {
        fprintf(stderr, "failed to open output file: b354756.txt\n");
        return 1;
    }

    fprintf(out, "0 1\n");
    fflush(out);

    mpz_t ans;
    mpz_init(ans);

    for (int n = 1; n <= 38; ++n) {
        uint64_t *row_masks = (uint64_t *)calloc((size_t)n, sizeof(uint64_t));
        int16_t *coldeg = (int16_t *)calloc((size_t)n, sizeof(int16_t));
        if (!row_masks || !coldeg) {
            fprintf(stderr, "allocation failed\n");
            mpz_clear(ans);
            fclose(out);
            return 1;
        }
        build_row_masks_and_coldeg(n, row_masks, coldeg);

        double bits_upper = log2_bregman_upper_bound(row_masks, n);
        int pcnt = choose_prime_count_for_bits(bits_upper);
        if (pcnt < 0) {
            fprintf(stderr, "not enough primes for n=%d\n", n);
            free(row_masks);
            free(coldeg);
            mpz_clear(ans);
            fclose(out);
            return 1;
        }

        uint64_t residues[MAX_PRIMES];
        for (int i = 0; i < pcnt; ++i) {
            residues[i] = permanent_mod_glynn_parallel(row_masks, coldeg, n, PRIMES[i]);
        }
        free(row_masks);
        free(coldeg);

        crt_reconstruct(residues, pcnt, ans);
        char *s = mpz_get_str(NULL, 10, ans);
        if (!s) {
            fprintf(stderr, "mpz_get_str failed\n");
            mpz_clear(ans);
            fclose(out);
            return 1;
        }
        fprintf(out, "%d %s\n", n, s);
        free(s);
        fflush(out);
    }

    mpz_clear(ans);
    fclose(out);
    return 0;
}
