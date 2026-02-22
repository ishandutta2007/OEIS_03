#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <gmp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr int MAX_N = 38;
constexpr int MAX_PRIMES = 6;

constexpr uint64_t OLD_PRIMES[6] = {
    2305843009213693951ULL,
    2305843009213693921ULL,
    2305843009213693907ULL,
    2305843009213693723ULL,
    2305843009213693693ULL,
    2305843009213693669ULL
};

inline uint64_t add_mod(uint64_t a, uint64_t b, uint64_t p) {
    uint64_t s = a + b;
    return (s >= p) ? (s - p) : s;
}

inline uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t p) {
    return (a >= b) ? (a - b) : (a + (p - b));
}

inline uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t p) {
    return static_cast<uint64_t>((__uint128_t)a * (__uint128_t)b % p);
}

uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t p) {
    uint64_t r = 1;
    while (e) {
        if (e & 1ULL) r = mul_mod(r, a, p);
        a = mul_mod(a, a, p);
        e >>= 1ULL;
    }
    return r;
}

inline uint64_t mod_inv_prime(uint64_t a, uint64_t p) {
    return pow_mod(a, p - 2, p);
}

inline uint64_t gcd_u64(uint64_t a, uint64_t b) {
    while (b) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

inline uint64_t lcm_u64(uint64_t a, uint64_t b) {
    return (a / gcd_u64(a, b)) * b;
}

inline uint64_t i64_to_mod(int64_t x, uint64_t p) {
    if (x >= 0) return static_cast<uint64_t>(x) % p;
    uint64_t t = static_cast<uint64_t>(-x) % p;
    return (t == 0) ? 0 : (p - t);
}

bool is_old_prime(uint64_t p) {
    for (uint64_t q : OLD_PRIMES) {
        if (p == q) return true;
    }
    return false;
}

bool is_prime_u64(uint64_t n) {
    if (n < 2) return false;
    static const uint64_t small[] = {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL};
    for (uint64_t p : small) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    uint64_t d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    static const uint64_t bases[] = {
        2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL
    };
    for (uint64_t a0 : bases) {
        uint64_t a = a0 % n;
        if (a == 0) continue;
        uint64_t x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
                comp = false;
                break;
            }
        }
        if (comp) return false;
    }
    return true;
}

std::vector<uint64_t> generate_primes(int need) {
    std::vector<uint64_t> ps;
    ps.reserve(need);
    uint64_t cand = ((1ULL << 61) - 1ULL) - 2ULL;
    if ((cand & 1ULL) == 0ULL) --cand;

    while (static_cast<int>(ps.size()) < need) {
        if (!is_old_prime(cand) && is_prime_u64(cand)) ps.push_back(cand);
        cand -= 2ULL;
    }
    return ps;
}

void build_row_masks_and_coldeg(int n, std::vector<uint64_t>& row_masks, std::vector<int16_t>& coldeg) {
    row_masks.assign(n, 0ULL);
    coldeg.assign(n, 0);
    for (int i = 0; i < n; ++i) {
        uint64_t mask = 0ULL;
        const uint64_t u = static_cast<uint64_t>(i + 1);
        for (int j = 0; j < n; ++j) {
            const uint64_t v = static_cast<uint64_t>(j + 1);
            if (lcm_u64(u, v) <= static_cast<uint64_t>(n)) {
                mask |= (1ULL << j);
                ++coldeg[j];
            }
        }
        row_masks[i] = mask;
    }
}

double log2_factorial(int n) {
    static bool inited = false;
    static double pref[65];
    if (!inited) {
        pref[0] = 0.0;
        for (int i = 1; i <= 64; ++i) pref[i] = pref[i - 1] + std::log2(static_cast<double>(i));
        inited = true;
    }
    if (n <= 64) return pref[n];
    double s = pref[64];
    for (int i = 65; i <= n; ++i) s += std::log2(static_cast<double>(i));
    return s;
}

double log2_bregman_upper_bound(const std::vector<uint64_t>& row_masks) {
    double bits = 0.0;
    for (uint64_t m : row_masks) {
        int r = __builtin_popcountll(m);
        if (r == 0) return 0.0;
        bits += log2_factorial(r) / static_cast<double>(r);
    }
    return bits;
}

int choose_prime_count_for_bits(double bits_upper) {
    const int need = static_cast<int>(std::ceil(bits_upper)) + 16;
    int got = 0;
    for (int i = 0; i < MAX_PRIMES; ++i) {
        got += 61;
        if (got >= need) return i + 1;
    }
    return -1;
}

uint64_t permanent_mod_glynn_parallel(const std::vector<uint64_t>& row_masks,
                                      const std::vector<int16_t>& coldeg,
                                      int n, uint64_t p) {
    if (n == 0) return 1ULL;
    if (n == 1) return (row_masks[0] & 1ULL) ? 1ULL : 0ULL;

    const int vars = n - 1;
    const uint64_t total = 1ULL << vars;
    const int offset = n;
    const int span = 2 * n + 1;

    std::vector<uint64_t> colsum_mod(span);
    for (int v = -n; v <= n; ++v) colsum_mod[v + offset] = i64_to_mod(v, p);

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    int chunks = nthreads * 8;
    if (static_cast<uint64_t>(chunks) > total) chunks = static_cast<int>(total);
    if (chunks < 1) chunks = 1;

    std::vector<uint64_t> chunk_acc(chunks, 0ULL);

#pragma omp parallel for schedule(static) if(chunks > 1)
    for (int c = 0; c < chunks; ++c) {
        const uint64_t start = (static_cast<uint64_t>(c) * total) / static_cast<uint64_t>(chunks);
        const uint64_t end = (static_cast<uint64_t>(c + 1) * total) / static_cast<uint64_t>(chunks);
        if (start >= end) {
            chunk_acc[c] = 0ULL;
            continue;
        }

        std::vector<int8_t> delta(n, 1);
        std::vector<int16_t> colsum = coldeg;

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
        uint64_t acc = 0ULL;

        for (uint64_t t = start; t < end; ++t) {
            uint64_t prod = 1ULL;
            for (int j = 0; j < n; ++j) {
                prod = mul_mod(prod, colsum_mod[colsum[j] + offset], p);
                if (prod == 0ULL) break;
            }

            if (!odd) {
                acc = add_mod(acc, prod, p);
            } else {
                acc = sub_mod(acc, prod, p);
            }

            if (t + 1 < end) {
                int b = __builtin_ctzll(t + 1);
                int row = b + 1;
                delta[row] = static_cast<int8_t>(-delta[row]);
                int16_t change = static_cast<int16_t>(2 * static_cast<int>(delta[row]));
                uint64_t m = row_masks[row];
                while (m) {
                    int j = __builtin_ctzll(m);
                    colsum[j] = static_cast<int16_t>(colsum[j] + change);
                    m &= (m - 1);
                }
                odd ^= 1;
            }
        }

        chunk_acc[c] = acc;
    }

    uint64_t sum = 0ULL;
    for (int c = 0; c < chunks; ++c) sum = add_mod(sum, chunk_acc[c], p);
    uint64_t scale = mod_inv_prime(pow_mod(2ULL, static_cast<uint64_t>(n - 1), p), p);
    return mul_mod(sum, scale, p);
}

void crt_reconstruct(const std::vector<uint64_t>& residues,
                     const std::vector<uint64_t>& primes,
                     mpz_t out) {
    mpz_t x, mod;
    mpz_init_set_ui(x, 0);
    mpz_init_set_ui(mod, 1);

    for (size_t i = 0; i < residues.size(); ++i) {
        uint64_t p = primes[i];
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

}  // namespace

int main() {
    std::FILE* out = std::fopen("b354756_1.txt", "w");
    if (!out) {
        std::fprintf(stderr, "failed to open output file: b354756_1.txt\n");
        return 1;
    }

    std::fprintf(out, "0 1\n");
    std::fflush(out);

    mpz_t ans;
    mpz_init(ans);

    const std::vector<uint64_t> prime_pool = generate_primes(MAX_PRIMES);

    std::vector<uint64_t> row_masks;
    std::vector<int16_t> coldeg;

    for (int n = 1; n <= MAX_N; ++n) {
        build_row_masks_and_coldeg(n, row_masks, coldeg);
        double bits_upper = log2_bregman_upper_bound(row_masks);
        int pcnt = choose_prime_count_for_bits(bits_upper);
        if (pcnt < 0) {
            std::fprintf(stderr, "not enough primes for n=%d\n", n);
            mpz_clear(ans);
            std::fclose(out);
            return 1;
        }

        std::vector<uint64_t> residues(pcnt);
        for (int i = 0; i < pcnt; ++i) {
            residues[i] = permanent_mod_glynn_parallel(row_masks, coldeg, n, prime_pool[i]);
        }

        std::vector<uint64_t> primes(prime_pool.begin(), prime_pool.begin() + pcnt);
        crt_reconstruct(residues, primes, ans);
        char* s = mpz_get_str(nullptr, 10, ans);
        if (!s) {
            std::fprintf(stderr, "mpz_get_str failed\n");
            mpz_clear(ans);
            std::fclose(out);
            return 1;
        }
        std::fprintf(out, "%d %s\n", n, s);
        std::free(s);
        std::fflush(out);
    }

    mpz_clear(ans);
    std::fclose(out);
    return 0;
}
