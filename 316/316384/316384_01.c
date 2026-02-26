#include <stdio.h>
#include <string.h>

#define MAXLEN 80

// a.reverse == f(a) を検証（早期打ち切りあり）
static int check(const int *a, int alen) {
    int cur[MAXLEN], nxt[MAXLEN];
    int curlen = alen;
    memcpy(cur, a, alen * sizeof(int));

    for (int fidx = 0; ; fidx++) {
        // 全ゼロ判定
        int all_zero = 1;
        for (int i = 0; i < curlen; i++) {
            if (cur[i] != 0) { all_zero = 0; break; }
        }
        if (all_zero) return fidx == alen;  // f の長さが a と一致するか
        if (fidx >= alen) return 0;         // f が a より長い

        int s = 0, nxtlen = 0;
        for (int i = 0; i < curlen; i++) {
            if (cur[i] == 1) {
                nxt[nxtlen++] = 0;
                s += 1;
            } else if (cur[i] > 1) {
                nxt[nxtlen++] = cur[i] - 2;
                s += 2;
            }
            // cur[i] == 0 は無視（次の配列に含めない）
        }

        // f_ary[fidx] == a[alen-1-fidx] でなければ不一致
        if (s != a[alen - 1 - fidx]) return 0;

        memcpy(cur, nxt, nxtlen * sizeof(int));
        curlen = nxtlen;
    }
}

static long long count;
static int seq[MAXLEN];

// g(n, m, depth): seq[0..depth-1] が現在の部分列, m は現在の和
static void g(int n, int m, int depth) {
    if (m >= n) {
        if (m == n && check(seq, depth)) count++;
        return;
    }
    int last = depth > 0 ? seq[depth - 1] : 0;
    // 次の要素の最小値: 直前が奇数なら +1, 偶数なら +0（A006950 と同じ制約）
    int s = last > 0 ? last + (last & 1) : 1;
    int rem = n - m;
    for (int x = s; x <= rem; x++) {
        seq[depth] = x;
        g(n, m + x, depth + 1);
    }
}

static long long A316384(int n) {
    if (n == 0) return 1;  // 空列 [] は条件を満たす
    count = 0;
    g(n, 0, 0);
    return count;
}

int main() {
    for (int i = 0; i <= 150; i++) {
        printf("%d %lld\n", i, A316384(i));
    }
    return 0;
}
