#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

// n=1000 の場合、1001x1001のテーブルを使用
#define MAX_N 1001

int main() {
    int max_n = 1000;
    
    // dp[y][v] : 現在のステップ x において、高さ y、残り容量 v である状態の数
    // 巨大な数になるため mpz_t (GMPの整数型) を使用
    mpz_t (*dp)[MAX_N] = malloc(sizeof(mpz_t[MAX_N][MAX_N]));
    mpz_t (*next_dp)[MAX_N] = malloc(sizeof(mpz_t[MAX_N][MAX_N]));
    mpz_t *results = malloc(sizeof(mpz_t[MAX_N]));

    for (int i = 0; i < MAX_N; i++) {
        mpz_init(results[i]);
        for (int j = 0; j < MAX_N; j++) {
            mpz_init(dp[i][j]);
            mpz_init(next_dp[i][j]);
        }
    }

    // ベースケース: x=0, y=0, twoar=0 のとき 1
    mpz_set_ui(dp[0][0], 1);
    mpz_add(results[0], results[0], dp[0][0]);

    // x を 1 から max_n まで回す
    for (int x = 1; x <= max_n; x++) {
        // 次のステップのテーブルをクリア
        for (int y = 0; y < MAX_N; y++) {
            for (int v = 0; v < MAX_N; v++) mpz_set_ui(next_dp[y][v], 0);
        }

        for (int y = 0; y <= x; y++) {
            for (int v = 0; v <= max_n; v++) {
                if (mpz_cmp_ui(dp[y][v], 0) == 0) continue;

                // Mapleのロジックを逆方向に適用（x-1 から x へ）
                // 1. y=0 の場合 (Mapleの y=0 かつ x%2==0 の条件に相当)
                if (y == 0 && x % 2 == 0) {
                    // Motzk(x, 0, v) = Motzk(x-1, 1, v-1) の関係
                    // ここでは「逆」に、dp[1][v] が next_dp[0][v+1] に寄与すると考える
                }

                // --- 順方向のDP遷移 ---
                // y, v の状態から、次の x+1 ステップでどの y', v' に行けるか
            }
        }
        
        // ※ 実際の遷移ロジックは非常に複雑なため、
        // Mapleの「Motzk(x, y, v) = ...」をそのままループで回す形式が確実です
        for (int y = 0; y <= x; y++) {
            for (int v = 0; v <= max_n; v++) {
                // y == 0 のとき
                if (y == 0) {
                    if (x % 2 == 0 && v >= 1) 
                        mpz_set(next_dp[0][v], dp[1][v-1]);
                } 
                // y > 0 かつ y%2 == x%2 のとき
                else if (y % 2 == x % 2) {
                    // y+1, y, y-1 からの遷移をまとめる
                    if (y + 1 < MAX_N && v >= 2*y + 1)
                        mpz_add(next_dp[y][v], next_dp[y][v], dp[y+1][v-2*y-1]);
                    if (v >= 2*y)
                        mpz_add(next_dp[y][v], next_dp[y][v], dp[y][v-2*y]);
                    if (y - 1 >= 0 && v >= 2*y - 1)
                        mpz_add(next_dp[y][v], next_dp[y][v], dp[y-1][v-2*y+1]);
                } 
                // それ以外
                else {
                    if (v >= 2*y)
                        mpz_set(next_dp[y][v], dp[y][v-2*y]);
                }
                
                // A316585(v) は全x, yについての Motzk(x, y, v) の和
                mpz_add(results[v], results[v], next_dp[y][v]);
            }
        }

        // dp テーブルの更新 (スワップ)
        mpz_t (*temp)[MAX_N] = dp;
        dp = next_dp;
        next_dp = temp;
    }

    // 結果出力
    FILE *fp = fopen("b316585.txt", "w");
    for (int n = 0; n <= max_n; n++) {
        fprintf(fp, "%d ", n);
        mpz_out_str(fp, 10, results[n]);
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 0;
}