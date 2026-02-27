#!/usr/bin/env python3
"""
平らな面に三角形を積み上げるときの「左右対称な積み方」の個数 a(n) を数える。

考え方（実装で使うモデル）:
- 積み方を small Schroeder path で表す。
- n は「面積」に対応し、a(n) は面積 n の対称パス数になる。
- 左右対称パスは
    1) L + mirror(reverse(L))
  または
    2) L + F + mirror(reverse(L))   (中央に水平辺 F を1本入れる)
  の形に分解できる。

この性質を使い、まず「左半分 L」のみを DP で数えてから、
上の 1) 2) を合計して a(n) を得る。
"""

from __future__ import annotations

import argparse
from collections import defaultdict


def symmetric_stack_counts(n_max: int) -> list[int]:
    # dp[area2][h] = 左半分 L の個数
    #   area2: 面積の 2 倍値（整数で管理しやすいので 2 倍して持つ）
    #   h    : L の終点の高さ
    dp: list[defaultdict[int, int]] = [defaultdict(int) for _ in range(n_max + 1)]
    # 空パス（面積0, 高さ0）を 1 通りとして開始
    dp[0][0] = 1

    for area2 in range(n_max + 1):
        for h, ways in list(dp[area2].items()):
            # U: 右上がり 1 ステップ
            # 面積2倍の増分は 2*h + 1
            na = area2 + (2 * h + 1)
            if na <= n_max:
                dp[na][h + 1] += ways

            if h > 0:
                # D: 右下がり 1 ステップ（高さ 0 未満には行けない）
                # 面積2倍の増分は 2*h - 1
                na = area2 + (2 * h - 1)
                if na <= n_max:
                    dp[na][h - 1] += ways

                # F: 水平 2 ステップ（small Schroeder の規則で高さ>0 のみ）
                # 面積2倍の増分は 4*h
                na = area2 + 4 * h
                if na <= n_max:
                    dp[na][h] += ways

    ans = [0] * (n_max + 1)

    for n in range(n_max + 1):
        total = 0

        # ケース1: 中央に F を入れない対称形
        # full = L + mirror(reverse(L))
        # 全体面積 n はそのまま area2(L)=n に一致
        total += sum(dp[n].values())

        # ケース2: 中央に F を1本入れる対称形
        # full = L + F + mirror(reverse(L))
        # 全体面積 n = area2(L) + 2*h になるので、その条件で加算
        for a in range(n + 1):
            need = n - a
            for h, ways in dp[a].items():
                if h > 0 and 2 * h == need:
                    total += ways

        ans[n] = total

    return ans


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "n",
        type=int,
        nargs="?",
        default=5000,
        help="a(0) から a(n) までを出力 (default: 5000)",
    )
    args = parser.parse_args()

    seq = symmetric_stack_counts(args.n)
    for i, v in enumerate(seq):
        print(i, v)


if __name__ == "__main__":
    main()
