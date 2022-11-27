# Using graphillion

from graphillion import GraphSet
import graphillion.tutorial as tl

def A333758(n, k):
    universe = tl.grid(n - 1, k - 1)
    GraphSet.set_universe(universe)
    cycles = GraphSet.cycles()
    points = [i for i in range(1, k * n + 1) if i % k < 2 or ((i - 1) // k + 1) % n < 2]
    for i in points:
        cycles = cycles.including(i)
    return cycles.len()

def A358697(n):
    return A333758(6, n)

print([A358697(n) for n in range(2, 14)])
for n in range(2, 18):
    print(A358697(n))