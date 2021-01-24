# Using graphillion
from graphillion import GraphSet
import graphillion.tutorial as tl

def A116469(n, k):
    if n == 1 or k == 1: return 1
    universe = tl.grid(n-1, k-1)
    GraphSet.set_universe(universe)
    spanning_trees = GraphSet.trees(is_spanning=True)
    return spanning_trees.len()

def A334002(n):
    return A116469(n, 7)

# print([A334002(n) for n in range(1, 15)])

cnt = 1
for i in range(1, 201):
    print(str(cnt) + ' ' + str(A334002(i)))
    cnt += 1

