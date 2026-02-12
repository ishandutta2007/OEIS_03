# Using graphillion
from graphillion import GraphSet
import graphillion.tutorial as tl

def make_cylinder_graph(n):
    universe = list(tl.grid(3, n - 1))
    for i in range(1, n + 1):
        u = i
        v = 3 * n + i
        universe.append((u, v))
    return universe

def A390791(n):
    universe = make_cylinder_graph(2 * n + 2)
    GraphSet.set_universe(universe)
    start, goal = 1, 3 * (2 * n + 2)
    paths = GraphSet.paths(start, goal, is_hamilton=True)
    return paths.len()

print([A390791(n) for n in range(0, 18)])
# for n in range(1, 43):
#    print(A(n))