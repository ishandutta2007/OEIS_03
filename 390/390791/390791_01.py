# Using graphillion
from graphillion import GraphSet
def make_CnXPk(n, k):
    grids = []
    for i in range(1, k + 1):
        for j in range(1, n):
            grids.append((i + (j - 1) * k, i + j * k))
        grids.append((i + (n - 1) * k, i))
    for i in range(1, k * n, k):
        for j in range(1, k):
            grids.append((i + j - 1, i + j))
    return grids

def A390791(n):
    universe = make_CnXPk(4, 2 * n + 2)
    GraphSet.set_universe(universe)
    start, goal = 1, 3 * (2 * n + 2)
    paths = GraphSet.paths(start, goal, is_hamilton=True)
    return paths.len()

# print([A390791(n) for n in range(0, 18)])
for n in range(0, 200 + 1):
   print(n, A390791(n))

# def A(start, goal, n, k):
#     universe = make_CnXPk(n, k)
#     GraphSet.set_universe(universe)
#     paths = GraphSet.paths(start, goal, is_hamilton=True)
#     return paths.len()
# def B(n, k):
#     m = k * n
#     s = 0
#     for i in range(1, m):
#         for j in range(i + 1, m + 1):
#             s += A(i, j, n, k)
#     return s
# def A003752(n):
#     return B(4, n)
# print([A003752(n) for n in range(1, 8)])
# def A390791_2(n):
#     return A(1, 3 * (2 * n + 2), 4, 2 * n + 2)
# print([A390791_2(n) for n in range(1, 8)])