a=[1, 12, -45, 1404, -33048, 684288, -5847309, -440129376, 30809872071, -952939532952, -1846906652457];
f(x) = sum(n=1, #a, a[n]*x^n);
print("A(A(x)): ", Vec(f(f(x)) + x*O(x^#a)))
print("A(A(A(x))): ", Vec(f(f(f(x))) + x*O(x^#a)))

for(n=0, 15, print1(n^2 * 9^(n-1), ", "))



