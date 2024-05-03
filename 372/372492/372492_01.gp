a=[1, 2, 0, 4, -8, -8, 288, -1712, -1888, 105472, -288576, -10404800, 84940672, 1454871936, -24372060160, -255228956416, 8232158755328, 49829958005760];
f(x) = sum(n=1, #a, a[n]*x^n);
print("B(x): ", Vec(f(f(x)) + x*O(x^#a)))

\\ Let B(x) = A(A(x)). B(B(x)) = F(x).
b=[1, 4, 8, 16, 32, 0, 256, 768, -14848, 51200, 1077248, -12132352, -89997312, 2637987840, 5445910528, -693655371776, 1905039638528, 226021104680960];
f(x) = sum(n=1, #b, b[n]*x^n);
print("F(x): ", Vec(f(f(x)) + x*O(x^#b)))

for(n=0, 15, print1(n * 4^(n-1), ", "))

\\ B(x) = G(2*x)/2, where G(x) is the g.f. for A309509.
print("B(x): ", b)

c=[1, 2, 2, 2, 2, 0, 4, 6, -58, 100, 1052, -5924, -21972, 322020, 332392, -21168682, 29068598];
g(x) = sum(n=1, #c, c[n]*x^n);
Vec(g(2*x)/2 + x*O(x^#c))

