M=32;

my(A=1, n=M); for(i=1, n, A=1/(1-x+x*O(x^n))^4-x^4*A^2); Vec(A)   