seq(n) = my(A=1); for(i=1, n, A=exp(x*(1 + x +x*O(x^n))/A ) ); Vec(serlaplace(A));         
seq(18)
