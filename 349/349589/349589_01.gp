my(A=1, n=22); for(i=1, n, A=exp((1-exp(-x*A))/(A+x*O(x^n)))); Vec(serlaplace(A))