my(A=1, n=21); for(i=1, n, A=(1-x*A+x*O(x^n))^(-log(1-x*A+x*O(x^n))^2/(6*A))); Vec(serlaplace(A)) 