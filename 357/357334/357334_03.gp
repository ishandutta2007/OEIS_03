my(A=0, n=17); for(i=1, n, A=-log(1-x+x*O(x^n))*exp(3*A+x*O(x^n))); concat(0, Vec(serlaplace(A)))