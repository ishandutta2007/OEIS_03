seq(n) = my(A=1); for(i=1, n, A=exp( 2*x*(1 + x +x*O(x^n)) * A^(1/2) ) ); Vec(serlaplace(A));      
seq(35)    
