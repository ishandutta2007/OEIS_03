a(n) = n! * sum(k=0, n, k^(n-k) * binomial(3*k+1,k)/( (3*k+1)*(n-k)! ) );                        
for(n=0, 17, print1(a(n),", "))        

a(n) = n!*sum(k=0, n, k^(n-k)*binomial(3*k, k)/((2*k+1)*(n-k)!));                    
for(n=0, 17, print1(a(n),", "))      
