\\ a(0) = 0; a(n) = (n+1) * n * a(n-1) + 1.
a(n) = if(n==0, 0, (n+1) * n * a(n-1) + 1);
for(n=0, 16, print1(a(n), ", "));  
