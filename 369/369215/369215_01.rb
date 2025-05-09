def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

# a(n) = (1/(n+1)) * Sum_{k=0..n} binomial(n+k,k) * binomial(4*n+2*k+2,n-k).
def A(n)
  (0..n).inject(0){|s, i| s + ncr(n + i, i) * ncr(4 * n + 2 * i + 2, n - i)} / (n + 1)
end

n = 1000
(0..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}