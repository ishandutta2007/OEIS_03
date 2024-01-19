def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

# a(n) = (1/(n+1)) * Sum_{k=0..n} (-1)^k * binomial(3*n+k+2,k) * binomial(7*n-k+5,n-k).
def A(n)
  (0..n).inject(0){|s, i| s + (-1) ** (i % 2) * ncr(3 * n + i + 2, i) * ncr(7 * n - i + 5, n - i)} / (n + 1)
end

n = 1000
(0..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
