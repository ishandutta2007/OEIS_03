def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

# a(n) = Sum_{k = 0..n} 4^(n-k)*binomial(n,k)*binomial(n-1,k)*binomial(2*k,k).
def A(n)
  (0..n).inject(0){|s, i| s + 4 ** (n - i) * ncr(n, i) * ncr(n - 1, i) * ncr(2 * i, i)}
end

n = 1000
(0..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
