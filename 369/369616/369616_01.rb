def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

# a(n) = sum(k=0, n, binomial(n, k)*binomial(3*k+1, k)/(k+1));
def A(n)
  (0..n).inject(0){|s, i| s + ncr(n, i) * ncr(3 * i + 1, i) / (i + 1)}
end

n = 1000
(0..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
