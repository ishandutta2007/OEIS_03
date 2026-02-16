def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

# a(n) = Sum_{k=0..n} (-1)^(n-k) * binomial(n,k)^2 * binomial(n+k,k).
def a(n)
  (0..n).inject(0){|s, k| s + (-1)**(n-k) * ncr(n, k)**2 * ncr(n + k, k)}
end

n = 1000
(0..n).each{|i|
  j = a(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}