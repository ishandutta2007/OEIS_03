def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

def A(k, n)
  (0..n / 2).inject(0){|s, i| s + ncr(n - i, i) * ncr(k * n + 1, n - i)} / (k * n + 1)
end

n = 1000
(0..n).each{|i|
  j = A(4, i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}