def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

def A(k, m, n)
  return 1 if n == 0
  (0..(n - 1) / 2).inject(0){|s, i| s + k ** (n - i) * ncr(n, i) * ncr(m * n - i, n - 1 - 2 * i)} / n
end

n = 1000
(0..n).each{|i|
  j = A(4, 3, i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}

