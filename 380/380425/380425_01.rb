def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

def A(n)
  return 1 if n == 0
  2 * (1..n).inject(0){|s, i| s + i ** (n - i) * (n + 2) ** (i - 1) * ncr(n, i)}
end

n = 1000
(0..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}