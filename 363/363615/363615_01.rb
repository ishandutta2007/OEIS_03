def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

def A(k, m, n)
  s = 0
  (1..n).each{|i|
    s += (-1) ** i * ncr(i + k, m) if n % i == 0
  }
  s
end

n = 10000
(1..n).each{|i|
  j = -A(-1, 2, i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}