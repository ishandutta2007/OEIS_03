def phi(n)
  (1..n).count{|i| i.gcd(n) == 1}
end

def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

def A(k, n)
  s = 0
  (1..n).each{|i|
    s += phi(n / i) * ncr(i + k - 1, k) if n % i == 0
  }
  s
end

n = 1000
(1..n).each{|i|
  j = A(i - 1, i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
