def f(n)
  return 1 if n < 2
  (1..n).inject(:*)
end

def phi(n)
  (1..n).count{|i| i.gcd(n) == 1}
end

def A(n)
  s = 0
  (1..n).each{|i|
    s += phi(n / i) * f(i) if n % i == 0
  }
  s
end

n = 500
(0..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
