def f(n)
  return 1 if n < 2
  (1..n).inject(:*)
end

def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

def A(n)
  (0..n - 1).inject(0){|s, i| s + ncr(n - 1, i) / f(i + n).to_r}
end

n = 500
(1..n).each{|i|
  j = A(i).denominator
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}

