def A(n)
  (1..n).inject(0){|s, i| s + (n / i) ** i}
end

n = 5000
(1..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
