def A(k, n)
  s = 0
  (1..n).each{|i|
    j = i * i
    break if j > n
    s += i ** (2 * k) if n % j == 0
  }
  s
end

n = 10000
(1..n).each{|i|
  j = A(8, i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}

