def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

def A(n)
  ary = (0..n).map{|i| ncr(n, i)}
  (0..n).inject(0){|s, i| s + (-1) ** (i % 2) * (0..i).inject(0){|t, j| t + ary[j]} ** 3}
end

n = 1000
(0..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
