def ncr(n, r)
  return 1 if r == 0
  (n - r + 1..n).inject(:*) / (1..r).inject(:*)
end

def A(k, n)
  ary = [1]
  (1..n).each{|i| 
    ary << (0..i - k).inject(ncr(i + k - 1, k - 1)){|s, j| s - ary[j] * ary[-k - j]}
  }
  ary
end

n = 1100
m = 1000
ary = A(5, n)
(0..m).each{|i|
  j = ary[i]
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}