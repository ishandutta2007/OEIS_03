require 'prime'

def f(n)
  return 1 if n < 2
  (1..n).inject(:*)
end

def A003557(n)
  s = 1
  n.prime_division.each{|i| s *= i[0] ** (i[1] - 1)}
  s
end

def A(a, k, n)
  [1] + (1..n).map{|i| (k / a) ** i * (0..i - 1).inject(1){|s, j| s * (a * j - 1)} / f(i)}
end

n = 20
a = A003557(16)
ary = A(a, 16, n)
(0..n).each{|i|
  j = ary[i]
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
