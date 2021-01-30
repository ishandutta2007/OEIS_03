require 'prime'

def power(a, n)
  return 1 if n == 0
  k = power(a, n >> 1)
  k *= k
  return k if n & 1 == 0
  return k * a
end

def sigma(x, i)
  sum = 1
  pq = i.prime_division
  pq.each{|a, n| sum *= (power(a, (n + 1) * x) - 1) / (power(a, x) - 1)}
  sum
end

def A(n)
  a = [0] + (1..n).map{|i| sigma(1, i * n)}
  ary = [1]
  (1..n).each{|i|
    ary << (1..i).inject(0){|s, j| s + a[j] * ary[-j]} / i
  }
  ary[-1]
end

n = 1000
(0..n).each{|i|
  j = A(i)
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
