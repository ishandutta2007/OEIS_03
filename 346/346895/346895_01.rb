def f(n)
  return 1 if n < 2
  (1..n).inject(:*)
end

def A(k, n)
  m = f(k)
  a = [1]
  ary = [1]
  (1..n).each{|i|
    a << 0
    b = [0]
    (0..i - 1).each{|j|
      b[j + 1] = a[j] + (j + 1) * a[j + 1]
    }
    a = b
    ary << (0..i / k).inject(0){|s, j| s + f(k * j) * a[k * j] / (m ** j).to_r}.to_i
  }
  ary
end

n = 500
ary = A(4, n)
(0..n).each{|i|
  j = ary[i]
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
