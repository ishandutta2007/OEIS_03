def f(n)
  return 1 if n < 2
  (1..n).inject(:*)
end

# m���ȉ������o��
def mul(f_ary, b_ary, m)
  s1, s2 = f_ary.size, b_ary.size
  ary = Array.new(s1 + s2 - 1, 0)
  (0..s1 - 1).each{|i|
    (0..s2 - 1).each{|j|
      ary[i + j] += f_ary[i] * b_ary[j]
    }
  }
  ary[0..m]
end

def A(k, n)
  a = [0] + (1..n).map{|i| i ** k / f(i).to_r}
  b = a
  ary = a.clone
  (2..n).each{|i|
    b = mul(a, b, n)
    i1 = (-1) ** (i % 2) / i.to_r
    (1..n).each{|j|
      ary[j] -= i1 * b[j]
    }
  }
  [0] + (1..n).map{|i| (ary[i] * f(i)).to_i}
end

def A308484(n)
  a = []
  (0..n).each{|i| a << A(i, n - i + 1)}
  ary = []
  (1..n).each{|i|
    (1..i).each{|j|
      ary << a[i - j][j]
    }
  }
  ary
end

n = 140
ary = A308484(n)
(1..ary.size).each{|i|
  j = ary[i - 1]
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
