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

def q(n)
  ary = [1]
  (1..n).each{|i|
    if i * (i + 1) / 2 <= n
      b_ary = Array.new(i * (i + 1) / 2 + 1, 0)
      b_ary[0], b_ary[-1] = 1, 1
      ary = mul(ary, b_ary, n)
    end
  }
  ary + [0] * (n + 1 - ary.size)
end

def r(n)
  ary = [1]
  (1..n).each{|i|
    if i * (i + 1) / 2 <= n
      b_ary = Array.new(i * (i + 1) / 2 + 1, 0)
      b_ary[0], b_ary[-1] = 1, -1
      ary = mul(ary, b_ary, n)
    end
  }
  ary
end

def I(ary, n)
  a = [1]
  (0..n - 1).each{|i| a << -(0..i).inject(0){|s, j| s + ary[1 + i - j] * a[j]}}
  a
end

n = 10100
ary0 = q(n)
ary1 = r(n)
ary2 = I(ary0, n)
ary = mul(ary2, ary1, n)
(0..10000).each{|i|
  j = ary[i]
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
