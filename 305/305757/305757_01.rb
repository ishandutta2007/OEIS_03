require 'prime'

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

# m���ȉ������o��
def power(ary, n, m)
  return [1] if n == 0
  k = power(ary, n >> 1, m)
  k = mul(k, k, m)
  return k if n & 1 == 0
  return mul(k, ary, m)
end

def A001158(n)
  s = 0
  (1..n).each{|i| s += i * i * i if n % i == 0}
  s
end

def A004009(n)
  a = [1] + (1..n).map{|i| 240 * A001158(i)}
end

def E4_3(n)
  ary = A004009(n)
  power(ary, 3, n)
end

def s(n)
  s = 0
  (1..n).each{|i| s += i if n % i == 0}
  s
end

def A(k, n)
  ary = [1]
  a = [0] + (1..n).map{|i| s(i)}
  (1..n).each{|i| ary << (1..i).inject(0){|s, j| s - k * a[j] * ary[-j]} / i}
  ary
end

def f(n)
  a = E4_3(n)
  b = A(-24, n)
  mul(a, b, n)
end

def A008683(n)
  ary = n.prime_division
  return (-1) ** (ary.size % 2) if ary.all?{|i| i[1] == 1}
  0
end

# ary[0] = 1
def inverse_Euler_transform(ary, n)
  c_ary = [0]
  (1..n).each{|i| c_ary << (1..i - 1).inject(i * ary[i]){|s, j| s - ary[j] * c_ary[-j]}}
  m_ary = [0] + (1..n).map{|i| A008683(i)}
  a = [0]
  (1..n).each{|i|
    s = 0
    (1..i).each{|j|
      s += m_ary[i / j] * c_ary[j] if i % j == 0
    }
    if s % i == 0
      a << s / i
    else
      a << s / i.to_r
    end
  }
  a
end

n = 11
a = f(n)
a[1] = 24
ary = inverse_Euler_transform(a, n)
(1..n).each{|i|
  j = ary[i]
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
