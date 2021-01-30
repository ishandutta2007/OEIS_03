# m次以下を取り出す
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

def s(f_ary, g_ary, n)
  s = 0
  (1..n).each{|i| s += i * f_ary[i] * g_ary[i] ** (n / i) if n % i == 0}
  s
end

def A(f_ary, g_ary, n)
  ary = [1]
  a = [0] + (1..n).map{|i| s(f_ary, g_ary, i)}
  (1..n).each{|i| ary << (1..i).inject(0){|s, j| s + a[j] * ary[-j]} / i}
  ary
end

def S(n)
  a = [1, 1]
  (2..n).each{|i|
    s = 0
    (1..i - 1).each{|j|
      s += a[j] * (-1) ** ((i / j) % 2) if i % j == 0
    }
    a << s
  }
  a
end

def B(n)
  ary1 = S(n)
  ary2 = Array.new(n + 1, 1)
  A(ary1, ary2, n)
end

n = 100
p ary = B(n)
a = Array.new(n + 1, 0)
a[0] = 1
(1..n).each{|i|
  b = Array.new(n + 1, 0)
  (0..n / i).each{|j| b[i * j] = ary[j]}
  a = mul(a, b, n)
}
b = Array.new(n + 1, 0)
(0..n / 2).each{|i| b[2 * i] = a[i]}
c = mul(b, b, n)
d = Array.new(n + 1, 1)
p e = mul(c, d, n)
p a == e
