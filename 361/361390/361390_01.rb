def carryless_mul(f_ary, b_ary)
  s1, s2 = f_ary.size, b_ary.size
  ary = Array.new(s1 + s2 - 1, 0)
  (0..s1 - 1).each{|i|
    (0..s2 - 1).each{|j|
      ary[i + j] += f_ary[i] * b_ary[j]
    }
  }
  ary.map{|i| i % 10}
end

def carryless_power(ary, n)
  return [1] if n == 0
  k = carryless_power(ary, n >> 1)
  k = carryless_mul(k, k)
  return k if n & 1 == 0
  return carryless_mul(k, ary)
end

def A(n, k)
  carryless_power(n.to_s.split('').map(&:to_i), k).join.to_i
end

n = 139
cnt = 0
(0..n).each{|i|
  (0..i).each{|j|
    k = A(j, i - j)
    break if k.to_s.size > 1000
    print cnt
    print ' '
    puts k
    cnt +=1
  }
}
