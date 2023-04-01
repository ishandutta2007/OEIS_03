require 'prime'

def A(n)
  return 1 if n == 1
  s = 1
  n.prime_division.each{|i, j|
    if i == 2
      if j == 1
        s *= -5
      elsif j == 2
        s *= 4
      else
        s *= 0
      end
    else
      if j == 1
        s *= -i * i
      else
        s *= 0
      end
    end
  }
  s
end

n = 10000
(1..n).each{|i| 
  print i
  print ' '
  puts A(i)
}

