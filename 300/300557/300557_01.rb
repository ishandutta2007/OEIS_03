def A300557(n)
  a = '1'
  [1] + (2..n).map{|i| a = (a.to_s + i.to_s(2) * i).to_i}
end
n = 100
ary = A300557(n)
(1..n).each{|i|
  j = ary[i - 1]
  break if j.to_s.size > 1000
  print i
  print ' '
  puts j
}
