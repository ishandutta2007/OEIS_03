def A(n)
  ary = [1]
  (1..n).each{|i|
    ary += ary.uniq.sort.map{|i| ary.count(i)}
  }
  ary
end

n = 100
ary = A(n)
p a = (1..ary.size).select{|i| ary[i - 1] == 6}
# (1..1000).each{|i|
#   print i
#   print ' '
#   puts a[i - 1]
# }