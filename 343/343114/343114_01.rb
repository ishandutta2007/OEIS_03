n = 1000
(1..n).each{|i|
  print i
  print ' '
  puts (1..i).inject(0){|s, j| s + j.gcd(i ** j)}
}
