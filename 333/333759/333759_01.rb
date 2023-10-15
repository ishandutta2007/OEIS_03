def f(n)
  # x == 0 && y == 1を除く
  (1..n * n - 1).select{|i| (i % n + 1) % n < 2 || (i / n + 1) % n < 2} - [n]
end

def search(x, y, n, used)
  return 0 if x < 0 || n <= x || y < 0 || n <= y || used[x + y * n]
  return 1 if x == 0 && y == 1 && f(n).all?{|i| used[i] == true}
  cnt = 0
  used[x + y * n] = true
  @move.each{|mo|
    cnt += search(x + mo[0], y + mo[1], n, used)
  }
  used[x + y * n] = false
  cnt
end

def A(n)
  return 1 if n < 3
  @move = [[1, 0], [-1, 0], [0, 1], [0, -1]]
  used = Array.new(n * n, false)
  search(0, 0, n, used)
end

def A333759(n)
  (2..n).map{|i| p A(i)}
end

p A333759(6)
