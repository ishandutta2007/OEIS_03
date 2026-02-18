# 和因子はmin以上max以下
def partition(n, min, max)
  return [[]] if n == 0
  [max, n].min.downto(min).flat_map{|i| partition(n - i, min, i).map{|rest| [i, *rest]}}
end

def A(n)
  cnt = 0
  partition(n, 1, n).each{|ary|
    cnt += 1 if 2 * ary.count(&:odd?) == 3 * ary.count(&:even?)
  }
  cnt
end

p (0..50).map{|i| A(i)}

