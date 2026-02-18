# 和因子はmin以上max以下
def partition(n, min, max)
  return [[]] if n == 0
  [max, n].min.downto(min).flat_map{|i| partition(n - i, min, i).map{|rest| [i, *rest]}}
end

def A(n, k)
  cnt = 0
  partition(n, 1, n).each{|ary|
    cnt += 1 if ary.count(&:odd?) == k * ary.count(&:even?)
  }
  cnt
end

p (0..25).map{|i| A(2*i, 4)}

