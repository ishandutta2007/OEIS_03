# 和因子はmin以上max以下
def partition(n, min, max)
  return [[]] if n == 0
  [max, n].min.downto(min).flat_map{|i| partition(n - i, min, i).map{|rest| [i, *rest]}}
end

def A(n, k)
  return 1 if n == 0
  cnt = 0
  partition(n, 1, n).each{|ary|
    m = ary.size
    cnt += 1 if ary[-1] >= k * m
    # cnt += 1 if 5 * ary[0] >= m
  }
  cnt
end 

# p (1..70).map{|i| A(i, 1)}
# p (0..70).map{|i| A(i, 2)}
# p (1..70).map{|i| A(i, 3)}
# p (1..70).map{|i| A(i, 4)}
p (0..70).map{|i| A(i, 5)}

# [1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 9, 10, 12, 14, 17, 19, 23, 26, 31, 35, 41, 46, 54, 61, 70, 79, 91, 102, 117, 131, 149, 167, 189, 211, 239, 266, 299, 333, 374, 415, 465, 515, 575, 637, 709, 783, 871, 961, 1065, 1174, 1299, 1429, 1579, 1735, 1913, 2100, 2311, 2533, 2785, 3049, 3345, 3659, 4010, 4380, 4794, 5231, 5717, 6233, 6804]
# [0, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 8, 8, 10, 11, 13, 14, 17, 18, 21, 23, 26, 28, 32, 34, 39, 42, 47, 51, 58, 62, 70, 76, 85, 92, 103, 111, 124, 134, 148, 160, 177, 190, 210, 226, 248, 267, 293, 315, 345, 371, 405, 436, 476, 511, 557, 599, 651, 700, 760, 816, 885, 950, 1028]
# [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 10, 11, 12, 14, 15, 17, 19, 21, 23, 26, 28, 31, 34, 37, 40, 44, 47, 51, 55, 59, 63, 69, 73, 79, 85, 92, 98, 107, 114, 124, 133, 144, 154, 168, 179, 194, 208, 225, 240, 260, 277, 299, 319, 343]
# [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 13, 13, 15, 16, 18, 19, 22, 23, 26, 28, 31, 33, 37, 39, 43, 46, 50, 53, 58, 61, 66, 70, 75, 79, 85, 89, 95, 100, 107, 112, 120, 126, 135, 142, 152]
# [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 16, 17, 19, 20, 22, 24, 26, 28, 31, 33, 36, 39, 42, 45, 49, 52, 56, 60, 64, 68, 73, 77, 82, 87, 92]

