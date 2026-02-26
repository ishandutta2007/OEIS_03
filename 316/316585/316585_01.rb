def a316585(max_n)
  # 結果を格納する配列（n=0..max_n）
  ans = Array.new(max_n + 1, 0)
  
  # dp[y][v] : 現在のステップxにおいて、高さy、累積値(twoar)がvであるパスの数
  dp = Array.new(max_n + 1){Array.new(max_n + 1, 0)}

  # 初期条件: x=0, y=0, twoar=0 のとき Motzk=1
  dp[0][0] = 1
  ans[0] += 1

  # x を 1 から max_n まで回す
  (1..max_n).each{|x|
    new_dp = Array.new(max_n + 1){Array.new(max_n + 1, 0)}
    
    # xステップ目における各高さyと値vの状態を計算
    (0..x).each{|y|
      (0..max_n).each{|v|
        # Mapleの再帰ロジックをループに変換
        if y == 0
          if x.even? && y + 1 <= max_n && v >= 1
            new_dp[0][v] = dp[1][v - 1]
          end
        elsif y % 2 == x % 2
          # y+1, y, y-1 からの遷移
          sum = 0
          sum += dp[y + 1][v - (2 * y + 1)] if y + 1 <= max_n && v >= (2 * y + 1)
          sum += dp[y][v - (2 * y)]         if v >= (2 * y)
          sum += dp[y - 1][v - (2 * y - 1)] if v >= (2 * y - 1)
          new_dp[y][v] = sum
        else
          # y からの遷移のみ
          new_dp[y][v] = dp[y][v - (2 * y)] if v >= (2 * y)
        end
        
        # A316585(v) は全x, yにおける総和
        ans[v] += new_dp[y][v] if new_dp[y][v] > 0
      }
    }
    
    dp = new_dp
  } 

  ans
end

# 実行と出力
n = 1000
ary = a316585(n)
(0..n).each{|i|
  print i
  print ' '
  puts ary[i]
}