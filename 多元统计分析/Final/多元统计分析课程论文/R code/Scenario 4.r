# 创建一个 48×48 的矩阵，初始化为 0
Beta_tens <- matrix(0, 48, 48)

# 圆的中心和半径
center_x <- 18
center_y <- 18
radius <- sqrt(0.10 * 48 * 48 / pi)

# 填充圆形区域为 1
for (i in 1:48) {
    for (j in 1:48) {
        # 计算 (i, j) 到圆心的距离
        if ((i - center_x)^2 + (j - center_y)^2 <= radius^2) {
            Beta_tens[i, j] <- 1
        }
    }
}
library(fields) # 用于 image.plot
image.plot(Beta_tens,
    col = gray.colors(25, start = 1, end = 0), axes = FALSE,
    main = "True Tensor Coefficient"
)
