# 参数设置
set.seed(123)
N <- 500
p <- c(48, 48)

# 生成图像数据张量 X
X <- array(rnorm(N * prod(p)), dim = c(N, p))

# 可视化前3张图像
library(fields)
for (i in 1:3) {
    image.plot(X[i, , ],
        col = gray.colors(25, start = 1, end = 0),
        main = paste("Sample Image", i),
        axes = FALSE
    )
    Sys.sleep(1)
}
