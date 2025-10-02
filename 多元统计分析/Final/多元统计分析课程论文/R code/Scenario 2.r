library(fields)

# ==== 张量维度与秩 ====
p <- c(48, 48)
R <- 3

# ==== 手动设置边缘向量 ====

a1 <- rep(0, 48)
a1[38:48] <- 1
b1 <- rep(0, 48)
b1[5:12] <- 1


a2 <- rep(0, 48)
a2[38:48] <- 1
b2 <- rep(0, 48)
b2[33:43] <- 1


a3 <- rep(0, 48)
a3[5:12] <- 1
b3 <- rep(0, 48)
b3[c(5:12, 33:43)] <- 1

# ==== 构造 PARAFAC 张量 ====
Beta_tens <- outer(a1, b1) + outer(a2, b2) + outer(a3, b3)

# ==== 可视化 ====
image.plot(Beta_tens, col = gray.colors(25, start = 1, end = 0), axes = FALSE)
mtext(text = seq(10, 50, 10), side = 2, line = 0.3, at = seq(10, 50, 10) / 48, las = 1, cex = 0.8)
mtext(text = seq(10, 50, 10), side = 1, line = 0.3, at = seq(10, 50, 10) / 48, las = 2, cex = 0.8)
