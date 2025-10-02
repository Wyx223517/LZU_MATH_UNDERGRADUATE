library(fields)

# ==== 参数 ====
p <- c(48, 48)
R <- 3

# ==== 构造 rank-R PARAFAC 张量 ====
generate_parafac_tensor <- function(p, R, prob = 0.2, size = 2) {
    A_list <- list()
    B_list <- list()

    for (r in 1:R) {
        a_r <- rbinom(p[1], size = size, prob = prob)
        b_r <- rbinom(p[2], size = size, prob = prob)
        A_list[[r]] <- a_r
        B_list[[r]] <- b_r
    }

    tensor <- matrix(0, nrow = p[1], ncol = p[2])
    for (r in 1:R) {
        tensor <- tensor + outer(A_list[[r]], B_list[[r]])
    }

    # 标准化为最大值为1
    tensor <- tensor / max(tensor)
    return(tensor)
}

# ==== 生成张量 ====
set.seed(123)
Beta_tens <- generate_parafac_tensor(p = c(48, 48), R = 3)

# ==== 保存图片 ====
png("Beta_tensor_Scenario1.png", width = 800, height = 800, res = 120)
image.plot(Beta_tens, col = gray.colors(25, start = 1, end = 0), axes = FALSE)
mtext(text = seq(10, 50, 10), side = 2, line = 0.3, at = seq(10, 50, 10) / 48, las = 1, cex = 0.8)
mtext(text = seq(10, 50, 10), side = 1, line = 0.3, at = seq(10, 50, 10) / 48, las = 2, cex = 0.8)
dev.off()
