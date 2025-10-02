# ==== 加载包 ====
library(LiblineaR)
library(fields)
library(imager)
library(magrittr)
library(dplyr)
library(abind)

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

load_images <- function(image_dir, target_size = c(48, 48)) {
    images <- list()
    labels <- c()

    for (label in c("yes", "no")) {
        folder_path <- file.path(image_dir, label)
        label_value <- ifelse(label == "yes", 1, -1)

        image_files <- list.files(folder_path, full.names = TRUE)

        for (file in image_files) {
            img <- load.image(file)

            # 如果是RGB彩色图，则转为灰度图
            if (spectrum(img) == 3) {
                img <- grayscale(img)
            }

            img_resized <- resize(img, target_size[1], target_size[2])
            img_matrix <- as.array(img_resized)

            # 转为 (height, width, 1) 格式
            img_array <- array(img_matrix[, , 1, 1], dim = c(target_size[2], target_size[1], 1))

            images[[length(images) + 1]] <- img_array
            labels <- c(labels, label_value)
        }
    }

    # 去掉通道维度，转为 (N, 48, 48)
    images_tensor <- abind(images, along = 0)[, , , 1]

    return(list(images = images_tensor, labels = labels))
}

# # ==== Scenario1 ====
# Beta_tens <- generate_parafac_tensor(p = c(48, 48), R = 3)

# # ==== Scenario 2 ====
# a1 <- rep(0, 48)
# a1[38:48] <- 1
# b1 <- rep(0, 48)
# b1[5:12] <- 1

# a2 <- rep(0, 48)
# a2[38:48] <- 1
# b2 <- rep(0, 48)
# b2[33:43] <- 1

# a3 <- rep(0, 48)
# a3[5:12] <- 1
# b3 <- rep(0, 48)
# b3[c(5:12, 33:43)] <- 1

# Beta_tens <- outer(a1, b1) + outer(a2, b2) + outer(a3, b3)

# # ==== Scenario 3 ====
# Beta_tens <- matrix(0, 48, 48)
# for (i in 15:40) {
#     for (j in 10:35) {
#         Beta_tens[i, j] <- 1
#     }
# }

# ==== Scenario 4 ====
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

# # ==== 脑肿瘤数据 ====
# # 设置路径
# train_dir <- "E:/BRAIN_TUMOR/BRAIN_TUMOR/train"
# val_dir <- "E:/BRAIN_TUMOR/BRAIN_TUMOR/val"

# # 加载数据
# train_data <- load_images(train_dir)
# val_data <- load_images(val_dir)

# # 合并数据
# all_images <- abind(train_data$images, val_data$images, along = 1)
# all_labels <- c(train_data$labels, val_data$labels)

# N <- dim(all_images)[1] # 总样本数
# idx <- sample(1:N) # 打乱索引
# split_point <- floor(0.7 * N)

# train_idx <- idx[1:split_point]
# val_idx <- idx[(split_point + 1):N]

# x.train <- all_images[train_idx, , ]
# y.train <- all_labels[train_idx]

# x.test <- all_images[val_idx, , ]
# y.test <- all_labels[val_idx]

# # 归一化
# x.train <- x.train / 255
# x.test <- x.test / 255

# 生成 2D 图像和二元相应变量
N <- 500
p <- c(48, 48)
rank <- 3
X <- array(rnorm(N * prod(p)), dim = c(N, p))
Y <- sapply(1:N, function(x) sum(X[x, , ] * Beta_tens, na.rm = T))
hist(Y)

# SVM Loss
Ylabel <- rep(0, N)
Ylabel[Y >= 0] <- 1
Ylabel[Y < 0] <- -1

# # Lasso loss
# p <- 1 / (1 + exp(-Y))
# Ylabel <- ifelse(p > 0.5, 1, -1)

# ==== 划分训练测试集 ====
train_index <- sort(sample(1:N, 0.7 * N, replace = FALSE))
x.train <- X[train_index, , ]
y.train <- Ylabel[train_index]
x.test <- X[-train_index, , ]
y.test <- Ylabel[-train_index]

# ==== 向量化 ====
x.train.mat <- t(apply(x.train, 1, c))
x.test.mat <- t(apply(x.test, 1, c))

# ==== 训练 L1-regularized L2-loss SVM ====
# type = 5 => L1-regularized L2-loss SVM (classification)
model <- LiblineaR(data = x.train.mat, target = y.train, type = 5, cost = 1, bias = TRUE, verbose = FALSE)

# ==== 提取系数并还原张量 ====
coef_vector <- model$W[-length(model$W)]
tensor_est <- matrix(coef_vector, nrow = 48, ncol = 48)

# ==== RMSE 和相关系数 ====
rmse_val <- sqrt(mean((c(Beta_tens) - c(tensor_est))^2))
cor_val <- cor(c(Beta_tens), c(tensor_est))
cat(sprintf("张量系数估计的 RMSE 为: %.6f\n", rmse_val))
cat(sprintf("张量系数估计与真实值的相关系数为: %.6f\n", cor_val))

# ==== 可视化 ====
image.plot(tensor_est,
    col = gray.colors(25, start = 1, end = 0), axes = F,
    main = "L1-SVM Estimated Tensor Coefficients"
)
mtext(text = seq(10, 50, 10), side = 2, line = 0.3, at = seq(10, 50, 10) / 48, las = 1, cex = 0.8)
mtext(text = seq(10, 50, 10), side = 1, line = 0.3, at = seq(10, 50, 10) / 48, las = 2, cex = 0.8)

# ==== 预测 ====
bias <- model$W[length(model$W)]
pred <- x.test.mat %*% coef_vector + bias
clust.test <- ifelse(pred > 0, 1, -1)

# ==== 误分类率 ====
missclassrate <- 1 - sum(clust.test == y.test) / length(y.test)
cat(sprintf("测试集误分类率为: %.4f\n", missclassrate))

# ==== F1-score ====
TP <- sum(clust.test == 1 & y.test == 1)
FP <- sum(clust.test == 1 & y.test == -1)
FN <- sum(clust.test == -1 & y.test == 1)
f1score <- TP / (TP + (FP + FN) / 2)
cat(sprintf("F1-score 为: %.4f\n", f1score))
