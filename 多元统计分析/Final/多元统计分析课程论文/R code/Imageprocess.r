library(imager)
library(magrittr)
library(dplyr)
library(abind)

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

# 设置路径
train_dir <- "E:/BRAIN_TUMOR/BRAIN_TUMOR/train"
val_dir <- "E:/BRAIN_TUMOR/BRAIN_TUMOR/val"

# 加载数据
train_data <- load_images(train_dir)
val_data <- load_images(val_dir)

# 拆分图像与标签
X_train <- train_data$images
Y_train <- train_data$labels

X_val <- val_data$images
Y_val <- val_data$labels

# 归一化
X_train <- X_train / 255
X_val <- X_val / 255

# 查看维度
cat("训练集维度：", dim(X_train), "\n")
cat("验证集维度：", dim(X_val), "\n")

# 可视化一张图像
image(t(X_train[200, , ]),
    col = gray.colors(256),
    useRaster = TRUE, main = paste("Label:", Y_train[1])
)
