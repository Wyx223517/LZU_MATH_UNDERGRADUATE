library(e1071)
library(ROCR)
library(MASS)

setwd("E:/719465160/课程作业/多元统计分析/HW03/判别分析/Rdata")

## 读取Heart数据
Heart <- read.csv(
    "Heart.csv",
    header = TRUE, row.names = 1,
    stringsAsFactors = TRUE
)
d <- na.omit(Heart)
set.seed(1)
train_id <- sample(nrow(d), size = round(0.7 * nrow(d)))
train <- rep(FALSE, nrow(d))
train[train_id] <- TRUE
test <- (!train)
d[["AHD"]] <- factor(d[["AHD"]], levels = c("No", "Yes"))

## 定义错判率函数
classifier.error <- function(truth, pred) {
    tab1 <- table(truth, pred)
    err <- 1 - sum(diag(tab1)) / sum(c(tab1))
    err
}

## 进行线性判别法LDA
res.ld <- lda(AHD ~ ., data = d[train, ], prior = rep(1 / 2, 2))
fit.ld <- predict(res.ld)$class
tab1 <- table(truth = d[train, "AHD"], fitted = fit.ld)
tab1
cat(
    "LDA错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)

## 支持向量判别法，其中取 cost=1
res.svc <- svm(
    AHD ~ .,
    data = d[train, ],
    kernel = "linear", cost = 1, scale = TRUE
)
fit.svc <- predict(res.svc)
summary(res.svc)

tab1 <- table(truth = d[train, "AHD"], fitted = fit.svc)
tab1
cat(
    "SVC错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)
## tune()函数选取较好的调节参数
set.seed(101)
res.tune <- tune(
    svm, AHD ~ .,
    data = d[train, ], kernel = "linear", scale = TRUE,
    ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100, 1000))
)
summary(res.tune)
summary(res.tune$best.model)

## 将两种方法训练的模型用于测试集
pred.ld <- predict(res.ld, d[test, ])$class
tab1 <- table(truth = d[test, "AHD"], predict = pred.ld)
tab1

cat(
    "LDA错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)

pred.svc <- predict(res.tune$best.model, newdata = d[test, ])
tab1 <- table(truth = d[test, "AHD"], predict = pred.svc)
tab1

cat(
    "SVC错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)

## ROC曲线分析
png(filename = "roc_curve_training.png", width = 800, height = 600)
pred.ld <- predict(res.ld,
    newdata = d[train, ],
    decision.values = TRUE
)
decf.ld <- pred.ld$posterior[, "Yes"]
predobj <- prediction(decf.ld, d[train, "AHD"],
    label.ordering = c("No", "Yes")
)
perf <- performance(predobj, "tpr", "fpr")
plot(perf, main = "ROC Curve (Training Data)", col = "cyan")

pred.svc <- predict(
    res.tune$best.model,
    newdata = d[train, ],
    decision.values = TRUE
)
decf.svc <- attributes(pred.svc)$decision.values
predobj <- prediction(
    decf.svc, d[train, "AHD"],
    label.ordering = c("No", "Yes")
)
perf <- performance(predobj, "tpr", "fpr")
plot(perf, add = TRUE, col = "blue")
legend("bottomright",
    lty = 1, col = c("cyan", "blue"),
    legend = c("LDA", "SVC")
)
dev.off()

## 多项式核SVM
res.svm1 <- svm(AHD ~ .,
    data = d[train, ], kernel = "polynomial",
    order = 2, cost = 0.1, scale = TRUE
)
fit.svm1 <- predict(res.svm1)
summary(res.svm1)
tab1 <- table(truth = d[train, "AHD"], fitted = fit.svm1)
tab1
cat(
    "2阶多项式核SVM错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)

## 寻找最优参数
set.seed(101)
res.tune2 <- tune(
    svm, AHD ~ .,
    data = d[train, ], kernel = "polynomial",
    order = 2, scale = TRUE,
    ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100, 1000))
)
summary(res.tune2)
fit.svm2 <- predict(res.tune2$best.model)
tab1 <- table(truth = d[train, "AHD"], fitted = fit.svm2)
tab1
cat(
    "2阶多项式核最优参数SVM错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)

## 用于测试集
pred.svm2 <- predict(res.tune2$best.model, d[test, ])
tab1 <- table(truth = d[test, "AHD"], predict = pred.svm2)
tab1
cat(
    "2阶多项式核最优参数SVM测试集错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)

## 径向核SVM
res.svm3 <- svm(
    AHD ~ .,
    data = d[train, ], kernel = "radial",
    gamma = 0.1, cost = 0.1, scale = TRUE
)
fit.svm3 <- predict(res.svm3)
summary(res.svm3)
tab1 <- table(truth = d[train, "AHD"], fitted = fit.svm3)
tab1
cat(
    "径向核(gamma=0.1, cost=0.1)SVM错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)

## 寻找最优参数
set.seed(101)
res.tune4 <- tune(
    svm, AHD ~ .,
    data = d[train, ], kernel = "radial",
    scale = TRUE,
    ranges = list(
        cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100, 1000),
        gamma = c(0.1, 0.01, 0.001)
    )
)
summary(res.tune4)
fit.svm4 <- predict(res.tune4$best.model)
tab1 <- table(truth = d[train, "AHD"], fitted = fit.svm4)
tab1
cat(
    "径向核最优参数SVM错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)

## 用于测试集
pred.svm4 <- predict(res.tune4$best.model, d[test, ])
tab1 <- table(truth = d[test, "AHD"], predict = pred.svm2)
tab1
cat(
    "径向核最优参数SVM测试集错判率:",
    round((tab1[1, 2] + tab1[2, 1]) / sum(c(tab1)), 2), "\n"
)
