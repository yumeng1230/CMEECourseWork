########################## 
# 2. 线性(多项式)回归对比 
##########################

library(dplyr)
library(broom)
library(MuMIn)  # 用于计算 AICc

# 读取数据
datalinear <- read.csv("../data/modified_growth_data.csv")

# 函数：对某个数据子集拟合二次 & 三次多项式，并输出结果
fit_models <- function(df) {
  if (nrow(df) < 8) return(NULL)  # 过滤过小数据集
  
  # 二次多项式
  model1 <- lm(logPopBio ~ Time + I(Time^2), data = df)
  summary1 <- glance(model1)
  coef1 <- tidy(model1) %>% filter(term == "Time") %>% pull(estimate)
  
  # 三次多项式
  model2 <- lm(logPopBio ~ Time + I(Time^2) + I(Time^3), data = df)
  summary2 <- glance(model2)
  coef2 <- tidy(model2) %>% filter(term == "Time") %>% pull(estimate)
  
  tibble(
    ID = unique(df$ID),
    Model1_AIC = AIC(model1),
    Model1_AICc = AICc(model1),
    Model1_BIC = BIC(model1),
    Model1_R2 = summary1$r.squared,
    Model1_B1 = coef1,
    
    Model2_AIC = AIC(model2),
    Model2_AICc = AICc(model2),
    Model2_BIC = BIC(model2),
    Model2_R2 = summary2$r.squared,
    Model2_B1 = coef2,
    
    Best_Model = ifelse(AICc(model1) < AICc(model2), "Quadratic", "Cubic")
  )
}

# 对每个 ID 分组拟合
results_linear <- datalinear %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(fit_models) %>%
  bind_rows()

# 仅选择包含 AICc 和 BIC 的列
results_selected <- results_linear %>%
  dplyr::select(ID, Model1_AICc, Model1_BIC, Model2_AICc, Model2_BIC)

# 保存完整结果
write.csv(results_linear, "../results/model_comparison.csv", row.names = FALSE)
message("Full model comparison results saved to '../results/model_comparison.csv'")

# 仅保存 AICc 和 BIC 结果
write.csv(results_selected, "../results/model_aicc_bic_comparison.csv", row.names = FALSE)
message("AICc and BIC comparison saved to '../results/model_aicc_bic_comparison.csv'")

