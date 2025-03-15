# 新版本：即使部分模型拟合失败也绘图，并将结果保存到一个新的文件夹
source("main.R")
# 定义新函数，允许部分模型拟合成功
fit_and_select_best_model_6_partial <- function(df_sub, n_starts = 50) {
  
  # 1) Quadratic
  fit_quad <- tryCatch(
    lm(logPopBio ~ Time + I(Time^2), data = df_sub),
    error = function(e) NULL
  )
  # 2) Cubic
  fit_cubic <- tryCatch(
    lm(logPopBio ~ Time + I(Time^2) + I(Time^3), data = df_sub),
    error = function(e) NULL
  )
  # 3) Logistic
  fit_logistic <- tryCatch(
    multi_start_fit_one_model(df_sub, "Logistic", n_starts),
    error = function(e) NULL
  )
  # 4) Gompertz
  fit_gompertz <- tryCatch(
    multi_start_fit_one_model(df_sub, "Gompertz", n_starts),
    error = function(e) NULL
  )
  # 5) Baranyi
  fit_baranyi <- tryCatch(
    multi_start_fit_one_model(df_sub, "Baranyi", n_starts),
    error = function(e) NULL
  )
  # 6) Three-phase linear
  fit_threephase <- tryCatch(
    multi_start_fit_three_phase_model(df_sub, n_starts),
    error = function(e) NULL
  )
  
  # 将所有模型拟合结果存入列表
  model_list <- list(
    Quadratic  = fit_quad,
    Cubic      = fit_cubic,
    Logistic   = fit_logistic,
    Gompertz   = fit_gompertz,
    Baranyi    = fit_baranyi,
    ThreePhase = fit_threephase
  )
  
  # 只保留拟合成功的模型（非 NULL）
  model_list <- model_list[!sapply(model_list, is.null)]
  
  # 如果没有任何模型拟合成功，则返回 NULL
  if(length(model_list) == 0) {
    return(NULL)
  }
  
  # 计算所有拟合成功模型的 AICc 值，选择 AICc 最小的模型
  aicc_values <- sapply(model_list, AICc)  # 使用 MuMIn::AICc
  best_model_name <- names(aicc_values)[which.min(aicc_values)]
  
  # 生成预测网格
  t_min <- min(df_sub$Time, na.rm = TRUE)
  t_max <- max(df_sub$Time, na.rm = TRUE)
  new_time <- seq(t_min, t_max, length.out = 200)
  
  # 对每个拟合成功的模型分别生成预测数据
  preds_list <- lapply(names(model_list), function(mod_name) {
    model_fit <- model_list[[mod_name]]
    data.frame(
      Time = new_time,
      logPopBio = predict(model_fit, newdata = data.frame(Time = new_time)),
      Model = mod_name
    )
  })
  
  preds_all <- do.call(rbind, preds_list)
  
  # 返回包含预测数据和最佳模型名称的列表
  return(list(
    predictions = preds_all,
    best_model  = best_model_name,
    models = model_list
  ))
}

# 设置输出结果的新文件夹路径
output_dir <- "../results_partial/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 对所有 ID 进行循环处理
for(this_id in all_ids) {
  # 取出该 ID 对应的数据子集
  df_sub <- growth_data %>% filter(ID == this_id)
  
  # 跳过数据点少于 10 的子集
  if(nrow(df_sub) < 8) {
    next
  }
  
  # 使用新函数进行 6 模型拟合，允许部分模型失败
  fit_result <- fit_and_select_best_model_6_partial(df_sub, n_starts = 150)
  
  # 如果所有模型都拟合失败，仍然绘制原始数据散点图
  if(is.null(fit_result)) {
    p <- ggplot(df_sub, aes(x = Time, y = logPopBio)) +
      geom_point() +
      labs(
        title = paste0("ID = ", this_id, " (No model fitted successfully)"),
        x = "Time",
        y = "log(PopBio)"
      ) +
      theme_minimal()
    
    png_filename <- paste0(output_dir, "6models_ID_", this_id, ".png")
    png(png_filename, width = 2000, height = 1500, res = 200)
    print(p)
    dev.off()
    
    message("输出PNG (仅原始数据): ", png_filename)
    next
  }
  
  # 获取预测数据和最佳模型名称
  pred_df   <- fit_result$predictions
  best_name <- fit_result$best_model
  
  # 绘制图像：原始散点图 + 所有拟合成功模型的预测曲线
  p <- ggplot() + 
    geom_point(data = df_sub, aes(x = Time, y = logPopBio)) + 
    geom_line(data = pred_df, aes(x = Time, y = logPopBio, color = Model), size = 0.8) + 
    labs(
      title = paste0("ID = ", this_id, ", Best Model (AICc) = ", best_name),
      x = "Time",
      y = "log(PopBio)"
    ) + 
    theme_minimal() + 
    scale_color_discrete(name = "Model") + 
    theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold")
    )
  
  # 保存图像为 PNG 文件
  png_filename <- paste0(output_dir, "6models_ID_", this_id, ".png")
  png(png_filename, width = 2000, height = 1500, res = 200)
  print(p)
  dev.off()  
  
  message("输出PNG: ", png_filename, " ; Best model = ", best_name) 
}

message("全部完成: 现在已有包括三相线性模型在内的 6 模型对比结果及图像保存在 ", output_dir)
