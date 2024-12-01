# This code is part of the final project for 2024 Fall STAT 689. 
# You are free to use, modify, and distribute this code with proper attribution.
# For questions or feedback, please contact: chaoan@tamu.edu


# Load necessary packages
if (!requireNamespace("nloptr")) install.packages("nloptr")
library(nloptr)
if (!requireNamespace("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# Parameter settings
set.seed(123)
n_genes <- 10000       # Number of genes
prop_effect <- 0.01    # Proportion of truly significant genes
alpha <- 0.05         # FDR threshold
reward <- 0.05        

# Construct hypotheses (H₀ and H₁)
true_effect <- sample(c(TRUE, FALSE), n_genes, replace = TRUE, prob = c(prop_effect, 1 - prop_effect))

# Simulated data
# - For H₀ (null hypothesis), generate test statistics that follow a standard normal distribution
# - For H₁ (alternative hypothesis), generate test statistics with a shift in the mean (e.g., mean = 2)
z_scores <- ifelse(true_effect, rnorm(n_genes, mean = 5, sd = 1), rnorm(n_genes, mean = 0, sd = 1))

# Convert z-scores to p-values
p_values <- 2 * (1 - pnorm(abs(z_scores)))

# Use the Benjamini-Hochberg method to adjust p-values
#adjusted_p_values <- p.adjust(p_values, method = "BH")

# Reject score function
RejectScore <- function(false_positive,true_positive) {
  RS = true_positive * reward - false_positive * (1 - reward)
  return(RS)
}

# All methods and parameters
methods <- c("BH", "BC", "FBC", "FBH", "ST")
function_types <- c("log", "sqrt", "linear", "square", "exp")
lambdas <- seq(0, 0.0006, by = 0.0001)

# Define the custom function phi
Flexfunction <- function(function_type, x) {
  if (function_type == "log") {
    return(log(x + 1) / log(2))  # log(x+1)/(log(2))
  } else if (function_type == "sqrt") {
    return(sqrt(x))  # sqrt(x)
  } else if (function_type == "linear") {
    return(x)  # x
  } else if (function_type == "square") {
    return(x^2)  # x^2
  } else if (function_type == "exp") {
    return((exp(x) - 1) / (exp(1) - 1))  # (e^x-1)/(e-1)
  } else {
    stop("Invalid type")
  }
}

# Set range of x values
x <- seq(0, 1, length.out = 100)

# Calculate the values of each function
y1 <- sapply(x, function(x) Flexfunction("log", x))
y2 <- sapply(x, function(x) Flexfunction("sqrt", x))
y3 <- sapply(x, function(x) Flexfunction("linear", x))
y4 <- sapply(x, function(x) Flexfunction("square", x))
y5 <- sapply(x, function(x) Flexfunction("exp", x))

# Plot the functions
plot(x, y1, type = "l", col = "red", lwd = 2, ylim = c(0, 1),
     xlab = "x", ylab = "y", main = "Flexfunction Visualization")
lines(x, y2, col = "blue", lwd = 2)
lines(x, y3, col = "green", lwd = 2)
lines(x, y4, col = "purple", lwd = 2)
lines(x, y5, col = "orange", lwd = 2)
# Add a legend
legend("bottomright", 
       legend = c("log(x+1)/log 2", "sqrt(x)", "x", "x^2", "(e^x-1)/(e-1)"),
       col = c("red", "blue", "green", "purple", "orange"), 
       lwd = 2, 
       cex = 0.7)

# Calculate R_(t)
# Modified Reject function
Reject <- function(t, method, function_type = NULL) {
  if (method == "BC" || method == "BH"||method == "ST") {
    R <- sum(p_values <= t)
  } else if (method == "FBC" || method == "FBH") {
    phi <- Flexfunction(function_type, p_values)
    R <- sum(phi <= t)
  }
  return(R)
}

# Null proportion function
Npropo <- function(lambda) {
  Np = (1 + n_genes - Reject(lambda, "ST")) / ((1 - lambda) * n_genes)
  return(Np)
}

# Calculate m_(t)
m_t <- function(t, method, function_type = NULL, lambda = NULL) {
  if (method == "BH") {
    m = n_genes * t
  } else if (method == "BC") {
    m = 1 + sum(p_values >= 1 - t)
  } else if (method == "FBH") {
    target_function <- function(p) {
      return(Flexfunction(function_type, p) - t)
    }
    ur = uniroot(target_function, interval = c(0, 1))
    m = n_genes * ur$root
  } else if (method == "FBC") {
    phi <- Flexfunction(function_type, 1 - p_values)
    m = 1 + sum(phi <= t)
  } else if (method == "ST") {
    m = n_genes * Npropo(lambda) * t
  }
  return(m)
}

# Calculate threshold T
threshold <- function(alpha, method, function_type = NULL, lambda = NULL) {
  # Define interval D
  if (method == "ST") {
    D <- lambda
  } else if (method == "FBH" || method == "BH") {
    D <- 1
  } else {
    D <- 0.5
  }
  
  # Define constraint function
  resctrict_function <- function(t) {
    method_t <- m_t(t, method, function_type = function_type, lambda = lambda)
    reject_t <- Reject(t, method, function_type = function_type)
    return(method_t / max(1, reject_t) - alpha)
  }
  
  # Define objective function
  target_function <- function(t) {
    if (t < 0 || t > D || is.na(resctrict_function(t)) || resctrict_function(t) >= 0) {
      return(1e6) # Penalty value
    } else {
      return(-t) # Maximize t
    }
  }
  
  # Use nloptr for optimization
  op <- nloptr::nloptr(
    x0 = 0, # Starting point
    eval_f = target_function,
    lb = 0, # Lower bound
    ub = D, # Upper bound
    opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-20)
  )
  
  # Check optimization result
  if (op$objective >= 1e6) {
    T <- 0 # If objective function value is too high, return 0
  } else {
    T <- -op$objective # Reverse to positive value
  }
  
  return(T)
}

# Output results table
results_base <- tibble(
  Gene_ID = 1:n_genes,
  True_Effect = true_effect,
  Z_Score = z_scores,
  P_Value = p_values
)
all_results <- list()
significants <- list()

# Calculate performance metrics
calculate_performance <- function(results, T) {
  true_positive <- sum(results$Significant & results$True_Effect)   # Correct discoveries
  false_positive <- sum(results$Significant & !results$True_Effect) # Incorrect discoveries
  false_negative <- sum(!results$Significant & results$True_Effect) # Undetected true signals
  true_negative <- sum(!results$Significant & !results$True_Effect) # Correct rejections
  
  # FDR (False Discovery Rate) = FP / (FP + TP)
  fdr <- false_positive / (max(0,false_positive + true_positive))
  find_power <- true_positive / (true_positive + false_negative) # Statistical power
  Reject_Rate <-(false_positive + true_positive)/n_genes
  RejectS=RejectScore(false_positive,true_positive)
  # Return performance metrics
  performance <- tibble(
    Metric = c("True Positive", "False Positive", "True Negative", "False Negative", "Threshold","Reject Rate","Reject Score", "FDR", "Finding Power"),
    Value = c(true_positive, false_positive, true_negative, false_negative, T,Reject_Rate, RejectS,fdr, find_power)
  )
  
  return(performance)
}

# Loop for method determination
for (method in methods) {
  if (method == "ST") {
    for (lambda in lambdas) {
      threshold_value <- threshold(alpha, method, lambda = lambda)
      significant <- p_values <= threshold_value
      significants[[paste(method, "lambda", lambda, sep = "_")]] <- significant
      results <- bind_cols(results_base, tibble(
        Threshold_T = threshold_value,
        Significant = significant))
      performance <- calculate_performance(results, threshold_value)
      all_results[[paste(method, "lambda", lambda, sep = "_")]] <- list(
        results = results,
        performance = performance
      )
    }
  } else if (method == "FBC" | method == "FBH") {
    for (function_type in function_types) {
      threshold_value <- threshold(alpha, method, function_type = function_type)
      significant <- p_values <= threshold_value
      significants[[paste(method, function_type, sep = "_")]] <- significant
      results <- bind_cols(results_base, tibble(
        Threshold_T = threshold_value,
        Significant = significant
      ))
      
      # Calculate performance metrics
      performance <- calculate_performance(results, threshold_value)
      
      # Save to results table
      all_results[[paste(method, function_type, sep = "_")]] <- list(
        results = results,
        performance = performance)
    }
  } else {
    threshold_value <- threshold(alpha, method)
    significant <- p_values <= threshold_value
    significants[[method]] <- significant
    results <- bind_cols(results_base, tibble(
      Threshold_T = threshold_value,
      Significant = significant
    ))
    
    # Calculate performance metrics
    performance <- calculate_performance(results, threshold_value)
    
    # Save to results table
    all_results[[method]] <- list(
      results = results,
      performance = performance
    )
  }
}

# Print all method results and performance metrics
for (result_key in names(all_results)) {
  cat("Method:", result_key, "\n")
  print(all_results[[result_key]]$results)
  print(all_results[[result_key]]$performance)
  cat("\n\n")
}

# Visualize the results
# Count the categorized results by method and add method name to the data frame
category_summary_all <- bind_rows(
  lapply(names(all_results), function(result_key) {
    all_results[[result_key]]$results %>%
      mutate(Method = result_key) %>%  # Use result_key as method name
      mutate(
        Category = case_when(
          Significant & True_Effect ~ "True Positive",
          Significant & !True_Effect ~ "False Positive",
          !Significant & True_Effect ~ "False Negative",
          !Significant & !True_Effect ~ "True Negative"
        )
      ) %>%
      count(Method, Category)  # Count by method and category
  })
)

# Calculate FDR and add it to the table
category_summary_table <- category_summary_all %>%
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) %>%
  mutate(
    FDR = `False Positive` / max(1,(`False Positive` + `True Positive`)),  # Calculate FDR
    Find_Power = `True Positive` / (`True Positive` + `False Negative`),  # Calculate Find Power
    Reject_Rate = (`False Positive` + `True Positive`) / n_genes,  # Calculate Reject Rate
    RejectS = RejectScore(`False Positive`, `True Positive`)  # Calculate Reject Score
  )

print(category_summary_table)

# Combine all results into a large data frame and add method column
visualization_data <- bind_rows(
  lapply(names(all_results), function(result_key) {
    all_results[[result_key]]$results %>%
      mutate(Method = result_key) %>%  # Use result_key as method name
      mutate(
        Category = case_when(
          Significant & True_Effect ~ "True Positive",
          Significant & !True_Effect ~ "False Positive",
          !Significant & True_Effect ~ "False Negative",
          !Significant & !True_Effect ~ "True Negative"
        )
      )
  })
)

# Create facet plot, facet by method
visualization_data %>%
  ggplot(aes(x = Gene_ID, y = P_Value, color = Category)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(aes(yintercept = Threshold_T), linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(
    values = c(
      "True Positive" = "blue",
      "False Positive" = "orange",
      "True Negative" = "gray",
      "False Negative" = "red"
    )
  ) +
  labs(
    title = "Adjusted P-Value Analysis with Result Categories",
    x = "Gene ID",
    y = "Adjusted P-Value",
    color = "Category"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ Method, scales = "free_y")  # Facet by method

# Filtered visualization
visualization_data %>%
  # Filter the part where P-Value is between 0 and Threshold_T
  filter(P_Value <= Threshold_T) %>%
  ggplot(aes(x = Gene_ID, y = P_Value, color = Category)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(aes(yintercept = Threshold_T), linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(
    values = c(
      "True Positive" = "blue",
      "False Positive" = "orange",
      "True Negative" = "gray",
      "False Negative" = "red"
    )
  ) +
  ggtitle("P-Value Distribution") +  # Add global title at the center of the plot
  labs(
    x = "Gene ID",
    y = "P-Value",
    color = "Category"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 12),  # x-axis title
    axis.title.y = element_text(size = 12),  # y-axis title
    strip.background = element_rect(fill = "lightgray", color = "black"),  # Facet label background
    strip.text = element_text(size = 12),  # Facet label font size
    plot.title = element_text(hjust = 0.5)  # Center the global title
  ) +
  facet_wrap(~ Method, scales = "free_y", labeller = label_value)  # Use Method for faceting

