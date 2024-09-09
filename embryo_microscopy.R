library(tidyverse)
library(viridis)

# Function to identify and remove outliers using IQR method
remove_outliers <- function(df, column) {
  Q1 <- quantile(df[[column]], 0.25)
  Q3 <- quantile(df[[column]], 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  df <- df %>% filter(df[[column]] >= lower_bound & df[[column]] <= upper_bound)
  return(df)
}

data <- read_csv("Statistics_v2.csv")

data$genotype <- as.factor(data$genotype)
data$dividing <- as.factor(data$dividing)

data <- data %>% filter(DAF != 3)


colour_mapping <- c("WT" = rgb(122/255,151/255,243/255), 
                    "dn4mt1-#6" = rgb(18/255,163/255,125/255), 
                    "dn4mt1-#21" = rgb(18/255,163/255,125/255))

# Apply outlier removal for each relevant column grouped by genotype
data_clean <- data %>%
  group_by(genotype) %>%
  do(remove_outliers(., "nuclei_number")) %>%
  ungroup()


# Summarize data to get the mean values and n (sample size) by genotype and DAF after cleaning
summary_data <- data_clean %>%
  group_by(genotype, DAF) %>%
  summarise(
    mean_nuclei_number = mean(nuclei_number),
    sem_nuclei_number = sd(nuclei_number) / sqrt(n()),
    n = n(),
    mean_longitudinal_diameter = mean(logitudinal_diameter_um),
    sem_longitudinal_diameter = sd(logitudinal_diameter_um) / sqrt(n()),
    mean_transverse_diameter = mean(transverse_diameter_um),
    sem_transverse_diameter = sd(transverse_diameter_um) / sqrt(n())
  )

# Summarize the total n for each genotype across all DAFs
genotype_n_summary <- data_clean %>%
  group_by(genotype) %>%
  summarise(
    total_n = n()
  )

# Create custom labels for the genotypes including the total sample size (total_n)
custom_labels <- genotype_n_summary %>%
  mutate(label = paste(genotype, "(n =", total_n, ")"))

# Create a named vector for the labels to use in the plot
label_vector <- setNames(custom_labels$label, custom_labels$genotype)



# Filter the data to focus on the early days (1-4 DAF)
early_data <- data %>% filter(DAF >= 1 & DAF <= 5)

early_data_clean <- early_data %>%
  group_by(genotype) %>%
  do(remove_outliers(., "nuclei_number")) %>%
  ungroup()

# Summarize data to get the mean values by genotype and DAF after cleaning
early_summary_data <- early_data_clean %>%
  group_by(genotype, DAF) %>%
  summarise(
    mean_nuclei_number = mean(nuclei_number, na.rm = TRUE),
    sem_nuclei_number = sd(nuclei_number, na.rm = TRUE) / sqrt(n()),
    mean_longitudinal_diameter = mean(logitudinal_diameter_um, na.rm = TRUE),
    sem_longitudinal_diameter = sd(logitudinal_diameter_um, na.rm = TRUE) / sqrt(n()),
    mean_transverse_diameter = mean(transverse_diameter_um, na.rm = TRUE),
    sem_transverse_diameter = sd(transverse_diameter_um, na.rm = TRUE) / sqrt(n())
  )




# Plot 1: Nuclei number across genotypes and DAF (Line plot with error bars)

ggplot(summary_data, aes(x = DAF, y = mean_nuclei_number, color = genotype)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_nuclei_number - sem_nuclei_number, ymax = mean_nuclei_number + sem_nuclei_number), width = 0.2) +
  scale_color_manual(values = colour_mapping, labels = label_vector) +
  labs(title = "Nucleus Number Across Genotypes and Developmental Stages",
       x = "Days After Fertilization (DAF)",
       y = "Mean Nucleus Number") +
  scale_x_continuous(breaks = c(unique(summary_data$DAF), 3)) + # Adding 3 to the x-axis
  theme_classic() +
  theme(legend.position = "right")


ggplot(early_summary_data, aes(x = DAF, y = mean_nuclei_number, color = genotype)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_nuclei_number - sem_nuclei_number, ymax = mean_nuclei_number + sem_nuclei_number), width = 0.2) +
  scale_color_manual(values = colour_mapping) +
  labs(title = "Nuclei Number Across Genotypes and Developmental Stages (DAF 1-5)",
       x = "Days After Fertilization (DAF)",
       y = "Mean Nuclei Number") +
  theme_minimal()



# Plot 2: Longitudinal diameter across genotypes and DAF
ggplot(summary_data, aes(x = DAF, y = mean_longitudinal_diameter, color = genotype)) +
  geom_line() +
  geom_point() +
  labs(title = "Longitudinal Diameter Across Genotypes and Developmental Stages",
       x = "Days After Fertilization (DAF)",
       y = "Mean Longitudinal Diameter (um)") +
  theme_minimal()

# Plot 3: Transverse diameter across genotypes and DAF
ggplot(summary_data, aes(x = DAF, y = mean_transverse_diameter, color = genotype)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = colour_mapping) +
  scale_x_continuous(breaks = c(unique(summary_data$DAF), 3)) +
  labs(title = "Transverse Diameter Across Genotypes and Developmental Stages",
       x = "Days After Fertilization (DAF)",
       y = "Mean Transverse Diameter (um)") +
  theme_classic()


# Boxplot comparisons

# Boxplot 1: Nuclei number by genotype, split by DAF
ggplot(data_clean, aes(x = genotype, y = nuclei_number, fill = genotype)) +
  geom_boxplot() +
  facet_wrap(~DAF) +
  scale_color_manual(values = colour_mapping) +
  labs(title = "Nuclei Number by Genotype Split by Developmental Stage",
       x = "Genotype",
       y = "Nuclei Number") +
  theme_minimal()


# Calculate the proportion of dividing cells by genotype and DAF
division_summary <- early_data %>%
  group_by(genotype, DAF) %>%
  summarise(
    total_cells = n(),
    dividing_cells = sum(dividing == "TRUE"),
    proportion_dividing = dividing_cells / total_cells
  )

ggplot(division_summary, aes(x = factor(DAF), y = proportion_dividing, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportion of Dividing Cells by Genotype and Developmental Stage (1-4 DAF)",
       x = "Days After Fertilization (DAF)",
       y = "Proportion of Dividing Cells") +
  theme_minimal()

# Line Plot with Points
ggplot(division_summary, aes(x = DAF, y = proportion_dividing, color = genotype)) +
  geom_line() +
  geom_point(size = 4) +
  scale_x_continuous(breaks = c(unique(summary_data$DAF), 3)) +
  scale_color_manual(values = colour_mapping) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportion of Dividing Cells by Genotype and Developmental Stage (1-5 DAF)",
       x = "Days After Fertilization (DAF)",
       y = "Proportion of Dividing Cells") +
  theme_classic()


ggplot(division_summary, aes(x = factor(DAF), y = genotype, fill = proportion_dividing)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue", labels = scales::percent_format()) +
  labs(title = "Heatmap of Dividing Cells by Genotype and Developmental Stage (1-4 DAF)",
       x = "Days After Fertilization (DAF)",
       y = "Genotype") +
  theme_minimal()


# Calculate circularity
data_clean <- data_clean %>%
  mutate(circularity = transverse_diameter_um / logitudinal_diameter_um)

# Summarize the data: calculate mean and standard error of the mean (SEM)
circularity_summary <- data_clean %>%
  group_by(genotype, DAF) %>%
  summarise(
    mean_circularity = mean(circularity, na.rm = TRUE),
    sem_circularity = sd(circularity, na.rm = TRUE) / sqrt(n())
  )

# Plot Circularity Over Time (DAF) with Error Bars
ggplot(circularity_summary, aes(x = DAF, y = mean_circularity, color = genotype)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = colour_mapping) +
  scale_x_continuous(breaks = c(unique(summary_data$DAF), 3)) +
  geom_errorbar(aes(ymin = mean_circularity - sem_circularity, ymax = mean_circularity + sem_circularity), width = 0.2) +
  labs(title = "Mean Circularity Over Time by Genotype",
       x = "Days After Fertilization (DAF)",
       y = "Mean Circularity (Transverse / Longitudinal)") +
  theme_classic()

# Create bar graph with error bars
ggplot(summary_data, aes(x = factor(DAF), y = mean_nuclei_number, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_nuclei_number - sem_nuclei_number, ymax = mean_nuclei_number + sem_nuclei_number), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(title = "Mean Number of Nuclei by Genotype and DAF",
       x = "Days After Fertilization (DAF)",
       y = "Mean Nuclei Number") +
  scale_fill_manual(values = colour_mapping) +
  theme_minimal() +
  theme(legend.position = "top")

# Create bar graph with error bars
ggplot(early_summary_data, aes(x = factor(DAF), y = mean_nuclei_number, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_nuclei_number - sem_nuclei_number, ymax = mean_nuclei_number + sem_nuclei_number), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(title = "Mean Number of Nuclei by Genotype and DAF",
       x = "Days After Fertilization (DAF)",
       y = "Mean Nuclei Number") +
  scale_fill_manual(values = colour_mapping) +
  theme_minimal() +
  theme(legend.position = "top")



