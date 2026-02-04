# This is to construct a linear regression model to show that the farther away the LoxP sites, the better the donor score predicted by SpliceAI.
library(tidyverse)
library(ggplot2)
sgRNA <- read.csv("Intron2sgRNA.csv")
head(sgRNA)
sgRNA <- sgRNA %>%
  rename(donor_score = Splice.AI.donor.score)
sgRNA <- sgRNA %>%
  rename(position = Position)
sgRNA <- sgRNA %>%
  mutate(position = as.numeric(position),
         donor_score = as.numeric(donor_score))
model <- lm(donor_score ~ position, data = sgRNA)
summary(model)
plot(model)
plot <- ggplot(sgRNA, aes(x = position, y = donor_score)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(
    data = sgRNA[sgRNA$position == 29658, ],
    aes(position, donor_score),
    color = "magenta",
    size = 3
  ) +
  annotate("text",
           x = -Inf,
           y = 0.325,
           label = "original sgRNA",
           hjust = -0.1,
           vjust = -0.7,
           color = "magenta") +
  geom_hline(yintercept = 0.515, color = "red", linetype = "dashed", alpha = 0.5) +
  annotate("text",
           x = -Inf,
           y = 0.515,
           label = "WT donor score",
           hjust = -0.1,
           vjust = -0.5,
           color = "red") +
  theme_minimal() +
  ggtitle("sgRNA Designs at Intron 2")
plot
