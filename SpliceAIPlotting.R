library(readr)
library(ggplot2)
library(tidyr)

#Set working directory
dir <- "/Users/rluo/Documents/TillotsonLab/SpliceAI/ATRX/SpliceAIPlotting"
setwd(dir)
list.files(pattern = "*.tsv")
df <- read_tsv("RL_Ex2FLEx_output_threshold0.25_new_cons.tsv")
#Rename Columns
colnames(df) <- c("Sequence.ID",
                  "Position",
                  "Nucleotide",
                  "Acceptor.Score",
                  "Donor.Score")
unique(df$Sequence.ID)
#Split it into separate tables based on unique sequence ID
df_wt <- subset(df, Sequence.ID == "Ex2WT_for_SpliceAI(New)")
df_reco <- subset(df, Sequence.ID == "Ex2FLEx_Recombination_for_SpliceAI(New)")

df_wt <- df_wt[,c(2,4,5)]
df_reco <- df_reco[,c(2,4,5)]
#Now plot
df_wt_long <- pivot_longer(df_wt, cols = c("Acceptor.Score", "Donor.Score"), 
                           names_to = "Type", values_to = "Probability")

# based on the sequence and the elements to preset
exons.wt <- data.frame(
  xmin = c(1785),
  xmax = c(1897),
  ymin = 0.15,
  ymax = 0.1,
  label = c("ex.2")
)

wt_plot <- ggplot() +
  geom_point(data = df_wt_long, aes(x = Position, y = Probability, color = Type), alpha = 0.6)+
  labs(title = "SpliceAI Predictions mATRX WT", x = "Position", y = "Splicing Probability") +
  scale_color_manual(values = c("red","blue"), labels = c("Acceptor", "Donor")) +  # Acceptor = blue, Donor = red
  ylim(c(0.1,1))+
  xlim(c(0,2603))+
  geom_hline(yintercept = 0.5, color = "grey20", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1785, color = "red", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1897, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_rect(data = exons.wt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_text(data = exons.wt, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), 
            color = "white", size = 2) +
  theme_classic()
wt_plot

#Now plot recombination one. I am not plotting the FLEx one because it lost its donor and acceptor completely which is what we are expecting.
df_reco_long <- pivot_longer(df_reco, cols = c("Acceptor.Score", "Donor.Score"), 
                             names_to = "Type", values_to = "Probability")

# Position of Ex2 has changed after recombination
exons.rec <- data.frame(
  xmin = c(1819),
  xmax = c(1931),
  ymin = 0.15,
  ymax = 0.1,
  label = c("ex.2")
)

rec_plot <- ggplot() +
  geom_point(data = df_reco_long, aes(x = Position, y = Probability, color = Type), alpha = 0.6)+
  labs(title = "SpliceAI Predictions mATRX Recombination", x = "Position", y = "Splicing Probability") +
  scale_color_manual(values = c("red","blue"), labels = c("Acceptor", "Donor")) +  # Acceptor = blue, Donor = red
  ylim(c(0.1,1))+
  xlim(c(0,2671))+
  geom_hline(yintercept = 0.5, color = "grey20", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1819, color = "red", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1931, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_rect(data = exons.rec, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_text(data = exons.rec, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), 
            color = "white", size = 2) +
  theme_classic()
rec_plot

library(gridExtra)
grid.arrange(wt_plot, rec_plot, ncol = 2)