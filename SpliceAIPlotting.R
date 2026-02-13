library(readr)
library(ggplot2)
library(tidyr)

#Set working directory
dir <- "/Users/rluo/Documents/TillotsonLab/SpliceAI/ATRX/SpliceAIPlotting"
setwd(dir)
list.files(pattern = "*.tsv")
df <- read_tsv("RL_Ex2FLEx_output_threshold0.50_long.tsv")
#Rename Columns
colnames(df) <- c("Sequence.ID",
                  "Position",
                  "Nucleotide",
                  "Acceptor.Score",
                  "Donor.Score")
unique(df$Sequence.ID)
#Split it into separate tables based on unique sequence ID
df_wt <- subset(df, Sequence.ID == "Ex2WT_for_SpliceAI_long")
df_FLEx <- subset(df, Sequence.ID == "Ex2FLEx_for_SpliceAI_long")
df_reco <- subset(df, Sequence.ID == "Ex2FLEx_Recombination_for_SpliceAI_long")

df_wt <- df_wt[,c(2,4,5)]
df_FLEx <- df_FLEx[,c(2,4,5)]
df_reco <- df_reco[,c(2,4,5)]
#Now plot
df_wt_long <- pivot_longer(df_wt, cols = c("Acceptor.Score", "Donor.Score"), 
                           names_to = "Type", values_to = "Probability")

# based on the sequence and the elements to preset
exon2.wt <- data.frame(
  xmin = c(1190),
  xmax = c(1302),
  ymin = 0.15,
  ymax = 0.1,
  label = c("ex.2")
)

exon3.wt <- data.frame(
  xmin = c(10211),
  xmax = c(10268),
  ymin = 0.15,
  ymax = 0.1,
  label = c("ex.3")
)

wt_plot <- ggplot() +
  geom_point(data = df_wt_long, aes(x = Position, y = Probability, color = Type), alpha = 0.6)+
  labs(title = "SpliceAI Predictions mATRX WT", x = "Position", y = "Splicing Probability") +
  scale_color_manual(values = c("red","blue"), labels = c("Acceptor", "Donor")) +  # Acceptor = blue, Donor = red
  ylim(c(0.1,1))+
  xlim(c(0,10560))+
  geom_hline(yintercept = 0.5, color = "grey20", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1190, color = "red", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1302, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 10211, color = "red", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 10268, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_rect(data = exon2.wt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_text(data = exon2.wt, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), 
            color = "white", size = 2) +
  geom_rect(data = exon3.wt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_text(data = exon3.wt, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), 
            color = "white", size = 2) +
  theme_classic()
wt_plot

df_FLEx_long <- pivot_longer(df_FLEx, cols = c("Acceptor.Score", "Donor.Score"),
                             names_to = "Type", values_to = "Probability")
exon2.flex <- data.frame(
  xmin = c(1351),
  xmax = c(1464),
  ymin = 0.15,
  ymax = 0.1,
  label = c("ex.2")
)

exon3.flex <- data.frame(
  xmin = c(10447),
  xmax = c(10504),
  ymin = 0.15,
  ymax = 0.1,
  label = c("ex.3")
)

flex_plot <- ggplot() +
  geom_point(data = df_FLEx_long, aes(x = Position, y = Probability, color = Type), alpha = 0.6)+
  labs(title = "SpliceAI Predictions mATRX FLEx", x = "Position", y = "Splicing Probability") +
  scale_color_manual(values = c("red","blue"), labels = c("Acceptor", "Donor")) +  # Acceptor = blue, Donor = red
  ylim(c(0.1,1))+
  xlim(c(0,10800))+
  geom_hline(yintercept = 0.5, color = "grey20", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1351, color = "red", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1464, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 10447, color = "red", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 10504, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_rect(data = exon2.flex, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_text(data = exon2.flex, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), 
            color = "white", size = 2) +
  geom_rect(data = exon3.flex, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_text(data = exon3.flex, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), 
            color = "white", size = 2) +
  theme_classic()
flex_plot

#Now plot recombination one. I am not plotting the FLEx one because it lost its donor and acceptor completely which is what we are expecting.
df_reco_long <- pivot_longer(df_reco, cols = c("Acceptor.Score", "Donor.Score"), 
                             names_to = "Type", values_to = "Probability")

# Position of Ex2 has changed after recombination
exon2.rec <- data.frame(
  xmin = c(1224),
  xmax = c(1336),
  ymin = 0.15,
  ymax = 0.1,
  label = c("ex.2")
)

exon3.rec <- data.frame(
  xmin = c(10279),
  xmax = c(10336),
  ymin = 0.15,
  ymax = 0.1,
  label = c("ex.3")
)

rec_plot <- ggplot() +
  geom_point(data = df_reco_long, aes(x = Position, y = Probability, color = Type), alpha = 0.6)+
  labs(title = "SpliceAI Predictions mATRX Recombination", x = "Position", y = "Splicing Probability") +
  scale_color_manual(values = c("red","blue"), labels = c("Acceptor", "Donor")) +  # Acceptor = blue, Donor = red
  ylim(c(0.1,1))+
  xlim(c(0,10700))+
  geom_hline(yintercept = 0.5, color = "grey20", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1224, color = "red", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 1336, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 10279, color = "red", linetype = "dashed", alpha = 0.2) +
  geom_vline(xintercept = 10336, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_rect(data = exon2.rec, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_text(data = exon2.rec, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), 
            color = "white", size = 2) +
  geom_rect(data = exon3.rec, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_text(data = exon3.rec, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), 
            color = "white", size = 2) +
  theme_classic()
rec_plot

library(gridExtra)
grid.arrange(wt_plot, flex_plot, rec_plot, ncol = 1)