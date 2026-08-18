# Install the package if it has not been installed yet
install.packages("rmcorr")

# Load the package
library(rmcorr)

# Read the data
data <- read.csv("C:/Users/Desktop/CMCSMR.csv")

# View the first few rows
# head(data)

# Repeated-measures correlation between SRBS_MRBD and CMC
rmc <- rmcorr(
  participant = subject,
  measure1 = SRBS_MRBD,
  measure2 = CMC,
  dataset = data
)

rmc

# Plot the repeated-measures correlation
plot(rmc)


library(dplyr)

set.seed(1)

# Function to calculate the within-subject correlation
rmcorr_r <- function(d) {
  
  d <- d |>
    group_by(subject) |>
    mutate(
      CMC_c = CMC - mean(CMC),
      SRBS_MRBD_c = SRBS_MRBD - mean(SRBS_MRBD)
    ) |>
    ungroup()
  
  cor(d$SRBS_MRBD_c, d$CMC_c)
}

# Observed within-subject correlation
r_obs <- rmcorr_r(data)

# Number of permutations
B <- 10000

# Within-subject permutation test:
# SRBS_MRBD values are shuffled separately within each subject
r_perm <- replicate(B, {
  
  dp <- data |>
    group_by(subject) |>
    mutate(
      SRBS_MRBD = sample(SRBS_MRBD)
    ) |>
    ungroup()
  
  rmcorr_r(dp)
})

# Two-sided permutation p-value
p_perm <- mean(
  abs(r_perm) >= abs(r_obs)
)

# Print the observed correlation and permutation p-value
c(
  r_obs = r_obs,
  p_perm = p_perm
)


# Refit the repeated-measures correlation with CMC on the x-axis
# and SRBS_MRBD on the y-axis
rm2 <- rmcorr(
  participant = subject,
  measure1 = CMC,
  measure2 = SRBS_MRBD,
  dataset = data
)

# Reduce plot margins
op <- par(
  mar = c(4, 4, 1, 1)
)

# Plot with specified axis limits
plot(
  rm2,
  xlim = c(0.03, 0.18),
  ylim = c(0, 140),
  xaxs = "i",
  yaxs = "i"
)

# Restore the original graphical parameters
par(op)