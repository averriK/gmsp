library(expm)
library(MASS)

# Example values
dt <- 0.01
K <- 1.0
C <- 0.1
s <- rnorm(1000)  # Random data

# Define A
A <- matrix(c(0, 1, -K, -C), nrow = 2, byrow = TRUE)

# Exponential of A*dt
Ae <- expm(A * dt)

# Calculate AeB
AeB <- (Ae - diag(2)) %*% ginv(A) %*% matrix(c(0, 1), nrow = 2)

# Initialize y
y <- matrix(0, nrow = 2, ncol = length(s))

# Update y over time
for (k in 2:length(s)) {
  y[, k] <- Ae %*% y[, k-1] + AeB * s[k]
}

# Calculate displ_max
displ_max <- max(abs(y[1, ]))  # Spectral relative displacement (m)

# Calculate PSA
omega_n <- 2 * pi / 0.5  # Example period
PSA <- displ_max * omega_n^2  # Pseudo spectral acceleration (m/s^2)

print(PSA)
