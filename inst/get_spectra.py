import numpy as np
from scipy.linalg import expm
from numpy.linalg import pinv

# Example values
dt = 0.01
K = 1.0
C = 0.1
data = np.random.rand(1000)

# Define A
A = np.array([[0, 1], [-K, -C]])

# Exponential of A*dt
Ae = expm(A * dt)

# Calculate AeB
AeB = (Ae - np.eye(2)) @ pinv(A) @ np.array([[0], [1]])

# Initialize y
y = np.zeros((2, len(data)))

# Update y over time
for k in range(1, len(data)):
    y[:, k] = Ae @ y[:, k-1] + AeB.flatten() * data[k]

# Calculate displ_max
displ_max = np.max(np.abs(y[0, :]))  # Spectral relative displacement (m)

# Calculate PSA
omega_n = 2 * np.pi / 0.5  # Example period
PSA = displ_max * omega_n ** 2  # Pseudo spectral acceleration (m/s^2)

print(PSA)
