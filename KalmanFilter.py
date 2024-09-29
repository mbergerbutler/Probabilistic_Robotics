import numpy as np

A = np.array([1])
B = np.array([1])
eta = np.array([0.01])




class Robot():
    def __init__(self, A, B, eta):
        self.A = A
        self.B = B
        self.eta = eta
        self.mu = 0
        self.sigma = 0
        
    def update(self, u, z):
        [self.mu, self.sigma] = kalman_filter(self, u, z)

def kalman_filter(mu_prev, sigma_prev, u, z, A, B, R, C, Q)
    mu_prediction = A * mu_prev + B * u
    sigma_prediction = A * sigma_prev * np.transpose(A) + R
    
    K = sigma_prediction * np.transpose(C) * (C * sigma_prediction * np.transpose(C) + Q)^-1
    mu = mu_prediction + K * (z - C * mu_prediction)
    sigma = (np.identity(A.shpae[0]) - K * C) * sigma_prediction
    return mu, sigma
    
## State transition function
# x_t = A_t * x_t-1 + B_t * u_t + eta_t

# - x and u must have same dimensions nx1 
# - A_t is a square matrix of size nXn where n is the dimension of the state vector
# - B_t is of size nxm with m being the dimension of the control vector u_t
# - eta_t is a Gaussian random vecotr that models the uncertainty introduced by the state transition.
#   It is of the same dimension as the state vecotr (nx1). It's mean is zero and covariance is denoted R_t
# - z_t = C_t*x_t + delta_t
# - C_t is matrix of size kxn where k is the dimension of the measurment vector z_t
# - delta_t is a multivariate Gaussian with zero mean and covariance Q_t

# Kalman_filter(mu_t-1, sigma_t-1, u_t, z_t)
    ## Prediction Step
    # mu_prediction = A_t * mu_t-1 + B_t * u_t
    # sigma_prediction = A_t * sigma_t-1 * transpose(A_t) + R_t
    ## Update Step
    # K_t = sigma_prediction * transpose(C_t) * (C_t * sigma_prediction * transpose(C_t) + Q_t)^-1 # Kalman Gain
    # mu_t = mu_prediction + K_t * (z_t - C_t * mu_prediction)
    # sigma_t = (I - K_t * C_t) * sigma_prediction
    # Return mu_t, sigma_t