# # def generate_full_factorial_combinations(values):
# #     """
# #     Generates a combination of parameters based on a full factorial design.
    
# #     :param values: list of values for each parameter.
# #     :return: Array of parameter combinations.
# #     """
# #     # Calculate the number of levels for each parameter
# #     levels = [len(vals) for vals in values]
    
# #     # Generate full factorial design matrix
# #     design = fullfact(levels)
    
# #     # Conversion of design matrices to actual parameter values
# #     scaled_design = np.array([[values[j][int(i)] for j, i in enumerate(row)] for row in design])

# #     return scaled_design

# # Define parameter
# Purkinje_edgelen_values = [1, 2, 3]
# Purkinje_nsplit_values = [1, 2, 3]
# EP_sigma_values = [0.1, 0.2, 0.3, 0.4, 0.5] # sigma11 sigma22 sigma33
# Fiber_alpha_beta_values = [-100, -101, -102, -103] #alpha_beta_values

# # values = [
# #     Purkinje_edgelen_values, 
# #     Purkinje_nsplit_values, 
# #     EP_sigma_values, 
# #     EP_sigma_values, 
# #     EP_sigma_values, 
# #     Fiber_alpha_beta_values, 
# #     Fiber_alpha_beta_values
# # ]

# combinations = [
#     (edgelen, nsplit, sigma11, sigma22, sigma33, alpha, beta)
#     for edgelen in Purkinje_edgelen_values
#     for nsplit in Purkinje_nsplit_values
#     for sigma11 in EP_sigma_values
#     for sigma22 in EP_sigma_values
#     for sigma33 in EP_sigma_values
#     for alpha in Fiber_alpha_beta_values
#     for beta in Fiber_alpha_beta_values
# ]

# print('combinations', combinations)
# print('nbr combination:', len(combinations))



# import pandas as pd

# df = pd.DataFrame(combinations, columns=["Purkinje_edgelen", "Purkinje_nsplit", "EP_sigma11", "EP_sigma22", "EP_sigma33", "fiber_alpha", "fiber_beta"])
# csv_file_path = r'C:\Users\xuhu\pyheart-lib\examples\simulator\parameter_combinations.csv'
# df.to_csv(csv_file_path, index=False)


import numpy as np
from pyDOE import lhs

def scale_lhs_samples(samples, ranges):
    """ Scale Latin hypercube samples(0-1) to a specified range """
    return np.array([samples[:, i] * (ranges[i][1] - ranges[i][0]) + ranges[i][0] for i in range(samples.shape[1])]).T

def generate_lhs_combinations(num_samples=5):
    """    
    param num_samples: number of samples for each parameter.
    """
    ranges = [(0.5, 2), (0.2, 1), (3, 10), (5, 10)]
    lhs_samples = lhs(len(ranges), samples=num_samples)
    scaled_samples = scale_lhs_samples(lhs_samples, ranges)

    return scaled_samples

combinations = generate_lhs_combinations()

print('combinations', combinations)
print('nbr combination:', len(combinations))

import pandas as pd
df = pd.DataFrame(combinations, columns=["Purkinje_edgelen", "SigmaX", "Ratio", "Ratio2"])
csv_file_path = r'C:\Users\xuhu\pyheart-lib\examples\simulator\parameter_combinations.csv'
df.to_csv(csv_file_path, index=False)