import numpy as np
import pandas as pd
from pyDOE import lhs


def scale_lhs_samples(samples, ranges):
    """Scale Latin hypercube samples(0-1) to a specified range"""
    return np.array(
        [
            samples[:, i] * (ranges[i][1] - ranges[i][0]) + ranges[i][0]
            for i in range(samples.shape[1])
        ]
    ).T


def generate_lhs_combinations(num_samples=150):
    """num_samples: number of samples for each parameter."""
    '''
    sigmaX values: (0.2, 2)
    ratio2 (in keywords.EmMat001): [1,10] (such that: sigmaPurkinje=ratio2 * sigmaX)
    '''
    ranges = [(0.2, 2), (1, 10)]

    lhs_samples = lhs(len(ranges), samples=num_samples)
    scaled_samples = scale_lhs_samples(lhs_samples, ranges)

    return scaled_samples

def plot_DoE_with_two_params():
    '''verification the algorithm of DoE points generated: normalization'''
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(figsize=(10, 6))

    # Plotting the data
    plt.scatter(df['sigmaX'], df['Ratio2'])
    plt.xlabel('sigmaX')
    plt.ylabel('Ratio2')
    plt.grid(True)
    plt.show()

combinations = generate_lhs_combinations()

print("combinations", combinations)
print("nbr combination:", len(combinations))

# df = pd.DataFrame(combinations, columns=["Purkinje_edgelen", "SigmaX", "Ratio", "Ratio2"])
df = pd.DataFrame(combinations, columns=["sigmaX", "Ratio2"])

csv_file_path = r"./two_parameter_combinations.csv"
csv_file_path = r"D:\xuhu\pyansys-heart\examples\simulator\two_parameter_combinations150.csv"
df.to_csv(csv_file_path, index=False)

plot_DoE_with_two_params()