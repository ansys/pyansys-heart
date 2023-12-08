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


def generate_lhs_combinations(num_samples=50):
    """num_samples: number of samples for each parameter."""
    ranges = [(0.5, 2), (0.2, 1), (3, 10), (5, 10)]
    lhs_samples = lhs(len(ranges), samples=num_samples)
    scaled_samples = scale_lhs_samples(lhs_samples, ranges)

    return scaled_samples


combinations = generate_lhs_combinations()

print("combinations", combinations)
print("nbr combination:", len(combinations))

df = pd.DataFrame(combinations, columns=["Purkinje_edgelen", "SigmaX", "Ratio", "Ratio2"])
csv_file_path = r".\parameter_combinations.csv"
df.to_csv(csv_file_path, index=False)