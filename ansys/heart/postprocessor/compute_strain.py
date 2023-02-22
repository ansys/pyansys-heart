"""Compute myocardial strain."""
from ansys.heart.preprocessor.models import HeartModel
import numpy as np


def compute_myocardial_strain(model: HeartModel, deformation_gradient, reference=None):
    """Compute left ventricle myocardial strain."""
    # model info
    e_l, e_r, e_c = model.compute_left_ventricle_element_cs()
    ele_id = np.where(~np.isnan(model.aha_ids))[0]

    strain = np.zeros((len(ele_id), 3))
    def_grad = deformation_gradient[ele_id]
    if reference is not None:
        pass

    # todo: vectorization
    for i_ele in range(len(ele_id)):
        if reference is not None:
            pass

        else:
            right_cauchy_green = np.matmul(
                def_grad[i_ele, :].reshape(3, 3),
                def_grad[i_ele, :].reshape(3, 3).T,
            )

        # Green Lagrangian strain: E = 0.5*(lambda**2-1)
        # lambda = sqrt(e*right_cauchy_green*e)
        strain[i_ele, 0] = 0.5 * (
            np.matmul(np.matmul(e_l[i_ele].T, right_cauchy_green), e_l[i_ele]) - 1
        )
        strain[i_ele, 1] = 0.5 * (
            np.matmul(np.matmul(e_r[i_ele].T, right_cauchy_green), e_r[i_ele]) - 1
        )
        strain[i_ele, 2] = 0.5 * (
            np.matmul(np.matmul(e_c[i_ele].T, right_cauchy_green), e_c[i_ele]) - 1
        )

    return strain


def compute_AHA17_segment_strain(model: HeartModel, element_strain):
    """Average elemental strain for AHA17 segments."""
    aha_strain = np.zeros((17, 3))

    # get aha17 label for left ventricle elements
    aha17_label = model.aha_ids[~np.isnan(model.aha_ids)]

    for i in range(1, 18):
        # get index in strain table
        indices = np.where(aha17_label == i)[0]
        # average
        aha_strain[i - 1] = np.mean(element_strain[indices, :], axis=0)
    return aha_strain
