# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""

Define materials
----------------
This example show you how to create a mechanical material and assign it to a heart part.
"""

###############################################################################
# Import material module
# ~~~~~~~~~~~~~~~~~~~~~~

from pathlib import Path

import matplotlib.pyplot as plt

from ansys.health.heart.settings.material.curve import (
    ActiveCurve,
    constant_ca2,
    kumaraswamy_active,
)
from ansys.health.heart.settings.material.ep_material import CellModel, EPMaterial
from ansys.health.heart.settings.material.material import (
    ACTIVE,
    ANISO,
    ISO,
    ActiveModel,
    Mat295,
    NeoHookean,
)

###############################################################################
# .. note::
#    Unit system used for heart modeling in LS-DYNA is ["MPa", "mm", "N", "ms", "g"]


###############################################################################
# Create a material
# ~~~~~~~~~~~~~~~~~

# Neo-Hookean material can be created as following
neo = NeoHookean(rho=0.001, c10=1, nu=0.499)
###############################################################################
## The recommended approach is to create a Neo-Hookean material by
# activating only the isotropic module in MAT295.
neo2 = Mat295(rho=0.001, iso=ISO(itype=1, beta=2, kappa=1, mu1=0.05, alpha1=2))

###############################################################################
# .. note::
#    Please refer to LS-DYNA manual for more details of MAT_295

# More steps to create MAT295 which is used for myocardium

# step 1: create an isotropic module
iso = ISO(k1=1, k2=1, kappa=100)

# step 2: create an anisotropoc moddule
fiber = ANISO.HGOFiber(k1=1, k2=1)
aniso1 = ANISO(fibers=[fiber])

# Create fiber with sheet, and their interactions
sheet = ANISO.HGOFiber(k1=1, k2=1)
aniso2 = ANISO(fibers=[fiber, sheet], k1fs=1, k2fs=1)

# step3: create the active module

# example 1:
# create active model 1
ac_model1 = ActiveModel.Model1()
# create Ca2+ curve
ac_curve1 = ActiveCurve(constant_ca2(tb=800, ca2ionm=ac_model1.ca2ionm), type="ca2", threshold=0.5)
# build active module
active = ACTIVE(model=ac_model1, ca2_curve=ac_curve1)

## Active model 1 needs a constant ca2ion
# but the curve needs to cross threshold at every start of heart beat

# You can plot Ca2+ with threshold
fig = active.ca2_curve.plot_time_vs_ca2()
plt.show()

# example 2
# create active model 3
ac_model3 = ActiveModel.Model3()
# create a stress curve and show
ac_curve3 = ActiveCurve(kumaraswamy_active(t_end=800), type="stress")
fig = ac_curve3.plot_time_vs_stress()
plt.show()

###############################################################################
# .. note::
#   With setting eta=0 is model 3, stress curve will be the active stress for all elements.
#   If eta!=0, this is idealized active stress when fiber stretch stays to 1.

# PyAnsys-Heart will convert the stress curve to Ca2+ curve (input of MAT_295)
fig = ac_curve3.plot_time_vs_ca2()
plt.show()

# build active module
active3 = ACTIVE(model=ac_model3, ca2_curve=ac_curve3)

###############################################################################
# Finally, MAT295 can be created with the above modules
iso_mat = Mat295(rho=1, iso=iso, aniso=None, active=None)
passive_mat = Mat295(rho=1, iso=iso, aniso=aniso1, active=None)
active_mat = Mat295(rho=1, iso=iso, aniso=aniso1, active=active)

###############################################################################
# EP materials can be created as follows
ep_mat_active = EPMaterial.Active(
    sigma_fiber=1, sigma_sheet=0.5, beta=140, cm=0.01, cell_model=CellModel.Tentusscher()
)
epinsulator = EPMaterial.Insulator()

# import pyvista as pv


# ##############################################################################
# .. note::
#    Ca2+ curve will be ignored if the simulation is coupled with electrophysiology
#
# ##############################################################################
# Assign material to a part
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign the materials to the heart model


###############################################################################
# Load a heart model
###############################################################################
# .. note::
#    You need to complete the full heart preprocessing example first.
import ansys.health.heart.examples as examples
import ansys.health.heart.models as models

heart_model_vtu, heart_model_partinfo, _ = examples.get_preprocessed_fullheart()
workdir = str(Path.home() / "pyansys-heart" / "Rodero2021")

# load a full heart model.
heartmodel: models.FullHeart = models.FullHeart(working_directory=workdir)
heartmodel.load_model_from_mesh(heart_model_vtu, heart_model_partinfo)

heartmodel.mesh.set_active_scalars("_volume-id")
heartmodel.mesh.plot()

# Print default materials
print(heartmodel.left_ventricle.meca_material)
print(heartmodel.left_ventricle.ep_material)

###############################################################################
# .. note::
#    If no material is set before writing k files, default material from ```settings```
# will be set.

# Assign the material we just created
heartmodel.left_ventricle.meca_material = active_mat
heartmodel.left_ventricle.ep_material = ep_mat_active

print(heartmodel.left_ventricle.meca_material)
print(heartmodel.left_ventricle.ep_material)
###############################################################################
# Create a new part and set material

# # A new part can be created by elements IDs
# ids = np.where(heartmodel.mesh.point_data_to_cell_data()["uvc_longitudinal"] > 0.9)[0]
# new_part: Part = heartmodel.create_part_by_ids(ids, "new_part")

# # Show the part
# plotter = heartmodel.plot_part(new_part)
