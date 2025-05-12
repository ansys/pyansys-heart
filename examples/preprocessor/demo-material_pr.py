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
This example shows how to create a mechanical material and assign it to a heart part.
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
#    The unit system for heart modeling in LS-DYNA is ``["MPa", "mm", "N", "ms", "g"]``.

###############################################################################
# Create a material
# ~~~~~~~~~~~~~~~~~
# Create a Neo-Hookean material as follows.
neo = NeoHookean(rho=0.001, c10=1, nu=0.499)
###############################################################################
# The recommended way to create a Neo-Hookean material is by
# activating only the isotropic module in MAT_295.
neo2 = Mat295(rho=0.001, iso=ISO(itype=1, beta=2, kappa=1, mu1=0.05, alpha1=2))

###############################################################################
# .. note::
#   For more information on MAT_295, see the `LS-DYNA manuals <https://lsdyna.ansys.com/manuals/>`_.

# Additional steps follow for creating MAT_295, which is used for myocardium.

# step 1: create an isotropic module
iso = ISO(k1=1, k2=1, kappa=100)

# step 2: create an anisotropic module
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

## Active model 1 must have a constant ca2ion,
# but the curve must cross the threshold at every start of the heart beat.

###############################################################################
# Plot Ca2+ with threshold
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Ca2+ with the threshold.
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
#   When eta=0 in model 3, the stress curve is the active stress for all elements.
#   If eta!=0, this is the idealized active stress when fiber stretch stays at 1.

# PyAnsys Heart converts the stress curve to Ca2+ curve (input of MAT_295)
fig = ac_curve3.plot_time_vs_ca2()
plt.show()

# build active module
active3 = ACTIVE(model=ac_model3, ca2_curve=ac_curve3)

###############################################################################
# Create MAT_295 with modules
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create MAT_295 with the preceding modules.
iso_mat = Mat295(rho=1, iso=iso, aniso=None, active=None)
passive_mat = Mat295(rho=1, iso=iso, aniso=aniso1, active=None)
active_mat = Mat295(rho=1, iso=iso, aniso=aniso1, active=active)

###############################################################################
# Create EP materials
# ~~~~~~~~~~~~~~~~~~~
ep_mat_active = EPMaterial.Active(
    sigma_fiber=1, sigma_sheet=0.5, beta=140, cm=0.01, cell_model=CellModel.Tentusscher()
)
epinsulator = EPMaterial.Insulator()

# import pyvista as pv


# ##############################################################################
# .. note::
#    The Ca2+ curve is ignored if the simulation is coupled with electrophysiology.

###############################################################################
# Assign materials to a part
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign the materials to the heart model.


###############################################################################
# Load a heart model
# ~~~~~~~~~~~~~~~~~~

###############################################################################
# .. note::
#    You must complete the full heart preprocessing example first.

import ansys.health.heart.examples as examples
import ansys.health.heart.models as models

heart_model_vtu, heart_model_partinfo, _ = examples.get_preprocessed_fullheart()
workdir = str(Path.home() / "pyansys-heart" / "Rodero2021")

# Load a full-heart model.
heartmodel: models.FullHeart = models.FullHeart(working_directory=workdir)
heartmodel.load_model_from_mesh(heart_model_vtu, heart_model_partinfo)

heartmodel.plot_mesh(show_edges=False)

# Print the default material. You should see that the material is empty.
print(heartmodel.left_ventricle.meca_material)
print(heartmodel.left_ventricle.ep_material)

###############################################################################
# .. note::
#    If no material is set before writing k files, the default material
#    from the ``settings`` object is used.

# Assign the material that you just created.
heartmodel.left_ventricle.meca_material = active_mat
heartmodel.left_ventricle.ep_material = ep_mat_active

# Print it. You should see the following:
# MAT295(rho=1, iso=ISO(itype=-3, beta=0.0, nu=0.499, k1=1, k2=1), aopt=2.0, aniso=ANISO(atype=-1, fibers=[ANISO.HGOFiber(k1=1, k2=1, a=0.0, b=1.0, _theta=0.0, _ftype=1, _fcid=0)], k1fs=None, k2fs=None, vec_a=(1.0, 0.0, 0.0), vec_d=(0.0, 1.0, 0.0), nf=1, intype=0), active=ActiveModel.Model1(t0=None, ca2ion=None, ca2ionm=4.35, n=2, taumax=0.125, stf=0.0, b=4.75, l0=1.58, l=1.85, dtmax=150, mr=1048.9, tr=-1629.0))  # noqa
print(heartmodel.left_ventricle.meca_material)
print(heartmodel.left_ventricle.ep_material)
###############################################################################
# Create a new part and set material

# # A new part can be created by elements IDs
# ids = np.where(heartmodel.mesh.point_data_to_cell_data()["uvc_longitudinal"] > 0.9)[0]
# new_part: Part = heartmodel.create_part_by_ids(ids, "new_part")

# # Show the part
# plotter = heartmodel.plot_part(new_part)
