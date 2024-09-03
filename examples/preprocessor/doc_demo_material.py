# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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
import os
from pathlib import Path

import matplotlib.pyplot as plt

from ansys.heart.simulator.settings.material.curve import (
    ActiveCurve,
    Kumaraswamy_active,
    constant_ca2,
)
from ansys.heart.simulator.settings.material.ep_material import CellModel, EPMaterial
from ansys.heart.simulator.settings.material.material import (
    ACTIVE,
    ANISO,
    ISO,
    MAT295,
    ActiveModel,
    NeoHookean,
)

# sphinx_gallery_start_ignore
docs_images_folder = Path(Path(__file__).resolve().parents[2], "doc", "source", "_static", "images")
# sphinx_gallery_end_ignore

###############################################################################
# .. note::
#    Unit system used for heart modeling in LS-DYNA is ["MPa", "mm", "N", "ms", "g"]


###############################################################################
# Create a material
# ~~~~~~~~~~~~~~~~~

## Neo-Hookean material can be created as following
neo = NeoHookean(rho=0.001, c10=1, nu=0.499)

###############################################################################
# .. note::
#    Please refer to LS-DYNA manual for more details of MAT_295

## More steps to create MAT295 which is used for myocardium

# step 1: create an isotropic module
iso = ISO(k1=1, k2=1, nu=0.499)

# step 2: create an anisotropoc moddule
fiber = ANISO.HGO_Fiber(k1=1, k2=1)
aniso1 = ANISO(fibers=[fiber])

# Create fiber with sheet, and their interactions
sheet = ANISO.HGO_Fiber(k1=1, k2=1)
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
###############################################################################
# .. image:: /_static/images/model1_ca2.png
#   :width: 400pt
#   :align: center

# sphinx_gallery_start_ignore
fig.savefig(os.path.join(docs_images_folder, "model1_ca2.png"))
# sphinx_gallery_end_ignore

# example 2
# create active model 3
ac_model3 = ActiveModel.Model3()
# create a stress curve and show
ac_curve3 = ActiveCurve(Kumaraswamy_active(t_end=800), type="stress")
fig = ac_curve3.plot_time_vs_stress()
plt.show()

###############################################################################
# .. image:: /_static/images/model3_stress.png
#   :width: 400pt
#   :align: center

# sphinx_gallery_start_ignore
fig.savefig(os.path.join(docs_images_folder, "model3_stress.png"))
# sphinx_gallery_end_ignore

###############################################################################
# .. note::
#   With setting eta=0 is model 3, stress curve will be the active stress for all elements.
#   If eta!=0, this is idealized active stress when fiber stretch stays to 1.

# PyAnsys-Heart will convert the stress curve to Ca2+ curve (input of MAT_295)
fig = ac_curve3.plot_time_vs_ca2()
plt.show()

###############################################################################
# .. image:: /_static/images/model3_ca2+.png
#   :width: 400pt
#   :align: center

# sphinx_gallery_start_ignore
fig.savefig(os.path.join(docs_images_folder, "model3_ca2+.png"))
# sphinx_gallery_end_ignore

# build active module
active3 = ACTIVE(model=ac_model3, ca2_curve=ac_curve3)

###############################################################################
## Finally, MAT295 can be created with the above modules
iso_mat = MAT295(rho=1, iso=iso, aniso=None, active=None)
passive_mat = MAT295(rho=1, iso=iso, aniso=aniso1, active=None)
active_mat = MAT295(rho=1, iso=iso, aniso=aniso1, active=active)

###############################################################################
## EP materials can be created as follows
ep_mat_active = EPMaterial.Active(
    sigma_fiber=1, sigma_sheet=0.5, beta=140, cm=0.01, cell_model=CellModel.Tentusscher()
)
epinsulator = EPMaterial.Insulator()

###############################################################################
# .. note::
#    Ca2+ curve will be ignored if the simulation is coupled with electrophysiology

###############################################################################
# Assign material to a part
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign the materials to the heart model

import numpy as np

from ansys.heart.preprocessor.mesh.objects import Part
import ansys.heart.preprocessor.models as models

###############################################################################
# Load a heart model

###############################################################################
# .. note::
#    run doc_preprocess_fullheart_rodero_01.py in the same folder if loading failed

workdir = Path(Path(__file__).resolve().parents[2], "downloads", "Rodero2021", "01", "FullHeart")
path_to_model = os.path.join(workdir, "heart_model.vtu")

# load four chamber heart model.
heartmodel: models.FullHeart = models.FullHeart(models.ModelInfo(work_directory=workdir))
heartmodel.load_model_from_mesh(path_to_model, path_to_model.replace(".vtu", ".partinfo.json"))

# Print default material and you should see
# Material is empty.
print(heartmodel.left_ventricle.meca_material)
print(heartmodel.left_ventricle.ep_material)

###############################################################################
# .. note::
#    If no material is set before writing k files, default material from ```settings``` will be set.

# Assign the material we just created
heartmodel.left_ventricle.meca_material = active_mat
heartmodel.left_ventricle.ep_material = ep_mat_active

# Print it, you should see
# MAT295(rho=1, iso=ISO(itype=-3, beta=0.0, nu=0.499, k1=1, k2=1), aopt=2.0, aniso=ANISO(atype=-1, fibers=[ANISO.HGO_Fiber(k1=1, k2=1, a=0.0, b=1.0, _theta=0.0, _ftype=1, _fcid=0)], k1fs=None, k2fs=None, vec_a=(1.0, 0.0, 0.0), vec_d=(0.0, 1.0, 0.0), nf=1, intype=0), active=ActiveModel.Model1(t0=None, ca2ion=None, ca2ionm=4.35, n=2, taumax=0.125, stf=0.0, b=4.75, l0=1.58, l=1.85, dtmax=150, mr=1048.9, tr=-1629.0))  # noqa
print(heartmodel.left_ventricle.meca_material)

print(heartmodel.left_ventricle.ep_material)
###############################################################################
# Create a new part and set material

# A new part can be created by elements IDs
ids = np.where(heartmodel.mesh.point_data_to_cell_data()["uvc_longitudinal"] > 0.9)[0]
new_part: Part = heartmodel.create_part_by_ids(ids, "new_part")

# Show the part
plotter = heartmodel.plot_part(new_part)

###############################################################################
# .. image:: /_static/images/show_a_part.png
#   :width: 400pt
#   :align: center

# sphinx_gallery_start_ignore
import pyvista as pv

plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(heartmodel.mesh, opacity=0.5, color="white")
part = heartmodel.mesh.extract_cells(new_part.element_ids)
plotter.add_mesh(part, opacity=0.95, color="red")
plotter.screenshot(os.path.join(docs_images_folder, "show_a_part.png"))
# sphinx_gallery_end_ignore

## set passive anisotropic material for it
new_part.fiber = True
new_part.active = False
new_part.meca_material = passive_mat
## and set it to an EP insulator
new_part.ep_material = epinsulator

print(new_part.meca_material)
print(new_part.ep_material)
