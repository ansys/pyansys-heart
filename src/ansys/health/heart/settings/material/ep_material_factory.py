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

"""Factory for creating EP material models."""

from typing import Literal

from pint import Quantity

from ansys.health.heart.settings.material.ep_material import (
    Active,
    ActiveBeam,
    EPSolverType,
    Passive,
)


def get_default_myocardium_material(
    ep_solver_type: EPSolverType | Literal["Monodomain", "Eikonal", "Reaction-Eikonal"],
) -> Active:
    """Get the default myocardium material for a solver type.

    Parameters
    ----------
    ep_solver_type : EPSolverType | str
        The type of EP solver to select appropriate defaults for.

    Returns
    -------
    Active
        The default myocardium EP material model.
    """
    if isinstance(ep_solver_type, str):
        ep_solver_type = EPSolverType(ep_solver_type)

    # import defaults depending on solver type.
    if ep_solver_type in (EPSolverType.REACTION_EIKONAL, EPSolverType.EIKONAL):
        from ansys.health.heart.settings.defaults.electrophysiology import (
            default_myocardium_material_eikonal as defaults,
        )
    if ep_solver_type == EPSolverType.MONODOMAIN:
        from ansys.health.heart.settings.defaults.electrophysiology import (
            default_myocardium_material_monodomain as defaults,
        )

    # Remove units from default Quantity values.
    defaults = {k: (v.m if isinstance(v, Quantity) else v) for k, v in defaults.items()}

    return Active(solver_type=ep_solver_type, **defaults)


def get_passive_material(
    ep_solver_type: EPSolverType | Literal["Monodomain", "Eikonal", "Reaction-Eikonal"],
) -> Passive:
    """Get the default passive material for a solver type.

    Parameters
    ----------
    ep_solver_type : EPSolverType | str
        The type of EP solver to select appropriate defaults for.

    Returns
    -------
    Passive
        The default passive EP material model.
    """
    if isinstance(ep_solver_type, str):
        ep_solver_type = EPSolverType(ep_solver_type)

    # import defaults depending on solver type.
    if ep_solver_type in (EPSolverType.REACTION_EIKONAL, EPSolverType.EIKONAL):
        from ansys.health.heart.settings.defaults.electrophysiology import (
            default_myocardium_material_eikonal as defaults,
        )
    if ep_solver_type == EPSolverType.MONODOMAIN:
        from ansys.health.heart.settings.defaults.electrophysiology import (
            default_myocardium_material_monodomain as defaults,
        )

    # Remove units from default Quantity values.
    defaults = {k: (v.m if isinstance(v, Quantity) else v) for k, v in defaults.items()}

    del defaults["sigma_sheet"]
    del defaults["sigma_sheet_normal"]

    return Passive(**defaults)


def get_default_conduction_system_material(
    ep_solver_type: EPSolverType | Literal["Monodomain", "Eikonal", "Reaction-Eikonal"],
) -> ActiveBeam:
    """Get the default conduction-system (beam) material for a solver type.

    Parameters
    ----------
    ep_solver_type : EPSolverType | str
        The type of EP solver to select appropriate defaults for.

    Returns
    -------
    ActiveBeam
        The default conduction system material.
    """
    if isinstance(ep_solver_type, str):
        ep_solver_type = EPSolverType(ep_solver_type)

    # import defaults depending on solver type.
    if ep_solver_type in (EPSolverType.REACTION_EIKONAL, EPSolverType.EIKONAL):
        from ansys.health.heart.settings.defaults.electrophysiology import (
            default_beam_material_eikonal as defaults,
        )
    elif ep_solver_type == EPSolverType.MONODOMAIN:
        from ansys.health.heart.settings.defaults.electrophysiology import (
            default_beam_material_monodomain as defaults,
        )

    # Remove units from default Quantity values.
    defaults = {k: (v.m if isinstance(v, Quantity) else v) for k, v in defaults.items()}

    return ActiveBeam(solver_type=ep_solver_type, **defaults)
