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

from __future__ import annotations

from typing import Literal

from pint import Quantity

from ansys.health.heart import LOG as LOGGER
import ansys.health.heart.models as models
from ansys.health.heart.parts import Artery
from ansys.health.heart.settings.material.ep_material import (
    Active,
    ActiveBeam,
    EPMaterialModel,
    EPSolverType,
    Insulator,
    Passive,
)


def get_default_myocardium_material(
    ep_solver_type: EPSolverType | Literal["Monodomain", "Eikonal", "ReactionEikonal"],
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

    return Active(**defaults)


def get_default_passive_material(
    ep_solver_type: EPSolverType | Literal["Monodomain", "Eikonal", "ReactionEikonal"],
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
    ep_solver_type: EPSolverType | Literal["Monodomain", "Eikonal", "ReactionEikonal"],
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

    return ActiveBeam(**defaults)


def assign_default_ep_materials(
    model: models.HeartModel
    | models.FourChamber
    | models.FullHeart
    | models.BiVentricle
    | models.LeftVentricle,
    solvertype: Literal["Monodomain", "Eikonal", "ReactionEikonal"] | EPSolverType,
) -> None:
    """Assign default EP materials to all parts of a model that do not have a material assigned.

    This function modifies the input model in-place by assigning appropriate
    EP materials to parts and conduction paths based on their characteristics
    and the specified solver type.

    Parameters
    ----------
    model : models.HeartModel | models.FourChamber | models.FullHeart | models.BiVentricle |
            models.LeftVentricle
        Heart model to assign materials to. Modified in-place.
    solvertype : str | ep_materials.EPSolverType
        EP solver type for material selection. Must be one of:
        "Monodomain", "Eikonal", or "ReactionEikonal".

    Raises
    ------
    ValueError
        If solvertype is not recognized by the material factory.

    Examples
    --------
    >>> import ansys.health.heart.models as models
    >>> from ansys.health.heart.settings.material.ep_material_factory import (
    ...     assign_default_ep_materials,
    ... )
    >>> from ansys.health.heart.settings.material.ep_material import EPSolverType
    >>> model = models.BiVentricle()
    >>> assign_default_ep_materials(model, EPSolverType.MONODOMAIN)
    >>> all_parts_have_materials = all(part.ep_material is not None for part in model.parts)
    >>> print(all_parts_have_materials)
    True

    Notes
    -----
    This function follows the PyAnsys Heart material assignment hierarchy:

    - Active parts receive myocardium materials based on solver type
    - Passive parts receive passive materials (no cell model associated with them)
    - Conduction paths receive specialized conduction system materials
    - Arteries and veins are treated as insulators
    """
    try:
        solvertype = EPSolverType(solvertype)
    except ValueError as e:
        raise ValueError(
            f"Unrecognized solver type: {solvertype}. Use one of {list(EPSolverType)}"
        ) from e

    for part in model.parts:
        if not isinstance(part.ep_material, EPMaterialModel) or part.ep_material is None:
            if part.active:
                part.ep_material = get_default_myocardium_material(solvertype)
                LOGGER.info(
                    f"Part {part.name} does not have an EP material assigned. "
                    "Assigning default active EP material."
                )
            elif isinstance(part, Artery):
                part.ep_material = Insulator()
                LOGGER.info(
                    f"Part {part.name} does not have an EP material assigned. "
                    "Assigning an insulator material."
                )
            else:
                LOGGER.info(
                    f"Part {part.name} does not have an EP material assigned. "
                    "Assigning default passive EP material."
                )
                part.ep_material = get_default_passive_material(solvertype)

    for conduction_path in model.conduction_paths:
        if conduction_path.ep_material is None:
            conduction_path.ep_material = get_default_conduction_system_material(solvertype)
            LOGGER.warning(
                f"Conduction path {conduction_path.name} does not have an EP material assigned. "
                "Assigning default EP material."
            )

    return
