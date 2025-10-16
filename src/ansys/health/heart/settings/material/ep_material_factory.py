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

    Raises
    ------
    ValueError
        If ep_solver_type is not recognized.
    RuntimeError
        If material creation fails.

    Examples
    --------
    >>> material = get_default_myocardium_material(EPSolverType.MONODOMAIN)
    >>> print(type(material))
    <class 'ansys.health.heart.settings.material.ep_material.Active'>
    """
    try:
        if isinstance(ep_solver_type, str):
            ep_solver_type = EPSolverType(ep_solver_type)

        # Import defaults depending on solver type
        if ep_solver_type in (EPSolverType.REACTION_EIKONAL, EPSolverType.EIKONAL):
            from ansys.health.heart.settings.defaults.electrophysiology import (
                default_myocardium_material_eikonal as defaults,
            )
        elif ep_solver_type == EPSolverType.MONODOMAIN:
            from ansys.health.heart.settings.defaults.electrophysiology import (
                default_myocardium_material_monodomain as defaults,
            )
        else:
            raise ValueError(f"Unsupported EP solver type: {ep_solver_type}")

        # Remove units from default Quantity values
        processed_defaults = {
            k: (v.m if isinstance(v, Quantity) else v) for k, v in defaults.items()
        }

        return Active(**processed_defaults)

    except Exception as e:
        error_msg = f"Failed to create myocardium material for {ep_solver_type}: {e}"
        LOGGER.error(error_msg)
        raise RuntimeError(error_msg) from e


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

    Raises
    ------
    ValueError
        If ep_solver_type is not recognized.
    RuntimeError
        If material creation fails.

    Examples
    --------
    >>> material = get_default_passive_material("Eikonal")
    >>> print(type(material))
    <class 'ansys.health.heart.settings.material.ep_material.Passive'>
    """
    try:
        if isinstance(ep_solver_type, str):
            ep_solver_type = EPSolverType(ep_solver_type)

        # Import defaults depending on solver type
        if ep_solver_type in (EPSolverType.REACTION_EIKONAL, EPSolverType.EIKONAL):
            from ansys.health.heart.settings.defaults.electrophysiology import (
                default_myocardium_material_eikonal as defaults,
            )
        elif ep_solver_type == EPSolverType.MONODOMAIN:
            from ansys.health.heart.settings.defaults.electrophysiology import (
                default_myocardium_material_monodomain as defaults,
            )
        else:
            raise ValueError(f"Unsupported EP solver type: {ep_solver_type}")

        # Remove units from default Quantity values
        processed_defaults = {
            k: (v.m if isinstance(v, Quantity) else v) for k, v in defaults.items()
        }

        # Remove sheet conductivities for passive materials using safe pop method
        processed_defaults.pop("sigma_sheet", None)
        processed_defaults.pop("sigma_sheet_normal", None)

        return Passive(**processed_defaults)

    except Exception as e:
        error_msg = f"Failed to create passive material for {ep_solver_type}: {e}"
        LOGGER.error(error_msg)
        raise RuntimeError(error_msg) from e


def get_default_conduction_system_material(
    ep_solver_type: (
        EPSolverType | Literal["Monodomain", "Eikonal", "ReactionEikonal"]
    ) = EPSolverType.MONODOMAIN,
) -> ActiveBeam:
    """Get the default conduction-system (beam) material for a solver type.

    Parameters
    ----------
    ep_solver_type : EPSolverType | str, default: EPSolverType.MONODOMAIN
        The type of EP solver to select appropriate defaults for.

    Returns
    -------
    ActiveBeam
        The default conduction system material.

    Raises
    ------
    ValueError
        If ep_solver_type is not recognized.
    RuntimeError
        If material creation fails.

    Examples
    --------
    >>> material = get_default_conduction_system_material(EPSolverType.MONODOMAIN)
    >>> print(type(material))
    <class 'ansys.health.heart.settings.material.ep_material.ActiveBeam'>
    """
    try:
        if isinstance(ep_solver_type, str):
            ep_solver_type = EPSolverType(ep_solver_type)

        # Import defaults depending on solver type
        if ep_solver_type in (EPSolverType.REACTION_EIKONAL, EPSolverType.EIKONAL):
            from ansys.health.heart.settings.defaults.electrophysiology import (
                default_beam_material_eikonal as defaults,
            )
        elif ep_solver_type == EPSolverType.MONODOMAIN:
            from ansys.health.heart.settings.defaults.electrophysiology import (
                default_beam_material_monodomain as defaults,
            )
        else:
            raise ValueError(f"Unsupported EP solver type: {ep_solver_type}")

        # Remove units from default Quantity values
        processed_defaults = {
            k: (v.m if isinstance(v, Quantity) else v) for k, v in defaults.items()
        }

        return ActiveBeam(**processed_defaults)

    except Exception as e:
        error_msg = f"Failed to create conduction system material for {ep_solver_type}: {e}"
        LOGGER.error(error_msg)
        raise RuntimeError(error_msg) from e


def assign_default_ep_materials(
    model: (
        models.HeartModel
        | models.FourChamber
        | models.FullHeart
        | models.BiVentricle
        | models.LeftVentricle
    ),
    solver_type: (
        EPSolverType | Literal["Monodomain", "Eikonal", "ReactionEikonal"]
    ) = EPSolverType.MONODOMAIN,
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
    solver_type : EPSolverType | str, default: EPSolverType.MONODOMAIN
        EP solver type for material selection. Must be one of:
        "Monodomain", "Eikonal", or "ReactionEikonal".

    Raises
    ------
    ValueError
        If solver_type is not recognized by the material factory.
    RuntimeError
        If material assignment fails for any component.

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

    The function will skip parts that already have valid EP materials assigned.
    """
    try:
        solver_type = EPSolverType(solver_type)
    except ValueError as e:
        valid_types = [solver.value for solver in EPSolverType]
        error_msg = (
            f"Unknown EP solver type: '{solver_type}'. Valid options: {', '.join(valid_types)}"
        )
        LOGGER.error(error_msg)
        raise ValueError(error_msg) from e

    assignments = 0

    # Assign materials to parts (modifies objects in-place)
    for part in model.parts:
        if not isinstance(part.ep_material, EPMaterialModel) or part.ep_material is None:
            try:
                if part.active:
                    part.ep_material = get_default_myocardium_material(solver_type)
                    LOGGER.info(
                        f"Assigned active EP material to part '{part.name}' "
                        f"for {solver_type.value} solver."
                    )
                elif isinstance(part, Artery):
                    part.ep_material = Insulator()
                    LOGGER.info(f"Assigned insulator EP material to artery part '{part.name}'.")
                else:
                    part.ep_material = get_default_passive_material(solver_type)
                    LOGGER.info(
                        f"Assigned passive EP material to part '{part.name}' "
                        f"for {solver_type.value} solver."
                    )

                assignments += 1

            except Exception as e:
                error_msg = f"Failed to assign EP material to part '{part.name}': {e}"
                LOGGER.error(error_msg)
                raise RuntimeError(error_msg) from e

    # Assign materials to conduction paths (modifies objects in-place)
    for conduction_path in model.conduction_paths:
        if conduction_path.ep_material is None:
            try:
                conduction_path.ep_material = get_default_conduction_system_material(solver_type)
                LOGGER.warning(
                    f"Conduction path '{conduction_path.name}' did not have an "
                    f"EP material assigned. Assigned default conduction system material "
                    f"for {solver_type.value} solver."
                )
                assignments += 1

            except Exception as e:
                error_msg = (
                    f"Failed to assign EP material to conduction path '{conduction_path.name}': {e}"
                )
                LOGGER.error(error_msg)
                raise RuntimeError(error_msg) from e

    LOGGER.info(
        f"Successfully assigned default EP materials to "
        f"{assignments}/{len(model.parts) + len(model.conduction_paths)} components "
        f"using {solver_type.value} solver."
    )
