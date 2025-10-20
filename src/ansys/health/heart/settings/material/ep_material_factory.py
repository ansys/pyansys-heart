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

"""Factory functions for creating electrophysiology material models.

This module provides factory functions for creating EP material instances with
appropriate default parameters for different cardiac tissue types and solver
configurations. The factory functions handle solver-specific parameter selection,
unit conversion from Quantity objects, and proper material instantiation using
the existing Pydantic-based material models.

The module supports three main EP material types:
- Active materials for excitable myocardium with cell models
- Passive materials for non-excitable cardiac tissue
- ActiveBeam materials for specialized conduction system elements

All factory functions are optimized for PyAnsys Heart simulation workflows and
follow cardiac electrophysiology modeling best practices.

Examples
--------
>>> from ansys.health.heart.settings.material.ep_material import EPSolverType
>>> from ansys.health.heart.settings.material.ep_material_factory import (
...     get_default_myocardium_material,
...     get_default_passive_material,
...     assign_default_ep_materials,
... )
>>> import ansys.health.heart.models as models

>>> # Create materials for different solver types
>>> active_mono = get_default_myocardium_material(EPSolverType.MONODOMAIN)
>>> active_eik = get_default_myocardium_material("Eikonal")
>>> passive = get_default_passive_material(EPSolverType.MONODOMAIN)

>>> # Assign materials to a complete model
>>> model = models.BiVentricle()
>>> assign_default_ep_materials(model, EPSolverType.MONODOMAIN)
"""

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
    """Create default active EP material for myocardium tissue.

    This function creates an Active EP material instance with solver-specific
    default parameters appropriate for myocardium tissue simulation. The material
    includes cell model properties and conductivity values optimized for the
    specified EP solver type.

    Parameters
    ----------
    ep_solver_type : EPSolverType or Literal["Monodomain", "Eikonal", "ReactionEikonal"]
        The electrophysiology solver type to configure material for.
        Determines conductivity units and parameter values.

    Returns
    -------
    Active
        Configured active EP material with solver-specific defaults and
        associated Tentusscher cell model.

    Raises
    ------
    ValueError
        If ep_solver_type is not one of the supported solver types.
    RuntimeError
        If material creation or parameter processing fails.

    Examples
    --------
    >>> from ansys.health.heart.settings.material.ep_material import EPSolverType
    >>> material = get_default_myocardium_material(EPSolverType.MONODOMAIN)
    >>> print(material.sigma_fiber)
    0.5
    >>> print(type(material.cell_model))
    <class 'ansys.health.heart.settings.material.cell_models.Tentusscher'>

    >>> # Using string literal
    >>> material_eik = get_default_myocardium_material("Eikonal")
    >>> print(material_eik.sigma_fiber)
    0.7
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
    """Create default passive EP material for non-excitable tissue.

    This function creates a Passive EP material instance with solver-specific
    default parameters appropriate for passive cardiac tissue. Passive materials
    do not have cell models or sheet conductivities, representing tissue that
    does not actively generate electrical activity.

    Parameters
    ----------
    ep_solver_type : EPSolverType or Literal["Monodomain", "Eikonal", "ReactionEikonal"]
        The electrophysiology solver type to configure material for.
        Determines conductivity units and parameter values.

    Returns
    -------
    Passive
        Configured passive EP material with solver-specific defaults.
        Sheet conductivities are automatically excluded for passive materials.

    Raises
    ------
    ValueError
        If ep_solver_type is not one of the supported solver types.
    RuntimeError
        If material creation or parameter processing fails.

    Examples
    --------
    >>> from ansys.health.heart.settings.material.ep_material import EPSolverType
    >>> material = get_default_passive_material(EPSolverType.EIKONAL)
    >>> print(material.sigma_fiber)
    0.7
    >>> print(hasattr(material, "cell_model"))
    False

    >>> # Using string literal
    >>> material_mono = get_default_passive_material("Monodomain")
    >>> print(material_mono.sigma_fiber)
    0.5
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
    """Create default active beam EP material for conduction system elements.

    This function creates an ActiveBeam EP material instance optimized for
    specialized conduction system components such as Purkinje fibers, bundle
    branches, and AV node pathways. These materials are isotropic (fiber-only
    conduction) and have enhanced conductivity properties.

    Parameters
    ----------
    ep_solver_type : EPSolverType or Literal["Monodomain", "Eikonal", "ReactionEikonal"], \
default: EPSolverType.MONODOMAIN
        The electrophysiology solver type to configure material for.
        Determines conductivity units and parameter values for beam elements.

    Returns
    -------
    ActiveBeam
        Configured active beam EP material with enhanced conduction properties.
        Sheet conductivities are automatically set to None for isotropic behavior.

    Raises
    ------
    ValueError
        If ep_solver_type is not one of the supported solver types.
    RuntimeError
        If material creation or parameter processing fails.

    Examples
    --------
    >>> from ansys.health.heart.settings.material.ep_material import EPSolverType
    >>> material = get_default_conduction_system_material(EPSolverType.MONODOMAIN)
    >>> print(material.sigma_fiber)
    1.0
    >>> print(material.sigma_sheet is None)
    True

    >>> # Default solver type
    >>> material_default = get_default_conduction_system_material()
    >>> print(material_default.sigma_fiber)
    1.0

    >>> # Using string literal for Eikonal
    >>> material_eik = get_default_conduction_system_material("Eikonal")
    >>> print(material_eik.sigma_fiber)
    1.0
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
    """Assign default EP materials to all model components without existing materials.

    This function performs in-place assignment of appropriate EP materials to all
    parts and conduction paths in a heart model that do not already have valid
    EP materials assigned. The material selection follows cardiac physiology
    principles and is optimized for the specified solver type.

    Parameters
    ----------
    model : HeartModel or FourChamber or FullHeart or BiVentricle or LeftVentricle
        Heart model instance to assign materials to. The model is modified in-place
        with new material assignments.
    solver_type : EPSolverType or Literal["Monodomain", "Eikonal", "ReactionEikonal"], \
default: EPSolverType.MONODOMAIN
        Electrophysiology solver type for material parameter selection.
        Determines conductivity units and optimization.

    Raises
    ------
    ValueError
        If solver_type is not one of the supported EP solver types.
    RuntimeError
        If material assignment fails for any model component.

    Examples
    --------
    >>> import ansys.health.heart.models as models
    >>> from ansys.health.heart.settings.material.ep_material import EPSolverType
    >>>
    >>> # Create model and assign materials
    >>> model = models.BiVentricle()
    >>> assign_default_ep_materials(model, EPSolverType.MONODOMAIN)
    >>>
    >>> # Verify all parts have materials
    >>> all_parts_assigned = all(part.ep_material is not None for part in model.parts)
    >>> print(all_parts_assigned)
    True
    >>>
    >>> # Check material types
    >>> active_parts = [p for p in model.parts if p.active]
    >>> passive_parts = [
    ...     p for p in model.parts
    ...     if not p.active and not isinstance(p, models.Artery)
    ... ]
    >>> arteries = [p for p in model.parts if isinstance(p, models.Artery)]
    >>>
    >>> print(f"Active parts: {len(active_parts)}")
    >>> print(f"Passive parts: {len(passive_parts)}")
    >>> print(f"Insulator parts: {len(arteries)}")

    Notes
    -----
    **Material Assignment Rules:**

    * **Active parts** (contractile myocardium): Receive Active materials with
      Tentusscher cell models and solver-specific conductivities
    * **Passive parts** (non-contractile tissue): Receive Passive materials
      without cell models or sheet conductivities
    * **Conduction paths** (Purkinje system, bundles): Receive ActiveBeam
      materials with enhanced fiber conductivity
    * **Arteries and veins**: Receive Insulator materials with zero conductivity

    **Solver-Specific Optimizations:**

    * **Monodomain**: Uses conductivity values in mS/mm units
    * **Eikonal/ReactionEikonal**: Uses velocity values in mm/ms units

    The function logs all material assignments and provides a summary of the
    assignment process. Parts with existing valid EP materials are skipped.
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

    assignments_parts = 0

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

                assignments_parts += 1

            except Exception as e:
                error_msg = f"Failed to assign EP material to part '{part.name}': {e}"
                LOGGER.error(error_msg)
                raise RuntimeError(error_msg) from e

    assignments_conduction = 0
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
                assignments_conduction += 1

            except Exception as e:
                error_msg = (
                    f"Failed to assign EP material to conduction path '{conduction_path.name}': {e}"
                )
                LOGGER.error(error_msg)
                raise RuntimeError(error_msg) from e

    if assignments_parts > 0:
        LOGGER.info(
            f"Successfully assigned default EP materials to "
            f"{assignments_parts}/{len(model.parts)} parts "
            f"using {solver_type.value} solver."
        )
    if assignments_conduction > 0:
        LOGGER.info(
            f"Successfully assigned default EP materials to "
            f"{assignments_conduction}/{len(model.conduction_paths)} conduction paths "
            f"using {solver_type.value} solver."
        )
