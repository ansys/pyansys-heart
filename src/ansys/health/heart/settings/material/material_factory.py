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
"""Factory functions for creating mechanical material models.

This module provides factory functions for creating mechanical material instances
with appropriate default parameters for cardiac tissue simulation. The factory
functions handle parameter extraction from Quantity objects, material component
assembly (isotropic, anisotropic, active), and proper material instantiation using
the existing Pydantic-based Mat295 material model.

The module supports mechanical material creation for:
- Myocardium tissue with isotropic, anisotropic, and active components
- Passive tissue with isotropic properties only
- Electro-mechanical coupling configurations for multi-physics simulations

All factory functions are optimized for PyAnsys Heart simulation workflows and
follow cardiac mechanics modeling best practices including proper fiber orientation,
hyperelastic constitutive models, and active stress generation.

Examples
--------
>>> import ansys.health.heart.models as models
>>> from ansys.health.heart.settings.material.material_factory import (
...     assign_default_mechanics_materials,
... )

>>> # Create model and assign mechanical materials
>>> model = models.BiVentricle()
>>> assign_default_mechanics_materials(model, ep_coupled=True)

>>> # Verify material assignment
>>> all_assigned = all(part.meca_material is not None for part in model.parts)
>>> print(all_assigned)
True

>>> # Check for electro-mechanical coupling
>>> active_parts = [p for p in model.parts if p.active and p.fiber]
>>> has_active_materials = all(p.meca_material.active is not None for p in active_parts)
>>> print(has_active_materials)
True
"""

from ansys.health.heart import LOG as LOGGER
import ansys.health.heart.models as models
from ansys.health.heart.settings.material.curve import constant_ca2
from ansys.health.heart.settings.material.material import (
    ACTIVE,
    ANISO,
    ISO,
    ActiveCurve,
    ActiveModel1,
    ActiveModel3,
    HGOFiber,
    Mat295,
    MechanicalMaterialModel,
)


def _default_myocardium_material(ep_coupled: bool = False) -> Mat295:
    """Create default Mat295 myocardium material with complete material components.

    This function creates a comprehensive myocardium material model including
    isotropic (HGO), anisotropic (fiber/sheet), and active stress components.
    The material is configured for cardiac tissue simulation with appropriate
    hyperelastic and active properties.

    Parameters
    ----------
    ep_coupled : bool, default: False
        Whether to configure the material for electro-mechanical coupling.
        If True, uses ActiveModel3 with calcium-dependent activation.
        If False, uses ActiveModel1 with predefined activation curve.

    Returns
    -------
    Mat295
        Complete myocardium material with isotropic, anisotropic, and active
        components configured according to cardiac mechanics standards.

    Raises
    ------
    ValueError
        If material settings are incomplete or contain invalid values.
    RuntimeError
        If material creation fails due to component assembly errors.

    Examples
    --------
    >>> # Standard myocardium material for mechanics-only simulation
    >>> material = _default_myocardium_material(ep_coupled=False)
    >>> print(material.active.ca2_curve is not None)
    True
    >>> print(isinstance(material.active.model, ActiveModel1))
    True

    >>> # EP-coupled material for multi-physics simulation
    >>> ep_material = _default_myocardium_material(ep_coupled=True)
    >>> print(ep_material.active.ca2_curve is None)
    True
    >>> print(isinstance(ep_material.active.model, ActiveModel3))
    True

    >>> # Verify complete material structure
    >>> print(material.iso is not None)
    True
    >>> print(material.aniso is not None)
    True
    >>> print(material.active is not None)
    True

    Notes
    -----
    The material uses default mechanical properties from the mechanics defaults
    configuration. For electro-mechanical coupling, the active model switches
    from time-based activation (ActiveModel1) to calcium-dependent activation
    (ActiveModel3) to enable coupling with electrophysiology simulations.
    """
    try:
        from ansys.health.heart.settings.defaults.mechanics import material

        return _get_myocardium_material(material["myocardium"], ep_coupled=ep_coupled)

    except Exception as e:
        error_msg = f"Failed to create default myocardium material (EP coupled: {ep_coupled}): {e}"
        LOGGER.error(error_msg)
        raise RuntimeError(error_msg) from e


def _default_passive_material() -> Mat295:
    """Create default Mat295 passive material with isotropic properties only.

    This function creates a passive mechanical material model suitable for
    non-contractile cardiac tissue such as connective tissue or
    inactive myocardium regions. The material contains only isotropic
    hyperelastic properties without anisotropic or active components.

    Returns
    -------
    Mat295
        Passive material with isotropic Ogden hyperelastic model configured
        for non-contractile tissue.

    Raises
    ------
    ValueError
        If passive material settings are incomplete or contain invalid values.
    RuntimeError
        If material creation fails due to configuration errors.

    Examples
    --------
    >>> # Create passive material for non-contractile tissue
    >>> material = _default_passive_material()
    >>> print(material.active is None)
    True
    >>> print(material.aniso is None)
    True
    >>> print(material.iso is not None)
    True

    >>> # Verify isotropic configuration
    >>> print(material.iso.itype)  # Should be Ogden model type
    1
    >>> print(material.iso.mu1 > 0)  # Should have positive shear modulus
    True
    >>> print(material.iso.kappa > 0)  # Should have positive bulk modulus
    True

    Notes
    -----
    The passive material uses default mechanical properties from the mechanics
    defaults configuration. It employs an Ogden hyperelastic model suitable for
    soft tissue mechanics without contractile behavior. This material type is
    typically assigned to parts without fiber orientation or active properties.
    """
    try:
        from ansys.health.heart.settings.defaults.mechanics import material

        return _get_passive_material(material["passive"])

    except Exception as e:
        error_msg = f"Failed to create default passive material: {e}"
        LOGGER.error(error_msg)
        raise RuntimeError(error_msg) from e


def _get_myocardium_material(settings: dict, ep_coupled: bool = False) -> Mat295:
    """Construct Mat295 myocardium material from configuration settings.

    This function assembles a complete myocardium material model by extracting
    parameters from a settings dictionary and creating the isotropic, anisotropic,
    and active material components. It handles unit conversion from Quantity objects
    and configures the appropriate active model based on coupling requirements.

    Parameters
    ----------
    settings : dict
        Material configuration dictionary containing nested sections:
        - "isotropic": HGO model parameters (rho, kappa, k1, k2)
        - "anisotropic": Fiber/sheet parameters (k1f, k2f, k1s, k2s, k1fs, k2fs)
        - "active": Contractility parameters (taumax, beat_time, ss, sn)
    ep_coupled : bool, default: False
        Whether to configure for electro-mechanical coupling.
        Determines active model selection (ActiveModel1 vs ActiveModel3).

    Returns
    -------
    Mat295
        Assembled myocardium material with all components configured according
        to the provided settings and coupling requirements.

    Raises
    ------
    ValueError
        If required settings sections are missing or material construction fails.
    KeyError
        If required material properties are missing from settings.

    Examples
    --------
    >>> # Example settings structure
    >>> settings = {
    ...     "isotropic": {
    ...         "rho": Quantity(1.0, "g/mm^3"),
    ...         "kappa": Quantity(1000.0, "MPa"),
    ...         "k1": Quantity(0.5, "MPa"),
    ...         "k2": Quantity(0.7, "dimensionless"),
    ...     },
    ...     "anisotropic": {"k1f": Quantity(15.0, "MPa"), "k2f": Quantity(10.0, "dimensionless")},
    ...     "active": {
    ...         "taumax": Quantity(120.0, "kPa"),
    ...         "beat_time": Quantity(800.0, "ms"),
    ...         "ss": 0.0,
    ...         "sn": 0.0,
    ...     },
    ... }
    >>> material = _get_myocardium_material(settings, ep_coupled=True)
    >>> print(material.active is not None)
    True
    >>> print(isinstance(material.active.model, ActiveModel3))
    True

    Notes
    -----
    **Material Component Assembly:**

    * **Isotropic**: HGO hyperelastic model with bulk and shear response
    * **Anisotropic**: Fiber and optional sheet directions with HGO fibers
    * **Active**: Time-based (ActiveModel1) or calcium-based (ActiveModel3) stress

    **Electro-Mechanical Coupling:**

    * **EP Coupled (True)**: Uses ActiveModel3 with calcium threshold activation
    * **Mechanics Only (False)**: Uses ActiveModel1 with predefined time curve

    The function extracts magnitude values from Quantity objects using the `.m`
    attribute and validates that all required material sections are present.
    """
    # Validate required settings
    required_sections = ["isotropic", "anisotropic", "active"]
    missing_sections = [section for section in required_sections if section not in settings]
    if missing_sections:
        raise ValueError(
            f"Incomplete myocardium material settings. "
            f"Missing sections: {', '.join(missing_sections)}"
        )

    try:
        rho = settings["isotropic"]["rho"].m

        iso = ISO(
            kappa=settings["isotropic"]["kappa"].m,
            k1=settings["isotropic"]["k1"].m,
            k2=settings["isotropic"]["k2"].m,
            beta=2,
        )

        fibers = [
            HGOFiber(k1=settings["anisotropic"]["k1f"].m, k2=settings["anisotropic"]["k2f"].m)
        ]

        if "k1s" in settings["anisotropic"]:
            sheet = HGOFiber(
                k1=settings["anisotropic"]["k1s"].m, k2=settings["anisotropic"]["k2s"].m
            )
            fibers.append(sheet)

        if "k1fs" in settings["anisotropic"]:
            k1fs, k2fs = settings["anisotropic"]["k1fs"].m, settings["anisotropic"]["k2fs"].m
        else:
            k1fs, k2fs = None, None
        aniso = ANISO(fibers=fibers, k1fs=k1fs, k2fs=k2fs)

        max = settings["active"]["taumax"].m
        bt = settings["active"]["beat_time"].m
        ss = settings["active"]["ss"]
        sn = settings["active"]["sn"]

        if not ep_coupled:
            ac_mdoel = ActiveModel1(taumax=max)  # use default field in Model1 except taumax
            curve = ActiveCurve(func=constant_ca2(tb=bt), threshold=0.1, type="ca2")
            active = ACTIVE(
                ss=ss,
                sn=sn,
                model=ac_mdoel,
                ca2_curve=curve,
            )
        else:
            ac_mdoel = ActiveModel3(
                ca2ion50=0.001,
                n=2,
                f=0.0,
                l=1.9,
                eta=1.45,
                sigmax=max,  # MPa
            )

            active = ACTIVE(
                ss=ss,
                sn=sn,
                acthr=0.0002,
                model=ac_mdoel,
                ca2_curve=None,
            )

        return Mat295(rho=rho, iso=iso, aniso=aniso, active=active)

    except KeyError as e:
        raise ValueError(f"Missing required material property: {e}") from e
    except Exception as e:
        raise ValueError(f"Failed to construct myocardium material: {e}") from e


def _get_passive_material(passive_settings: dict) -> Mat295:
    """Construct passive Mat295 material from configuration settings.

    This function creates a passive mechanical material model by extracting
    parameters from a settings dictionary. The resulting material contains
    only isotropic hyperelastic properties suitable for non-contractile tissue.

    Parameters
    ----------
    passive_settings : dict
        Passive material configuration dictionary containing:
        - "rho": Material density as Quantity object
        - "itype": Hyperelastic model type (Ogden vs HGO)
        - "kappa": Bulk modulus as Quantity object
        - "mu1": Shear modulus as Quantity object (for Ogden model)
        - "alpha1": Ogden exponent parameter

    Returns
    -------
    Mat295
        Passive material with isotropic Ogden hyperelastic model configured
        according to the provided settings.

    Raises
    ------
    KeyError
        If required passive material properties are missing from settings.
    ValueError
        If material construction fails due to invalid parameter values.

    Examples
    --------
    >>> # Example passive settings structure
    >>> passive_settings = {
    ...     "rho": Quantity(1.0, "g/mm^3"),
    ...     "itype": 1,  # Ogden model
    ...     "kappa": Quantity(1000.0, "MPa"),
    ...     "mu1": Quantity(10.0, "MPa"),
    ...     "alpha1": 18.5,
    ... }
    >>> material = _get_passive_material(passive_settings)
    >>> print(material.active is None)
    True
    >>> print(material.aniso is None)
    True
    >>> print(material.iso.itype)
    1

    Notes
    -----
    **Passive Material Configuration:**

    * **Density**: Extracted from Quantity object using `.m` attribute
    * **Hyperelastic Model**: Typically Ogden model (itype=1) for soft tissue
    * **Bulk Response**: Controlled by kappa parameter for near-incompressibility
    * **Shear Response**: Controlled by mu1 and alpha1 for nonlinear elasticity

    The function uses the existing ISO Pydantic model for validation and
    parameter checking. No anisotropic or active components are included.
    """
    try:
        passive = Mat295(
            rho=passive_settings["rho"].m,
            iso=ISO(
                itype=passive_settings["itype"],
                beta=2,
                kappa=passive_settings["kappa"].m,
                mu1=passive_settings["mu1"].m,
                alpha1=passive_settings["alpha1"],
            ),
        )
        return passive

    except KeyError as e:
        raise KeyError(f"Missing required passive material property: {e}") from e
    except Exception as e:
        raise ValueError(f"Failed to construct passive material: {e}") from e


def assign_default_mechanics_materials(
    model: models.HeartModel
    | models.FullHeart
    | models.BiVentricle
    | models.LeftVentricle
    | models.FourChamber,
    ep_coupled: bool = False,
) -> None:
    """Assign default mechanical materials to all model parts without existing materials.

    This function performs in-place assignment of appropriate mechanical materials
    to all parts in a heart model that do not already have valid mechanical materials.
    The material selection is based on part characteristics (fiber orientation, active
    properties).

    Parameters
    ----------
    model : HeartModel or FullHeart or BiVentricle or LeftVentricle or FourChamber
        Heart model instance to assign materials to. The model is modified in-place
        with new mechanical material assignments.
    ep_coupled : bool, default: False
        Whether to configure materials for electro-mechanical coupling.
        Affects active model selection and activation mechanisms.

    Raises
    ------
    ValueError
        If material creation fails for any part due to configuration errors.
    RuntimeError
        If material assignment fails due to part validation or access errors.

    Examples
    --------
    >>> import ansys.health.heart.models as models
    >>> from ansys.health.heart.settings.material.material_factory import (
    ...     assign_default_mechanics_materials,
    ... )
    >>> # Create model and assign mechanical materials
    >>> model = models.BiVentricle()
    >>> assign_default_mechanics_materials(model, ep_coupled=True)
    >>> # Verify all parts have materials assigned
    >>> all_assigned = all(part.meca_material is not None for part in model.parts)
    >>> print(all_assigned)
    True
    >>> # Check material types based on part characteristics
    >>> fiber_parts = [p for p in model.parts if p.fiber]
    >>> non_fiber_parts = [p for p in model.parts if not p.fiber]
    >>> # Verify fiber parts have complete myocardium materials
    >>> fiber_materials_complete = all(
    ...     p.meca_material.iso is not None and p.meca_material.aniso is not None
    ...     for p in fiber_parts
    ... )
    >>> print(fiber_materials_complete)
    True
    >>> # Verify non-fiber parts have passive materials only
    >>> non_fiber_passive = all(
    ...     p.meca_material.aniso is None and p.meca_material.active is None
    ...     for p in non_fiber_parts
    ... )
    >>> print(non_fiber_passive)
    True

    Notes
    -----
    **Material Assignment Rules:**

    * **Parts with fiber orientation**: Receive complete myocardium materials with
      isotropic, anisotropic, and potentially active components
    * **Parts without fiber orientation**: Receive passive materials with isotropic
      properties only
    * **Active vs passive parts**: Active parts retain active material components,
      passive parts have active components set to None
    * **Existing materials**: Parts with valid mechanical materials are skipped

    **Material Component Configuration:**

    * **Isotropic**: HGO or Ogden hyperelastic models for bulk tissue response
    * **Anisotropic**: Fiber and sheet directions with directional stiffening
    * **Active**: Time-based or calcium-based contractile stress generation

    **Electro-Mechanical Coupling:**

    * **EP Coupled (True)**: Active materials use calcium-dependent activation
      for coupling with electrophysiology simulations
    * **Mechanics Only (False)**: Active materials use predefined activation curves
      for standalone mechanical analysis

    The function logs all material assignments and provides a summary of the
    assignment process including the number of parts processed.
    """
    assignments = 0

    for part in model.parts:
        if (
            not isinstance(part.meca_material, MechanicalMaterialModel)
            or part.meca_material is None
        ):
            try:
                if part.fiber:
                    part.meca_material = _default_myocardium_material(ep_coupled=ep_coupled)
                    LOGGER.info(
                        f"Assigned myocardium mechanical material to part '{part.name}' "
                        f"(EP coupled: {ep_coupled})."
                    )

                    # Disable the active module for passive parts
                    if not part.active:
                        part.meca_material.active = None
                        LOGGER.debug(f"Disabled active component for passive part '{part.name}'.")
                else:
                    part.meca_material = _default_passive_material()
                    LOGGER.info(f"Assigned passive mechanical material to part '{part.name}'.")

                assignments += 1

            except Exception as e:
                error_msg = f"Failed to assign mechanical material to part '{part.name}': {e}"
                LOGGER.error(error_msg)
                raise RuntimeError(error_msg) from e

    LOGGER.info(
        f"Successfully assigned default mechanical materials to "
        f"{assignments}/{len(model.parts)} parts "
    )

    return
