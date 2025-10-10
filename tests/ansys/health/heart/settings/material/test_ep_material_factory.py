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

import importlib.util as importlib_util
import sys

from ansys.health.heart.settings.material import ep_material_factory as factory
from ansys.health.heart.settings.material.ep_material import EPSolverType


def _make_dummy_defaults():
    """Create a dummy defaults module with eikonal and monodomain entries."""
    module_name = "ansys.health.heart.settings.defaults.electrophysiology"
    spec = importlib_util.spec_from_loader(module_name, loader=None)
    dummy = importlib_util.module_from_spec(spec)

    dummy.default_myocardium_material_eikonal = {
        "sigma_fiber": 0.7,
        "sigma_sheet": 0.2,
    }
    dummy.default_myocardium_material_monodomain = {
        "sigma_fiber": 0.5,
        "sigma_sheet": 0.1,
    }
    dummy.default_beam_material_eikonal = {"sigma": 1.0}
    dummy.default_beam_material_monodomain = {"sigma": 1.0}

    return module_name, dummy


# Mocks ActiveBeam and Active EP Material classes.
class MockActive:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


class MockActiveBeam:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


def test_get_default_myocardium_material_selects_correct_defaults(monkeypatch):
    module_name, dummy = _make_dummy_defaults()
    monkeypatch.setitem(sys.modules, module_name, dummy)

    # Capture calls to Active
    monkeypatch.setattr(factory, "Active", MockActive)

    # Eikonal via enum
    res = factory.get_default_myocardium_material(EPSolverType.EIKONAL)
    assert isinstance(res, MockActive)
    assert res.kwargs == dummy.default_myocardium_material_eikonal

    # Reaction-Eikonal should behave like Eikonal
    res = factory.get_default_myocardium_material(EPSolverType.REACTION_EIKONAL)
    assert isinstance(res, MockActive)
    assert res.kwargs == dummy.default_myocardium_material_eikonal

    # Monodomain via string input
    res = factory.get_default_myocardium_material("Monodomain")
    assert isinstance(res, MockActive)
    assert res.kwargs == dummy.default_myocardium_material_monodomain


def test_get_default_conduction_system_material_selects_correct_defaults(monkeypatch):
    module_name, dummy = _make_dummy_defaults()
    monkeypatch.setitem(sys.modules, module_name, dummy)

    monkeypatch.setattr(factory, "ActiveBeam", MockActiveBeam)

    # Eikonal via enum
    res = factory.get_default_conduction_system_material(EPSolverType.EIKONAL)
    assert isinstance(res, MockActiveBeam)
    assert res.kwargs == dummy.default_beam_material_eikonal

    # Monodomain via string
    res = factory.get_default_conduction_system_material("Monodomain")
    assert isinstance(res, MockActiveBeam)
    assert res.kwargs == dummy.default_beam_material_monodomain
