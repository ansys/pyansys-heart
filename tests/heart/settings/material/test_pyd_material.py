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

"""Test module for pyd_material.py."""

from pydantic import ValidationError
import pytest

from ansys.health.heart.settings.material._pyd_material import (
    ACTIVE,
    ANISO,
    ISO,
    ActiveModel1,
    ActiveModel3,
    HGOFiber,
    Mat295,
    MechanicalMaterialModel,
    NeoHookean,
)
from ansys.health.heart.settings.material.curve import ActiveCurve, constant_ca2


class TestISO:
    """Test ISO class."""

    def test_iso_hgo_model_valid(self):
        """Test valid HGO model parameters."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        assert iso.itype == -3
        assert iso.k1 == 1.0
        assert iso.k2 == 2.0

    def test_iso_ogden_model_valid(self):
        """Test valid Ogden model parameters."""
        iso = ISO(itype=-1, mu1=1.0, alpha1=2.0)
        assert iso.itype == -1
        assert iso.mu1 == 1.0
        assert iso.alpha1 == 2.0

    def test_iso_incompatible_itype_with_hgo(self):
        """Test incompatible itype with HGO parameters."""
        with pytest.raises(ValueError) as exc_info:
            ISO(itype=-1, k1=1.0, k2=2.0)
        assert "both mu1 and alpha1 must be specified" in str(exc_info.value)

    def test_iso_incompatible_itype_with_ogden(self):
        """Test incompatible itype with Ogden parameters."""
        with pytest.raises(ValueError) as exc_info:
            ISO(itype=-3, mu1=1.0, alpha1=2.0)
        assert "both k1 and k2 must be specified" in str(exc_info.value)

    def test_iso_invalid_input(self):
        """Test invalid input without proper parameters."""
        with pytest.raises(ValueError) as exc_info:
            ISO(itype=-3)
        assert "both k1 and k2 must be specified" in str(exc_info.value)

    def test_iso_kappa_updates_nu(self):
        """Test that kappa updates Poisson's ratio."""
        iso_hgo = ISO(itype=-3, k1=1.0, k2=2.0, kappa=100.0)
        expected_nu = (3 * 100.0 - 2 * 1.0) / (6 * 100.0 + 2 * 1.0)
        assert abs(iso_hgo.nu - expected_nu) < 1e-10

        iso_ogden = ISO(itype=-1, mu1=1.0, alpha1=2.0, kappa=100.0)
        expected_nu = (3 * 100.0 - 2 * 1.0) / (6 * 100.0 + 2 * 1.0)
        assert abs(iso_ogden.nu - expected_nu) < 1e-10

    def test_iso_default_values(self):
        """Test default values."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        assert iso.beta == 0.0
        assert iso.nu == 0.499

    def test_iso_serialization(self):
        """Test serialization and deserialization."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0, beta=0.1)
        data = iso.model_dump()
        iso_restored = ISO(**data)
        assert iso_restored.itype == iso.itype
        assert iso_restored.k1 == iso.k1
        assert iso_restored.k2 == iso.k2
        assert iso_restored.beta == iso.beta


class TestHGOFiber:
    """Test HGOFiber class."""

    def test_hgo_fiber_defaults(self):
        """Test default values."""
        fiber = HGOFiber()
        assert fiber.k1 is None
        assert fiber.k2 is None
        assert fiber.a == 0.0
        assert fiber.b == 1.0
        assert fiber._theta is None
        assert fiber._ftype == 1
        assert fiber._fcid == 0

    def test_hgo_fiber_custom_values(self):
        """Test custom values."""
        fiber = HGOFiber(k1=10.0, k2=20.0, a=0.5, b=0.8)
        assert fiber.k1 == 10.0
        assert fiber.k2 == 20.0
        assert fiber.a == 0.5
        assert fiber.b == 0.8

    def test_hgo_fiber_serialization(self):
        """Test serialization excludes private fields."""
        fiber = HGOFiber(k1=10.0, k2=20.0, _theta=45.0)
        data = fiber.model_dump()
        # Private fields should be excluded
        assert "_theta" not in data
        assert "_ftype" not in data
        assert "_fcid" not in data
        # Public fields should be included
        assert data["k1"] == 10.0
        assert data["k2"] == 20.0


class TestANISO:
    """Test ANISO class."""

    def test_aniso_default_single_fiber(self):
        """Test default behavior with single fiber."""
        aniso = ANISO()
        assert len(aniso.fibers) == 1
        assert aniso._nf == 1
        assert aniso._intype == 0
        assert aniso.fibers[0]._theta == 0.0

    def test_aniso_two_fibers(self):
        """Test with two fibers."""
        fiber1 = HGOFiber(k1=1.0, k2=2.0)
        fiber2 = HGOFiber(k1=3.0, k2=4.0)
        aniso = ANISO(fibers=[fiber1, fiber2])
        assert len(aniso.fibers) == 2
        assert aniso._nf == 2
        assert aniso._intype == 0
        assert aniso.fibers[0]._theta == 0.0
        assert aniso.fibers[1]._theta == 90.0

    def test_aniso_two_fibers_with_interaction(self):
        """Test two fibers with interaction terms."""
        fiber1 = HGOFiber(k1=1.0, k2=2.0)
        fiber2 = HGOFiber(k1=3.0, k2=4.0)
        aniso = ANISO(fibers=[fiber1, fiber2], k1fs=5.0, k2fs=6.0)
        assert aniso._nf == 2
        assert aniso._intype == 1

    def test_aniso_single_fiber_with_interaction_error(self):
        """Test error when single fiber has interaction terms."""
        fiber1 = HGOFiber(k1=1.0, k2=2.0)
        with pytest.raises(ValidationError) as exc_info:
            ANISO(fibers=[fiber1], k1fs=5.0, k2fs=6.0)
        assert "One fiber cannot have an interaction term" in str(exc_info.value)

    def test_aniso_invalid_fiber_count(self):
        """Test error with invalid number of fibers."""
        fiber1 = HGOFiber(k1=1.0, k2=2.0)
        fiber2 = HGOFiber(k1=3.0, k2=4.0)
        fiber3 = HGOFiber(k1=5.0, k2=6.0)
        with pytest.raises(ValidationError) as exc_info:
            ANISO(fibers=[fiber1, fiber2, fiber3])
        assert "No. of fibers must be 1 or 2" in str(exc_info.value)

    def test_aniso_repr_includes_computed_fields(self):
        """Test that __repr__ includes computed fields."""
        aniso = ANISO()
        repr_str = repr(aniso)
        assert "nf=1" in repr_str
        assert "intype=0" in repr_str

    def test_aniso_default_values(self):
        """Test default values."""
        aniso = ANISO()
        assert aniso.atype == -1
        assert aniso.k1fs is None
        assert aniso.k2fs is None
        assert aniso.vec_a == (1.0, 0.0, 0.0)
        assert aniso.vec_d == (0.0, 1.0, 0.0)


class TestActiveModel:
    """Test ActiveModel nested classes."""

    def test_active_model1_defaults(self):
        """Test ActiveModel.Model1 default values."""
        model1 = ActiveModel1()
        assert model1.t0 is None
        assert model1.ca2ion is None
        assert model1.ca2ionm == 4.35
        assert model1.n == 2
        assert model1.taumax == 0.125
        assert model1.stf == 0.0
        assert model1.b == 4.75
        assert model1.l0 == 1.58
        assert model1.l == 1.85
        assert model1.dtmax == 150
        assert model1.mr == 1048.9
        assert model1.tr == -1629.0

    def test_active_model1_custom_values(self):
        """Test ActiveModel.Model1 with custom values."""
        model1 = ActiveModel1(t0=1.0, ca2ion=2.0, n=3)
        assert model1.t0 == 1.0
        assert model1.ca2ion == 2.0
        assert model1.n == 3

    def test_active_model3_defaults(self):
        """Test ActiveModel.Model3 default values."""
        model3 = ActiveModel3()
        assert model3.t0 is None
        assert model3.ca2ion50 == 1.0
        assert model3.n == 1.0
        assert model3.f == 0.0
        assert model3.l == 1.0
        assert model3.eta == 0.0
        assert model3.sigmax is None

    def test_active_model3_custom_values(self):
        """Test ActiveModel.Model3 with custom values."""
        model3 = ActiveModel3(t0=2.0, ca2ion50=1.5, sigmax=100.0)
        assert model3.t0 == 2.0
        assert model3.ca2ion50 == 1.5
        assert model3.sigmax == 100.0


class TestACTIVE:
    """Test ACTIVE class."""

    def test_active_defaults(self):
        """Test default values."""
        active = ACTIVE()
        assert active.acid is None
        assert active._actype == 1  # Deduced from Model1
        assert active.acthr == 0.1  # From ca2_curve threshold
        assert active.acdir == 1
        assert active.sf == 1.0
        assert active.ss == 0.0
        assert active.sn == 0.0
        assert isinstance(active.model, ActiveModel1)
        assert active.ca2_curve is not None

    def test_active_with_model1(self):
        """Test ACTIVE with Model1."""
        model1 = ActiveModel1(t0=1.0)
        active = ACTIVE(model=model1)
        assert active._actype == 1
        assert isinstance(active.model, ActiveModel1)
        assert active.model.t0 == 1.0

    def test_active_with_model3(self):
        """Test ACTIVE with Model3."""
        model3 = ActiveModel3(t0=2.0)
        active = ACTIVE(model=model3)
        assert active._actype == 3
        assert isinstance(active.model, ActiveModel3)
        assert active.model.t0 == 2.0

    def test_active_threshold_from_curve(self):
        """Test threshold is taken from ca2_curve."""
        curve = ActiveCurve(constant_ca2(), threshold=0.5, type="ca2")
        active = ACTIVE(ca2_curve=curve)
        assert active.acthr == 0.5

    def test_active_custom_values(self):
        """Test custom values."""
        active = ACTIVE(acdir=2, sf=0.8, ss=0.1, sn=0.05)
        assert active.acdir == 2
        assert active.sf == 0.8
        assert active.ss == 0.1
        assert active.sn == 0.05


class TestMat295:
    """Test Mat295 class."""

    def test_mat295_minimal(self):
        """Test Mat295 with minimal required parameters."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        mat = Mat295(rho=1000.0, iso=iso)
        assert mat.rho == 1000.0
        assert mat.iso == iso
        assert mat.aopt == 2.0
        assert mat.aniso is None
        assert mat.active is None

    def test_mat295_with_aniso(self):
        """Test Mat295 with anisotropic module."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        aniso = ANISO()
        mat = Mat295(rho=1000.0, iso=iso, aniso=aniso)
        assert mat.aniso == aniso

    def test_mat295_with_active(self):
        """Test Mat295 with active module."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        active = ACTIVE()
        mat = Mat295(rho=1000.0, iso=iso, active=active)
        assert mat.active == active

    def test_mat295_complete(self):
        """Test Mat295 with all modules."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        aniso = ANISO()
        active = ACTIVE()
        mat = Mat295(rho=1000.0, iso=iso, aniso=aniso, active=active, aopt=3.0)
        assert mat.rho == 1000.0
        assert mat.iso == iso
        assert mat.aniso == aniso
        assert mat.active == active
        assert mat.aopt == 3.0

    def test_mat295_serialization(self):
        """Test Mat295 serialization."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        mat = Mat295(rho=1000.0, iso=iso)
        data = mat.model_dump()
        mat_restored = Mat295(**data)
        assert mat_restored.rho == mat.rho
        assert mat_restored.iso.itype == mat.iso.itype
        assert mat_restored.aopt == mat.aopt


class TestNeoHookean:
    """Test NeoHookean class (deprecated)."""

    def test_neohookean_basic(self):
        """Test basic NeoHookean creation."""
        neo = NeoHookean(rho=1000.0, c10=1.0)
        assert neo.rho == 1000.0
        assert neo.c10 == 1.0
        assert neo.kappa is None
        assert neo.nu is None

    def test_neohookean_with_kappa(self):
        """Test NeoHookean with bulk modulus."""
        neo = NeoHookean(rho=1000.0, c10=1.0, kappa=100.0)
        mu = neo.c10 * 2
        expected_nu = (3 * 100.0 - 2 * mu) / (6 * 100.0 + 2 * mu)
        assert abs(neo.nu - expected_nu) < 1e-10

    def test_neohookean_with_nu(self):
        """Test NeoHookean with explicit Poisson's ratio."""
        neo = NeoHookean(rho=1000.0, c10=1.0, nu=0.45)
        assert neo.nu == 0.45


class TestMechanicalMaterialModel:
    """Test MechanicalMaterialModel base class."""

    def test_dummy_material(self):
        """Test DummyMaterial."""
        dummy = MechanicalMaterialModel.DummyMaterial()
        assert repr(dummy) == "Material is empty."

    def test_inheritance(self):
        """Test that all material classes inherit from MechanicalMaterialModel."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        mat295 = Mat295(rho=1000.0, iso=iso)
        neo = NeoHookean(rho=1000.0, c10=1.0)

        assert isinstance(mat295, MechanicalMaterialModel)
        assert isinstance(neo, MechanicalMaterialModel)


class TestIntegration:
    """Integration tests combining multiple classes."""

    def test_complete_material_setup(self):
        """Test creating a complete material with all modules."""
        # Create ISO module
        iso = ISO(itype=-3, k1=10.0, k2=20.0, beta=0.1)

        # Create fibers for ANISO module
        fiber1 = HGOFiber(k1=5.0, k2=15.0, a=0.2)
        fiber2 = HGOFiber(k1=7.0, k2=17.0, a=0.3)
        aniso = ANISO(fibers=[fiber1, fiber2], k1fs=2.0, k2fs=3.0)

        # Create ACTIVE module with Model3
        model3 = ActiveModel3(ca2ion50=1.5, sigmax=150.0)
        curve = ActiveCurve(constant_ca2(), threshold=0.2, type="ca2")
        active = ACTIVE(model=model3, ca2_curve=curve)

        # Create complete Mat295
        mat = Mat295(rho=1050.0, iso=iso, aniso=aniso, active=active)

        # Verify all components
        assert mat.rho == 1050.0
        assert mat.iso.k1 == 10.0
        assert mat.aniso._nf == 2
        assert mat.aniso._intype == 1
        assert mat.active._actype == 3
        assert mat.active.acthr == 0.2

    @pytest.mark.xfail(reason="Serialization of the ACTIVE model not fully functional.")
    def test_serialization_roundtrip(self):
        """Test complete serialization and deserialization."""
        # Create complex material
        iso = ISO(itype=-3, k1=10.0, k2=20.0)
        aniso = ANISO()
        active = ACTIVE()
        mat = Mat295(rho=1000.0, iso=iso, aniso=aniso, active=active)

        # Serialize to dict
        data = mat.model_dump()

        # Deserialize back
        mat_restored = Mat295(**data)

        # Verify everything matches
        assert mat_restored.rho == mat.rho
        assert mat_restored.iso.k1 == mat.iso.k1
        assert mat_restored.iso.k2 == mat.iso.k2
        assert mat_restored.aniso._nf == mat.aniso._nf
        assert mat_restored.active._actype == mat.active._actype

    def test_json_serialization(self):
        """Test JSON serialization."""
        iso = ISO(itype=-3, k1=1.0, k2=2.0)
        mat = Mat295(rho=1000.0, iso=iso)

        # Test JSON serialization
        json_str = mat.model_dump_json()
        assert isinstance(json_str, str)
        assert "rho" in json_str
        assert "1000.0" in json_str

        # Test JSON deserialization
        mat_from_json = Mat295.model_validate_json(json_str)
        assert mat_from_json.rho == mat.rho
        assert mat_from_json.iso.k1 == mat.iso.k1
