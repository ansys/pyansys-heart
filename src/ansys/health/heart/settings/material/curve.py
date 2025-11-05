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

"""Module for active stress curve."""

from typing import Literal, Tuple

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_serializer,
    field_validator,
    model_validator,
)

from ansys.health.heart import LOG as LOGGER


def strocchi_active(t_end=800, t_act=0) -> tuple[np.ndarray, np.ndarray]:
    """
    Active stress in doi.org/10.1371/journal.pone.0235145.

    T_peak is described in MAT_295

    Parameters
    ----------
    t_end : int, default: 800
        Heart beat period.
    t_act : int, default: 0
        Start time.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (time, stress) array
    """
    # parameters used in Strocchi in ms
    tau_r = 130 * 800 / t_end
    tau_d = 100 * 800 / t_end
    tau_dur = 550 * 800 / t_end
    tau_emd = 0.0  # EM coupling delay
    # t_act = 0.0  # activation time from Eikonel model

    def _stress():
        # Active tension
        t = np.linspace(0, t_end, 101)
        active_stress = np.zeros(t.shape)
        ts = t - t_act - tau_emd
        for i, tt in enumerate(ts):
            if 0 < tt < tau_dur:
                active_stress[i] = np.tanh(tt / tau_r) ** 2 * np.tanh((tau_dur - tt) / tau_d) ** 2
        return (t, active_stress)

    return _stress()


def kumaraswamy_active(t_end=1000) -> tuple[np.ndarray, np.ndarray]:
    """
    Active stress in  Gaëtan Desrues doi.org/10.1007/978-3-030-78710-3_43.

    T_peak is described in MAT295

    Parameters
    ----------
    t_end : int, default: 1000
        Heart beat duration.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (timen,stress) array
    """
    apd90 = 250 * t_end / 1000  # action potential duration
    time_repolarization = 750 * t_end / 1000  #  repolarization time

    time = np.linspace(0, t_end, 101)
    stress = np.zeros(time.shape)

    def _kumaraswamy(a, b, x):
        return 1 - (1 - x**a) ** b

    for i, t in enumerate(time):
        if t < apd90:
            stress[i] = _kumaraswamy(2, 1.5, t / apd90)
        elif t < time_repolarization:
            stress[i] = -_kumaraswamy(2, 3, (t - apd90) / (time_repolarization - apd90)) + 1
    return (time, stress)


def constant_ca2(tb: float = 800, ca2ionm: float = 4.35) -> tuple[np.ndarray, np.ndarray]:
    """Constant ca2 curve for Active model 1.

    Parameters
    ----------
    tb : float, default: 800
        Heart beat period.
    ca2ionm : float, default: 4.35
        Amplitude, which equals ca2ionm in MAT_295.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (time, stress) array
    """
    t = np.linspace(0, tb, 101)
    v = np.ones((101)) * ca2ionm
    # set to 0 so across threshold at start of beat
    v[0:1] = 0
    # set to below threshold at the last 95% of period
    # eg. active stress will be disabled at 760ms with tb of 800ms
    v[-5:] = 0
    return (t, v)


class ActiveCurve(BaseModel):
    """Pydantic-backed ActiveCurve."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    func: Tuple[np.ndarray, np.ndarray] = None
    type: Literal["stress", "ca2"] = "ca2"
    threshold: float = 0.5e-6
    n_beat: int = 5

    # Derived values. exclude these from
    # json serialization.
    time: np.ndarray | None = Field(default=None, exclude=True)
    t_beat: float | None = Field(default=None, exclude=True)
    ca2: np.ndarray | None = Field(default=None, exclude=True)
    stress: np.ndarray | None = Field(default=None, exclude=True)

    @field_validator("func", mode="before")
    def _func_validator(cls, v):  # noqa: N805
        """Accept lists/tuples or numpy arrays and return tuple[np.ndarray, np.ndarray]."""
        if v is None:
            raise ValueError("func must be provided as (time, values) arrays")

        # Expect a sequence of length 2
        if not (isinstance(v, (list, tuple)) and len(v) == 2):
            raise ValueError("func must be a tuple/list of (time, values)")

        t, y = v
        t_arr = np.asarray(t)
        y_arr = np.asarray(y)

        if t_arr.ndim != 1 or y_arr.ndim != 1:
            raise ValueError("func arrays must be 1-dimensional")
        if t_arr.shape != y_arr.shape:
            raise ValueError("func arrays must have the same shape")
        if t_arr.size == 0:
            raise ValueError("func arrays must not be empty")
        if np.any(np.diff(t_arr) <= 0):
            raise ValueError("func time array must be strictly increasing")

        return (t_arr, y_arr)

    @model_validator(mode="after")
    def _post_init(self):
        # preserve public API names used by callers
        self.time = self.func[0]
        self.t_beat = float(self.time[-1])

        if self.type == "stress":
            # reset threshold as current implementation does
            self.threshold = 0.5e-6
            self.stress = self.func[1]
            self.ca2 = self._stress_to_ca2(self.stress)
        else:
            self.ca2 = self.func[1]
            self.stress = None

        # run the same checks
        self._check_threshold()
        return self

    # optional: serialize numpy arrays to lists for model_dump / JSON
    @field_serializer("func")
    def _serialize_func(self, func: tuple[np.ndarray, np.ndarray], info):
        if isinstance(func[0], np.ndarray) and isinstance(func[1], np.ndarray):
            return (func[0].tolist(), func[1].tolist())
        else:
            LOGGER.error("Failed to serialize func")
            return None

    def _check_threshold(self):
        if np.max(self.ca2) < self.threshold or np.min(self.ca2) > self.threshold:
            raise ValueError("Threshold must cross ca2+ curve at least once")

    def _stress_to_ca2(self, stress: np.ndarray) -> np.ndarray:
        if np.min(stress) < 0 or np.max(stress) > 1.0:
            raise ValueError("Stress curve must be between 0-1.")
        ca2 = 1 / (1 - 0.999 * stress) - 1
        ca2[0] = 0.0
        ca2[1:] += 2 * self.threshold
        return ca2

    def _repeat(self, curve: Tuple[np.ndarray, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
        t = np.copy(curve[0])
        v = np.copy(curve[1])
        for ii in range(1, self.n_beat):
            t = np.append(t, curve[0][1:] + ii * self.t_beat)
            v = np.append(v, curve[1][1:])
        return (t, v)

    @property
    def dyna_input(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return LS-DYNA input arrays."""
        return self._repeat((self.time, self.ca2))

    def plot_time_vs_ca2(self):
        """Plot time vs ca2."""
        import matplotlib.pyplot as plt

        t, v = self.dyna_input
        fig, ax = plt.subplots()
        ax.plot(t, v, label="ca2")
        ax.axhline(self.threshold, color="r", linestyle="--", label="threshold")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Ca2+")
        ax.set_title("Active Ca2+ Curve")
        ax.legend()
        return fig

    def plot_time_vs_stress(self):
        """Plot time vs stress."""
        if self.type != "stress":
            raise ValueError("Curve type is not 'stress', cannot plot stress.")

        import matplotlib.pyplot as plt

        t, v = self._repeat((self.time, self.stress))
        fig, ax = plt.subplots()
        ax.plot(t, v, label="stress")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Stress (normalized)")
        ax.set_title("Active Stress Curve")
        ax.legend()
        return fig
