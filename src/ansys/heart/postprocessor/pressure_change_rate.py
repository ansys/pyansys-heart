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

"""Compute rate of pressure change."""

from ansys.heart.preprocessor.models import HeartModel
from ansys.heart.postprocessor.dpf_utils import ICVoutReader
import numpy as np
from matplotlib import pyplot as plt

class PressureChangeRate:
    
    def __init__(self, model: HeartModel, binout_file):
        self.model = model
        self.icvout = ICVoutReader(binout_file)
        
        self.t = self.icvout.get_time()
        self.q = self.icvout.get_flowrate(1)
        self.p = self.icvout.get_pressure(1)*7500.6
        self.dp = np.gradient(self.p, self.t*1e-3, edge_order=2)
        Q_ratio = np.max(abs(self.q))/np.min(abs(self.q))
        tolerance = Q_ratio/1000
        q0 = np.max(abs(self.q))/tolerance
        t0 = int(np.argwhere(abs(self.q)< q0)[0])
        self.dp_max = np.max(self.dp[t0:]) 
    
        
    def plot_dpdt(self):
        fig = plt.figure()
        plt.subplot(211)
        plt.plot(self.t, self.p)
        plt.ylabel('Pressure (mmHg)')
        plt.title('Pressure')
        plt.subplot(212)
        plt.plot(self.t, self.dp)
        plt.xlabel('Time (ms)')
        plt.ylabel('dP/dt (mmHg/s)')
        plt.title('dP/dt')
        plt.plot(self.t[np.argwhere(self.dp==self.dp_max)], self.dp_max, 'r*', label='(dP/dt)max')
        plt.legend()
        return fig
    