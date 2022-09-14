import matplotlib.pyplot as plt
import numpy as np


class EDPVR:
    """
    ref: Klotz, Stefan, Marc L. Dickstein, and Daniel Burkhoff. Nature protocols 2.9 (2007): 2152-2158.
    """

    # human constant
    An = 27.78  # mmHg
    Bn = 2.76  # mmHg

    def __init__(self, vm, pm):
        self.vm = vm
        self.pm = pm

        self.v0 = self.vm * (0.6 - 0.006 * self.pm)
        self.v30 = self.v0 + (self.vm - self.v0) / (self.pm / self.An) ** (1 / self.Bn)

        if self.pm <= 22:
            self.Beta = np.log10(self.pm / 30) / np.log10(self.vm / self.v30)
            self.Alpha = 30 / self.v30 ** self.Beta
        else:
            v15 = 0.8 * (self.v30 - self.v0) + self.v0
            self.Beta = np.log10(self.pm / 15) / np.log10(self.vm / v15)
            self.Alpha = self.pm / self.vm ** self.Beta

    def get_constants(self):
        return self.Alpha, self.Beta

    def get_pressure(self, volume):
        return self.Alpha * volume ** self.Beta

    def get_volume(self, pressure):
        volume = np.zeros(pressure.shape)
        for i, p in enumerate(pressure):
            volume[i] = (p / self.Alpha) ** (1 / self.Beta)
            # handle singular issue in Klotz curve
            if volume[i] <= self.v0:
                volume[i] = self.v0
        return volume

    def plot_EDPVR(self):

        vv = np.linspace(0, 1.1 * self.vm, num=101)
        pp = self.get_pressure(vv)

        plt.plot(vv, pp, label="Klotz curve")
        plt.plot(self.v0, 0, "o", label="V0")
        plt.plot(self.vm, self.pm, "o", label="Vm,Pm")
        plt.title("EDVPR", fontsize=14)
        plt.xlabel("Volume (mL)", fontsize=14)
        plt.ylabel("Pressure (mmHg)", fontsize=14)
        plt.legend()
        plt.show()
        return


if __name__ == "__main__":
    # healthy baseline
    v_ed = 1.752e02  # mL
    p_ed = 15  # mmHg

    # v_ed = 120  # mL
    # p_ed = 7  # mmHg
    edpvr = EDPVR(v_ed, p_ed)

    edpvr.plot_EDPVR()
