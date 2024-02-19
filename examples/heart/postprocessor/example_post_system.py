""" For information only, not yet stable
    require previously launched simulation
    change paths accordingly"""

import matplotlib.pyplot as plt

from ansys.heart.postprocessor.SystemModelPost import SystemModelPost

if __name__ == "__main__":
    base_dir = r"my_base_directory"
    result = SystemModelPost(base_dir)
    # result._check_output()

    fig = result.plot_pv_loop()
    fig = result.plot_pressure_flow_volume(result.lv)
    fig = result.plot_pressure_flow_volume(result.rv_system)
    plt.show()
