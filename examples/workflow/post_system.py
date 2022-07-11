import numpy as np
from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
import matplotlib.pyplot as plt

if __name__ == "__main__":
    p_ed = np.array([2, 0.5333])
    v_ed = np.array([172.320, 234.799])

    base_dir = (
        r"\\LYOTECSPARE4\wye\pyheartlib_models\h01_bv_corase\bv_closed_impose_filling2"
    )
    result = SystemModelPost(base_dir, p_ed, v_ed, closed_loop=True)

    fig = result.plot_PV(last_loop=True)
    # fig.savefig('PV.png')

    fig = result.check_total_volume(plot_all=True)

    fig = result.plot_pressure_flow_volume("lv", ignore_filling=True, last_loop=True)
    fig = result.check_prefilling("lv")
    fig = result.check_output("lv")

    if result.type == "BV":
        result.plot_pressure_flow_volume("rv", ignore_filling=True, last_loop=True)
        result.check_prefilling("rv")
        result.check_output("rv")
    plt.show()
