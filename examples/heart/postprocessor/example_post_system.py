from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
import matplotlib.pyplot as plt

if __name__ == "__main__":
    base_dir = r"D:\Heart20\healthy20\health03_BV_2mm\simulation\main-mechanics"
    result = SystemModelPost(base_dir)
    # result._check_output()

    fig = result.plot_pv_loop()
    fig = result.plot_pressure_flow_volume(result.lv)
    fig = result.plot_pressure_flow_volume(result.rv)
    plt.show()
