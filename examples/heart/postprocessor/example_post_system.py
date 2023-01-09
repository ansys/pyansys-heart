from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
import matplotlib.pyplot as plt

if __name__ == "__main__":
    base_dir = r"full"
    result = SystemModelPost(base_dir, closed_loop=False)

    fig = result.plot_pv_loop(t_start=3000, t_end=5000)
    # fig.savefig('PV.png')

    if result.closed_loop:
        fig = result._check_total_volume(plot_all=True)

    fig = result.plot_pressure_flow_volume("lv")
    # fig = result.check_prefilling("lv")
    fig = result._check_output("lv")

    if result.type == "BV":
        result.plot_pressure_flow_volume("rv")
        # result.check_prefilling("rv")
        result._check_output("rv")
    plt.show()
