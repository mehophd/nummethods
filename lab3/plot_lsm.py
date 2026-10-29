import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

for path in ("/usr/share/fonts/truetype/chinese/NotoSansSC-Regular.ttf",
             "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
    try:
        fm.fontManager.addfont(path)
    except Exception:
        pass
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Noto Sans SC"]
plt.rcParams["axes.unicode_minus"] = False


def read_data(path):
    table_x, table_y, grid_x, grid_f1, grid_f2 = [], [], [], [], []
    block = 0
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "Блок 1" in line:
                    block = 1
                elif "Блок 2" in line:
                    block = 2
                continue
            vals = [float(p) for p in line.split()]
            if block == 1 and len(vals) >= 2:
                table_x.append(vals[0])
                table_y.append(vals[1])
            elif block == 2 and len(vals) >= 3:
                grid_x.append(vals[0])
                grid_f1.append(vals[1])
                grid_f2.append(vals[2])
    return (np.array(table_x), np.array(table_y),
            np.array(grid_x), np.array(grid_f1), np.array(grid_f2))


def main():
    data_file = sys.argv[1] if len(sys.argv) > 1 else "lsm_plot_data.txt"
    out_png = sys.argv[2] if len(sys.argv) > 2 else "lsm_graph.png"

    tx, ty, gx, gf1, gf2 = read_data(data_file)

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=150, constrained_layout=True)

    xs = np.linspace(-0.85, 0.95, 400)
    ax.plot(xs, np.arccos(xs), ":", color="grey", lw=1.3,
            label="$y=\\arccos x$ — приближаемая функция")

    ax.plot(gx, gf1, "-", color="#1f77b4", lw=2,
            label="$F_1(x)$ — многочлен 1-й степени (МНК)")
    ax.plot(gx, gf2, "--", color="#d62728", lw=2,
            label="$F_2(x)$ — многочлен 2-й степени (МНК)")
    ax.plot(tx, ty, "o", color="black", ms=7, zorder=5,
            label="Табличные данные")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Задание 3.3. Приближение табличной функции по методу\n"
                 "наименьших квадратов (МНК)")
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend(loc="upper right")

    fig.savefig(out_png)
    print("График сохранён в файл:", out_png)


if __name__ == "__main__":
    main()
