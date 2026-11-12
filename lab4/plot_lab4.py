import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

for path in (
    "/usr/share/fonts/truetype/chinese/NotoSansSC-Regular.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
):
    try:
        fm.fontManager.addfont(path)
    except Exception:
        pass
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Noto Sans SC"]
plt.rcParams["axes.unicode_minus"] = False


def read_columns(path, ncols):
    """Читает таблицу из ncols вещественных столбцов (строки с # - заголовки)."""
    cols = [[] for _ in range(ncols)]
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vals = [float(p) for p in line.split()]
            if len(vals) >= ncols:
                for i in range(ncols):
                    cols[i].append(vals[i])
    return [np.array(c) for c in cols]


def plot_cauchy(data_dir):
    x, yex, yeu, yrk, yad = read_columns(
        os.path.join(data_dir, "cauchy_plot_data.txt"), 5
    )

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=150, constrained_layout=True)

    ax.plot(
        x,
        yex,
        "-",
        color="grey",
        lw=1.8,
        label="$y=(x^2+1/x)\\,e^{\\sqrt{x}}$ — точное решение",
    )
    ax.plot(x, yeu, "o--", color="#1f77b4", lw=1.4, ms=5, label="Метод Эйлера (4.2)")
    ax.plot(
        x, yrk, "s:", color="#d62728", lw=1.6, ms=5, label="Метод Рунге-Кутты 4 (4.16)"
    )
    ax.plot(x, yad, "^-.", color="#2ca02c", lw=1.4, ms=5, label="Метод Адамса 4 (4.25)")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        "Задание 4.1. Задача Коши для ОДУ 2-го порядка\n"
        "$y''-\\frac{1}{\\sqrt{x}}\\,y'+\\frac{x+\\sqrt{x}-8}{4x^2}\\,"
        "y=0$,  $y(1)=2e$,  $y'(1)=2e$"
    )
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend(loc="upper left", fontsize=9)
    fig.savefig(os.path.join(data_dir, "cauchy_graph.png"))
    plt.close(fig)
    print("График сохранён:", os.path.join(data_dir, "cauchy_graph.png"))


def plot_bvp(data_dir):
    x, yex, ysh, yfd = read_columns(os.path.join(data_dir, "bvp_plot_data.txt"), 4)

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=150, constrained_layout=True)

    ax.plot(x, yex, "-", color="grey", lw=1.8, label="$y=e^x(x^2+1)$ — точное решение")
    ax.plot(
        x,
        ysh,
        "o--",
        color="#d62728",
        lw=1.5,
        ms=5,
        label="Метод стрельбы (РК4, 4-й порядок)",
    )
    ax.plot(
        x,
        yfd,
        "s:",
        color="#1f77b4",
        lw=1.5,
        ms=5,
        label="Конечно-разностный метод (2-й порядок)",
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        "Задание 4.2. Краевая задача для ОДУ 2-го порядка\n"
        "$x\\,y''-(2x+1)\\,y'+(x+1)\\,y=0$,  "
        "$y'(0)=1$,  $y'(1)-2y(1)=0$"
    )
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend(loc="upper left", fontsize=9)
    fig.savefig(os.path.join(data_dir, "bvp_graph.png"))
    plt.close(fig)
    print("График сохранён:", os.path.join(data_dir, "bvp_graph.png"))


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    plot_cauchy(data_dir)
    plot_bvp(data_dir)


if __name__ == "__main__":
    main()
