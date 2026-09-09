"""
Графики сравнения V–V моделирования.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np

from vv_modeling import DiffusionMethod, ModelingBatch, ModelingRun, _resolve_particles

DiffusionFilter = DiffusionMethod | Literal["all"]

METHOD_LABELS: dict[str, str] = {
    "none": "без дифф.",
    "analytic_factor": "аналит. фактор",
    "method_of_lines": "метод прямых",
}
METHOD_LINESTYLES: dict[str, str] = {
    "none": "-",
    "analytic_factor": "--",
    "method_of_lines": ":",
}
DIFFUSION_METHODS: tuple[DiffusionMethod, ...] = (
    "none",
    "analytic_factor",
    "method_of_lines",
)

DEFAULT_XLIM_BY_V: dict[int, tuple[float, float]] = {
    0: (0, 5),
    1: (0, 5),
    2: (0, 3),
    3: (0, 3),
    4: (0, 10),
    5: (0, 10),
}

DEFAULT_RUN_STYLES: dict[str, dict] = {
    "regression": {"color": "red", "label": None},
    "regression_without_diffusion": {"color": "blue", "label": None},
    "fhofer": {"color": "purple", "label": "FHO-FR fhovv"},
}

_AUTO_RUN_COLORS = ("red", "blue", "purple", "tab:orange", "green", "brown", "teal", "crimson")


def xlim_by_t_max_us(
    t_max_us: Sequence[float] | dict[int, float],
    *,
    t_min_us: float = 0.0,
) -> dict[int, tuple[float, float]]:
    if isinstance(t_max_us, dict):
        return {int(v): (t_min_us, float(t_end)) for v, t_end in t_max_us.items()}
    return {v: (t_min_us, float(t_end)) for v, t_end in enumerate(t_max_us)}


def _err_for_plot(err_spec, n_pts: int, axis_hint: str = ""):
    if err_spec is None:
        return None
    arr = np.atleast_1d(np.asarray(err_spec, dtype=float)).ravel()
    if arr.size == 1:
        v = float(arr[0])
        return None if v == 0.0 else v
    if arr.size != n_pts:
        raise ValueError(
            f"погрешность{axis_hint}: нужна длина 1 или {n_pts}, сейчас {arr.size}"
        )
    return None if np.all(arr == 0) else arr


def _methods_to_plot(diffusion_filter: DiffusionFilter) -> tuple[str, ...]:
    if diffusion_filter == "all":
        return DIFFUSION_METHODS
    return (diffusion_filter,)


def run_display_name(run: ModelingRun, method: str | None = None) -> str:
    title = (run.config.name or run.name.split("::")[0]).strip()
    return f"{title}"


def _legend_sim(
    run: ModelingRun,
    style: str,
    method: str,
    *,
    legend_by_run: dict[tuple[str, str, str], str],
) -> str:
    key_run = (run.name, style, method)
    if key_run in legend_by_run:
        return legend_by_run[key_run]
    return run_display_name(run, method)


@dataclass
class PopulationPlotConfig:
    """Настройки сетки графиков популяций по v."""

    diffusion_per_run: bool = True
    diffusion_filter: DiffusionFilter = "all"
    max_v_plots: int = 6
    ncols: int = 3
    figsize: tuple[float, float] = (16, 10)
    fontsize: float = 10

    plot_styles: dict[str, dict] = field(
        default_factory=lambda: {"default": {"color": "red"}}
    )
    run_styles: dict[str, dict] = field(default_factory=lambda: dict(DEFAULT_RUN_STYLES))
    xlim_by_v: dict[int, tuple[float, float]] = field(
        default_factory=lambda: dict(DEFAULT_XLIM_BY_V)
    )
    t_max_us_by_v: Sequence[float] | dict[int, float] | None = None

    legend_by_run: dict[tuple[str, str, str], str] = field(default_factory=dict)
    legend_exp_dots: str = "Экспериментальные данные"
    legend_exp_fho: str = "FHO-FR (статья)"
    legend_exp_csv: str = "эксп. {style} (v{v}.csv)"
    exp_csv_colors: dict[str, str] = field(
        default_factory=lambda: {"solid": "red", "dashed": "blue"}
    )

    show_experiment: bool = True
    show_article_fho: bool = True
    show_exp_csv_styles: bool = True

    def __post_init__(self) -> None:
        if self.t_max_us_by_v is None:
            return
        xlim = dict(self.xlim_by_v)
        xlim.update(xlim_by_t_max_us(self.t_max_us_by_v))
        self.xlim_by_v = xlim


def completed_run_names(
    batch: ModelingBatch,
    runs: tuple[str, ...] | None = None,
) -> tuple[str, ...]:
    return tuple(r.name for r in _completed_runs(batch, runs))


def completed_run_labels(
    batch: ModelingBatch,
    runs: tuple[str, ...] | None = None,
) -> tuple[str, ...]:
    return tuple(run_display_name(r) for r in _completed_runs(batch, runs))


def _completed_runs(batch: ModelingBatch, runs: tuple[str, ...] | None) -> list[ModelingRun]:
    names = runs if runs is not None else batch.names()
    out = []
    for name in names:
        if name not in batch:
            continue
        run = batch[name]
        if run.suite is not None:
            out.append(run)
    return out


def _run_color(
    run: ModelingRun,
    style: str,
    plot_cfg: PopulationPlotConfig,
    run_index: int,
) -> str:
    rs = plot_cfg.run_styles.get(run.name, {})
    if rs.get("color") is not None:
        return rs["color"]
    return _AUTO_RUN_COLORS[run_index % len(_AUTO_RUN_COLORS)]


def _methods_for_run(run: ModelingRun, plot_cfg: PopulationPlotConfig) -> tuple[str, ...]:
    if plot_cfg.diffusion_per_run:
        return (run.config.diffusion_method,)
    return _methods_to_plot(plot_cfg.diffusion_filter)


def _plot_simulation_layers(
    ax,
    v: int,
    runs: list[ModelingRun],
    plot_cfg: PopulationPlotConfig,
) -> None:
    for run_index, run in enumerate(runs):
        sims = run.simulations
        methods = _methods_for_run(run, plot_cfg)
        for style in plot_cfg.plot_styles:
            if style not in sims:
                continue
            color = _run_color(run, style, plot_cfg, run_index)
            for method in methods:
                if method not in sims[style]:
                    continue
                entry = sims[style][method]
                t = entry["sol"].t * 1e6
                frac = entry["fractional"][v, :]
                ax.plot(
                    t,
                    frac,
                    color=color,
                    linestyle=METHOD_LINESTYLES[method],
                    linewidth=1.5,
                    label=_legend_sim(
                        run,
                        style,
                        method,
                        legend_by_run=plot_cfg.legend_by_run,
                    ),
                )


def _plot_experiment_layers(
    ax,
    v: int,
    run: ModelingRun,
    plot_cfg: PopulationPlotConfig,
) -> None:
    exp_data = run.exp_data
    exp_mod = run.exp_mod_data
    if exp_data is None:
        return

    if plot_cfg.show_experiment and "dots" in exp_data:
        exp_dots = exp_data["dots"][v] if v < len(exp_data["dots"]) else None
        if exp_dots is not None:
            xv, yv, yev, xev = exp_dots
            n_iv = len(xv)
            ax.errorbar(
                xv,
                yv,
                xerr=_err_for_plot(xev, n_iv, " x"),
                yerr=_err_for_plot(yev, n_iv, " y"),
                fmt="s",
                markersize=5.5,
                capsize=3,
                color="tab:orange",
                ecolor="tab:orange",
                elinewidth=1,
                alpha=0.88,
                zorder=5,
                label=plot_cfg.legend_exp_dots.format(v=v),
            )

    if plot_cfg.show_article_fho and exp_mod and "fho" in exp_mod:
        fho = exp_mod["fho"][v] if v < len(exp_mod["fho"]) else None
        if fho is not None and len(fho) > 0:
            t_mod = fho.iloc[:, 0].values.astype(float)
            pop_mod = fho.iloc[:, 1].values.astype(float)
            lbl = plot_cfg.legend_exp_fho.format(v=v)
            ax.plot(t_mod, pop_mod, color="green", alpha=0.5, linewidth=1, zorder=3, label=lbl)

    if plot_cfg.show_exp_csv_styles:
        for style, series in exp_data.items():
            if style in ("dots", "fho") or not isinstance(series, list):
                continue
            if v >= len(series):
                continue
            df = series[v]
            if df is None or len(df) == 0:
                continue
            t_csv = df.iloc[:, 0].values.astype(float)
            pop_csv = df.iloc[:, 1].values.astype(float)
            color = plot_cfg.exp_csv_colors.get(style, "gray")
            lbl = plot_cfg.legend_exp_csv.format(style=style, v=v)
            ax.plot(
                t_csv,
                pop_csv,
                "o-",
                ms=4,
                color=color,
                alpha=0.75,
                label=lbl,
                zorder=4,
            )


def plot_population_grid(
    batch: ModelingBatch,
    *,
    runs: tuple[str, ...] | None = None,
    exp_from: str | None = None,
    plot_cfg: PopulationPlotConfig | None = None,
    ax: np.ndarray | None = None,
) -> tuple[plt.Figure, np.ndarray]:
    """
    Сетка nrows×ncols: доля населённости по v (по умолчанию 2×3 при ncols=3).

    runs — какие прогоны batch рисовать (None = все завершённые).
    exp_from — имя прогона для эксп. точек (None = активный).
    """
    plot_cfg = plot_cfg or PopulationPlotConfig()
    completed = _completed_runs(batch, runs)
    if not completed:
        raise RuntimeError("Нет завершённых прогонов в batch — сначала batch.run_all()")

    ref_run = batch[exp_from] if exp_from else batch.active
    v_max = ref_run.config.v_max
    n_plots = min(plot_cfg.max_v_plots, v_max + 1)

    ncols = plot_cfg.ncols
    nrows = int(np.ceil(n_plots / ncols))
    if ax is None:
        fig, axes = plt.subplots(nrows, ncols, figsize=plot_cfg.figsize)
    else:
        fig = ax.flat[0].figure
        axes = ax

    axes_flat = np.atleast_1d(axes).flatten()
    for idx in range(len(axes_flat)):
        if idx >= n_plots:
            axes_flat[idx].set_visible(False)
            continue

        v = idx
        ax_v = axes_flat[idx]
        _plot_simulation_layers(ax_v, v, completed, plot_cfg)
        _plot_experiment_layers(ax_v, v, ref_run, plot_cfg)

        fs = plot_cfg.fontsize
        ax_v.set_xlabel("Время, мкс", fontsize=fs)
        ax_v.set_ylabel("Доля населённости", fontsize=fs)
        ax_v.set_title(f"v = {v}", fontsize=fs)
        ax_v.tick_params(labelsize=fs)
        ax_v.grid(True, alpha=0.3)
        if v in plot_cfg.xlim_by_v:
            ax_v.set_xlim(*plot_cfg.xlim_by_v[v])
        ax_v.legend(fontsize=fs, loc="best")

    fig.tight_layout()
    return fig, axes_flat


def print_population_summary(
    batch: ModelingBatch,
    run_name: str | None = None,
    *,
    diffusion_method: DiffusionFilter | None = None,
    style: str = "solid",
) -> None:
    run = batch[run_name] if run_name else batch.active
    method = diffusion_method if diffusion_method else run.config.diffusion_method
    if method == "all":
        method = "none"

    sims = run.simulations
    if style not in sims or method not in sims[style]:
        raise KeyError(f"Нет результатов для style={style!r}, method={method!r}")

    sol = sims[style][method]["sol"]
    frac = sims[style][method]["fractional"]
    mult = run.config.rate_multiplier

    print(f"\n=== {run_display_name(run, method)} | kinetics={run.config.kinetics} ===")
    print(f"    k_VV × {mult:g}")
    time_indices = [0, len(sol.t) // 10, len(sol.t) // 5, len(sol.t) // 2, -1]
    v_show = min(6, run.config.v_max + 1)
    for ti in time_indices:
        print(f"  t = {sol.t[ti] * 1e6:.2f} мкс:")
        for v in range(v_show):
            print(f"    v={v}: {frac[v, ti]:.6f}")


def print_population_summaries(
    batch: ModelingBatch,
    runs: tuple[str, ...] | None = None,
    *,
    diffusion_method: DiffusionFilter | None = None,
    style: str = "solid",
) -> None:
    for name in completed_run_names(batch, runs):
        run = batch[name]
        method = (
            diffusion_method
            if diffusion_method and diffusion_method != "all"
            else run.config.diffusion_method
        )
        print_population_summary(
            batch,
            run_name=name,
            diffusion_method=method,
            style=style,
        )


# --- сравнение k_VV ---

def plot_k_vv_comparison(
    run: ModelingRun,
    *,
    article_x: Sequence[float] | np.ndarray | None = None,
    article_y: Sequence[float] | np.ndarray | None = None,
    article_label: str = "статья",
    v_range: tuple[int, int] = (1, 10),
    T: float | None = None,
    show_mm: bool = True,
    show_fhofer: bool = True,
    show_regression: bool = True,
    figsize: tuple[float, float] = (7.5, 4.8),
) -> tuple[plt.Figure, plt.Axes]:
    """
    График сравнения k_VV: статья vs k_vv_mm vs FHO-FR analytical vs регрессия.

    Parameters
    ----------
    run : ModelingRun
        Прогон с настроенной кинетикой.
    article_x, article_y : array-like, optional
        Данные из статьи ( x —振动 quantum number, y — k_VV).
    article_label : str
        Подпись для данных статьи на графике.
    v_range : tuple[int, int]
        Диапазон v для расчёта (start, stop+1).
    T : float, optional
        Температура для расчёта. По умолчанию из run.config.T.
    show_mm : bool
        Показать k_vv_mm (численное интегрирование).
    show_fhofer : bool
        Показать k_vv_fhofer (аналитическая формула).
    show_regression : bool
        Показать регрессионные коэффициенты.
    figsize : tuple
        Размер фигуры.
    """
    import os
    from functools import partial
    from concurrent.futures import ProcessPoolExecutor

    from k_VV import k_vv_fhofer
    from k_VV import k_lookup

    cfg = run.config
    m1, m2 = _resolve_particles(cfg.particle_a, cfg.particle_b)
    T_calc = T if T is not None else cfg.T
    v_int = np.arange(v_range[0], v_range[1])

    fig, ax = plt.subplots(figsize=figsize)

    if article_x is not None and article_y is not None:
        ax.plot(article_x, article_y, "s-", ms=7, label=article_label)

    if show_mm:
        from k_VV import k_vv_mm_for_v

        _compute_k_vv = partial(k_vv_mm_for_v, m1.name, m2.name, T=T_calc)

        max_workers = min(len(v_int), os.cpu_count() or 8)
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            k_from_mm = np.array(list(ex.map(_compute_k_vv, v_int.tolist(), chunksize=1)))
        k_mm_plot = np.where(k_from_mm > 0, k_from_mm, np.nan)
        ax.plot(v_int, k_mm_plot, "o-", ms=6, label=f"k_vv_mm, T = {T_calc} K")

    if show_fhofer:
        k_fhofer = np.array([
            k_vv_fhofer(m1, m2, i1=1, f1=0, i2=v - 1, f2=v, T=T_calc)
            for v in v_int
        ])
        ax.plot(v_int, k_fhofer, "o-", label="FHO-FR (аналит.)")

    if show_regression and cfg.coeffs:
        k_reg = np.array([
            k_lookup(cfg.coeffs, i1=1, f1=0, i2=v - 1, f2=v, T=T_calc)
            for v in v_int
            if k_lookup(cfg.coeffs, i1=1, f1=0, i2=v - 1, f2=v, T=T_calc) > 0
        ])
        v_reg = v_int[:len(k_reg)]
        if len(k_reg) > 0:
            ax.plot(v_reg, k_reg, "^-", ms=6, label="регрессия")

    ax.set_yscale("log")
    ax.set_xlabel(r"$v$ ($i_1=v$, $f_1=v-1$, $i_2=0$, $f_2=1$)")
    ax.set_ylabel(r"$k_{\mathrm{VV}}$ (см³/с)")
    ax.set_title(f"{cfg.particle_a}–{cfg.particle_b} V–V: k_VV при T = {T_calc} K")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()
    plt.tight_layout()
    return fig, ax
