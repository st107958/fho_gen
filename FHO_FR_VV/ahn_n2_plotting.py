"""
Графики сравнения V–V моделирования N₂–N₂ (Ahn & Adamovich).
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np

from ahn_n2_modeling import DiffusionMethod, ModelingBatch, ModelingRun

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

# цвет и шаблон подписи для прогонов batch (имя → параметры)
DEFAULT_RUN_STYLES: dict[str, dict] = {
    "regression": {"color": "red", "label": None},
    "regression_without_diffusion": {"color": "blue", "label": None},
    "fhofer": {"color": "purple", "label": "FHO-FR fhovv"},
    "fhofer_scaled": {"color": "tab:orange", "label": "FHO-FR fhovv ×{fortran_scale:.1g}"},
}

# цвета для прогонов без явного run_styles
_AUTO_RUN_COLORS = ("red", "blue", "purple", "tab:orange", "green", "brown", "teal", "crimson")


def xlim_by_t_max_us(
    t_max_us: Sequence[float] | dict[int, float],
    *,
    t_min_us: float = 0.0,
) -> dict[int, tuple[float, float]]:
    """Пределы оси X (мкс) по уровням v: (0, t_max) для каждого v."""
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
    """Подпись: имя конфига + метод диффузии."""
    title = (run.config.name or run.name.split("::")[0]).strip()
    m = method if method is not None else run.config.diffusion_method
    return f"{title} ({METHOD_LABELS[m]})"


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

    # True: для каждого прогона — только config.diffusion_method
    diffusion_per_run: bool = True
    # Используется при diffusion_per_run=False
    diffusion_filter: DiffusionFilter = "all"
    max_v_plots: int = 6
    ncols: int = 3
    figsize: tuple[float, float] = (16, 10)

    plot_styles: dict[str, dict] = field(
        default_factory=lambda: {"solid": {"color": "red"}}
    )
    run_styles: dict[str, dict] = field(default_factory=lambda: dict(DEFAULT_RUN_STYLES))
    xlim_by_v: dict[int, tuple[float, float]] = field(
        default_factory=lambda: dict(DEFAULT_XLIM_BY_V)
    )
    # Удобная альтернатива xlim_by_v: макс. время на графике (мкс) для v=0,1,…
    # Пример: t_max_us_by_v=(10, 10, 6, 6, 6, 5) или {0: 10, 4: 6}
    t_max_us_by_v: Sequence[float] | dict[int, float] | None = None

    legend_by_run: dict[tuple[str, str, str], str] = field(default_factory=dict)
    legend_exp_dots: str = "эксперимент (exp760v{v}.csv)"
    legend_exp_fho: str = "статья FHO-FR (fho760v{v}.csv)"
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
    """Внутренние ключи прогонов с готовыми результатами."""
    return tuple(r.name for r in _completed_runs(batch, runs))


def completed_run_labels(
    batch: ModelingBatch,
    runs: tuple[str, ...] | None = None,
) -> tuple[str, ...]:
    """Подписи для легенды/лога: config.name + метод диффузии."""
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
            ax.scatter(
                t_mod,
                pop_mod,
                marker="o",
                s=30,
                color="green",
                alpha=0.7,
                label=lbl,
                zorder=4,
            )
            ax.plot(t_mod, pop_mod, color="green", alpha=0.5, linewidth=1, zorder=3)

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
    Сетка 2×3 (или nrows×ncols): доля населённости по v.

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

        ax_v.set_xlabel("Время, мкс")
        ax_v.set_ylabel("Доля населённости")
        ax_v.set_title(f"v = {v}")
        ax_v.grid(True, alpha=0.3)
        if v in plot_cfg.xlim_by_v:
            ax_v.set_xlim(*plot_cfg.xlim_by_v[v])
        ax_v.legend(fontsize=6, loc="best")

    fig.tight_layout()
    return fig, axes_flat


def print_population_summary(
    batch: ModelingBatch,
    run_name: str | None = None,
    *,
    diffusion_method: DiffusionFilter | None = None,
    style: str = "solid",
) -> None:
    """Табличный вывод долей для выбранного прогона и метода диффузии."""
    run = batch[run_name] if run_name else batch.active
    method = diffusion_method if diffusion_method else run.config.diffusion_method
    if method == "all":
        method = "none"

    sims = run.simulations
    if style not in sims or method not in sims[style]:
        raise KeyError(f"Нет результатов для style={style!r}, method={method!r}")

    sol = sims[style][method]["sol"]
    frac = sims[style][method]["fractional"]
    mult = run.config.rate_multipliers.get(style, run.config.rate_multiplier)

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
    """Сводка по всем завершённым прогонам (те же, что на графике при runs=None)."""
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
