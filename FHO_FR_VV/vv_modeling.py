"""
Моделирование V–V кинетики с диффузией.

Конфигурация — ModelingConfig; пакет прогонов — ModelingBatch; расчёт — run_simulation_suite.
"""
from __future__ import annotations

import copy
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from collections.abc import Sequence
from typing import Literal

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from k_VV import k_vv_fhofer

# (t_µs, population, yerr, xerr) — xerr/yerr None или скаляр/массив
ExpDotTuple = tuple[
    Sequence[float] | np.ndarray,
    Sequence[float] | np.ndarray,
    Sequence[float] | np.ndarray | float | None,
    Sequence[float] | np.ndarray | float | None,
]

KineticsSource = Literal["regression", "fhofer"]
DiffusionMethod = Literal["none", "analytic_factor", "method_of_lines"]
MethodsArg = tuple[DiffusionMethod, ...] | Literal["auto"]

DEFAULT_EXP_YE_BY_V: dict[int, float] = {
    0: 0.010025380710659837,
    1: 0.010173160173160167,
    2: 0.004945770065075926,
    3: 0.0014414784394250514,
    4: 0.0015192743764172335,
    5: 0.0014965986394557824,
}


# --- реестр частиц ---

_PARTICLE_REGISTRY: dict[str, object] | None = None


def _resolve_particles(particle_a: str, particle_b: str):
    """Разрешение имён частиц в объекты Particle."""
    global _PARTICLE_REGISTRY
    if _PARTICLE_REGISTRY is None:
        from particles import N2, O2, CO
        _PARTICLE_REGISTRY = {"N2": N2, "O2": O2, "CO": CO}
    pa = _PARTICLE_REGISTRY.get(particle_a)
    pb = _PARTICLE_REGISTRY.get(particle_b)
    if pa is None:
        raise ValueError(f"Неизвестная частица: {particle_a!r}")
    if pb is None:
        raise ValueError(f"Неизвестная частица: {particle_b!r}")
    return pa, pb


# --- коэффициенты регрессии (из k_VV/k_VV_regr.py) ---
from k_VV import k_lookup, load_coeffs, load_coeffs_from_path, coeffs_dict_from_regression_csv


# --- V–V правая часть ---

def _k_transition(
    kinetics: KineticsSource,
    coeffs: dict,
    i1: int,
    f1: int,
    i2: int,
    f2: int,
    T: float,
    rate_multiplier: float,
    m1,
    m2,
) -> float:
    if kinetics == "regression":
        return rate_multiplier * k_lookup(coeffs, i1, f1, i2, f2, T)
    return rate_multiplier * k_vv_fhofer(m1, m2, i1=i1, f1=f1, i2=i2, f2=f2, T=T)


def R_VV(
    v: int,
    N: np.ndarray,
    v_max: int,
    T: float,
    *,
    kinetics: KineticsSource = "regression",
    coeffs: dict | None = None,
    rate_multiplier: float = 1.0,
    m1=None,
    m2=None,
) -> float:
    R = 0.0
    for v_ in range(v_max + 1):
        k1 = k2 = k3 = k4 = 0.0
        if v + 1 <= v_max and v_ + 1 <= v_max:
            k1 = _k_transition(
                kinetics, coeffs, v_, v_ + 1, v + 1, v, T, rate_multiplier, m1, m2
            )
        if v - 1 >= 0 and v_ - 1 >= 0:
            k2 = _k_transition(
                kinetics, coeffs, v_, v_ - 1, v - 1, v, T, rate_multiplier, m1, m2
            )
        if v - 1 >= 0 and v_ + 1 <= v_max:
            k3 = _k_transition(
                kinetics, coeffs, v_, v_ + 1, v, v - 1, T, rate_multiplier, m1, m2
            )
        if v + 1 <= v_max and v_ - 1 >= 0:
            k4 = _k_transition(
                kinetics, coeffs, v_, v_ - 1, v, v + 1, T, rate_multiplier, m1, m2
            )
        Nvp1 = N[v + 1] if v + 1 <= v_max else 0.0
        Nvm1 = N[v - 1] if v - 1 >= 0 else 0.0
        R += (k1 * Nvp1 + k2 * Nvm1 - (k3 + k4) * N[v]) * N[v_]
    return R


@dataclass
class VVRateMatrices:
    """Предвычисленные k(V–V) при фиксированных T и kinetics (nv × nv)."""

    k1: np.ndarray
    k2: np.ndarray
    k3: np.ndarray
    k4: np.ndarray


def build_vv_rate_matrices(
    v_max: int,
    *,
    kinetics: KineticsSource,
    coeffs: dict,
    T: float,
    rate_multiplier: float,
    m1=None,
    m2=None,
) -> VVRateMatrices:
    """Один раз перед method_of_lines — в RHS только умножения."""
    nv = v_max + 1
    k1 = np.zeros((nv, nv))
    k2 = np.zeros((nv, nv))
    k3 = np.zeros((nv, nv))
    k4 = np.zeros((nv, nv))
    for v in range(nv):
        for v_ in range(nv):
            if v + 1 <= v_max and v_ + 1 <= v_max:
                k1[v, v_] = _k_transition(
                    kinetics, coeffs, v_, v_ + 1, v + 1, v, T, rate_multiplier, m1, m2
                )
            if v - 1 >= 0 and v_ - 1 >= 0:
                k2[v, v_] = _k_transition(
                    kinetics, coeffs, v_, v_ - 1, v - 1, v, T, rate_multiplier, m1, m2
                )
            if v - 1 >= 0 and v_ + 1 <= v_max:
                k3[v, v_] = _k_transition(
                    kinetics, coeffs, v_, v_ + 1, v, v - 1, T, rate_multiplier, m1, m2
                )
            if v + 1 <= v_max and v_ - 1 >= 0:
                k4[v, v_] = _k_transition(
                    kinetics, coeffs, v_, v_ - 1, v, v + 1, T, rate_multiplier, m1, m2
                )
    return VVRateMatrices(k1=k1, k2=k2, k3=k3, k4=k4)


def R_VV_from_matrices(N: np.ndarray, rates: VVRateMatrices) -> np.ndarray:
    """
    Тот же R_VV, что цикл по v_, но без повторных k_lookup.
    N: (nv,) или (nv, n_r) → (nv,) или (nv, n_r).
    """
    if N.ndim == 1:
        return _R_VV_column_fast(N, rates)
    return _R_VV_grid_fast(N, rates)


def _R_VV_column_fast(N: np.ndarray, rates: VVRateMatrices) -> np.ndarray:
    nv = N.shape[0]
    Nvp1 = np.zeros(nv)
    Nvm1 = np.zeros(nv)
    if nv > 1:
        Nvp1[:-1] = N[1:]
        Nvm1[1:] = N[:-1]
    k34 = rates.k3 + rates.k4
    R = np.zeros(nv)
    for v_ in range(nv):
        coef = (
            rates.k1[:, v_] * Nvp1
            + rates.k2[:, v_] * Nvm1
            - k34[:, v_] * N
        )
        R += coef * N[v_]
    return R


def _R_VV_grid_fast(N: np.ndarray, rates: VVRateMatrices) -> np.ndarray:
    nv, _ = N.shape
    Nvp1 = np.zeros_like(N)
    Nvm1 = np.zeros_like(N)
    if nv > 1:
        Nvp1[:-1, :] = N[1:, :]
        Nvm1[1:, :] = N[:-1, :]
    k34 = rates.k3 + rates.k4
    R = np.zeros_like(N)
    for v_ in range(nv):
        coef = (
            rates.k1[:, v_][:, np.newaxis] * Nvp1
            + rates.k2[:, v_][:, np.newaxis] * Nvm1
            - k34[:, v_][:, np.newaxis] * N
        )
        R += coef * N[v_, np.newaxis, :]
    return R


# --- загрузка данных ---

def read_data_csv(base_dir: Path, filename: str, sep: str = ";") -> pd.DataFrame:
    return pd.read_csv(base_dir / filename, sep=sep, header=None)


def load_experimental_data(
    base_dir: Path,
    v_max: int,
    *,
    prefix: str = "exp",
    ye_by_v: dict[int, float] | None = None,
) -> dict[str, list]:
    exp_data: dict[str, list] = {"dots": []}
    ye_map = ye_by_v if ye_by_v is not None else DEFAULT_EXP_YE_BY_V

    for v in range(v_max + 1):
        try:
            df = read_data_csv(base_dir, f"{prefix}v{v}.csv")
            x_data = df.iloc[:, 0].values.astype(float)
            y_data = df.iloc[:, 1].values.astype(float)
            ye_data = (
                np.full_like(y_data, ye_map[v])
                if v in ye_map
                else y_data * 0.02
            )
            exp_data["dots"].append((x_data, y_data, ye_data, None))
        except FileNotFoundError:
            print(f"Предупреждение: {prefix}v{v}.csv не найден")
            exp_data["dots"].append(None)
    return exp_data


def _as_float_array(x: Sequence[float] | np.ndarray) -> np.ndarray:
    return np.asarray(x, dtype=float).ravel()


def load_experimental_data_inline(
    v_max: int,
    dots_by_v: Sequence[ExpDotTuple | None],
) -> dict[str, list]:
    """Точки эксперимента из ноутбука: список (x, y, ye, xe) по v."""
    dots: list = []
    for v in range(v_max + 1):
        if v >= len(dots_by_v) or dots_by_v[v] is None:
            dots.append(None)
            continue
        xv, yv, yev, xev = dots_by_v[v]
        ye_arr = _as_float_array(yev) if yev is not None else yv * 0.02
        xe_arr = None if xev is None else _as_float_array(xev)
        dots.append((_as_float_array(xv), _as_float_array(yv), ye_arr, xe_arr))
    return {"dots": dots}


def load_experimental_data_csv_styles(
    base_dir: Path,
    v_max: int,
    styles: Sequence[str],
    *,
    filename_fmt: str = "v{v}{style}.csv",
) -> dict[str, list]:
    """CSV solid/dashed: v0solid.csv, v1dashed.csv, …"""
    out: dict[str, list] = {style: [] for style in styles}
    for style in styles:
        for v in range(v_max + 1):
            fname = filename_fmt.format(v=v, style=style)
            try:
                out[style].append(read_data_csv(base_dir, fname))
            except FileNotFoundError:
                print(f"Предупреждение: {fname} не найден")
                out[style].append(None)
    return out


def load_modeling_data(
    base_dir: Path,
    v_max: int,
    *,
    prefix: str = "fho",
) -> dict[str, list]:
    exp_data: dict[str, list] = {"fho": []}
    for v in range(v_max + 1):
        try:
            exp_data["fho"].append(read_data_csv(base_dir, f"{prefix}v{v}.csv"))
        except FileNotFoundError:
            print(f"Предупреждение: {prefix}v{v}.csv не найден")
            exp_data["fho"].append(None)
    return exp_data


# --- пространственная диффузия (method of lines) ---

def laplacian_cyl(f: np.ndarray, r: np.ndarray, dr: float) -> np.ndarray:
    """Лапласиан в цилиндре; f — 1D (n_r) или 2D (..., n_r) по последней оси."""
    if f.ndim == 1:
        return _laplacian_cyl_1d(f, r, dr)
    return _laplacian_cyl_nd(f, r, dr)


def _laplacian_cyl_1d(f: np.ndarray, r: np.ndarray, dr: float) -> np.ndarray:
    lap = np.zeros_like(f)
    ri = r[1:-1]
    r_plus = ri + 0.5 * dr
    r_minus = ri - 0.5 * dr
    lap[1:-1] = (1.0 / ri) * (
        (r_plus * (f[2:] - f[1:-1]) - r_minus * (f[1:-1] - f[:-2])) / dr**2
    )
    lap[0] = 4.0 * (f[1] - f[0]) / dr**2
    lap[-1] = 0.0
    return lap


def _laplacian_cyl_nd(f: np.ndarray, r: np.ndarray, dr: float) -> np.ndarray:
    lap = np.zeros_like(f)
    ri = r[1:-1]
    r_plus = ri + 0.5 * dr
    r_minus = ri - 0.5 * dr
    lap[..., 1:-1] = (1.0 / ri) * (
        (
            r_plus * (f[..., 2:] - f[..., 1:-1])
            - r_minus * (f[..., 1:-1] - f[..., :-2])
        )
        / dr**2
    )
    lap[..., 0] = 4.0 * (f[..., 1] - f[..., 0]) / dr**2
    lap[..., -1] = 0.0
    return lap


class _MethodOfLinesRHS:
    """Кэшированный RHS для solve_ivp (method of lines)."""

    __slots__ = ("nv", "n_r", "r", "dr", "D", "rates", "rate_multiplier")

    def __init__(
        self,
        cfg: "ModelingConfig",
        rate_multiplier: float,
        n_r: int,
        r: np.ndarray,
        dr: float,
    ) -> None:
        self.nv = cfg.v_max + 1
        self.n_r = n_r
        self.r = r
        self.dr = dr
        self.D = cfg.D
        self.rate_multiplier = rate_multiplier
        m1, m2 = _resolve_particles(cfg.particle_a, cfg.particle_b)
        self.rates = build_vv_rate_matrices(
            cfg.v_max,
            kinetics=cfg.kinetics,
            coeffs=cfg.coeffs,
            T=cfg.T,
            rate_multiplier=rate_multiplier,
            m1=m1,
            m2=m2,
        )

    def __call__(self, t: float, state_flat: np.ndarray) -> np.ndarray:
        state = state_flat.reshape(self.nv, self.n_r)
        dstate = np.zeros_like(state)
        lap_v = laplacian_cyl(state[1:, :], self.r, self.dr)
        sum_lapl = lap_v.sum(axis=0)
        dstate[:] = R_VV_from_matrices(state, self.rates)
        dstate[1:, :] += self.D * lap_v
        dstate[0, :] -= self.D * sum_lapl
        return dstate.ravel()


def spatial_average(state: np.ndarray, r: np.ndarray, sigma_probe: float) -> np.ndarray:
    w = r * np.exp(-2.0 * r**2 / sigma_probe**2)
    w_norm = w / np.trapezoid(w, r)
    return np.trapezoid(state * w_norm[np.newaxis, :], r, axis=1)


# --- конфигурация ---

@dataclass
class ModelingConfig:
    """Параметры моделирования V–V кинетики с диффузией."""

    name: str = ""
    base_dir: Path = field(default_factory=lambda: Path(".").resolve())
    particle_a: str = ""
    particle_b: str = ""

    # кинетика
    kinetics: KineticsSource = "regression"
    rate_multiplier: float = 1.1
    regression_csv: str | None = None  # путь к CSV с коэффициентами (или имя файла)

    # физика
    T: float = 300.0
    v_max: int = 6
    n_total: float = 2.45e19
    f_v0: float = 0.63
    f_v1: float = 0.37

    # интегрирование
    t_end_s: float = 12e-6
    n_time: int = 600

    # диффузия
    diffusion_method: DiffusionMethod = "none"
    sigma_pump_cm: float = 46e-4
    sigma_probe_cm: float = 46e-4
    d_ref: float = 0.25
    p_torr: float = 760.0

    # сетка method of lines
    r_max_cm: float = 200e-4
    n_r: int = 80

    experiment_key: tuple[str, str, float] | None = None  # (m1, m2, p_torr) → registry

    exp_ye_by_v: dict[int, float] = field(default_factory=lambda: dict(DEFAULT_EXP_YE_BY_V))
    exp_dots_inline: Sequence[ExpDotTuple | None] | None = None
    exp_csv_styles: tuple[str, ...] | None = None

    _coeffs: dict | None = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        self.base_dir = Path(self.base_dir).resolve()

    @property
    def coeffs(self) -> dict:
        if self._coeffs is None:
            if self.regression_csv is not None:
                self._coeffs = load_coeffs_from_path(self.regression_csv)
            else:
                self._coeffs = load_coeffs(self.base_dir, self.particle_a, self.particle_b)
        return self._coeffs

    @property
    def initial_fractions(self) -> np.ndarray:
        fr = np.zeros(self.v_max + 1)
        fr[0] = self.f_v0
        fr[1] = self.f_v1
        return fr

    @property
    def t_span(self) -> tuple[float, float]:
        return (0.0, self.t_end_s)

    @property
    def t_eval(self) -> np.ndarray:
        return np.linspace(0.0, self.t_end_s, self.n_time)

    @property
    def D(self) -> float:
        return self.d_ref * (760.0 / self.p_torr)

    def build_initial_state_method_of_lines(self, n_r: int | None = None, r_max: float | None = None):
        n_r = n_r if n_r is not None else self.n_r
        r_max = r_max if r_max is not None else self.r_max_cm
        dr = r_max / (n_r - 1)
        r = np.linspace(0.0, r_max, n_r)
        profile = np.exp(-2.0 * r**2 / self.sigma_pump_cm**2)
        w_probe = r * np.exp(-2.0 * r**2 / self.sigma_probe_cm**2)
        integral_total = np.trapezoid(self.n_total * w_probe, r)
        integral_profile = np.trapezoid(profile * w_probe, r)
        A = (self.initial_fractions[1] * integral_total) / integral_profile
        n1_init = A * profile
        n0_init = self.n_total - n1_init
        state = np.zeros((self.v_max + 1, n_r))
        state[0, :] = n0_init
        state[1, :] = n1_init
        return state, r, dr

    def rhs_time(self, t: float, N: np.ndarray, rate_multiplier: float) -> np.ndarray:
        dNdt = np.zeros_like(N)
        m1, m2 = _resolve_particles(self.particle_a, self.particle_b)
        for v in range(self.v_max + 1):
            dNdt[v] = R_VV(
                v,
                N,
                self.v_max,
                self.T,
                kinetics=self.kinetics,
                coeffs=self.coeffs,
                rate_multiplier=rate_multiplier,
                m1=m1,
                m2=m2,
            )
        return dNdt

    def solve(
        self,
        rate_multiplier: float,
        method: DiffusionMethod = "none",
    ) -> tuple[object, np.ndarray]:
        if method in ("none", "analytic_factor"):
            N_init = self.initial_fractions * self.n_total
            sol = solve_ivp(
                lambda t, N: self.rhs_time(t, N, rate_multiplier),
                self.t_span,
                N_init,
                t_eval=self.t_eval,
                method="Radau",
            )
            if not sol.success:
                raise RuntimeError(f"solve_ivp failed: {sol.message}")
            if method == "analytic_factor":
                R = 1.0 / np.sqrt(
                    1.0
                    + 4.0
                    * self.D
                    * sol.t
                    / (self.sigma_pump_cm**2 + self.sigma_probe_cm**2)
                )
                fractional = np.zeros_like(sol.y)
                fractional[1:, :] = (
                    sol.y[1:, :] * R[np.newaxis, :] / self.n_total
                )
                fractional[0, :] = 1.0 - fractional[1:, :].sum(axis=0)
            else:
                fractional = sol.y / np.sum(sol.y, axis=0)
            return sol, fractional

        if method == "method_of_lines":
            state_init, r, dr = self.build_initial_state_method_of_lines()
            n_r = len(r)
            rhs = _MethodOfLinesRHS(self, rate_multiplier, n_r, r, dr)
            sol = solve_ivp(
                rhs,
                self.t_span,
                state_init.flatten(),
                t_eval=self.t_eval,
                method="Radau",
                rtol=1e-9,
                atol=1e6,
            )
            if not sol.success:
                raise RuntimeError(f"solve_ivp failed: {sol.message}")
            n_t_avg = np.zeros((self.v_max + 1, len(sol.t)))
            for idx in range(len(sol.t)):
                state_t = sol.y[:, idx].reshape((self.v_max + 1, n_r))
                n_t_avg[:, idx] = spatial_average(state_t, r, self.sigma_probe_cm)
            return sol, n_t_avg / self.n_total

        raise ValueError(f"Неизвестный метод диффузии: {method!r}")

    def _get_registry_entry(self) -> dict | None:
        """Возвращает запись из experiments/registry для experiment_key."""
        if self.experiment_key is None:
            return None
        from experiments.registry import get_experiment
        m1, m2, p = self.experiment_key
        return get_experiment(m1, m2, p)

    def load_experimental_data(self) -> dict:
        out: dict = {}

        # 1) Если задан experiment_key — загружаем из registry
        reg = self._get_registry_entry()
        if reg is not None:
            # inline точки из Python-модуля
            if "dots_inline_module" in reg:
                import importlib
                mod = importlib.import_module(reg["dots_inline_module"])
                dots = getattr(mod, reg["dots_inline_name"])
                out.update(load_experimental_data_inline(self.v_max, dots))
            # CSV стили (solid, dashed, …)
            if "styles" in reg:
                styles = reg["styles"]
                if "dots_inline_csv_dir" in reg:
                    csv_dir = Path(reg["dots_inline_csv_dir"])
                elif "dots_inline_module" in reg:
                    from experiments.registry import BASE_DIR
                    parts = reg["dots_inline_module"].split(".")
                    csv_dir = BASE_DIR / "/".join(parts[1:-1])
                else:
                    csv_dir = self.base_dir
                out.update(
                    load_experimental_data_csv_styles(
                        csv_dir, self.v_max, styles,
                        filename_fmt=reg.get("dots_inline_csv_fmt", "v{v}{style}.csv"),
                    )
                )
            # CSV dots по prefix
            if "dots_prefix" in reg:
                prefix = reg["dots_prefix"]
                out.update(
                    load_experimental_data(self.base_dir, self.v_max, prefix=prefix)
                )

        # 2) Fallback: legacy поля
        if not out:
            if self.exp_dots_inline is not None:
                out.update(
                    load_experimental_data_inline(self.v_max, self.exp_dots_inline)
                )
            if self.exp_csv_styles:
                out.update(
                    load_experimental_data_csv_styles(
                        self.base_dir, self.v_max, self.exp_csv_styles
                    )
                )
        if out:
            return out
        return load_experimental_data(
            self.base_dir, self.v_max, ye_by_v=self.exp_ye_by_v
        )

    def load_modeling_data(self) -> dict:
        # Если задан experiment_key и есть fho_prefix — используем registry
        reg = self._get_registry_entry()
        if reg is not None and "fho_prefix" in reg:
            prefix = reg["fho_prefix"]
            return load_modeling_data(self.base_dir, self.v_max, prefix=prefix)
        return load_modeling_data(self.base_dir, self.v_max)

    def with_kinetics(self, kinetics: KineticsSource) -> "ModelingConfig":
        """Копия конфига с другим источником k_VV."""
        cfg = copy.copy(self)
        cfg.kinetics = kinetics
        cfg._coeffs = self._coeffs
        return cfg


def run_simulation_suite(
    cfg: ModelingConfig,
    *,
    kinetics_sources: tuple[KineticsSource, ...] | None = None,
    methods: tuple[DiffusionMethod, ...] = ("none", "analytic_factor", "method_of_lines"),
    styles: tuple[str, ...] | None = None,
) -> dict[KineticsSource, dict[str, dict[DiffusionMethod, dict]]]:
    """
    Запуск набора симуляций.

    Возвращает results[kinetics][style][method] = {"sol", "fractional"}.
    Если kinetics_sources не задан — используется cfg.kinetics.
    """
    if kinetics_sources is None:
        kinetics_sources = (cfg.kinetics,)
    styles = styles if styles is not None else ("default",)
    out: dict = {}
    for kin in kinetics_sources:
        kin_cfg = cfg.with_kinetics(kin)
        out[kin] = {}
        for style in styles:
            out[kin][style] = {}
            for method in methods:
                sol, frac = kin_cfg.solve(cfg.rate_multiplier, method=method)
                out[kin][style][method] = {"sol": sol, "fractional": frac}
    return out


SuiteResult = dict[KineticsSource, dict[str, dict[DiffusionMethod, dict]]]


@dataclass
class ModelingRun:
    """Один именованный прогон с конфигом и результатами."""

    name: str
    config: ModelingConfig
    suite: SuiteResult | None = None
    exp_data: dict | None = None
    exp_mod_data: dict | None = None

    def run(
        self,
        *,
        kinetics_sources: tuple[KineticsSource, ...] | None = None,
        methods: tuple[DiffusionMethod, ...] = (
            "none",
            "analytic_factor",
            "method_of_lines",
        ),
        styles: tuple[str, ...] | None = None,
        load_data: bool = True,
    ) -> "ModelingRun":
        self.suite = run_simulation_suite(
            self.config,
            kinetics_sources=kinetics_sources,
            methods=methods,
            styles=styles,
        )
        if load_data:
            self.exp_data = self.config.load_experimental_data()
            self.exp_mod_data = self.config.load_modeling_data()
        return self

    @property
    def simulations(self) -> dict[str, dict[DiffusionMethod, dict]]:
        if self.suite is None:
            raise RuntimeError(f"Сначала вызовите batch.run('{self.name}') или run.run()")
        return self.suite[self.config.kinetics]

    def simulations_for(self, kinetics: KineticsSource) -> dict[str, dict[DiffusionMethod, dict]]:
        if self.suite is None:
            raise RuntimeError(f"Сначала вызовите batch.run('{self.name}')")
        return self.suite[kinetics]

    def expose_notebook_globals(
        self,
        kinetics: KineticsSource | None = None,
    ) -> dict[str, object]:
        kin = kinetics or self.config.kinetics
        cfg = self.config
        if self.suite is None:
            raise RuntimeError(f"Сначала вызовите batch.run('{self.name}')")
        return {
            "cfg": cfg,
            "METHOD": cfg.diffusion_method,
            "V_MAX": cfg.v_max,
            "T": cfg.T,
            "N_TOTAL": cfg.n_total,
            "initial_fractions": cfg.initial_fractions,
            "t_span": cfg.t_span,
            "t_eval": cfg.t_eval,
            "sigma_pump_cm": cfg.sigma_pump_cm,
            "sigma_probe_cm": cfg.sigma_probe_cm,
            "D": cfg.D,
            "COEFFS": cfg.coeffs,
            "BASE_DIR": cfg.base_dir,
            "suite": self.suite,
            "simulations": self.suite[kin],
            "exp_data": self.exp_data,
            "exp_mod_data": self.exp_mod_data,
        }


ModelingRunJob = tuple[
    str,
    ModelingConfig,
    tuple[DiffusionMethod, ...],
    tuple[str, ...] | None,
    bool,
    tuple[KineticsSource, ...] | None,
]


def _execute_modeling_run_job(job: ModelingRunJob) -> tuple[str, SuiteResult, dict | None, dict | None]:
    """Воркер для ProcessPoolExecutor (должен быть функцией верхнего уровня)."""
    name, config, methods, styles, load_data, kinetics_sources = job
    run = ModelingRun(name=name, config=config)
    run.run(
        kinetics_sources=kinetics_sources,
        methods=methods,
        styles=styles,
        load_data=load_data,
    )
    return name, run.suite, run.exp_data, run.exp_mod_data


def _make_run_key(cfg: ModelingConfig, index: int, used: set[str]) -> str:
    label = (cfg.name or f"run_{index}").strip()
    key = f"{label}::{cfg.diffusion_method}"
    if key in used:
        key = f"{key}::{index}"
    used.add(key)
    return key


class ModelingBatch:
    """
    Несколько моделирований с разными ModelingConfig.

    Пример::

        batch = ModelingBatch([
            ModelingConfig(name="baseline", diffusion_method="analytic_factor", ...),
            ModelingConfig(name="variant", diffusion_method="none", ...),
        ])
        batch.run_all()
    """

    def __init__(
        self,
        configs: Sequence[ModelingConfig] | dict[str, ModelingConfig] | None = None,
        /,
        **kwargs: ModelingConfig,
    ) -> None:
        self._runs: dict[str, ModelingRun] = {}
        self._active: str | None = None

        if isinstance(configs, dict):
            for name, cfg in configs.items():
                self.add(name, cfg)
            return

        if isinstance(configs, (list, tuple)):
            used: set[str] = set()
            for index, cfg in enumerate(configs):
                self.add(_make_run_key(cfg, index, used), cfg)
            return

        for name, cfg in kwargs.items():
            self.add(name, cfg)

    def add(
        self,
        name: str,
        config: ModelingConfig | None = None,
        /,
        **kwargs: object,
    ) -> ModelingRun:
        if config is None:
            config = ModelingConfig(**kwargs)  # type: ignore[arg-type]
        if not config.name:
            config.name = name.split("::")[0]
        self._runs[name] = ModelingRun(name=name, config=config)
        return self._runs[name]

    def _resolve_methods(
        self,
        name: str,
        methods: MethodsArg,
    ) -> tuple[DiffusionMethod, ...]:
        if methods == "auto":
            return (self._runs[name].config.diffusion_method,)
        return methods

    def run(
        self,
        name: str,
        *,
        kinetics_sources: tuple[KineticsSource, ...] | None = None,
        methods: MethodsArg = "auto",
        styles: tuple[str, ...] | None = None,
        load_data: bool = True,
    ) -> ModelingRun:
        return self._runs[name].run(
            kinetics_sources=kinetics_sources,
            methods=self._resolve_methods(name, methods),
            styles=styles,
            load_data=load_data,
        )

    def run_all(
        self,
        *,
        methods: MethodsArg = "auto",
        styles: tuple[str, ...] | None = None,
        load_data: bool = True,
        only: tuple[str, ...] | None = None,
        parallel: bool = True,
        max_workers: int | None = None,
    ) -> "ModelingBatch":
        names = only if only is not None else tuple(self._runs.keys())
        if not names:
            return self

        if not parallel or len(names) == 1:
            for name in names:
                self.run(name, methods=methods, styles=styles, load_data=load_data)
            return self

        n_workers = max_workers if max_workers is not None else min(
            len(names), os.cpu_count() or 1
        )
        jobs: list[ModelingRunJob] = [
            (
                name,
                self._runs[name].config,
                self._resolve_methods(name, methods),
                styles,
                load_data,
                None,
            )
            for name in names
        ]

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_execute_modeling_run_job, job): job[0] for job in jobs
            }
            for future in as_completed(futures):
                name, suite, exp_data, exp_mod_data = future.result()
                run = self._runs[name]
                run.suite = suite
                run.exp_data = exp_data
                run.exp_mod_data = exp_mod_data
                print(f"  готово: {name}")

        return self

    def set_active(self, name: str) -> "ModelingBatch":
        if name not in self._runs:
            raise KeyError(f"Нет прогона {name!r}; доступны: {self.names()}")
        self._active = name
        return self

    @property
    def active(self) -> ModelingRun:
        if self._active is None:
            if not self._runs:
                raise RuntimeError("Batch пуст — добавьте конфиг через add() или ModelingBatch(...)")
            self._active = next(iter(self._runs))
        return self._runs[self._active]

    def __getitem__(self, name: str) -> ModelingRun:
        return self._runs[name]

    def __contains__(self, name: str) -> bool:
        return name in self._runs

    def names(self) -> tuple[str, ...]:
        return tuple(self._runs.keys())

    def run_at(self, index: int = 0) -> ModelingRun:
        return self[self.names()[index]]

    def expose_notebook_globals(self, name: str | None = None) -> dict[str, object]:
        run = self._runs[name] if name is not None else self.active
        return run.expose_notebook_globals()
