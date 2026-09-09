"""
Экспериментальные данные.

Использование:
    from experiments import get_experiment
    exp = get_experiment("N2", "N2", 760)
"""
from experiments.registry import get_experiment, experiment_exists, REGISTRY

__all__ = ["get_experiment", "experiment_exists", "REGISTRY"]
