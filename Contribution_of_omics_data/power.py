"""Анализ мощности: какой прирост AUC можно было обнаружить при имеющихся выборках.

Используется перестановочное (нулевое) распределение прироста из run_increment.py: это прирост,
который получается при полностью случайном блоке на тех же пациентах с той же процедурой.
Порог значимости — его 95-й перцентиль (односторонний α = 0,05). Минимальный обнаруживаемый прирост
при мощности 80% ≈ порог + 0,84 × SD нулевого распределения (разброс при эффекте принят таким же).
Запуск: python power.py -> results/power.csv
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "results"
d = pd.DataFrame(json.load(open(OUT / "increment_results.json"))["results"])
d["null_q95"] = d.null_deltas.apply(lambda v: float(np.percentile(v, 95)))
d["null_sd"] = d.null_deltas.apply(lambda v: float(np.std(v)))
d["mde"] = d.null_q95 + 0.84 * d.null_sd
d[["endpoint", "block", "n", "n_pos", "n_neg", "delta", "null_q95", "null_sd", "mde"]].to_csv(OUT / "power.csv", index=False)
print(d.groupby("block")[["n", "mde"]].median().sort_values("mde").round(3).to_string())
print("\nПо всем парам: медиана", round(d.mde.median(), 3), "диапазон", round(d.mde.min(), 3), "–", round(d.mde.max(), 3))
print("Доля пар, где +0,05 обнаружим:", round((d.mde <= 0.05).mean(), 2))
