"""Прирост AUC от каждого блока данных к базовой модели M0 для каждой конечной точки.

Запуск:  python run_increment.py            (полный расчёт: 20 повторов CV, 100 перестановок)
         python run_increment.py --quick    (быстрая проверка: 3 повтора, 10 перестановок)
Результат: results/increment_results.csv и results/increment_results.json
"""
import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from data import load_all, load_drugs, COMPOSITIONAL
from model import evaluate_increment

OUT = Path(__file__).resolve().parent / "results"

# Базовые признаки M0: то, что есть у врача без омиксов (клиника + рутинные анализы).
LABS = ["log_calprotectin", "log_crp", "rbc", "hb", "wbc", "lbp"]
CLINIC_IBD = ["age", "sex_male", "is_CD", "extraintestinal"]

MIN_N, MIN_CLASS = 30, 8  # меньше — оценка прироста бессмысленна, пару пропускаем


def endpoints(clin):
    """Конечные точки: (исход 0/1, базовые признаки, блоки-исключения, описание)."""
    ibd = clin.diagnosis.isin(["CD", "UC"])
    e = {}
    e["CD_vs_UC"] = (
        (clin.diagnosis == "CD").astype(int)[ibd],
        ["age", "sex_male", "extraintestinal"] + LABS, [],
        "Дифференциальный диагноз: болезнь Крона (1) против язвенного колита (0)")
    e["IBD_vs_Control"] = (
        ibd.astype(int),
        ["age", "sex_male"] + LABS, [],
        "Риск: ВЗК (1) против контроля (0)")
    e["Exacerbation"] = (
        clin.exacerbation[ibd].dropna().astype(int),
        CLINIC_IBD + LABS, [],
        "Клиническая активность у ВЗК: обострение (1) против ремиссии (0)")
    act = clin.endo_activity[ibd].dropna()
    e["Endoscopic_A2plus"] = (
        (act >= 2).astype(int),
        CLINIC_IBD + LABS, [],
        "Эндоскопическая активность у ВЗК: A2–A3 (1) против A0–A1 (0)")
    atr = clin.ileal_atrophy_score[ibd].dropna()
    e["Ileal_atrophy"] = (
        (atr >= 1).astype(int),
        CLINIC_IBD + LABS, ["Гистология"],  # исход взят из гистологии — блок исключён
        "Атрофия слизистой подвздошной кишки у ВЗК (гистология ТОПК ≥1)")
    return e


def prepare_clin(clin):
    c = clin.copy()
    c["log_calprotectin"] = np.log1p(c.calprotectin)
    c["log_crp"] = np.log1p(c.crp)
    c["is_CD"] = (c.diagnosis == "CD").astype(float)
    return c


def build_jobs(clin, blocks):
    jobs, skipped = [], []
    for ep, (y, base_cols, excluded, _) in endpoints(clin).items():
        for bname, B in blocks.items():
            if bname in excluded:
                skipped.append((ep, bname, "блок является источником исхода"))
                continue
            ids = y.index.intersection(B.index)
            yy = y.loc[ids]
            n_min = int(min(yy.sum(), len(yy) - yy.sum()))
            if len(ids) < MIN_N or n_min < MIN_CLASS:
                skipped.append((ep, bname, f"мало данных: n={len(ids)}, меньший класс={n_min}"))
                continue
            jobs.append(dict(endpoint=ep, block_name=bname,
                             Xbase=clin.loc[ids, base_cols].values,
                             Xblock=B.loc[ids].values, y=yy.values,
                             compositional=bname in COMPOSITIONAL))
    return jobs, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=-1)
    a = ap.parse_args()
    cfg = dict(n_repeats=3, n_perm=10, perm_repeats=1, n_boot=300) if a.quick else \
        dict(n_repeats=20, n_perm=100, perm_repeats=3, n_boot=2000)

    clin, blocks = load_all()
    clin = prepare_clin(clin)
    jobs, skipped = build_jobs(clin, blocks)
    print(f"Пар «конечная точка × блок»: {len(jobs)}, пропущено: {len(skipped)}")
    for s in skipped:
        print("  пропуск:", *s)

    res = Parallel(n_jobs=a.jobs, verbose=5)(
        delayed(evaluate_increment)(**j, **cfg, seed=0) for j in jobs)

    OUT.mkdir(exist_ok=True)
    rows = [asdict(r) for r in res]
    meta = {ep: d for ep, (_, _, _, d) in endpoints(clin).items()}
    json.dump({"config": cfg, "endpoints": meta, "skipped": skipped, "results": rows},
              open(OUT / "increment_results.json", "w"), ensure_ascii=False, indent=1)
    df = pd.DataFrame(rows).drop(columns=["delta_repeats", "null_deltas"])
    for c in ["auc_m0_ci", "auc_m1_ci", "delta_ci"]:
        df[c] = df[c].apply(lambda t: f"[{t[0]:.3f}; {t[1]:.3f}]")
    df.to_csv(OUT / "increment_results.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(df[["endpoint", "block", "n", "n_pos", "n_neg", "auc_m0", "auc_m1", "delta",
                  "delta_ci", "p_perm"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
