"""Связь отдельных признаков каждого блока с клиническими конечными точками (задача 2).

Для каждой пары «конечная точка × блок» и каждого признака блока:
    признак (стандартизованный) ~ исход (0/1) + возраст + пол [+ тип ВЗК]  — основная модель
    признак ~ исход + возраст + пол [+ тип ВЗК] + терапия на визите 1      — чувствительность
Тип ВЗК (CD/UC) добавляется в задачах, где есть оба диагноза (обострение, эндоскопия, атрофия,
биологическая терапия). В задаче «биологическая/JAK-терапия» сама эта терапия в поправку не входит.
Эффект = разница средних между группами в стандартных отклонениях (как d Коэна, с поправками).
Поправка Бенджамини–Хохберга внутри каждой пары «конечная точка × блок».

Преобразования: доли (16S, метаболомы, липидом) — фильтр распространённости ≥20%, замена нулей
на половину минимума, CLR; КЦЖК — log1p; остальные блоки как есть.

Запуск: python associations.py  ->  results/associations.csv, results/associations_summary.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

from data import load_all, load_drugs, COMPOSITIONAL
from run_increment import endpoints, prepare_clin

OUT = Path(__file__).resolve().parent / "results"
DRUGS = ["препараты 5-АСК", "гормоны (ГКС)", "тиопурины / иммуносупрессоры",
         "генно-инженерная конструкция", "антибиотики"]
SKIP_BLOCKS = {"NOD2"}  # уже входит в блок Генетика_WES


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]; o = np.argsort(pv); m = len(pv)
    qq = pv[o] * m / np.arange(1, m + 1)
    qq = np.minimum.accumulate(qq[::-1])[::-1].clip(max=1)
    tmp = np.empty(m); tmp[o] = qq
    q[ok] = tmp
    return q


def transform_block(B, name, min_prev=0.2):
    B = B.astype(float)
    if name in COMPOSITIONAL:
        X = B.fillna(0.0)
        X = X.loc[:, (X > 0).mean() >= min_prev]
        floor = X.where(X > 0).min() / 2
        L = np.log(X.where(X > 0, floor, axis=1))
        return L.sub(L.mean(axis=1), axis=0)
    if name == "КЦЖК_кал":
        return np.log1p(B.clip(lower=0))
    return B


def fit_one(x, y, C):
    ok = x.notna() & y.notna() & C.notna().all(axis=1)
    x, y, C = x[ok], y[ok], C[ok]
    binary = set(y.unique()) <= {0, 1}
    if ok.sum() < 20 or x.std() == 0 or y.nunique() < 2 or (binary and min(y.sum(), (1 - y).sum()) < 5):
        return np.nan, np.nan, int(ok.sum())
    xs = (x - x.mean()) / x.std()
    Xd = sm.add_constant(pd.concat([y.rename("y"), C], axis=1).astype(float), has_constant="add")
    Xd = Xd.loc[:, Xd.std() > 0].assign(const=1.0)
    r = sm.OLS(xs.values, Xd.values).fit()
    i = list(Xd.columns).index("y")
    return float(r.params[i]), float(r.pvalues[i]), int(ok.sum())


def main():
    clin, blocks = load_all()
    clin = prepare_clin(clin)
    drugs = load_drugs().reindex(clin.index)[DRUGS]
    rows = []
    for ep, (y, _, excluded, _) in endpoints(clin).items():
        for bname, B in blocks.items():
            if bname in excluded or bname in SKIP_BLOCKS:
                continue
            ids = y.index.intersection(B.index)
            if len(ids) < 30 or min(y[ids].sum(), len(ids) - y[ids].sum()) < 8:
                continue
            Z = transform_block(B.loc[ids], bname)
            Z = Z.loc[:, Z.notna().mean() >= 0.6]
            cov_cols = ["age", "sex_male"]
            if ep in ("Exacerbation", "Endoscopic_A2plus", "Ileal_atrophy", "Advanced_therapy"):
                cov_cols.append("is_CD")
            cov = clin.loc[ids, cov_cols]
            dr = drugs.loc[ids]
            if ep == "Advanced_therapy":
                dr = dr.drop(columns=["генно-инженерная конструкция"])
            cov_t = pd.concat([cov, dr], axis=1)
            # для «ВЗК против контроля» поправка на терапию не делается: у контролей терапии нет,
            # и поправка смешала бы её с диагнозом
            for f in Z.columns:
                b, p, n = fit_one(Z[f], y[ids].astype(float), cov)
                bt, pt, _ = fit_one(Z[f], y[ids].astype(float), cov_t) if ep != "IBD_vs_Control" else (np.nan, np.nan, 0)
                rows.append(dict(endpoint=ep, block=bname, feature=f, n=n, effect_sd=b, p=p,
                                 effect_sd_adj_therapy=bt, p_adj_therapy=pt))
    df = pd.DataFrame(rows)
    df["q"] = df.groupby(["endpoint", "block"]).p.transform(lambda s: pd.Series(bh(s.values), index=s.index))
    df["q_adj_therapy"] = df.groupby(["endpoint", "block"]).p_adj_therapy.transform(
        lambda s: pd.Series(bh(s.values), index=s.index))
    df = df.sort_values(["endpoint", "block", "p"])
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "associations.csv", index=False)

    summ = df.groupby(["endpoint", "block"]).agg(
        n=("n", "max"), features=("feature", "size"),
        sig_q10=("q", lambda s: int((s < 0.1).sum())),
        sig_q10_adj_therapy=("q_adj_therapy", lambda s: int((s < 0.1).sum())),
        min_q=("q", "min")).reset_index()
    summ.to_csv(OUT / "associations_summary.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(summ.to_string(index=False))
        top = df[df.q < 0.1].sort_values("p")
        print(f"\nПризнаков с q < 0,1: {len(top)}")
        print(top[["endpoint", "block", "feature", "n", "effect_sd", "p", "q", "q_adj_therapy"]]
              .head(80).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
