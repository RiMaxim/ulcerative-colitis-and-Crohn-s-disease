"""Генотип × омиксы: связаны ли варианты FUT2, NOD2, ATG16L1 с микробиомом, КЦЖК и метаболомом кала.

Гипотезы заданы заранее (литература):
  FUT2 non-secretor (ALT/ALT)  — нет фукозилированных гликанов в слизи → меньше Bifidobacterium,
                                  ниже разнообразие;
  NOD2 (≥1 аллель риска)       — хуже контроль бактерий у стенки → сдвиг состава, меньше продуцентов бутирата;
  ATG16L1 T300A (0/1/2)        — нарушение аутофагии в клетках Панета.
Модель: признак ~ генотип + диагноз (CD/UC/контроль) + возраст + пол (МНК, признак стандартизован).
16S агрегирован до родов (меньше тестов, устойчивее), CLR. Поправка BH внутри «генотип × блок».
Состав в целом (PERMANOVA, расстояние Эйчисона) считается в R: genotype_omics.R.

Запуск: python genotype_omics.py -> results/genotype_omics.csv, results/genotype_omics_input.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd

from data import load_all, load_omics
from associations import bh, fit_one, transform_block

OUT = Path(__file__).resolve().parent / "results"


def genus_table(rel):
    g = rel.T.groupby(rel.columns.str.split().str[0].str.strip("[]")).sum().T
    return g


def clr(X, min_prev=0.2):
    X = X.loc[:, (X > 0).mean() >= min_prev]
    floor = X.where(X > 0).min() / 2
    L = np.log(X.where(X > 0, floor, axis=1))
    return L.sub(L.mean(axis=1), axis=0)


def main():
    clin, blocks = load_all()
    wes = blocks["Генетика_WES"]
    geno = pd.DataFrame({
        "FUT2_nonsecretor": wes.FUT2_nonsecretor,
        "NOD2_carrier": (wes.NOD2_risk_alleles >= 1).astype(float).where(wes.NOD2_risk_alleles.notna()),
        "ATG16L1_T300A": wes.ATG16L1_T300A,
    })
    meta = clin[["diagnosis", "age", "sex_male"]].join(geno, how="inner")
    cov = pd.get_dummies(meta.diagnosis, drop_first=True, dtype=float).join(meta[["age", "sex_male"]])

    rel = load_omics("16S.xlsx").fillna(0)
    gen = genus_table(rel)
    p = rel.where(rel > 0)
    alpha = pd.DataFrame({
        "Shannon": -(rel * np.log(p)).sum(axis=1),
        "Число видов": (rel > 0).sum(axis=1),
        "Bifidobacterium (CLR)": None,
    })
    G = clr(gen)
    alpha["Bifidobacterium (CLR)"] = G["Bifidobacterium"] if "Bifidobacterium" in G else np.nan

    tests = {
        "16S: разнообразие": alpha,
        "16S: роды": G,
        "КЦЖК_кал": transform_block(blocks["КЦЖК_кал"], "КЦЖК_кал"),
        "Метаболом_кал": transform_block(blocks["Метаболом_кал"], "Метаболом_кал"),
    }
    rows = []
    for g in geno.columns:
        for tname, T in tests.items():
            ids = T.index.intersection(meta.index)
            for f in T.columns:
                x = T.loc[ids, f].astype(float)
                gg = meta.loc[ids, g]
                if gg.nunique() < 2 or (gg > 0).sum() < 4:
                    continue
                b, pv, n = fit_one(x, gg, cov.loc[ids])
                rows.append(dict(genotype=g, block=tname, feature=f, n=n,
                                 carriers=int((gg.loc[x.notna() & gg.notna()] > 0).sum()),
                                 effect_sd=b, p=pv))
    df = pd.DataFrame(rows)
    df["q"] = df.groupby(["genotype", "block"]).p.transform(lambda s: pd.Series(bh(s.values), index=s.index))
    df = df.sort_values("p")
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "genotype_omics.csv", index=False)

    # вход для PERMANOVA в R: CLR родов + метаданные
    ids = G.index.intersection(meta.index)
    G.loc[ids].join(meta.loc[ids]).join(alpha.loc[ids, ["Shannon"]]).to_csv(OUT / "genotype_omics_input.csv")

    with pd.option_context("display.width", 200, "display.max_rows", 100):
        print("Носители (у кого есть 16S):", {g: int((meta.loc[ids, g] > 0).sum()) for g in geno.columns}, "n =", len(ids))
        print(df[df.block == "16S: разнообразие"].round(3).to_string(index=False))
        print("\nЛучшие 25 по p:")
        print(df.head(25).round(4).to_string(index=False))
        print("\nq < 0,1:", int((df.q < 0.1).sum()))


if __name__ == "__main__":
    main()
