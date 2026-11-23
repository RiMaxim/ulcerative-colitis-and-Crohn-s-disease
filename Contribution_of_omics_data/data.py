"""Загрузка и очистка таблиц когорты ВЗК (только визит 1, одна строка на пациента).

Все исходные xlsx лежат в папке на уровень выше (Desktop/2). Текст "ND" = пропуск.
Функция load_all() возвращает:
  clin    — клиника, рутинные анализы и конечные точки (DataFrame, индекс Patient);
  blocks  — словарь {имя блока: DataFrame признаков, индекс Patient}; в блок входят
            только пациенты, у которых этот блок реально измерен.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
VISIT = 1


def _num(s):
    """Число из ячейки: 'ND' -> NaN, запятая -> точка, ' 1:160' -> 160."""
    s = s.astype(str).str.strip().str.replace(",", ".", regex=False)
    s = s.str.replace(r"^1:", "", regex=True)
    return pd.to_numeric(s, errors="coerce")


def _short(txt):
    """Короткое имя столбца из многострочного заголовка Inventory."""
    return re.sub(r"\s+", " ", str(txt)).strip()


def load_inventory():
    raw = pd.read_excel(ROOT / "Inventory.xlsx", header=None)
    data = raw.iloc[3:].reset_index(drop=True)
    data = data[data[0].notna()]
    data = data[data[1].astype(float) == VISIT].copy()
    data.index = data[0].astype(int)
    data.index.name = "Patient"

    clin = pd.DataFrame(index=data.index)
    clin["diagnosis"] = data[2].astype(str).str.strip()
    clin["sex_male"] = (data[3].astype(str).str.strip() == "M").astype(float)
    clin["age"] = _num(data[4])
    clin["extraintestinal"] = _num(data[5])
    clin["exacerbation"] = _num(data[6])
    endo = data[7].astype(str).str.strip()
    clin["endo_location"] = endo.str.extract(r"^([LE]\d)")[0]
    clin["endo_activity"] = pd.to_numeric(endo.str.extract(r"A(\d)")[0], errors="coerce")

    labs = {8: "calprotectin", 9: "crp", 10: "rbc", 11: "hb", 12: "wbc", 13: "lbp"}
    for c, name in labs.items():
        clin[name] = _num(data[c])

    # Серология: ASCA — отн. ед/мл, ANCA — титр 1:N -> log2(N)
    ser = pd.DataFrame(index=data.index)
    ser["ASCA_IgG"] = np.log1p(_num(data[14]))
    ser["ASCA_IgA"] = np.log1p(_num(data[15]))
    ser["ANCA_IgG_log2titer"] = np.log2(_num(data[16]))
    ser["ANCA_IgA_log2titer"] = np.log2(_num(data[17]))

    # Гистология: ТОПК (подвздошная) — столбцы 18–28, КСК (восходящая ободочная) — 29–39, баллы 0–3
    names = [_short(raw.iloc[1, c]).replace("KCK", "КСК") for c in range(18, 40)]
    hist = pd.DataFrame({n: _num(data[c]) for n, c in zip(names, range(18, 40))}, index=data.index)
    clin["ileal_atrophy_score"] = hist[names[0]]

    # Фекальные КЦЖК и желчные кислоты (FSCA), столбцы 40–59
    fs_names = []
    for c in range(40, 60):
        top = raw.iloc[1, c] if pd.notna(raw.iloc[1, c]) and c not in (48, 52) else raw.iloc[2, c]
        prefix = "abs " if 48 <= c <= 51 else ("rel " if 52 <= c <= 55 else "")
        fs_names.append(prefix + _short(top).split(" рз")[0].split(" 5.")[0].split(" 1.")[0].split(" 0.")[0][:60])
    scfa = pd.DataFrame({n: _num(data[c]) for n, c in zip(fs_names, range(40, 60))}, index=data.index)

    blocks = {}
    for name, df, min_frac in [("Серология", ser, 0.5), ("Гистология", hist, 0.5), ("КЦЖК_кал", scfa, 0.5)]:
        keep = df.notna().mean(axis=1) >= min_frac
        blocks[name] = df[keep]
    return clin, blocks


def load_omics(fname):
    """Композиционные омиксы (строки нормированы на сумму 1). Пропущенные образцы отбрасываются."""
    d = pd.read_excel(ROOT / fname, na_values="ND")
    d = d[d["Visit number"] == VISIT].set_index("Patient").drop(columns="Visit number")
    d = d.apply(pd.to_numeric, errors="coerce")
    d = d[d.notna().mean(axis=1) >= 0.5]
    d.index = d.index.astype(int)
    return d


def load_drugs():
    d = pd.read_excel(ROOT / "Drugs.xlsx", na_values="ND")
    d = d[d["Visit number"] == VISIT].set_index("Patient").drop(columns="Visit number")
    d.index = d.index.astype(int)
    return d.apply(pd.to_numeric, errors="coerce")


OMICS = {
    "16S": "16S.xlsx",
    "Метаболом_кал": "Fecal metabolome.xlsx",
    "Метаболом_сыворотка": "Serum Metabolome.xlsx",
    "Липидом": "Lipidome.xlsx",
}
COMPOSITIONAL = set(OMICS)


def load_all():
    clin, blocks = load_inventory()
    for name, f in OMICS.items():
        blocks[name] = load_omics(f)
    return clin, blocks


if __name__ == "__main__":
    clin, blocks = load_all()
    print(clin.head())
    print(clin.diagnosis.value_counts())
    for k, v in blocks.items():
        print(f"{k:22s} {v.shape[0]:4d} пациентов × {v.shape[1]:4d} признаков")
