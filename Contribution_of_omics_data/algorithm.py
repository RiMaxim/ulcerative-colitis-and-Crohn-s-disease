"""Итоговый алгоритм (задача 1): базовая модель на клинике и рутинных анализах.

Для каждой конечной точки, на всех пациентах, у кого она определена:
  - AUC базовой модели вне выборки (20 × 5-кратная CV), 95% ДИ (бутстреп);
  - сравнение с одним кальпротектином (сколько даёт модель сверх него); если AUC одного
    кальпротектина вне выборки < 0,5 (сигнала нет, оценка «перевёрнута»), сравнение не делается;
  - AUC — по усреднённым по повторам вероятностям, как и 95% ДИ (так же, как в run_increment.py);
  - калибровка: наклон и смещение (идеал 1 и 0), точки по квинтилям риска для рисунка;
  - коэффициенты: отношение шансов на 1 стандартное отклонение признака с 95% ДИ
    (логистическая регрессия на всех пациентах, пропуски — медианой).

Запуск: python algorithm.py -> results/algorithm_models.csv, algorithm_coefs.csv, algorithm_calibration.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.special import logit

from data import load_all
from model import BaseTransformer
from run_increment import endpoints, prepare_clin

OUT = Path(__file__).resolve().parent / "results"
LABELS = {"age": "Возраст", "sex_male": "Мужской пол", "extraintestinal": "Внекишечные проявления",
          "is_CD": "Болезнь Крона (а не ЯК)", "log_calprotectin": "Кальпротектин (log)",
          "log_crp": "СРБ (log)", "rbc": "Эритроциты", "hb": "Гемоглобин", "wbc": "Лейкоциты", "lbp": "LBP"}


def cv_probs(X, y, n_rep=20, k=5):
    P = []
    for r in range(n_rep):
        p = np.zeros(len(y))
        for tr, te in StratifiedKFold(k, shuffle=True, random_state=r).split(X, y):
            t = BaseTransformer().fit(X[tr])
            m = LogisticRegression(C=1.0, max_iter=5000).fit(t.transform(X[tr]), y[tr])
            p[te] = m.predict_proba(t.transform(X[te]))[:, 1]
        P.append(p)
    return np.mean(P, 0), [roc_auc_score(y, p) for p in P]


def boot_ci(y, p, f, n=2000, seed=0):
    rng = np.random.default_rng(seed); v = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            v.append(f(y[i], p[i]) if not isinstance(p, tuple) else None)
    return np.percentile(v, [2.5, 97.5])


def calibration(y, p):
    lp = logit(np.clip(p, 1e-6, 1 - 1e-6))
    slope = sm.Logit(y, sm.add_constant(lp)).fit(disp=0).params[1]
    citl = sm.GLM(y, np.ones((len(y), 1)), family=sm.families.Binomial(), offset=lp).fit().params[0]
    return float(slope), float(citl)


def main():
    clin = prepare_clin(load_all()[0])
    mods, coefs, cal = [], [], []
    for ep, (y, base, _, desc) in endpoints(clin).items():
        X = clin.loc[y.index, base].values.astype(float)
        yy = y.values.astype(int)
        p, aucs = cv_probs(X, yy)
        pc, aucs_c = cv_probs(clin.loc[y.index, ["log_calprotectin"]].values.astype(float), yy)
        rng = np.random.default_rng(0); a, ac, d = [], [], []
        for _ in range(2000):
            i = rng.integers(0, len(yy), len(yy))
            if yy[i].min() == yy[i].max():
                continue
            u, uc = roc_auc_score(yy[i], p[i]), roc_auc_score(yy[i], pc[i])
            a.append(u); ac.append(uc); d.append(u - uc)
        slope, citl = calibration(yy, p)
        mods.append(dict(endpoint=ep, description=desc, n=len(yy), n_pos=int(yy.sum()),
                         auc=float(roc_auc_score(yy, p)), auc_lo=np.percentile(a, 2.5), auc_hi=np.percentile(a, 97.5),
                         auc_calpro=float(roc_auc_score(yy, pc)), auc_calpro_lo=np.percentile(ac, 2.5),
                         auc_calpro_hi=np.percentile(ac, 97.5),
                         delta_vs_calpro=float(roc_auc_score(yy, p) - roc_auc_score(yy, pc)),
                         delta_lo=np.percentile(d, 2.5), delta_hi=np.percentile(d, 97.5),
                         calib_slope=slope, calib_intercept=citl, brier=brier_score_loss(yy, p)))
        # коэффициенты на всех пациентах
        t = BaseTransformer().fit(X)
        Z = pd.DataFrame(t.transform(X), columns=base)
        Z = Z.loc[:, Z.std() > 0]
        r = sm.Logit(yy, sm.add_constant(Z)).fit(disp=0, maxiter=200)
        ci = r.conf_int()
        for f in Z.columns:
            coefs.append(dict(endpoint=ep, feature=f, label=LABELS.get(f, f), OR=np.exp(r.params[f]),
                              OR_lo=np.exp(ci.loc[f, 0]), OR_hi=np.exp(ci.loc[f, 1]), p=r.pvalues[f]))
        # точки калибровки по квинтилям
        q = pd.qcut(p, 5, labels=False, duplicates="drop")
        for g in np.unique(q):
            cal.append(dict(endpoint=ep, group=int(g), pred=float(p[q == g].mean()),
                            obs=float(yy[q == g].mean()), n=int((q == g).sum())))
    OUT.mkdir(exist_ok=True)
    M = pd.DataFrame(mods)
    bad = M.auc_calpro < 0.5
    M.loc[bad, ["delta_vs_calpro", "delta_lo", "delta_hi"]] = np.nan
    M.to_csv(OUT / "algorithm_models.csv", index=False)
    C = pd.DataFrame(coefs); C.to_csv(OUT / "algorithm_coefs.csv", index=False)
    pd.DataFrame(cal).to_csv(OUT / "algorithm_calibration.csv", index=False)
    with pd.option_context("display.width", 220):
        print(M.drop(columns="description").round(3).to_string(index=False))
        print(C[C.p < 0.05].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
