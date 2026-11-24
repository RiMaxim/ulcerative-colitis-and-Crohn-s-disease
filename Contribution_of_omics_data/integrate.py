"""Интеграция всех блоков сразу (задача 1): клиника + шкалы всех блоков в одной модели.

У каждого пациента свой набор блоков (все блоки сразу есть лишь у 39 человек), поэтому:
  - внутри каждого фолда для каждого блока обучается своя шкала на тех обучающих пациентах,
    у кого блок есть (как в model.py: L2-логистическая регрессия, C по внутренней CV);
  - шкала стандартизуется; у пациента без блока она равна 0 («блок ничего не говорит»);
  - итоговая модель — логистическая регрессия на [базовые признаки + шкалы всех блоков]
    (шкалы обучающих пациентов — out-of-fold, как в стекинге model.py).
Сравнение с базовой моделью на всех пациентах задачи; 20 × 5 CV, парный бутстреп 95% ДИ.
Для «ВЗК против контроля» берутся только блоки, измеренные у контролей (серология, 16S, генетика, HLA),
иначе само наличие анализа выдавало бы диагноз.

Запуск: python integrate.py -> results/integration_results.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

from data import load_all, COMPOSITIONAL
from model import BaseTransformer, BlockTransformer, _fit_block_model, _block_oof_scores
from run_increment import endpoints, prepare_clin

OUT = Path(__file__).resolve().parent / "results"
CONTROL_BLOCKS = {"Серология", "16S", "Генетика_WES", "HLA"}


def one_run(Xb, blocks, y, seed, k=5):
    n = len(y)
    p0, p1 = np.zeros(n), np.zeros(n)
    for f, (tr, te) in enumerate(StratifiedKFold(k, shuffle=True, random_state=seed).split(Xb, y)):
        bt = BaseTransformer().fit(Xb[tr])
        B_tr, B_te = bt.transform(Xb[tr]), bt.transform(Xb[te])
        p0[te] = LogisticRegression(C=1.0, max_iter=5000).fit(B_tr, y[tr]).predict_proba(B_te)[:, 1]
        S_tr, S_te = [], []
        for name, (Z, has) in blocks.items():
            s_tr, s_te = np.zeros(len(tr)), np.zeros(len(te))
            htr, hte = has[tr], has[te]
            ytr = y[tr][htr]
            if htr.sum() >= 20 and np.bincount(ytr, minlength=2).min() >= 4:
                blk = BlockTransformer(compositional=name in COMPOSITIONAL).fit(Z[tr][htr])
                Ztr, Zte = blk.transform(Z[tr][htr]), blk.transform(Z[te][hte])
                mb = _fit_block_model(Ztr, ytr, seed * 100 + f)
                C = float(np.atleast_1d(mb.C_)[0])
                o = _block_oof_scores(Ztr, ytr, C, seed * 100 + f)
                mu, sd = o.mean(), o.std() or 1.0
                s_tr[htr] = (o - mu) / sd
                if hte.any():
                    s_te[hte] = (mb.decision_function(Zte) - mu) / sd
            S_tr.append(s_tr); S_te.append(s_te)
        M_tr = np.column_stack([B_tr] + S_tr)
        M_te = np.column_stack([B_te] + S_te)
        p1[te] = LogisticRegression(C=1.0, max_iter=5000).fit(M_tr, y[tr]).predict_proba(M_te)[:, 1]
    return p0, p1


def evaluate(ep, y, Xb, blocks, n_rep=20, n_boot=2000):
    runs = [one_run(Xb, blocks, y, seed=r) for r in range(n_rep)]
    d = [roc_auc_score(y, p1) - roc_auc_score(y, p0) for p0, p1 in runs]
    P0 = np.mean([r[0] for r in runs], 0); P1 = np.mean([r[1] for r in runs], 0)
    rng = np.random.default_rng(0); bd = []
    for _ in range(n_boot):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            bd.append(roc_auc_score(y[i], P1[i]) - roc_auc_score(y[i], P0[i]))
    return dict(endpoint=ep, n=len(y), n_pos=int(y.sum()), blocks=", ".join(blocks),
                auc_m0=roc_auc_score(y, P0), auc_all=roc_auc_score(y, P1),
                delta=float(np.mean(d)), delta_lo=np.percentile(bd, 2.5), delta_hi=np.percentile(bd, 97.5),
                n_with_any_omics=None)


def main():
    clin, B = load_all()
    clin = prepare_clin(clin)
    B = {k: v for k, v in B.items() if k != "NOD2"}  # NOD2 уже входит в блок Генетика_WES
    jobs = []
    for ep, (y, base, excluded, _) in endpoints(clin).items():
        use = [b for b in B if b not in excluded and (ep != "IBD_vs_Control" or b in CONTROL_BLOCKS)]
        ids = y.index
        blocks = {}
        for b in use:
            Z = B[b].reindex(ids)
            has = Z.notna().any(axis=1).values
            blocks[b] = (Z.values.astype(float), has)
        jobs.append((ep, y.values.astype(int), clin.loc[ids, base].values.astype(float), blocks))
    res = Parallel(n_jobs=-1, verbose=0)(delayed(evaluate)(*j) for j in jobs)
    df = pd.DataFrame(res).drop(columns="n_with_any_omics")
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "integration_results.csv", index=False)
    print(df.drop(columns="blocks").round(3).to_string(index=False))


if __name__ == "__main__":
    main()
