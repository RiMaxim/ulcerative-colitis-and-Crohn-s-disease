"""Оценка прироста от блока данных: AUC(M0 + блок) − AUC(M0) на одних и тех же пациентах.

Схема (всё внутри фолдов кросс-валидации, тестовые пациенты не видны ни на одном шаге обучения):
  1. M0 — логистическая регрессия на базовых признаках (клиника + рутинные анализы).
  2. Блок → одна «шкала блока»: L2-логистическая регрессия на признаках блока,
     сила регуляризации подбирается внутренней CV по AUC.
  3. M0 + блок — логистическая регрессия на [базовые признаки, шкала блока] (стекинг).
     Шкала блока для обучающих пациентов берётся out-of-fold (внутренняя CV),
     иначе она была бы переобучена и стекинг завысил бы прирост.
  4. Повторённая стратифицированная k-fold CV (одинаковые разбиения для M0 и M0+блок).
  5. ΔAUC: среднее по повторам; 95% ДИ — парный бутстреп пациентов по усреднённым
     out-of-fold вероятностям; p — перестановочный тест (строки блока перемешиваются
     между пациентами, M0 и исход не трогаются: так проверяется именно прирост блока).
"""
from dataclasses import dataclass, field
import warnings
import numpy as np
from scipy.special import expit
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, brier_score_loss

C_GRID = np.logspace(-3, 1, 9)
warnings.filterwarnings("ignore", category=FutureWarning)  # смена API LogisticRegressionCV в sklearn 1.10


# ---------------------------------------------------------------- предобработка
class BaseTransformer(BaseEstimator, TransformerMixin):
    """Медианная импутация (по обучающей выборке) + стандартизация."""

    def fit(self, X, y=None):
        X = np.asarray(X, float)
        self.med_ = np.nanmedian(X, axis=0)
        self.med_ = np.where(np.isnan(self.med_), 0.0, self.med_)
        Xi = np.where(np.isnan(X), self.med_, X)
        self.mu_, self.sd_ = Xi.mean(0), Xi.std(0)
        self.sd_[self.sd_ == 0] = 1.0
        return self

    def transform(self, X):
        X = np.asarray(X, float)
        Xi = np.where(np.isnan(X), self.med_, X)
        return (Xi - self.mu_) / self.sd_


class BlockTransformer(BaseEstimator, TransformerMixin):
    """Предобработка признаков блока, параметры оцениваются только на обучающей выборке.

    compositional=True (16S, метаболомы, липидом — доли, сумма по строке = 1):
        фильтр распространённости (ненулевой у ≥ min_prev обучающих образцов),
        замена нулей на половину минимального положительного значения признака,
        CLR-преобразование, стандартизация.
    compositional=False (серология, гистология, КЦЖК):
        удаление признаков с > max_missing пропусков и нулевой дисперсией,
        медианная импутация, стандартизация.
    """

    def __init__(self, compositional=False, min_prev=0.2, max_missing=0.4):
        self.compositional = compositional
        self.min_prev = min_prev
        self.max_missing = max_missing

    def fit(self, X, y=None):
        X = np.asarray(X, float)
        if self.compositional:
            X = np.nan_to_num(X, nan=0.0)
            prev = (X > 0).mean(0)
            keep = prev >= self.min_prev
            if keep.sum() < 5:  # в очень разреженном блоке берём 5 самых распространённых
                keep = np.zeros_like(keep)
                keep[np.argsort(-prev)[:5]] = True
            self.keep_ = keep
            Xk = X[:, keep]
            pos = np.where(Xk > 0, Xk, np.inf).min(0)
            self.floor_ = np.where(np.isfinite(pos), pos / 2, 1e-9)
            Z = self._clr(Xk)
        else:
            miss = np.isnan(X).mean(0)
            keep = miss <= self.max_missing
            self.keep_ = keep
            Xk = X[:, keep]
            self.med_ = np.nanmedian(Xk, 0)
            Z = np.where(np.isnan(Xk), self.med_, Xk)
            var_ok = Z.std(0) > 0
            self.keep_[np.where(keep)[0][~var_ok]] = False
            self.med_ = self.med_[var_ok]
            Z = Z[:, var_ok]
        # все признаки постоянны в этой обучающей выборке (например, нет ни одного носителя
        # редкого аллеля) — блок даёт постоянную шкалу, т.е. «информации нет»
        self.empty_ = Z.shape[1] == 0
        self.mu_, self.sd_ = Z.mean(0), Z.std(0)
        self.sd_[self.sd_ == 0] = 1.0
        return self

    def _clr(self, Xk):
        L = np.log(np.maximum(Xk, self.floor_))
        return L - L.mean(1, keepdims=True)

    def transform(self, X):
        X = np.asarray(X, float)
        Xk = X[:, self.keep_]
        if self.compositional:
            Z = self._clr(np.nan_to_num(Xk, nan=0.0))
        else:
            Z = np.where(np.isnan(Xk), self.med_, Xk)
        if self.empty_:
            return np.zeros((X.shape[0], 1))
        return (Z - self.mu_) / self.sd_


# ---------------------------------------------------------------- модели
def _fit_block_model(Z, y, seed):
    """L2-логистическая регрессия; C подбирается внутренней CV по AUC."""
    n_min = np.bincount(y).min()
    cv = StratifiedKFold(min(5, n_min), shuffle=True, random_state=seed)
    m = LogisticRegressionCV(Cs=C_GRID, cv=cv, scoring="roc_auc", max_iter=5000)
    m.fit(Z, y)
    return m


def _block_oof_scores(Z, y, C, seed):
    """Out-of-fold шкала блока на обучающих пациентах (для стекинга)."""
    n_min = np.bincount(y).min()
    cv = StratifiedKFold(min(5, n_min), shuffle=True, random_state=seed + 7)
    s = np.zeros(len(y))
    for tr, te in cv.split(Z, y):
        m = LogisticRegression(C=C, max_iter=5000).fit(Z[tr], y[tr])
        s[te] = m.decision_function(Z[te])
    return s


def one_cv_run(Xb, Zblock, y, compositional, seed, k=5, base_C=1.0):
    """Одна повторность k-fold CV. Возвращает out-of-fold вероятности M0, блока и M0+блок."""
    n = len(y)
    p0, pb, p1 = np.zeros(n), np.zeros(n), np.zeros(n)
    cv = StratifiedKFold(k, shuffle=True, random_state=seed)
    for f, (tr, te) in enumerate(cv.split(Xb, y)):
        bt = BaseTransformer().fit(Xb[tr])
        B_tr, B_te = bt.transform(Xb[tr]), bt.transform(Xb[te])
        m0 = LogisticRegression(C=base_C, max_iter=5000).fit(B_tr, y[tr])
        p0[te] = m0.predict_proba(B_te)[:, 1]

        blk = BlockTransformer(compositional=compositional).fit(Zblock[tr])
        Z_tr, Z_te = blk.transform(Zblock[tr]), blk.transform(Zblock[te])
        mb = _fit_block_model(Z_tr, y[tr], seed * 100 + f)
        C = float(np.atleast_1d(mb.C_)[0])
        s_te = mb.decision_function(Z_te)
        pb[te] = expit(s_te)
        s_tr = _block_oof_scores(Z_tr, y[tr], C, seed * 100 + f)

        sd = s_tr.std() or 1.0
        S_tr = np.column_stack([B_tr, (s_tr - s_tr.mean()) / sd])
        S_te = np.column_stack([B_te, (s_te - s_tr.mean()) / sd])
        m1 = LogisticRegression(C=base_C, max_iter=5000).fit(S_tr, y[tr])
        p1[te] = m1.predict_proba(S_te)[:, 1]
    return p0, pb, p1


@dataclass
class IncrementResult:
    endpoint: str
    block: str
    n: int
    n_pos: int
    n_neg: int
    n_features: int
    auc_m0: float
    auc_m0_ci: tuple
    auc_block: float
    auc_m1: float
    auc_m1_ci: tuple
    delta: float
    delta_ci: tuple
    p_perm: float
    brier_m0: float
    brier_m1: float
    delta_repeats: list = field(repr=False, default_factory=list)
    null_deltas: list = field(repr=False, default_factory=list)


def _bootstrap(y, p0, p1, n_boot, rng):
    n = len(y)
    a0, a1, d = [], [], []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if y[i].min() == y[i].max():
            continue
        u0, u1 = roc_auc_score(y[i], p0[i]), roc_auc_score(y[i], p1[i])
        a0.append(u0); a1.append(u1); d.append(u1 - u0)
    q = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    return q(a0), q(a1), q(d)


def evaluate_increment(endpoint, block_name, Xbase, Xblock, y, compositional,
                       n_repeats=20, n_perm=100, perm_repeats=3, n_boot=2000, seed=0):
    """Прирост блока для одной конечной точки. Xbase, Xblock, y — на одних и тех же пациентах."""
    Xb = np.asarray(Xbase, float)
    Zb = np.asarray(Xblock, float)
    y = np.asarray(y, int)

    P0, PB, P1, deltas = [], [], [], []
    for r in range(n_repeats):
        p0, pb, p1 = one_cv_run(Xb, Zb, y, compositional, seed=seed + r)
        P0.append(p0); PB.append(pb); P1.append(p1)
        deltas.append(roc_auc_score(y, p1) - roc_auc_score(y, p0))
    P0, PB, P1 = map(lambda a: np.mean(a, 0), (P0, PB, P1))

    rng = np.random.default_rng(seed)
    ci0, ci1, cid = _bootstrap(y, P0, P1, n_boot, rng)

    obs = float(np.mean(deltas))
    null = []
    for b in range(n_perm):
        perm = rng.permutation(len(y))
        d = []
        for r in range(perm_repeats):
            p0, _, p1 = one_cv_run(Xb, Zb[perm], y, compositional, seed=10_000 + b * 10 + r)
            d.append(roc_auc_score(y, p1) - roc_auc_score(y, p0))
        null.append(float(np.mean(d)))
    p = (1 + sum(v >= obs for v in null)) / (1 + len(null)) if n_perm else np.nan

    n_feat = BlockTransformer(compositional=compositional).fit(Zb).keep_.sum()
    return IncrementResult(
        endpoint=endpoint, block=block_name, n=len(y), n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
        n_features=int(n_feat),
        auc_m0=float(roc_auc_score(y, P0)), auc_m0_ci=ci0,
        auc_block=float(roc_auc_score(y, PB)),
        auc_m1=float(roc_auc_score(y, P1)), auc_m1_ci=ci1,
        delta=obs, delta_ci=cid, p_perm=float(p),
        brier_m0=float(brier_score_loss(y, P0)), brier_m1=float(brier_score_loss(y, P1)),
        delta_repeats=[float(x) for x in deltas], null_deltas=null,
    )
