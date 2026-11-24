"""HTML-отчёт по приросту: python report.py -> results/increment_report.html"""
import base64
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

OUT = Path(__file__).resolve().parent / "results"
R = json.load(open(OUT / "increment_results.json"))
df = pd.DataFrame(R["results"])
# Поправка Бенджамини–Хохберга на все пары «конечная точка × блок»
_p = df.p_perm.values; _o = _p.argsort(); _m = len(_p)
_q = (_p[_o] * _m / (pd.Series(range(1, _m + 1)).values)).clip(max=1)
_q = pd.Series(_q[::-1]).cummin().values[::-1]
df["q_fdr"] = 0.0; df.loc[df.index[_o], "q_fdr"] = _q
EP_ORDER = ["CD_vs_UC", "IBD_vs_Control", "Exacerbation", "Endoscopic_A2plus", "Ileal_atrophy", "CD_ileal", "UC_extensive", "Advanced_therapy"]


def forest_png():
    eps = [e for e in EP_ORDER if e in set(df.endpoint)]
    fig, axes = plt.subplots(len(eps), 1, figsize=(8, 1.0 + 0.42 * len(df)), sharex=True,
                             gridspec_kw={"height_ratios": [max(1, (df.endpoint == e).sum()) for e in eps]})
    for ax, ep in zip(axes, eps):
        d = df[df.endpoint == ep].sort_values("delta")
        y = range(len(d))
        lo = [c[0] for c in d.delta_ci]; hi = [c[1] for c in d.delta_ci]
        col = ["#1a7f37" if (p < 0.05 and l > 0) else "#57606a" for p, l in zip(d.p_perm, lo)]
        ax.hlines(y, lo, hi, color=col, lw=2)
        ax.scatter(d.delta, y, color=col, zorder=3, s=28)
        ax.axvline(0, color="#999", lw=1, ls="--")
        ax.set_yticks(list(y)); ax.set_yticklabels([f"{b}  (n={n})" for b, n in zip(d.block, d.n)], fontsize=9)
        ax.set_title(R["endpoints"][ep][:70], fontsize=9.5, loc="left")
    axes[-1].set_xlabel("ΔAUC = AUC(M0 + блок) − AUC(M0), точка и 95% ДИ")
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=130); plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def table(ep):
    d = df[df.endpoint == ep].sort_values("delta", ascending=False)
    rows = []
    for _, r in d.iterrows():
        sig = r.p_perm < 0.05 and r.delta_ci[0] > 0
        rows.append(
            f"<tr{' class=sig' if sig else ''}><td>{r.block}</td><td>{r.n} ({r.n_pos}/{r.n_neg})</td>"
            f"<td>{r.n_features}</td>"
            f"<td>{r.auc_m0:.3f} <small>[{r.auc_m0_ci[0]:.2f}; {r.auc_m0_ci[1]:.2f}]</small></td>"
            f"<td>{r.auc_block:.3f}</td>"
            f"<td>{r.auc_m1:.3f} <small>[{r.auc_m1_ci[0]:.2f}; {r.auc_m1_ci[1]:.2f}]</small></td>"
            f"<td><b>{r.delta:+.3f}</b> <small>[{r.delta_ci[0]:+.3f}; {r.delta_ci[1]:+.3f}]</small></td>"
            f"<td>{r.p_perm:.3f}</td><td>{r.q_fdr:.2f}</td>"
            f"<td>{r.brier_m0:.3f} → {r.brier_m1:.3f}</td></tr>")
    return ("<table><tr><th>Блок</th><th>n (1/0)</th><th>Признаков</th><th>AUC M0</th>"
            "<th>AUC только блок</th><th>AUC M0+блок</th><th>ΔAUC [95% ДИ]</th><th>p (перест.)</th><th>q (FDR)</th>"
            "<th>Brier M0 → M0+блок</th></tr>" + "".join(rows) + "</table>")


cfg = R["config"]
html = f"""<!doctype html><html lang=ru><head><meta charset=utf-8><title>Прирост омиксов</title>
<style>body{{font-family:-apple-system,Segoe UI,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#1f2328}}
table{{border-collapse:collapse;width:100%;margin:8px 0 24px;font-size:14px}}td,th{{border-bottom:1px solid #d0d7de;padding:5px 8px;text-align:left}}
th{{background:#f6f8fa}}tr.sig td{{background:#dafbe1}}small{{color:#57606a}}code{{background:#f6f8fa;padding:1px 4px}}</style></head><body>
<h1>Прирост AUC от омиксных блоков к базовой модели</h1>
<p>Визит 1, одна строка на пациента. Генетика: 18 известных вариантов WES (10 признаков), NOD2 отдельно и HLA.
Повторов CV: {cfg['n_repeats']} × 5 фолдов; перестановок: {cfg['n_perm']}; бутстреп: {cfg['n_boot']}.</p>
<h2>Как читать</h2>
<ul><li><b>M0</b> — клиника и рутинные анализы: возраст, пол, кальпротектин, СРБ, эритроциты, Hb, лейкоциты, LBP;
для задач внутри ВЗК ещё тип ВЗК и внекишечные проявления.</li>
<li><b>ΔAUC</b> — насколько AUC растёт, если к M0 добавить блок. Считается на одних и тех же пациентах (у кого есть блок),
поэтому AUC M0 в разных строках немного различается.</li>
<li><b>95% ДИ</b> — парный бутстреп пациентов. <b>p</b> — перестановочный тест: блок «отвязывается» от пациентов, и проверяется,
бывает ли случайно прирост не меньше наблюдаемого.</li>
<li>Зелёным отмечены строки с p &lt; 0.05 и ДИ выше нуля. q — p-value с поправкой Бенджамини–Хохберга на все пары: устойчивым считаем прирост с q &lt; 0.1. Прирост около нуля или отрицательный значит, что блок ничего не добавляет к M0
(отрицательный — цена лишней сложности на малой выборке).</li></ul>
<img style="max-width:100%" src="data:image/png;base64,{forest_png()}">
""" + "".join(f"<h2>{R['endpoints'][ep]}</h2>{table(ep)}" for ep in EP_ORDER if ep in set(df.endpoint)) + \
    "<h2>Пропущенные пары</h2><ul>" + "".join(f"<li>{a} × {b}: {c}</li>" for a, b, c in R["skipped"]) + "</ul></body></html>"

(OUT / "increment_report.html").write_text(html, encoding="utf-8")
print("saved", OUT / "increment_report.html")
