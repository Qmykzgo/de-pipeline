#!/usr/bin/env python3
"""
Ultra-Lean Differential Expression Pipeline (depipeline1)
Dependencies: pandas, numpy, matplotlib (Zero extra 3rd-party libraries)
- PCA: Pure NumPy Singular Value Decomposition (np.linalg.svd)
- DE Engine: Welch's t-test on log2(CPM) + Continued Fraction Regularized Incomplete Beta (pure Python/NumPy)
- FDR Correction: Benjamini-Hochberg procedure (pure NumPy)
- Enrichment: Enrichr REST API via standard library urllib.request
"""
from __future__ import annotations
import argparse
import base64
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import sys
import time
from datetime import datetime, timezone
import urllib.request
import urllib.parse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DISCLAIMER = "This pipeline is a research and educational tool — not for clinical diagnostic use."
ENRICHR_LIBRARIES = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "WikiPathways_2023_Human",
]
ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"
PALETTE = {
    "neutral": "#6B7280",
    "up": "#DC2626",
    "down": "#2563EB",
    "accent": "#7C3AED",
    "highlight": "#F59E0B",
}

# ── Math & Statistics (Pure Python / NumPy) ───────────────────────────────────
def _betacf(a: float, b: float, x: float, max_iter: int = 150, eps: float = 3e-15) -> float:
    """Continued fraction for regularized incomplete beta using modified Lentz method."""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        # Even step
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        # Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        del_h = d * c
        h *= del_h
        if abs(del_h - 1.0) < eps:
            break
    return h


def _betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta function I_x(a, b) in pure Python."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    # Use symmetry relation if x > (a+1)/(a+b+2) for fast convergence
    if x > (a + 1.0) / (a + b + 2.0):
        return 1.0 - _betainc(b, a, 1.0 - x)
    log_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1.0 - x) - log_beta) / a
    return front * _betacf(a, b, x)


def _t_pvalue(t_stat: float, dof: float) -> float:
    """Compute two-tailed p-value for Student's / Welch's t distribution."""
    if math.isnan(t_stat) or math.isnan(dof) or dof <= 0:
        return 1.0
    t2 = t_stat * t_stat
    x = dof / (dof + t2)
    # Two-tailed p-value is I_x(dof/2, 1/2)
    p = _betainc(dof / 2.0, 0.5, x)
    return max(0.0, min(1.0, float(p)))


def _bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg step-up procedure for FDR adjustment."""
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = np.ones(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        v = ranked[i] * n / (i + 1)
        prev = min(prev, v)
        adj[i] = prev
    out = np.empty(n)
    out[order] = np.clip(adj, 0.0, 1.0)
    return out


# ── I/O & Parsing ─────────────────────────────────────────────────────────────
def _sep(path: Path) -> str:
    return "\t" if path.suffix.lower() == ".tsv" else ","


def load_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=_sep(path))
    if df.shape[1] < 3:
        raise ValueError("Counts file must contain 1 gene column and at least 2 sample columns")
    counts = df.set_index(df.columns[0]).apply(pd.to_numeric, errors="coerce")
    if counts.isna().any().any():
        raise ValueError("Non-numeric or missing values found in count matrix")
    if (counts < 0).any().any():
        raise ValueError("Negative values found in count matrix (raw counts required)")
    return counts.rename_axis("gene").rename(str, axis=1)


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep=_sep(path))
    if "sample_id" not in meta.columns:
        raise ValueError("Metadata table must contain a 'sample_id' column")
    meta["sample_id"] = meta["sample_id"].astype(str)
    return meta.set_index("sample_id")


def parse_formula(f: str) -> list[str]:
    f = f.strip()
    if not f.startswith("~"):
        raise ValueError("Design formula must start with '~' (e.g. ~ condition or ~ batch + condition)")
    terms = [t.strip() for t in f[1:].split("+") if t.strip()]
    if not terms:
        raise ValueError("Formula needs at least one design term")
    return terms


def parse_contrast(c: str) -> tuple[str, str, str]:
    parts = [p.strip() for p in c.split(",")]
    if len(parts) != 3 or any(not p for p in parts):
        raise ValueError("Contrast must be in format: factor,numerator,denominator (e.g. condition,treated,control)")
    return parts[0], parts[1], parts[2]


def validate(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    terms: list[str],
    factor: str,
    num: str,
    den: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = [s for s in counts.columns if s not in meta.index]
    if missing:
        raise ValueError(f"Samples in counts missing from metadata: {missing[:5]}")
    meta = meta.loc[counts.columns].copy()
    for t in terms:
        if t not in meta.columns:
            raise ValueError(f"Design term {t!r} not found in metadata columns")
    if factor not in meta.columns:
        raise ValueError(f"Contrast factor {factor!r} not found in metadata columns")
    g = meta[factor].astype(str)
    if num not in set(g):
        raise ValueError(f"Numerator level {num!r} not found in factor {factor!r}")
    if den not in set(g):
        raise ValueError(f"Denominator level {den!r} not found in factor {factor!r}")
    if int((g == num).sum()) < 2 or int((g == den).sum()) < 2:
        raise ValueError("Each contrast group requires at least 2 replicates for statistical testing")
    return counts, meta


# ── Quality Control & Normalization ───────────────────────────────────────────
def compute_qc(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    det = (counts > 0).sum(axis=0)
    med = counts.median(axis=0)
    return pd.DataFrame({
        "sample_id": counts.columns.tolist(),
        "library_size": lib.values,
        "detected_genes": det.values,
        "median_count": med.values,
    })


def filter_low(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 2) -> pd.DataFrame:
    mask = (counts >= min_count).sum(axis=1) >= min_samples
    out = counts.loc[mask].copy()
    if out.shape[0] < 2:
        raise ValueError("Fewer than 2 genes retained after low-count filtering")
    return out


def norm_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Compute log2(CPM + 1) normalized matrix."""
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


# ── Pure NumPy SVD PCA ────────────────────────────────────────────────────────
def run_pca(norm: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Perform Principal Component Analysis via Singular Value Decomposition (np.linalg.svd).
    Mathematically identical to sklearn PCA with zero extra dependencies.
    """
    # X shape: (samples, genes)
    X = norm.T.values.astype(float)
    n_samples, n_genes = X.shape
    # Center genes
    X_centered = X - np.mean(X, axis=0)

    # SVD: X_centered = U @ diag(S) @ Vt
    U, S, _ = np.linalg.svd(X_centered, full_matrices=False)
    coords = U * S  # PC projections

    # Variance explained = S_i^2 / (n_samples - 1)
    dof = max(n_samples - 1, 1)
    variance = (S ** 2) / dof
    total_var = np.sum(variance) if np.sum(variance) > 0 else 1.0
    var_ratio = variance / total_var

    pc1 = coords[:, 0].tolist() if coords.shape[1] > 0 else [0.0] * n_samples
    pc2 = coords[:, 1].tolist() if coords.shape[1] > 1 else [0.0] * n_samples

    df = pd.DataFrame({
        "sample_id": norm.columns.tolist(),
        "PC1": pc1,
        "PC2": pc2,
    })
    return df, var_ratio


# ── Ultra-Lean Differential Expression Engine ─────────────────────────────────
def run_de(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    factor: str,
    num: str,
    den: str,
) -> pd.DataFrame:
    """
    Perform Welch's t-test on log2(CPM) with pure Python/NumPy p-values and FDR.
    """
    g = meta[factor].astype(str)
    ns = g[g == num].index.tolist()
    ds = g[g == den].index.tolist()

    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    log2cpm = np.log2(cpm + 1.0)

    xn = log2cpm[ns]
    xd = log2cpm[ds]

    mn = xn.mean(axis=1)
    md = xd.mean(axis=1)
    vn = xn.var(axis=1, ddof=1)
    vd = xd.var(axis=1, ddof=1)

    nn = float(len(ns))
    nd_ = float(len(ds))

    vn_n = vn / nn + 1e-12
    vd_n = vd / nd_ + 1e-12

    se = np.sqrt(vn_n + vd_n)
    t_stat = (mn - md) / se

    # Welch-Satterthwaite degrees of freedom
    dof = (vn_n + vd_n) ** 2 / (vn_n ** 2 / max(nn - 1, 1) + vd_n ** 2 / max(nd_ - 1, 1))
    dof = np.maximum(dof, 1.0)

    # Compute p-values via regularized incomplete beta function
    pvals = np.array([
        _t_pvalue(float(t_i), float(d_i))
        for t_i, d_i in zip(t_stat.values, dof.values)
    ])
    padj = _bh(pvals)
    base = cpm.mean(axis=1)

    res = pd.DataFrame({
        "gene": counts.index.tolist(),
        "baseMean": base.values,
        "log2FoldChange": (mn - md).values,
        "stat": t_stat.values,
        "pvalue": pvals,
        "padj": padj,
    }).sort_values("padj").reset_index(drop=True)

    return res


# ── Enrichr REST Client (Standard Library urllib.request) ───────────────────────
def run_enrichr(genes: list[str]) -> pd.DataFrame:
    """Query Enrichr REST API using Python's standard library urllib.request."""
    if not genes:
        return pd.DataFrame()

    url_add = f"{ENRICHR_BASE}/addList"
    boundary = "----WebKitFormBoundary7MA4YWxkTrZu0gW"
    gene_list_str = "\n".join(genes)

    body = (
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="list"\r\n\r\n'
        f"{gene_list_str}\r\n"
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="description"\r\n\r\n'
        f"DE_genes\r\n"
        f"--{boundary}--\r\n"
    ).encode("utf-8")

    req = urllib.request.Request(
        url_add,
        data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
        method="POST",
    )

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            res_json = json.loads(resp.read().decode("utf-8"))
            user_list_id = res_json.get("userListId")
    except Exception as e:
        print(f"[WARN] Enrichr API submit failed: {e}")
        return pd.DataFrame()

    if not user_list_id:
        return pd.DataFrame()

    rows = []
    for lib in ENRICHR_LIBRARIES:
        time.sleep(0.3)
        params = urllib.parse.urlencode({"userListId": user_list_id, "backgroundType": lib})
        url_enrich = f"{ENRICHR_BASE}/enrich?{params}"
        try:
            req_enrich = urllib.request.Request(url_enrich)
            with urllib.request.urlopen(req_enrich, timeout=30) as resp:
                enr_json = json.loads(resp.read().decode("utf-8"))
                entries = enr_json.get(lib, [])
                for e in entries:
                    if len(e) < 7:
                        continue
                    rows.append({
                        "library": lib,
                        "term": e[1],
                        "pvalue": float(e[2]),
                        "odds_ratio": float(e[3]),
                        "combined_score": float(e[4]),
                        "overlap_genes": ";".join(e[5]) if isinstance(e[5], list) else str(e[5]),
                        "adj_pvalue": float(e[6]),
                    })
        except Exception as e:
            print(f"[WARN] Enrichr {lib}: {e}")

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values("combined_score", ascending=False).reset_index(drop=True)


# ── Visualizations (Matplotlib) ───────────────────────────────────────────────
def _meta_with_samples(meta: pd.DataFrame) -> pd.DataFrame:
    m = meta.copy()
    m.index.name = "sample_id"
    return m.reset_index()


def plot_pca(pca_df: pd.DataFrame, meta: pd.DataFrame, factor: str, var_ratio: np.ndarray, out: Path):
    meta_r = _meta_with_samples(meta)
    pca2 = pca_df.copy()
    pca2["sample_id"] = pca2["sample_id"].astype(str)
    meta_r["sample_id"] = meta_r["sample_id"].astype(str)
    merged = pca2.merge(meta_r, on="sample_id", how="left")

    fig, ax = plt.subplots(figsize=(7, 5))
    groups = merged[factor].unique()
    colors = plt.cm.Set2(np.linspace(0, 1, len(groups)))

    for group, color in zip(groups, colors):
        sub = merged[merged[factor] == group]
        ax.scatter(sub["PC1"], sub["PC2"], label=str(group), color=color, s=80, alpha=0.9, edgecolors="white", linewidths=0.5)
        for _, row in sub.iterrows():
            ax.annotate(row["sample_id"], (row["PC1"], row["PC2"]), textcoords="offset points", xytext=(5, 3), fontsize=7, alpha=0.7)

    p1 = f"{var_ratio[0]*100:.1f}%" if len(var_ratio) > 0 else "?"
    p2 = f"{var_ratio[1]*100:.1f}%" if len(var_ratio) > 1 else "?"
    ax.set_xlabel(f"PC1 ({p1})", fontsize=11)
    ax.set_ylabel(f"PC2 ({p2})", fontsize=11)
    ax.set_title("PCA - Sample Separation (Pure SVD)", fontsize=13, fontweight="bold")
    ax.legend(title=factor, framealpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_volcano(de: pd.DataFrame, out: Path):
    df = de.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FoldChange", "padj"]).copy()
    y = -np.log10(df["padj"].clip(lower=1e-300))
    up = (df["padj"] < 0.05) & (df["log2FoldChange"] >= 1.0)
    down = (df["padj"] < 0.05) & (df["log2FoldChange"] <= -1.0)
    ns = ~(up | down)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.scatter(df.loc[ns, "log2FoldChange"], y[ns], s=10, alpha=0.45, color=PALETTE["neutral"], label="NS", rasterized=True)
    ax.scatter(df.loc[up, "log2FoldChange"], y[up], s=14, alpha=0.8, color=PALETTE["up"], label=f"Up ({up.sum()})")
    ax.scatter(df.loc[down, "log2FoldChange"], y[down], s=14, alpha=0.8, color=PALETTE["down"], label=f"Down ({down.sum()})")

    sig_df = df.loc[up | down]
    if not sig_df.empty:
        top5 = sig_df.nsmallest(5, "padj")
        for _, row in top5.iterrows():
            ax.annotate(
                row["gene"],
                (row["log2FoldChange"], -np.log10(row["padj"])),
                fontsize=7.5,
                ha="center",
                xytext=(0, 6),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
            )

    ax.axvline(1.0, ls="--", lw=0.8, color="gray")
    ax.axvline(-1.0, ls="--", lw=0.8, color="gray")
    ax.axhline(-np.log10(0.05), ls="--", lw=0.8, color="gray")
    ax.set_xlabel("log2 Fold Change", fontsize=11)
    ax.set_ylabel("-log10 adj.p", fontsize=11)
    ax.set_title("Volcano Plot", fontsize=13, fontweight="bold")
    ax.legend(framealpha=0.8, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_ma(de: pd.DataFrame, out: Path):
    df = de.replace([np.inf, -np.inf], np.nan).dropna(subset=["baseMean", "log2FoldChange", "padj"]).copy()
    x = np.log10(df["baseMean"].clip(lower=1e-6))
    y = df["log2FoldChange"]
    sig = df["padj"] < 0.05

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(x[~sig], y[~sig], s=8, alpha=0.4, color=PALETTE["neutral"], label="NS", rasterized=True)
    ax.scatter(x[sig], y[sig], s=10, alpha=0.8, color=PALETTE["accent"], label=f"Sig ({sig.sum()})")
    ax.axhline(0, ls="--", lw=0.8, color="gray")
    ax.set_xlabel("log10 baseMean", fontsize=11)
    ax.set_ylabel("log2 Fold Change", fontsize=11)
    ax.set_title("MA Plot", fontsize=13, fontweight="bold")
    ax.legend(framealpha=0.8, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_heatmap(norm: pd.DataFrame, de: pd.DataFrame, meta: pd.DataFrame, factor: str, out: Path, top_n: int = 30):
    sig_df = de.dropna(subset=["padj"])
    sig_df = sig_df[sig_df["padj"] < 0.05].nsmallest(top_n, "padj")
    genes = sig_df["gene"].tolist()
    genes = [g for g in genes if g in norm.index]
    if not genes:
        return
    mat = norm.loc[genes]
    std = mat.std(axis=1).replace(0, 1)
    z = mat.subtract(mat.mean(axis=1), axis=0).div(std, axis=0)

    meta_r = _meta_with_samples(meta)
    meta_r["sample_id"] = meta_r["sample_id"].astype(str)
    col_order = meta_r.sort_values(factor)["sample_id"].tolist()
    col_order = [c for c in col_order if c in z.columns]
    z = z[col_order]

    fig, ax = plt.subplots(figsize=(max(7, len(col_order) * 0.7), max(6, len(genes) * 0.35)))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_title(f"Top {len(genes)} DE Genes (Z-scored log2 CPM)", fontsize=12, fontweight="bold")
    plt.colorbar(im, ax=ax, shrink=0.6, label="Z-score")

    group_order = meta_r.set_index("sample_id").loc[col_order, factor].values
    for i in range(len(group_order) - 1):
        if group_order[i] != group_order[i + 1]:
            ax.axvline(i + 0.5, color="black", lw=1.5)

    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_enrichment_bubble(enrichment: pd.DataFrame, out: Path, top_n: int = 15):
    if enrichment.empty:
        return
    df = enrichment[enrichment["adj_pvalue"] < 0.05].nlargest(top_n, "combined_score")
    if df.empty:
        df = enrichment.nlargest(top_n, "combined_score")
    df = df.copy()
    df["term_short"] = df["term"].str[:55]
    df["n_overlap"] = df["overlap_genes"].apply(lambda x: len(x.split(";")) if x else 0)

    fig, ax = plt.subplots(figsize=(8, 0.45 * len(df) + 2))
    sc = ax.scatter(
        df["combined_score"],
        df["term_short"],
        s=df["n_overlap"] * 10 + 20,
        c=-np.log10(df["adj_pvalue"].clip(1e-50)),
        cmap="Oranges",
        alpha=0.85,
        edgecolors="gray",
        linewidths=0.5,
    )
    plt.colorbar(sc, ax=ax, label="-log10 adj.p")
    ax.set_xlabel("Combined Score (Enrichr)", fontsize=11)
    ax.set_title("Pathway Enrichment", fontsize=12, fontweight="bold")
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ── Reproducibility Bundle ────────────────────────────────────────────────────
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def write_repro(out_dir: Path, counts_path: Path, meta_path: Path, formula: str, contrast: str):
    repro = out_dir / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"python3 de_pipeline.py \\\n"
        f"  --counts {counts_path.resolve()} \\\n"
        f"  --metadata {meta_path.resolve()} \\\n"
        f"  --formula \"{formula}\" \\\n"
        f"  --contrast \"{contrast}\" \\\n"
        f"  --output {out_dir.resolve()}\n"
    )
    (repro / "commands.sh").write_text(cmd)
    env = (
        "name: depipeline1\n"
        "channels:\n"
        "  - conda-forge\n"
        "dependencies:\n"
        "  - python>=3.10\n"
        "  - pandas\n"
        "  - numpy\n"
        "  - matplotlib\n"
    )
    (repro / "environment.yml").write_text(env)
    cs = []
    for p in [counts_path, meta_path]:
        if p.exists():
            cs.append(f"{_sha256(p)}  {p.name}")
    for p in sorted((out_dir / "tables").glob("*.csv")):
        cs.append(f"{_sha256(p)}  tables/{p.name}")
    for p in sorted((out_dir / "figures").glob("*.png")):
        cs.append(f"{_sha256(p)}  figures/{p.name}")
    (repro / "checksums.sha256").write_text("\n".join(cs) + "\n")


# ── Reports (Markdown & Self-Contained HTML) ──────────────────────────────────
def write_md_report(out_dir: Path, n_samples: int, n_pre: int, n_post: int, formula: str, contrast: str, de: pd.DataFrame, enr: pd.DataFrame):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    sig = de[(de["padj"].notna()) & (de["padj"] < 0.05) & (de["log2FoldChange"].abs() >= 1.0)]
    n_up = int((sig["log2FoldChange"] > 0).sum())
    n_down = int((sig["log2FoldChange"] < 0).sum())
    rows = "\n".join(
        f"| `{r.gene}` | {r.log2FoldChange:+.3f} | {r.pvalue:.3e} | {r.padj:.3e} |"
        for r in de.head(10).itertuples()
    )
    report = f"""# Differential Expression Report (Ultra-Lean Pipeline)
**Generated**: {now}
**Samples**: {n_samples} | **Genes pre-filter**: {n_pre} | **Genes post-filter**: {n_post}
**Formula**: `{formula}` | **Contrast**: `{contrast}`
**Backend**: `Welch log2(CPM) + Pure NumPy SVD PCA` (Dependencies: `pandas`, `numpy`, `matplotlib`)
---
## 1. QC & PCA
QC table: `tables/qc_summary.csv`
![PCA](figures/pca.png)
---
## 2. Differential Expression
**Significant** (padj < 0.05 & |log2FC| >= 1.0): **{len(sig)}** ({n_up} up, {n_down} down)
Full results: `tables/de_results.csv`
### Top 10 Genes
| Gene | log2FC | p-value | padj |
|------|-------:|--------:|-----:|
{rows}
![Volcano](figures/volcano.png)
![MA Plot](figures/ma_plot.png)
![Heatmap](figures/heatmap.png)
---
## 3. Pathway Enrichment
"""
    if not enr.empty:
        report += "Results: `tables/enrichment_results.csv`\n![Enrichment](figures/enrichment_bubble.png)\n"
        top5 = enr[enr["adj_pvalue"] < 0.05].nlargest(5, "combined_score")
        if not top5.empty:
            report += "| Term | Library | adj.p | Score |\n|------|---------|------:|------:|\n"
            for r in top5.itertuples():
                report += f"| {r.term[:60]} | {r.library} | {r.adj_pvalue:.3e} | {r.combined_score:.1f} |\n"
    else:
        report += "_No enrichment results._\n"
    report += f"\n---\n## 4. Reproducibility\n`reproducibility/commands.sh` | `reproducibility/environment.yml` | `reproducibility/checksums.sha256`\n\n> {DISCLAIMER}\n"
    (out_dir / "report.md").write_text(report)


def write_html(out_dir: Path):
    md = (out_dir / "report.md").read_text()
    figs = out_dir / "figures"
    for fname in ["pca.png", "volcano.png", "ma_plot.png", "heatmap.png", "enrichment_bubble.png"]:
        p = figs / fname
        if p.exists():
            b64 = base64.b64encode(p.read_bytes()).decode()
            img = f'<img src="data:image/png;base64,{b64}" style="max-width:100%;border-radius:6px;margin:12px 0;">'
            md = re.sub(r"!\[[^\]]*\]\(figures/" + fname.replace(".", r"\.") + r"\)", img, md)
    lines = md.split("\n")
    out_l = []
    in_code = False
    in_table = False

    def mi(t: str) -> str:
        t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
        t = re.sub(r"`(.+?)`", r"<code>\1</code>", t)
        return t

    for line in lines:
        if line.startswith("```"):
            if in_code:
                out_l.append("</code></pre>")
                in_code = False
            else:
                out_l.append("<pre><code>")
                in_code = True
            continue
        if in_code:
            out_l.append(line.replace("<", "&lt;").replace(">", "&gt;"))
            continue
        if line.startswith("<img"):
            out_l.append(line)
            continue
        if line.strip() == "---":
            out_l.append("<hr>")
            continue
        matched = False
        for lvl in range(6, 0, -1):
            if line.startswith("#" * lvl + " "):
                out_l.append(f"<h{lvl}>{mi(line[lvl+1:])}</h{lvl}>")
                matched = True
                break
        if matched:
            continue
        if line.startswith("|"):
            if not in_table:
                out_l.append('<table border="1">')
                in_table = True
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if all(set(c) <= set("-: ") for c in cells):
                continue
            out_l.append("<tr>" + "".join(f"<td>{mi(c)}</td>" for c in cells) + "</tr>")
        else:
            if in_table:
                out_l.append("</table>")
                in_table = False
            if line.strip():
                p = re.sub(r"> (.+)", r"<blockquote>\1</blockquote>", mi(line))
                out_l.append(f"<p>{p}</p>")
            else:
                out_l.append("")
    if in_table:
        out_l.append("</table>")
    css = (
        "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:960px;margin:40px auto;padding:0 20px;color:#1f2937;line-height:1.6;background:#f9fafb}"
        "h1{color:#111827;border-bottom:3px solid #7C3AED;padding-bottom:8px}"
        "h2{color:#374151;border-bottom:1px solid #e5e7eb;padding-bottom:4px;margin-top:2em}"
        "table{border-collapse:collapse;width:100%;margin:1em 0;font-size:.9em}"
        "th,td{border:1px solid #e5e7eb;padding:8px 12px}"
        "th{background:#f3f4f6;font-weight:600}"
        "code{background:#f3f4f6;padding:2px 5px;border-radius:4px;font-size:.88em}"
        "pre code{display:block;padding:12px;overflow-x:auto}"
        "blockquote{border-left:4px solid #7C3AED;margin:1em 0;padding:4px 16px;color:#6b7280}"
        "hr{border:none;border-top:1px solid #e5e7eb;margin:2em 0}"
        "img{box-shadow:0 2px 8px rgba(0,0,0,.12);max-width:100%}"
        "p{margin:.6em 0}"
    )
    html = (
        f'<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">'
        f'<title>Ultra-Lean DE Report</title><style>{css}</style></head>'
        f"<body>\n" + "\n".join(out_l) + "</body></html>"
    )
    (out_dir / "report.html").write_text(html)


# ── Main Pipeline Execution ───────────────────────────────────────────────────
def run_pipeline(
    counts_path: Path,
    meta_path: Path,
    formula: str,
    contrast: str,
    out_dir: Path,
    min_count: int = 10,
    min_samples: int = 2,
    run_enrichment: bool = True,
) -> dict:
    if out_dir.exists() and any(out_dir.iterdir()):
        raise FileExistsError(f"Output directory {out_dir!r} is not empty — specify a new directory or remove existing files")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(exist_ok=True)
    (out_dir / "tables").mkdir(exist_ok=True)

    print("Step 1/7  Loading & validating inputs ...", flush=True)
    terms = parse_formula(formula)
    factor, num, den = parse_contrast(contrast)
    counts = load_counts(counts_path)
    meta = load_metadata(meta_path)
    counts, meta = validate(counts, meta, terms, factor, num, den)
    n_pre = counts.shape[0]
    print(f"  {counts.shape[1]} samples × {n_pre} genes loaded", flush=True)

    print("Step 2/7  Computing QC statistics ...", flush=True)
    qc = compute_qc(counts)
    qc.to_csv(out_dir / "tables" / "qc_summary.csv", index=False)
    print(f"  Library sizes: {qc['library_size'].min():,} – {qc['library_size'].max():,} reads", flush=True)

    print("Step 3/7  Filtering & CPM normalization ...", flush=True)
    filtered = filter_low(counts, min_count, min_samples)
    n_post = filtered.shape[0]
    norm = norm_cpm(filtered)
    norm.to_csv(out_dir / "tables" / "normalized_counts.csv")
    print(f"  Genes retained: {n_post}/{n_pre}", flush=True)

    print("Step 4/7  Performing PCA (Pure NumPy SVD) ...", flush=True)
    pca_df, var_ratio = run_pca(norm)
    plot_pca(pca_df, meta, factor, var_ratio, out_dir / "figures" / "pca.png")

    print(f"Step 5/7  Differential expression testing (Welch log2-CPM) ...", flush=True)
    de = run_de(filtered, meta, factor, num, den)
    de.to_csv(out_dir / "tables" / "de_results.csv", index=False)
    sig = de[(de["padj"].notna()) & (de["padj"] < 0.05) & (de["log2FoldChange"].abs() >= 1.0)]
    print(f"  Significant genes: {len(sig)} (FDR < 0.05, |log2FC| >= 1.0)", flush=True)

    print("Step 6/7  Generating publication figures ...", flush=True)
    plot_volcano(de, out_dir / "figures" / "volcano.png")
    plot_ma(de, out_dir / "figures" / "ma_plot.png")
    plot_heatmap(norm, de, meta, factor, out_dir / "figures" / "heatmap.png")

    enr = pd.DataFrame()
    if run_enrichment and not sig.empty:
        print("Step 7/7  Querying Enrichr pathways (urllib.request) ...", flush=True)
        enr = run_enrichr(sig["gene"].tolist())
        if not enr.empty:
            enr.to_csv(out_dir / "tables" / "enrichment_results.csv", index=False)
            plot_enrichment_bubble(enr, out_dir / "figures" / "enrichment_bubble.png")
            print(f"  {len(enr)} enriched pathway terms retrieved", flush=True)
    else:
        print("Step 7/7  Skipping enrichment (no significant genes or disabled)", flush=True)

    write_md_report(out_dir, counts.shape[1], n_pre, n_post, formula, contrast, de, enr)
    write_html(out_dir)
    write_repro(out_dir, counts_path, meta_path, formula, contrast)

    summary = {
        "samples": counts.shape[1],
        "genes_pre": n_pre,
        "genes_post": n_post,
        "formula": formula,
        "contrast": contrast,
        "backend_used": "welch_log2cpm_svd",
        "dependencies": ["pandas", "numpy", "matplotlib"],
        "n_significant": len(sig),
        "n_up": int((sig["log2FoldChange"] > 0).sum()),
        "n_down": int((sig["log2FoldChange"] < 0).sum()),
        "n_enriched_terms": len(enr),
        "output_dir": str(out_dir),
        "disclaimer": DISCLAIMER,
    }
    (out_dir / "result.json").write_text(json.dumps(summary, indent=2))
    print(f"\nPipeline finished successfully -> {out_dir}/", flush=True)
    print("  View report.html for the interactive visual report", flush=True)
    return summary


# ── CLI Interface ─────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description="Ultra-Lean Differential Expression Pipeline (pandas + numpy + matplotlib only)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--counts", help="Path to raw count matrix (CSV/TSV)")
    p.add_argument("--metadata", help="Path to sample metadata table (CSV/TSV)")
    p.add_argument("--formula", default="~ condition", help="Design formula (e.g. ~ condition or ~ batch + condition)")
    p.add_argument("--contrast", default="condition,treated,control", help="Contrast: factor,numerator,denominator")
    p.add_argument("--min-count", type=int, default=10, help="Minimum raw count threshold for filtering")
    p.add_argument("--min-samples", type=int, default=2, help="Minimum number of samples meeting min-count")
    p.add_argument("--output", required=True, help="Output directory path")
    p.add_argument("--demo", action="store_true", help="Run with bundled demo benchmark dataset")
    p.add_argument("--no-enrichment", action="store_true", help="Skip Enrichr pathway enrichment step")
    args = p.parse_args()

    here = Path(__file__).resolve().parent
    if args.demo:
        counts_path = here / "data" / "demo_counts.csv"
        meta_path = here / "data" / "demo_metadata.csv"
        formula = "~ batch + condition"
        contrast = "condition,treated,control"
        print("Running in DEMO mode with bundled benchmark dataset ...")
    else:
        if not args.counts or not args.metadata:
            p.error("Provide --counts and --metadata, or run with --demo")
        counts_path = Path(args.counts)
        meta_path = Path(args.metadata)
        formula = args.formula
        contrast = args.contrast

    result = run_pipeline(
        counts_path=counts_path,
        meta_path=meta_path,
        formula=formula,
        contrast=contrast,
        out_dir=Path(args.output),
        min_count=args.min_count,
        min_samples=args.min_samples,
        run_enrichment=not args.no_enrichment,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
