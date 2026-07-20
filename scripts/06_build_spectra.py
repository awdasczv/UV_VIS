"""
06_build_spectra.py
-------------------
TD-DFT 결과 -> CSV + 연속 스펙트럼 + 비교 그래프 (요구사항 11~14번).

산출물
  results/transitions.csv        전이 상세
  results/conformer_energies.csv 컨포머 에너지·볼츠만 가중치
  results/spectra_all.csv        파장별 연속 스펙트럼 (모든 시리즈)
  results/comparison.png         비교 그래프
  results/calc_settings.json     계산 설정 스냅샷

실행:  .\scripts\run.ps1 scripts\06_build_spectra.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import (CALCULATIONS, CONFORMERS, INPUTS, RESULTS,
                       gaussian_spectrum, load_checkpoint, save_checkpoint)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- dataviz 스킬의 검증된 카테고리 팔레트 (light mode, 슬롯 순서 고정) ---
PALETTE = ["#2a78d6", "#008300", "#e87ba4", "#eda100",
           "#1baf7a", "#eb6834", "#4a3aa7", "#e34948"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

TAUT_LABEL = {"enolA": "에놀 A (OH = tBu쪽)",
              "enolB": "에놀 B (OH = OMe쪽)",
              "diketo": "디케토"}


# ------------------------------------------------------------------ 로드
def load_results() -> list[dict]:
    p = CALCULATIONS / "02_tddft" / "all_results.json"
    d = load_checkpoint(p)
    if not d:
        print(f"TD-DFT 결과가 없습니다: {p}")
        return []
    return [r for r in d["results"] if r.get("ok")]


def load_config() -> dict:
    return json.loads((INPUTS / "calc_config.json").read_text(encoding="utf-8"))


# ------------------------------------------------------------------ CSV
def write_transitions_csv(results: list[dict]) -> pd.DataFrame:
    rows = []
    for r in results:
        for t in r["transitions"]:
            rows.append({
                "tautomer": r["tautomer"],
                "conformer_id": r["conf_id"],
                "level_id": r["level_id"],
                "functional": r["functional"],
                "basis": r["basis"],
                "solvent": r["solvent"],
                "geometry_source": r.get("geometry_source", ""),
                "rel_energy_kcalmol": r.get("rel_energy_kcalmol"),
                "boltzmann_weight": r.get("boltzmann_weight"),
                "state": t["state"],
                "excitation_energy_eV": t["energy_eV"],
                "wavelength_nm": t["wavelength_nm"],
                "oscillator_strength": t["osc_strength"],
                "major_orbital_transitions": t["orbital_transitions_str"],
            })
    df = pd.DataFrame(rows)
    out = RESULTS / "transitions.csv"
    df.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"  전이 상세 {len(df)} 행 -> {out.name}")
    return df


def write_conformer_csv() -> pd.DataFrame:
    rows = []
    for taut_dir in sorted(CONFORMERS.glob("*/ensemble.json")):
        ens = json.loads(taut_dir.read_text(encoding="utf-8"))
        sel_path = taut_dir.parent / "selected.json"
        selected = set()
        if sel_path.exists():
            selected = {s["conf_id"] for s in
                        json.loads(sel_path.read_text(encoding="utf-8"))["selected"]}
        for c in ens["conformers"]:
            rows.append({
                "tautomer": ens["tautomer"],
                "conformer_id": c["conf_id"],
                "method": ens["level"],
                "solvent_model": ens.get("solvent_xtb_alpb"),
                "energy_hartree": c["energy_hartree"],
                "rel_energy_kcalmol": c["rel_energy_kcalmol"],
                "boltzmann_weight_298K": c["boltzmann_weight"],
                "selected_for_dft": c["conf_id"] in selected,
            })
    df = pd.DataFrame(rows)
    out = RESULTS / "conformer_energies.csv"
    df.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"  컨포머 {len(df)} 개 -> {out.name}")
    return df


# ------------------------------------------------------------- 스펙트럼
def series_spectrum(recs: list[dict], grid: np.ndarray, fwhm: float,
                    weighted: bool) -> np.ndarray:
    """
    여러 컨포머 결과를 하나의 스펙트럼으로 합친다.
      weighted=True  -> 볼츠만 가중 앙상블 평균
      weighted=False -> 가장 안정한 컨포머 하나만
    """
    if not recs:
        return np.zeros_like(grid)
    if not weighted:
        recs = [min(recs, key=lambda r: r.get("rel_energy_kcalmol") or 0.0)]
        weights = [1.0]
    else:
        w = np.array([r.get("boltzmann_weight") or 1.0 for r in recs], dtype=float)
        weights = (w / w.sum()).tolist()

    eps = np.zeros_like(grid)
    for r, w in zip(recs, weights):
        lam = [t["wavelength_nm"] for t in r["transitions"]]
        f = [t["osc_strength"] for t in r["transitions"]]
        eps += gaussian_spectrum(lam, f, grid, fwhm_ev=fwhm, weight=w)
    return eps


def group(results: list[dict], **kw) -> list[dict]:
    return [r for r in results
            if all(r.get(k) == v for k, v in kw.items())]


# ------------------------------------------------------------------ 그래프
def annotate_peak(ax, grid, y, color, label, ymax_frac=1.0):
    """곡선의 최대점에 직접 라벨을 단다 (팔레트 대비 WARN 대응 = relief rule)."""
    if y.max() <= 0:
        return
    i = int(np.argmax(y))
    ax.annotate(f"{label}\n{grid[i]:.0f} nm",
                xy=(grid[i], y[i]), xytext=(6, 4), textcoords="offset points",
                fontsize=8, color=color, fontweight="bold",
                ha="left", va="bottom")


def style_axes(ax, title: str, xlabel: bool = True):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, fontsize=10, color=INK, loc="left", pad=8)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=MUTED, labelsize=8)
    if xlabel:
        ax.set_xlabel("파장 (nm)", fontsize=9, color=INK2)
    ax.set_ylabel("몰흡광계수 ε (L mol⁻¹ cm⁻¹)", fontsize=9, color=INK2)


def mark_experimental(ax, markers: dict):
    for name, lam in markers.items():
        ax.axvline(lam, color=MUTED, linestyle=(0, (4, 3)), linewidth=1.2, zorder=1)
        ax.annotate(f"실험 {lam:.0f} nm", xy=(lam, ax.get_ylim()[1]),
                    xytext=(3, -12), textcoords="offset points",
                    fontsize=7.5, color=INK2, rotation=90, va="top")


def make_figure(results: list[dict], cfg: dict, grid: np.ndarray,
                fwhm_list: list[float], out: Path) -> None:
    markers = cfg["experimental_markers_nm"]
    base_fwhm = fwhm_list[len(fwhm_list) // 2]
    levels = sorted({r["level_id"] for r in results})
    main_level = levels[0] if levels else None
    solvents = sorted({r["solvent"] for r in results})
    solv_pref = "ethanol" if "ethanol" in solvents else solvents[0]

    plt.rcParams["font.family"] = ["Malgun Gothic", "sans-serif"]
    plt.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor=SURFACE)
    fig.suptitle("아보벤존 TD-DFT UV–Vis 스펙트럼 비교",
                 fontsize=14, color=INK, x=0.02, ha="left", y=0.98)

    # --- 패널 A: 토토머 비교 (용매 적용, 앙상블 평균) ---
    ax = axes[0, 0]
    for i, taut in enumerate(["enolA", "enolB", "diketo"]):
        recs = group(results, tautomer=taut, level_id=main_level, solvent=solv_pref)
        if not recs:
            continue
        y = series_spectrum(recs, grid, base_fwhm, weighted=True)
        c = PALETTE[i]
        ax.plot(grid, y, color=c, linewidth=2.0, label=TAUT_LABEL[taut])
        annotate_peak(ax, grid, y, c, TAUT_LABEL[taut].split(" (")[0])
    style_axes(ax, f"A. 토토머 비교  ({main_level}, {solv_pref}, 앙상블 평균, "
                   f"FWHM {base_fwhm} eV)")
    ax.legend(fontsize=8, frameon=False, labelcolor=INK2)
    mark_experimental(ax, markers)

    # --- 패널 B: 용매 효과 (enolA) ---
    ax = axes[0, 1]
    for i, solv in enumerate(solvents):
        recs = group(results, tautomer="enolA", level_id=main_level, solvent=solv)
        if not recs:
            continue
        y = series_spectrum(recs, grid, base_fwhm, weighted=True)
        c = PALETTE[i]
        lbl = "용매 미적용 (기체상)" if solv == "none" else f"PCM ({solv})"
        ax.plot(grid, y, color=c, linewidth=2.0, label=lbl)
        annotate_peak(ax, grid, y, c, lbl)
    style_axes(ax, f"B. 용매 모델의 영향  (에놀 A, {main_level}, 앙상블 평균)")
    ax.legend(fontsize=8, frameon=False, labelcolor=INK2)
    mark_experimental(ax, markers)

    # --- 패널 C: 단일 최저구조 vs 앙상블 평균 ---
    ax = axes[1, 0]
    recs = group(results, tautomer="enolA", level_id=main_level, solvent=solv_pref)
    if recs:
        for i, (weighted, lbl) in enumerate([(False, "최저에너지 단일 구조"),
                                             (True, "컨포머 앙상블 평균")]):
            y = series_spectrum(recs, grid, base_fwhm, weighted=weighted)
            c = PALETTE[i]
            ax.plot(grid, y, color=c, linewidth=2.0,
                    linestyle="-" if weighted else (0, (5, 2)), label=lbl)
            annotate_peak(ax, grid, y, c, lbl)
    style_axes(ax, f"C. 컨포머 평균의 영향  (에놀 A, {main_level}, {solv_pref})")
    ax.legend(fontsize=8, frameon=False, labelcolor=INK2)
    mark_experimental(ax, markers)

    # --- 패널 D: 선폭(broadening) 비교 ---
    ax = axes[1, 1]
    # 선폭은 순서가 있는 값이므로 단일 색상 ramp(파랑 계열) 사용
    blue_ramp = ["#86b6ef", "#2a78d6", "#0d366b"]
    recs = group(results, tautomer="enolA", level_id=main_level, solvent=solv_pref)
    if recs:
        for i, fw in enumerate(fwhm_list):
            y = series_spectrum(recs, grid, fw, weighted=True)
            c = blue_ramp[min(i, len(blue_ramp) - 1)]
            ax.plot(grid, y, color=c, linewidth=2.0, label=f"FWHM {fw:.2f} eV")
            annotate_peak(ax, grid, y, c, f"{fw:.2f} eV")
    style_axes(ax, f"D. 가우시안 선폭 비교  (에놀 A, {main_level}, {solv_pref})")
    ax.legend(fontsize=8, frameon=False, labelcolor=INK2)
    mark_experimental(ax, markers)

    for ax in axes.ravel():
        ax.set_xlim(grid.min(), grid.max())

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"  그래프 -> {out.name}")


# ------------------------------------------------------------------ main
def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    cfg = load_config()
    results = load_results()
    if not results:
        return 1
    print(f"TD-DFT 성공 결과 {len(results)} 건")

    sp = cfg["spectrum"]
    grid = np.arange(sp["grid_nm"]["min"], sp["grid_nm"]["max"] + 1e-9,
                     sp["grid_nm"]["step"])
    fwhm_list = sp["broadening_fwhm_eV"]

    print("\nCSV 작성:")
    write_transitions_csv(results)
    write_conformer_csv()

    # 모든 시리즈의 연속 스펙트럼
    print("  연속 스펙트럼 계산 중...")
    spec = {"wavelength_nm": grid}
    for level in sorted({r["level_id"] for r in results}):
        for solv in sorted({r["solvent"] for r in results}):
            for taut in sorted({r["tautomer"] for r in results}):
                recs = group(results, tautomer=taut, level_id=level, solvent=solv)
                if not recs:
                    continue
                for fw in fwhm_list:
                    for weighted, wtag in [(True, "ensemble"), (False, "lowest")]:
                        key = f"{taut}|{level}|{solv}|fwhm{fw:.2f}|{wtag}"
                        spec[key] = series_spectrum(recs, grid, fw, weighted)
    df = pd.DataFrame(spec)
    out = RESULTS / "spectra_all.csv"
    df.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"  스펙트럼 {df.shape[1]-1} 시리즈 x {len(grid)} 점 -> {out.name}")

    make_figure(results, cfg, grid, fwhm_list, RESULTS / "comparison.png")

    save_checkpoint(RESULTS / "calc_settings.json", {
        "config": cfg,
        "n_tddft_results": len(results),
        "levels": sorted({r["level_id"] for r in results}),
        "solvents": sorted({r["solvent"] for r in results}),
        "tautomers": sorted({r["tautomer"] for r in results}),
        "spectrum_grid_nm": {"min": float(grid.min()), "max": float(grid.max()),
                             "step": sp["grid_nm"]["step"], "n_points": len(grid)},
        "broadening_fwhm_eV": fwhm_list,
    })
    print(f"  설정 스냅샷 -> calc_settings.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
