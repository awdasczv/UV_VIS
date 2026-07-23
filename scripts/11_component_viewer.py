"""
11_component_viewer.py
----------------------
여러 성분(분자)의 계산된 흡수 스펙트럼을 한 화면에 겹쳐 보는 인터랙티브 뷰어를 만든다.

왜 필요한가 (최종 목표)
  복합 선크림 제형의 흡광 스펙트럼 ~= 성분별 스펙트럼의 농도 가중합(+ML 잔차).
  그러려면 먼저 성분별 스펙트럼이 있어야 하고, "이들이 어떻게 더해지는가"를
  눈으로 봐야 한다. 이 뷰어는 성분 스펙트럼을 겹쳐 그리고, 가중치 슬라이더로
  Beer-Lambert 합(sum_i w_i * eps_i)을 실시간으로 보여준다.

무엇을 그리는가
  - 각 분자의 '주 흡수 종'(config.principal_species) 최고수준 TD-DFT 스펙트럼
  - 절대 몰흡광계수 eps(L mol^-1 cm^-1) 로 broadening (가중합이 물리적으로 의미)
  - UVB(280-315)/UVA(315-400) 밴드, 분자별 실험 lambda_max 마커
  - 가중치 슬라이더로 성분 합 곡선

분자 선택
  인자로 이름을 주거나(예: avobenzone dhhb), 생략하면 molecules/ 아래에서
  TD-DFT 결과가 있는 분자를 자동 검색한다. (UV_MOLECULE 과 무관하게 동작)

실행:  .\scripts\run.ps1 scripts\11_component_viewer.py [avobenzone dhhb ...]
산출:  results/component_viewer.html , results/component_spectra.csv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import MOLECULES, ROOT, gaussian_spectrum

GRID = np.arange(200.0, 500.0 + 1e-9, 1.0)
FWHM_EV = 0.30

# 성분 색상 팔레트 (dataviz 검증 팔레트)
PALETTE = ["#2a78d6", "#e8710a", "#1a9e5f", "#c0392b", "#7b48c8",
           "#0f9bd1", "#d64b9a", "#8a6d00"]

# 최고수준 선택 우선순위
BASIS_RANK = {"6-31+G(d)": 4, "def2-TZVP": 3, "def2-SVPD": 3,
              "def2-SVP": 2, "6-31G(d)": 1}
FUNC_RANK = {"B3LYP": 3, "PBE0": 2, "CAM-B3LYP": 1, "WB97X-D": 1}


def load(p: Path):
    try:
        return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:                                   # noqa: BLE001
        return None


def result_score(rec: dict) -> tuple:
    """최고수준 대표를 고르는 정렬 키 (클수록 우선)."""
    return (
        1 if rec.get("solvent") == "ethanol" else 0,
        1 if rec.get("geometry_source") == "DFT" else 0,
        BASIS_RANK.get(rec.get("basis"), 0),
        FUNC_RANK.get(rec.get("functional"), 0),
        float(rec.get("boltzmann_weight") or 0.0),
    )


def principal_species(cfg: dict, results: list[dict]) -> str | None:
    if cfg and cfg.get("principal_species"):
        return cfg["principal_species"]
    # 폴백: 가장 밝은 전이(f 최대)를 가진 종
    best, best_f = None, -1.0
    for r in results:
        b = r.get("brightest") or {}
        if float(b.get("osc_strength") or 0) > best_f:
            best_f = float(b.get("osc_strength") or 0)
            best = r.get("tautomer")
    return best


def collect_molecule(mol_dir: Path) -> dict | None:
    """한 분자의 주 흡수 종 최고수준 스펙트럼을 모은다."""
    cfg = load(mol_dir / "config.json") or {}
    exp = load(mol_dir / "experimental_reference.json")
    td_root = mol_dir / "calculations" / "02_tddft_orca"
    if not td_root.exists():
        return None

    results = [r for r in (load(p) for p in td_root.rglob("result.json"))
               if r and r.get("ok") and r.get("transitions")]
    if not results:
        return None

    species = principal_species(cfg, results)
    cand = [r for r in results if r.get("tautomer") == species]
    if not cand:
        return None
    best = max(cand, key=result_score)

    lam = [t["wavelength_nm"] for t in best["transitions"]]
    osc = [t["osc_strength"] for t in best["transitions"]]
    eps = gaussian_spectrum(lam, osc, GRID, fwhm_ev=FWHM_EV)
    b = best.get("brightest") or {}

    # 실험 lambda_max (스키마가 분자마다 조금 달라 관대하게 탐색)
    exp_nm, exp_conf = None, None
    if exp:
        for band in exp.get("experimental", {}).values():
            pt = (band or {}).get("primary_target") if isinstance(band, dict) else None
            if pt and pt.get("lambda_max_nm"):
                exp_nm = pt["lambda_max_nm"]
                exp_conf = pt.get("confidence") or pt.get("ref")
                break

    m = cfg.get("molecule", {})
    return {
        "molecule": mol_dir.name,
        "display_name": m.get("common_name", mol_dir.name),
        "species": species,
        "level": f"{best.get('functional')}/{best.get('basis')}",
        "solvent": best.get("solvent"),
        "geometry_source": best.get("geometry_source"),
        "calc_lambda_max_nm": round(b.get("wavelength_nm"), 1) if b.get("wavelength_nm") else None,
        "calc_osc": round(b.get("osc_strength"), 3) if b.get("osc_strength") else None,
        "exp_lambda_max_nm": exp_nm,
        "exp_confidence": exp_conf,
        "eps": [round(float(v), 1) for v in eps],
    }


def write_csv(series: list[dict], out: Path) -> None:
    cols = ["wavelength_nm"] + [s["molecule"] for s in series]
    lines = [",".join(cols)]
    for i, nm in enumerate(GRID):
        row = [f"{nm:.1f}"] + [f"{s['eps'][i]:.1f}" for s in series]
        lines.append(",".join(row))
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_html(series: list[dict]) -> str:
    data = {
        "grid": [float(x) for x in GRID],
        "series": [{k: s[k] for k in
                    ("molecule", "display_name", "species", "level", "solvent",
                     "geometry_source", "calc_lambda_max_nm", "calc_osc",
                     "exp_lambda_max_nm", "exp_confidence", "eps")}
                   for s in series],
        "palette": PALETTE,
    }
    return HTML_TEMPLATE.replace("/*__DATA__*/", json.dumps(data, ensure_ascii=False))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("molecules", nargs="*",
                    help="그릴 분자 이름들. 생략하면 결과가 있는 분자 자동검색.")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.molecules:
        dirs = [MOLECULES / m for m in args.molecules]
    else:
        dirs = sorted(p for p in MOLECULES.iterdir()
                      if (p / "calculations" / "02_tddft_orca").exists())

    series = []
    for d in dirs:
        s = collect_molecule(d)
        if s:
            series.append(s)
            print(f"  {s['molecule']:12s} {s['species']:8s} {s['level']:18s} "
                  f"{s['solvent']:8s} calc lambda_max={s['calc_lambda_max_nm']} nm "
                  f"(exp {s['exp_lambda_max_nm']})")
        else:
            print(f"  {d.name}: TD-DFT 결과 없음 (건너뜀)")

    if not series:
        print("그릴 성분이 없습니다. 먼저 05b_tddft_orca.py 를 실행하세요.")
        return 1

    out_dir = ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    html = out_dir / (args.out or "component_viewer.html")
    csv = out_dir / "component_spectra.csv"
    html.write_text(build_html(series), encoding="utf-8")
    write_csv(series, csv)
    print(f"\n뷰어 -> {html.relative_to(ROOT)}")
    print(f"CSV  -> {csv.relative_to(ROOT)}")
    print(f"성분 {len(series)} 개.")
    return 0


HTML_TEMPLATE = r"""<!doctype html>
<html lang="ko"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>선크림 성분 흡수 스펙트럼</title>
<style>
  :root{
    --bg:#f7f8fa; --panel:#ffffff; --ink:#1a1f2b; --muted:#5b6472;
    --grid:#e3e7ee; --line:#c7ccd6; --uva:#efe6fb; --uvb:#e6f0fb; --accent:#2a78d6;
  }
  @media (prefers-color-scheme: dark){
    :root{ --bg:#0f1319; --panel:#161b24; --ink:#e7ebf2; --muted:#95a0b2;
           --grid:#232a36; --line:#39424f; --uva:#241a3a; --uvb:#152437; --accent:#5aa0ee; }
  }
  :root[data-theme=light]{ --bg:#f7f8fa; --panel:#fff; --ink:#1a1f2b; --muted:#5b6472;
    --grid:#e3e7ee; --line:#c7ccd6; --uva:#efe6fb; --uvb:#e6f0fb; --accent:#2a78d6; }
  :root[data-theme=dark]{ --bg:#0f1319; --panel:#161b24; --ink:#e7ebf2; --muted:#95a0b2;
    --grid:#232a36; --line:#39424f; --uva:#241a3a; --uvb:#152437; --accent:#5aa0ee; }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--ink);
    font-family:-apple-system,"Segoe UI",Roboto,"Noto Sans KR",sans-serif;font-size:14px}
  .wrap{max-width:1080px;margin:0 auto;padding:22px 18px 60px}
  h1{font-size:20px;margin:0 0 2px;letter-spacing:-.01em}
  .sub{color:var(--muted);font-size:13px;margin-bottom:18px}
  .card{background:var(--panel);border:1px solid var(--grid);border-radius:12px;
    padding:16px;margin-bottom:16px}
  .chartwrap{position:relative}
  svg{width:100%;height:auto;display:block;touch-action:none}
  .ctl{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin-bottom:12px}
  .ctl label{display:flex;gap:6px;align-items:center;color:var(--muted);cursor:pointer}
  .pill{font-variant-numeric:tabular-nums}
  table{border-collapse:collapse;width:100%;font-size:13px}
  th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--grid)}
  th{color:var(--muted);font-weight:600}
  td.num{text-align:right;font-variant-numeric:tabular-nums}
  .sw{width:11px;height:11px;border-radius:3px;display:inline-block;vertical-align:middle;margin-right:7px}
  .comp{display:flex;flex-direction:column;gap:9px}
  .comprow{display:grid;grid-template-columns:150px 1fr 46px;gap:10px;align-items:center}
  .comprow input[type=range]{width:100%}
  .comprow .nm{white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  button{background:transparent;border:1px solid var(--line);color:var(--ink);
    border-radius:8px;padding:5px 11px;cursor:pointer;font-size:13px}
  button:hover{border-color:var(--accent)}
  .tip{position:absolute;pointer-events:none;background:var(--panel);border:1px solid var(--line);
    border-radius:8px;padding:8px 10px;font-size:12px;opacity:0;transition:opacity .08s;
    box-shadow:0 6px 20px rgba(0,0,0,.14);min-width:150px}
  .tip b{font-variant-numeric:tabular-nums}
  .muted{color:var(--muted)}
</style></head>
<body><div class="wrap">
  <div style="display:flex;justify-content:space-between;align-items:flex-start;gap:12px">
    <div>
      <h1>선크림 성분 흡수 스펙트럼</h1>
      <div class="sub">계산된 성분별 스펙트럼과 이들의 Beer-Lambert 가중합.
        최종 목표(복합 제형 예측)를 향한 성분 라이브러리.</div>
    </div>
    <button id="theme">🌓 테마</button>
  </div>

  <div class="card">
    <div class="ctl">
      <label><input type="checkbox" id="norm"> 정규화(피크=1)</label>
      <label><input type="checkbox" id="showsum" checked> 성분 합 표시</label>
      <label><input type="checkbox" id="showexp" checked> 실험 λmax 마커</label>
      <span class="muted pill" id="readout"></span>
    </div>
    <div class="chartwrap">
      <svg id="chart" viewBox="0 0 900 460" preserveAspectRatio="xMidYMid meet"></svg>
      <div class="tip" id="tip"></div>
    </div>
  </div>

  <div class="card">
    <div class="sub" style="margin:0 0 10px">성분 가중치 (Beer-Lambert 합 = Σ wᵢ·εᵢ).
      농도를 대신하는 상대 가중치다 — 실제 제형 비율을 넣어보며 합 스펙트럼을 확인.</div>
    <div class="comp" id="comp"></div>
  </div>

  <div class="card">
    <table><thead><tr>
      <th>성분</th><th>주 흡수 종</th><th>이론수준</th>
      <th class="num">계산 λmax</th><th class="num">실험 λmax</th><th class="num">오차</th>
    </tr></thead><tbody id="tbody"></tbody></table>
    <div class="sub" style="margin-top:10px">FWHM 0.30 eV 가우시안 broadening,
      절대 몰흡광계수 ε (L·mol⁻¹·cm⁻¹). 용매·구조는 각 성분 최고수준(에탄올/DFT 구조).</div>
  </div>
</div>
<script>
const DATA = /*__DATA__*/;
const G = DATA.grid, S = DATA.series, PAL = DATA.palette;
S.forEach((s,i)=>{ s.color = PAL[i%PAL.length]; s.weight = 1; s.on = true; });

const NM0=G[0], NM1=G[G.length-1];
const svg=document.getElementById('chart'), tip=document.getElementById('tip');
const W=900,H=460, mL=58,mR=16,mT=16,mB=42;
const px=nm=>mL+(nm-NM0)/(NM1-NM0)*(W-mL-mR);
const nmAt=x=>NM0+(x-mL)/(W-mL-mR)*(NM1-NM0);
let curMaxY=1;

function weightedSum(){ const y=G.map(()=>0);
  S.forEach(s=>{ if(!s.on) return; for(let i=0;i<G.length;i++) y[i]+=s.weight*s.eps[i]; });
  return y; }
function seriesY(s){ return s.eps; }

function draw(){
  const norm=document.getElementById('norm').checked;
  const showsum=document.getElementById('showsum').checked;
  const showexp=document.getElementById('showexp').checked;
  const active=S.filter(s=>s.on);
  let maxY=0;
  const disp=[];
  active.forEach(s=>{ let y=seriesY(s);
    if(norm){ const mx=Math.max(...y)||1; y=y.map(v=>v/mx); }
    disp.push({s,y}); maxY=Math.max(maxY,...y); });
  let sumY=null;
  if(showsum){ sumY=weightedSum();
    if(norm){ const mx=Math.max(...sumY)||1; sumY=sumY.map(v=>v/mx); }
    maxY=Math.max(maxY,...sumY); }
  maxY=maxY*1.08||1; curMaxY=maxY;
  const py=v=>H-mB-(v/maxY)*(H-mT-mB);

  let el='';
  // UVB / UVA 밴드
  el+=`<rect x="${px(280)}" y="${mT}" width="${px(315)-px(280)}" height="${H-mT-mB}" fill="var(--uvb)"/>`;
  el+=`<rect x="${px(315)}" y="${mT}" width="${px(400)-px(315)}" height="${H-mT-mB}" fill="var(--uva)"/>`;
  el+=`<text x="${(px(280)+px(315))/2}" y="${mT+13}" fill="var(--muted)" font-size="11" text-anchor="middle">UVB</text>`;
  el+=`<text x="${(px(315)+px(400))/2}" y="${mT+13}" fill="var(--muted)" font-size="11" text-anchor="middle">UVA</text>`;
  // x 그리드/눈금
  for(let nm=250;nm<=450;nm+=50){ el+=`<line x1="${px(nm)}" y1="${mT}" x2="${px(nm)}" y2="${H-mB}" stroke="var(--grid)"/>`;
    el+=`<text x="${px(nm)}" y="${H-mB+16}" fill="var(--muted)" font-size="11" text-anchor="middle">${nm}</text>`; }
  el+=`<text x="${(mL+W-mR)/2}" y="${H-6}" fill="var(--muted)" font-size="12" text-anchor="middle">파장 (nm)</text>`;
  // y 눈금
  for(let k=0;k<=4;k++){ const v=maxY*k/4, yy=py(v);
    el+=`<line x1="${mL}" y1="${yy}" x2="${W-mR}" y2="${yy}" stroke="var(--grid)"/>`;
    el+=`<text x="${mL-8}" y="${yy+4}" fill="var(--muted)" font-size="11" text-anchor="end">${norm?v.toFixed(2):Math.round(v)}</text>`; }
  el+=`<text transform="translate(15,${(mT+H-mB)/2}) rotate(-90)" fill="var(--muted)" font-size="12" text-anchor="middle">${norm?'정규화 흡광':'ε (L·mol⁻¹·cm⁻¹)'}</text>`;

  const path=y=>{ let d=''; for(let i=0;i<G.length;i++){ d+=(i?'L':'M')+px(G[i]).toFixed(1)+' '+py(y[i]).toFixed(1);} return d; };
  disp.forEach(({s,y})=>{ el+=`<path d="${path(y)}" fill="none" stroke="${s.color}" stroke-width="2"/>`; });
  if(sumY){ el+=`<path d="${path(sumY)}" fill="none" stroke="var(--ink)" stroke-width="2.4" stroke-dasharray="2 4"/>`; }
  // 실험 마커
  if(showexp){ active.forEach(s=>{ if(!s.exp_lambda_max_nm) return;
    const x=px(s.exp_lambda_max_nm);
    el+=`<line x1="${x}" y1="${mT}" x2="${x}" y2="${H-mB}" stroke="${s.color}" stroke-width="1.3" stroke-dasharray="4 3" opacity=".8"/>`;
    el+=`<text x="${x}" y="${mT+28}" fill="${s.color}" font-size="10" text-anchor="middle">${s.exp_lambda_max_nm}</text>`; }); }
  el+=`<line id="cross" x1="0" y1="${mT}" x2="0" y2="${H-mB}" stroke="var(--line)" opacity="0"/>`;
  svg.innerHTML=el;
}

function renderComp(){
  const c=document.getElementById('comp'); c.innerHTML='';
  S.forEach((s,i)=>{ const row=document.createElement('div'); row.className='comprow';
    row.innerHTML=`<div class="nm"><span class="sw" style="background:${s.color}"></span>
      <label style="display:inline"><input type="checkbox" ${s.on?'checked':''} data-on="${i}"> ${s.display_name}</label></div>
      <input type="range" min="0" max="2" step="0.05" value="${s.weight}" data-w="${i}">
      <span class="pill" data-wv="${i}">${s.weight.toFixed(2)}</span>`;
    c.appendChild(row); });
  c.querySelectorAll('[data-w]').forEach(r=>r.addEventListener('input',e=>{
    const i=+e.target.dataset.w; S[i].weight=+e.target.value;
    document.querySelector(`[data-wv="${i}"]`).textContent=S[i].weight.toFixed(2); draw(); }));
  c.querySelectorAll('[data-on]').forEach(r=>r.addEventListener('change',e=>{
    S[+e.target.dataset.on].on=e.target.checked; draw(); }));
}

function renderTable(){
  const tb=document.getElementById('tbody'); tb.innerHTML='';
  S.forEach(s=>{ const err=(s.calc_lambda_max_nm&&s.exp_lambda_max_nm)?
      (s.calc_lambda_max_nm-s.exp_lambda_max_nm):null;
    const tr=document.createElement('tr');
    tr.innerHTML=`<td><span class="sw" style="background:${s.color}"></span>${s.display_name}</td>
      <td>${s.species}</td><td class="muted">${s.level} · ${s.solvent} · ${s.geometry_source}</td>
      <td class="num">${s.calc_lambda_max_nm??'–'}</td>
      <td class="num">${s.exp_lambda_max_nm??'–'}${s.exp_confidence==='approximate'?'*':''}</td>
      <td class="num">${err==null?'–':(err>0?'+':'')+err.toFixed(1)}</td>`;
    tb.appendChild(tr); });
}

svg.addEventListener('pointermove',e=>{
  const r=svg.getBoundingClientRect(); const x=(e.clientX-r.left)/r.width*W;
  if(x<mL||x>W-mR){ tip.style.opacity=0; return; }
  const nm=nmAt(x); const idx=Math.round(nm-NM0);
  const cross=document.getElementById('cross'); if(cross){cross.setAttribute('x1',x);cross.setAttribute('x2',x);cross.setAttribute('opacity','.5');}
  const norm=document.getElementById('norm').checked;
  let rows=`<b>${nm.toFixed(0)} nm</b>`;
  S.filter(s=>s.on).forEach(s=>{ let v=s.eps[idx]; rows+=`<br><span class="sw" style="background:${s.color}"></span>${Math.round(v)}`; });
  tip.innerHTML=rows; tip.style.opacity=1;
  let tx=e.clientX-r.left+14; if(tx>r.width-140) tx=e.clientX-r.left-150;
  tip.style.left=tx+'px'; tip.style.top=(e.clientY-r.top+12)+'px';
  document.getElementById('readout').textContent=`${nm.toFixed(0)} nm 에서 읽는 중`;
});
svg.addEventListener('pointerleave',()=>{ tip.style.opacity=0;
  const cross=document.getElementById('cross'); if(cross)cross.setAttribute('opacity','0'); });

['norm','showsum','showexp'].forEach(id=>document.getElementById(id).addEventListener('change',draw));
document.getElementById('theme').addEventListener('click',()=>{
  const cur=document.documentElement.getAttribute('data-theme');
  const next=cur==='dark'?'light':(cur==='light'?'dark':(matchMedia('(prefers-color-scheme: dark)').matches?'light':'dark'));
  document.documentElement.setAttribute('data-theme',next); draw(); });

renderComp(); renderTable(); draw();
</script></body></html>
"""


if __name__ == "__main__":
    sys.exit(main())
