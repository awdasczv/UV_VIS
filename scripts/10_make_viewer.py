"""
10_make_viewer.py
-----------------
계산된 스펙트럼을 브라우저에서 인터랙티브하게 볼 수 있는 HTML 을 만든다.
결과: results/spectrum_viewer.html  (자체 완결 — 더블클릭하면 브라우저에서 열림)

무엇을 담나
  - 세 토토머(에놀 A/B, 디케토)의 최고 수준 스펙트럼
  - 비교 시리즈: 용매 유/무, xTB vs DFT 구조, B3LYP vs CAM-B3LYP
  - UVB(280~315) / UVA(315~400) 구간 음영
  - 실험 최대흡수파장 마커
  - 곡선 토글, 프리셋, 마우스 크로스헤어(파장별 흡광도 읽기)

데이터는 HTML 안에 JSON 으로 심는다(외부 파일 의존 없음).
다른 분자를 추가한 뒤 다시 실행하면 갱신된다.

실행:  .\scripts\run.ps1 scripts\10_make_viewer.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import CALCULATIONS, RESULTS, gaussian_spectrum, load_checkpoint

GRID = np.arange(200, 501, 1.0)     # 200~500 nm, 1 nm 간격
FWHM = 0.30                          # eV

# 보여줄 시리즈 정의: (id, 라벨, 그룹, 색슬롯, 토토머, 수준, 용매, 구조)
# 색슬롯: dataviz 검증 팔레트 순서
SERIES = [
    ("enolA_best",  "에놀 A (최종)",       "main", 0, "enolA",  "b3lyp_631+gd", "ethanol", "dftopt22"),
    ("enolB_best",  "에놀 B (최종)",       "main", 1, "enolB",  "b3lyp_631+gd", "ethanol", "dftopt22"),
    ("diketo_best", "디케토 (최종)",       "main", 2, "diketo", "b3lyp_631+gd", "ethanol", "dftopt22"),

    ("enolA_gas",   "에놀 A · 기체상",      "solvent", 3, "enolA", "b3lyp_631+gd", "none",    "dftopt22"),
    ("enolA_etoh",  "에놀 A · 에탄올",      "solvent", 0, "enolA", "b3lyp_631+gd", "ethanol", "dftopt22"),

    ("enolA_xtb",   "에놀 A · xTB 구조",    "geom", 5, "enolA", "b3lyp_631+gd", "ethanol", "xtb"),
    ("enolA_dft",   "에놀 A · DFT 구조",    "geom", 0, "enolA", "b3lyp_631+gd", "ethanol", "dftopt22"),

    ("enolA_b3lyp", "에놀 A · B3LYP",       "func", 0, "enolA", "b3lyp_def2svp",   "ethanol", "dftopt"),
    ("enolA_cam",   "에놀 A · CAM-B3LYP",   "func", 1, "enolA", "camb3lyp_def2svp","ethanol", "dftopt22"),
]

PALETTE = ["#2a78d6", "#008300", "#e87ba4", "#eda100",
           "#1baf7a", "#eb6834", "#4a3aa7", "#e34948"]


def collect() -> list[dict]:
    recs = []
    root = CALCULATIONS / "02_tddft_orca"
    for p in sorted(root.rglob("result.json")):
        d = load_checkpoint(p)
        if not d or not d.get("ok"):
            continue
        parts = p.relative_to(root).parts
        d["geom_label"] = parts[4] if len(parts) >= 6 else "xtb"
        recs.append(d)
    return recs


def find(recs, taut, level, solvent, geom):
    cand = [r for r in recs if r["tautomer"] == taut and r["level_id"] == level
            and r["solvent"] == solvent and r["geom_label"] == geom]
    if not cand:
        return None
    return min(cand, key=lambda r: r.get("rel_energy_kcalmol") or 0.0)


def main() -> int:
    recs = collect()
    series_out = []
    for sid, label, group, slot, taut, level, solv, geom in SERIES:
        r = find(recs, taut, level, solv, geom)
        if not r:
            print(f"  [건너뜀] {sid}: 데이터 없음 ({taut}/{level}/{solv}/{geom})")
            continue
        lam = [t["wavelength_nm"] for t in r["transitions"]]
        f = [t["osc_strength"] for t in r["transitions"]]
        eps = gaussian_spectrum(lam, f, GRID, fwhm_ev=FWHM)
        i = int(np.argmax(eps))
        # 스틱(막대) 스펙트럼: 진동자 세기 있는 전이만
        sticks = [[round(t["wavelength_nm"], 1), round(t["osc_strength"], 4)]
                  for t in r["transitions"] if t["osc_strength"] > 0.001]
        series_out.append({
            "id": sid, "label": label, "group": group,
            "color": PALETTE[slot % len(PALETTE)],
            "eps": [int(round(v)) for v in eps],
            "peakNm": round(float(GRID[i]), 1),
            "peakEps": int(round(float(eps[i]))),
            "brightNm": round(r["brightest"]["wavelength_nm"], 1),
            "f": round(r["brightest"]["osc_strength"], 3),
            "sticks": sticks,
        })

    data = {
        "grid": [int(x) for x in GRID],
        "fwhm_eV": FWHM,
        "series": series_out,
        "markers": [
            {"nm": 354.9, "label": "실험 에놀 355 nm"},
            {"nm": 265.0, "label": "실험 디케토 265 nm"},
        ],
        "bands": [
            {"from": 280, "to": 315, "label": "UVB", "color": "#eda100"},
            {"from": 315, "to": 400, "label": "UVA", "color": "#2a78d6"},
        ],
        "presets": {
            "토토머 비교": ["enolA_best", "enolB_best", "diketo_best"],
            "용매 효과 (에놀 A)": ["enolA_gas", "enolA_etoh"],
            "구조 효과 (에놀 A)": ["enolA_xtb", "enolA_dft"],
            "함수 효과 (에놀 A)": ["enolA_b3lyp", "enolA_cam"],
        },
    }

    html = HTML_TEMPLATE.replace("__DATA__", json.dumps(data, ensure_ascii=False))
    out = RESULTS / "spectrum_viewer.html"
    out.write_text(html, encoding="utf-8")
    print(f"뷰어 생성 -> {out}")
    print(f"  시리즈 {len(series_out)} 개, 격자 {len(GRID)} 점")
    print("  브라우저에서 열어보세요 (더블클릭 또는 start results\\spectrum_viewer.html)")
    return 0


HTML_TEMPLATE = r"""<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>아보벤존 UV-Vis 스펙트럼 뷰어</title>
<style>
  :root{
    --bg:#f7f8fa; --surface:#ffffff; --ink:#12151a; --ink2:#565d68;
    --muted:#8a919c; --grid:#e6e9ee; --axis:#c2c8d0; --line:#0b0b0b;
    --accent:#2a78d6; --shadow:0 1px 3px rgba(20,30,50,.08),0 8px 24px rgba(20,30,50,.06);
  }
  @media (prefers-color-scheme:dark){
    :root{ --bg:#0e1116; --surface:#161b22; --ink:#eef1f5; --ink2:#aab2bd;
      --muted:#6e7681; --grid:#232a33; --axis:#39414c; --line:#eef1f5;
      --shadow:0 1px 3px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.3); }
  }
  :root[data-theme="dark"]{ --bg:#0e1116; --surface:#161b22; --ink:#eef1f5; --ink2:#aab2bd;
    --muted:#6e7681; --grid:#232a33; --axis:#39414c; --line:#eef1f5;
    --shadow:0 1px 3px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.3); }
  :root[data-theme="light"]{ --bg:#f7f8fa; --surface:#ffffff; --ink:#12151a; --ink2:#565d68;
    --muted:#8a919c; --grid:#e6e9ee; --axis:#c2c8d0; --line:#0b0b0b; }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--ink);
    font-family:system-ui,-apple-system,"Segoe UI","Malgun Gothic",sans-serif;
    font-size:15px;line-height:1.5}
  .wrap{max-width:1080px;margin:0 auto;padding:28px 20px 60px}
  header h1{font-size:22px;margin:0 0 4px;letter-spacing:-.01em}
  header p{margin:0;color:var(--ink2);font-size:14px}
  .card{background:var(--surface);border-radius:14px;box-shadow:var(--shadow);
    padding:18px;margin-top:20px}
  .presets{display:flex;flex-wrap:wrap;gap:8px;margin:2px 0 4px}
  .presets button{font:inherit;font-size:13px;padding:7px 13px;border-radius:999px;
    border:1px solid var(--axis);background:transparent;color:var(--ink2);cursor:pointer;
    transition:all .15s}
  .presets button:hover{border-color:var(--accent);color:var(--accent)}
  .presets button.on{background:var(--accent);border-color:var(--accent);color:#fff}
  .chartbox{position:relative;margin-top:10px}
  svg{display:block;width:100%;height:auto;touch-action:none}
  .legend{display:flex;flex-wrap:wrap;gap:6px 16px;margin-top:14px}
  .legend label{display:inline-flex;align-items:center;gap:7px;font-size:13.5px;
    color:var(--ink2);cursor:pointer;user-select:none}
  .legend input{position:absolute;opacity:0;pointer-events:none}
  .swatch{width:22px;height:3px;border-radius:2px;flex:none;position:relative}
  .legend .off .swatch{opacity:.25}
  .legend .off span.txt{opacity:.4;text-decoration:line-through}
  .legend .dot{width:9px;height:9px;border-radius:50%;flex:none}
  .tip{position:absolute;pointer-events:none;background:var(--surface);
    border:1px solid var(--axis);border-radius:10px;box-shadow:var(--shadow);
    padding:9px 11px;font-size:12.5px;min-width:150px;opacity:0;transition:opacity .1s;
    font-variant-numeric:tabular-nums;z-index:5}
  .tip b{color:var(--ink)}
  .tip .row{display:flex;align-items:center;gap:7px;margin-top:3px;color:var(--ink2)}
  .tip .row .d{width:8px;height:8px;border-radius:50%;flex:none}
  .note{color:var(--muted);font-size:12.5px;margin-top:14px;line-height:1.6}
  .note code{background:var(--grid);padding:1px 5px;border-radius:5px;font-size:12px}
  .band-lab{font-size:11px;font-weight:600;letter-spacing:.03em}
  .toggle{float:right;font-size:12.5px;color:var(--ink2);cursor:pointer;
    border:1px solid var(--axis);border-radius:999px;padding:5px 11px;background:transparent}
</style>
</head>
<body>
<div class="wrap">
  <header>
    <button class="toggle" id="themeBtn">◐ 테마</button>
    <h1>아보벤존 UV–Vis 흡수 스펙트럼</h1>
    <p>분자 구조로부터 계산한 자외선 흡수 곡선. 곡선을 켜고/끄고, 마우스를 올려 값을 읽으세요.</p>
  </header>

  <div class="card">
    <div class="presets" id="presets"></div>
    <div class="chartbox" id="chartbox">
      <svg id="chart" viewBox="0 0 900 480" preserveAspectRatio="xMidYMid meet"></svg>
      <div class="tip" id="tip"></div>
    </div>
    <div class="legend" id="legend"></div>
    <div class="note">
      세로 점선 = 실험 최대흡수파장. 색 띠 = 자외선 차단 관심 구간
      (<b style="color:#c98500">UVB 280–315 nm</b>, <b style="color:#2a78d6">UVA 315–400 nm</b>).
      세로축은 몰흡광계수 ε(L mol⁻¹ cm⁻¹), 선폭 FWHM <code>0.30 eV</code> 가우시안.
      최종 수준: B3LYP/6-31+G(d) · DFT 최적화 구조 · CPCM(에탄올).
    </div>
  </div>
</div>

<script>
const DATA = __DATA__;
const NS="http://www.w3.org/2000/svg";
const svg=document.getElementById("chart");
const W=900,H=480,M={l:64,r:20,t:20,b:44};
const x0=M.l,x1=W-M.r,y0=H-M.b,y1=M.t;
const gmin=DATA.grid[0], gmax=DATA.grid[DATA.grid.length-1];
const visible=new Set(DATA.presets["토토머 비교"]);
let ymax=1;

function sx(nm){return x0+(nm-gmin)/(gmax-gmin)*(x1-x0);}
function sy(v){return y0-(v/ymax)*(y0-y1);}
function el(t,a){const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}

function computeYmax(){
  let m=1;
  for(const s of DATA.series){ if(!visible.has(s.id))continue;
    for(const v of s.eps) if(v>m)m=v; }
  ymax=m*1.08;
}
function niceStep(range,target){
  const raw=range/target, mag=Math.pow(10,Math.floor(Math.log10(raw)));
  const n=raw/mag; let s;
  if(n<1.5)s=1;else if(n<3)s=2;else if(n<7)s=5;else s=10;
  return s*mag;
}

function draw(){
  computeYmax();
  while(svg.firstChild)svg.removeChild(svg.firstChild);

  // 관심 구간 음영 (UVB/UVA)
  for(const b of DATA.bands){
    svg.appendChild(el("rect",{x:sx(b.from),y:y1,width:sx(b.to)-sx(b.from),
      height:y0-y1,fill:b.color,"fill-opacity":0.07}));
    const t=el("text",{x:(sx(b.from)+sx(b.to))/2,y:y1+14,"text-anchor":"middle",
      class:"band-lab",fill:b.color,"fill-opacity":0.9});
    t.textContent=b.label; svg.appendChild(t);
  }

  // 격자 + 축
  const xstep=50;
  for(let nm=Math.ceil(gmin/xstep)*xstep; nm<=gmax; nm+=xstep){
    svg.appendChild(el("line",{x1:sx(nm),y1:y1,x2:sx(nm),y2:y0,
      stroke:"var(--grid)","stroke-width":1}));
    const t=el("text",{x:sx(nm),y:y0+20,"text-anchor":"middle",
      "font-size":12,fill:"var(--muted)"}); t.textContent=nm; svg.appendChild(t);
  }
  const ystep=niceStep(ymax,5);
  for(let v=0; v<=ymax; v+=ystep){
    svg.appendChild(el("line",{x1:x0,y1:sy(v),x2:x1,y2:sy(v),
      stroke:"var(--grid)","stroke-width":1}));
    const t=el("text",{x:x0-8,y:sy(v)+4,"text-anchor":"end",
      "font-size":11,fill:"var(--muted)"});
    t.textContent=(v/1000).toFixed(0)+"k"; svg.appendChild(t);
  }
  // 축 라벨
  const yl=el("text",{x:16,y:(y0+y1)/2,"text-anchor":"middle","font-size":12,
    fill:"var(--ink2)",transform:`rotate(-90 16 ${(y0+y1)/2})`});
  yl.textContent="ε (L mol⁻¹ cm⁻¹)"; svg.appendChild(yl);
  const xl=el("text",{x:(x0+x1)/2,y:H-6,"text-anchor":"middle","font-size":12,
    fill:"var(--ink2)"}); xl.textContent="파장 (nm)"; svg.appendChild(xl);

  // 실험 마커
  for(const mk of DATA.markers){
    if(mk.nm<gmin||mk.nm>gmax)continue;
    svg.appendChild(el("line",{x1:sx(mk.nm),y1:y1,x2:sx(mk.nm),y2:y0,
      stroke:"var(--axis)","stroke-width":1.3,"stroke-dasharray":"5 4"}));
    const t=el("text",{x:sx(mk.nm)-4,y:y1+52,"text-anchor":"end","font-size":11,
      fill:"var(--ink2)",transform:`rotate(-90 ${sx(mk.nm)-4} ${y1+52})`});
    t.textContent=mk.label; svg.appendChild(t);
  }

  // 곡선
  for(const s of DATA.series){
    if(!visible.has(s.id))continue;
    let d="";
    for(let i=0;i<DATA.grid.length;i++){
      d+=(i?"L":"M")+sx(DATA.grid[i]).toFixed(1)+" "+sy(s.eps[i]).toFixed(1)+" ";
    }
    svg.appendChild(el("path",{d,fill:"none",stroke:s.color,"stroke-width":2.2,
      "stroke-linejoin":"round"}));
  }
  // 크로스헤어 그룹(위에)
  cross=el("g",{opacity:0}); svg.appendChild(cross);
  chLine=el("line",{y1,y2:y0,stroke:"var(--axis)","stroke-width":1}); cross.appendChild(chLine);
  chDots=[];
}
let cross,chLine,chDots=[];

function buildLegend(){
  const box=document.getElementById("legend"); box.innerHTML="";
  for(const s of DATA.series){
    const lab=document.createElement("label");
    lab.className=visible.has(s.id)?"":"off";
    const cb=document.createElement("input");
    cb.type="checkbox"; cb.checked=visible.has(s.id);
    cb.onchange=()=>{ if(cb.checked)visible.add(s.id);else visible.delete(s.id);
      lab.className=cb.checked?"":"off"; draw(); };
    const sw=document.createElement("span"); sw.className="swatch";
    sw.style.background=s.color;
    const tx=document.createElement("span"); tx.className="txt";
    tx.textContent=`${s.label} · ${s.brightNm}nm (f=${s.f})`;
    lab.append(cb,sw,tx); box.appendChild(lab);
  }
}

function buildPresets(){
  const box=document.getElementById("presets"); box.innerHTML="";
  for(const name in DATA.presets){
    const b=document.createElement("button"); b.textContent=name;
    b.onclick=()=>{ visible.clear(); DATA.presets[name].forEach(id=>visible.add(id));
      [...box.children].forEach(c=>c.classList.toggle("on",c===b));
      buildLegend(); draw(); };
    box.appendChild(b);
  }
  box.firstChild.classList.add("on");
}

// 마우스 크로스헤어 + 툴팁
const tip=document.getElementById("tip");
function onMove(ev){
  const r=svg.getBoundingClientRect();
  const px=(ev.clientX-r.left)/r.width*W;
  if(px<x0||px>x1){ cross.setAttribute("opacity",0); tip.style.opacity=0; return; }
  const nm=Math.round(gmin+(px-x0)/(x1-x0)*(gmax-gmin));
  const idx=nm-gmin;
  cross.setAttribute("opacity",1);
  chLine.setAttribute("x1",sx(nm)); chLine.setAttribute("x2",sx(nm));
  chDots.forEach(d=>d.remove()); chDots=[];
  let rows="";
  for(const s of DATA.series){
    if(!visible.has(s.id))continue;
    const v=s.eps[idx]||0;
    const dot=el("circle",{cx:sx(nm),cy:sy(v),r:4,fill:s.color,
      stroke:"var(--surface)","stroke-width":1.5});
    cross.appendChild(dot); chDots.push(dot);
    rows+=`<div class="row"><span class="d" style="background:${s.color}"></span>`
      +`${s.label}: <b style="color:var(--ink)">${v.toLocaleString()}</b></div>`;
  }
  tip.innerHTML=`<b>${nm} nm</b>${rows}`;
  const r2=svg.getBoundingClientRect();
  let tx=(sx(nm)/W)*r2.width+14, ty=ev.clientY-r.top-10;
  if(tx>r2.width-170)tx-=185;
  tip.style.left=tx+"px"; tip.style.top=ty+"px"; tip.style.opacity=1;
}
svg.addEventListener("mousemove",onMove);
svg.addEventListener("mouseleave",()=>{cross&&cross.setAttribute("opacity",0);tip.style.opacity=0;});

// 테마 토글
document.getElementById("themeBtn").onclick=()=>{
  const cur=document.documentElement.getAttribute("data-theme");
  const dark=cur? cur==="dark" : matchMedia("(prefers-color-scheme:dark)").matches;
  document.documentElement.setAttribute("data-theme",dark?"light":"dark");
  draw();
};

buildPresets(); buildLegend(); draw();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    sys.exit(main())
