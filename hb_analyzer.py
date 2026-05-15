#!/usr/bin/env python3
"""
PyMOL HB Interaction Analyzer - Python 核心分析模組
使用標準庫 + numpy 進行 PDB 解析與氫鍵計算
用法：python hb_analyzer.py <PDB_ID 或 .pdb 路徑> [--output report.html]
"""

import os
import sys
import json
import math
import argparse
import urllib.request
import urllib.error
from pathlib import Path

# ── 常數 ──────────────────────────────────────────────────────────────────────
HB_DA_CUTOFF   = 3.5   # Å  D…A heavy-atom distance
HB_SOFT_CUTOFF = 4.0   # Å  寬鬆接觸距離（for reporting）
CONTACT_CUTOFF = 5.0   # Å  近鄰殘基定義

# H-bond 捐體/受體的重原子元素
HB_ELEMENTS = {"N", "O", "S", "F"}

# 標準水分子 resname
WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD", "TIP", "SOL"}

# 胺基酸單字母碼
AA_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "HID":"H","HIE":"H","HIP":"H","CYX":"C","MSE":"M",
}

# 突變策略建議表
MUTATION_HINTS = {
    "SER": "Ser→Ala 移除 -OH，可破壞氫鍵；Ser→Thr 增大空間位阻",
    "THR": "Thr→Val 移除 -OH，消除氫鍵能力",
    "TYR": "Tyr→Phe 移除酚 -OH，去除氫鍵；Tyr→Trp 增大疏水接觸",
    "ASN": "Asn→Ala 移除醯胺基，破壞氫鍵網絡",
    "GLN": "Gln→Ala 縮短側鏈，消除氫鍵；Gln→Glu 引入負電荷",
    "ASP": "Asp→Asn 中和負電荷；Asp→Ala 完全移除羧基",
    "GLU": "Glu→Gln 中和負電荷；Glu→Ala 完全移除羧基",
    "HIS": "His→Ala 移除咪唑環，消除氫鍵與配位能力；His→Phe 保留芳香環但去除氫鍵",
    "ARG": "Arg→Lys 縮短側鏈；Arg→Ala 完全移除胍基正電荷",
    "LYS": "Lys→Arg 維持正電荷但改變幾何；Lys→Ala 移除正電荷",
    "CYS": "Cys→Ser 以 -OH 替換 -SH；Cys→Ala 移除極性基",
    "TRP": "Trp→Ala 移除吲哚環，影響氫鍵與堆疊交互作用",
    "MET": "Met→Leu 移除硫原子；Met→Ala 大幅縮減側鏈",
    "PHE": "Phe→Ala 移除苯環；Phe→Tyr 引入 -OH 增加氫鍵能力",
    "LEU": "Leu→Ala 縮短疏水側鏈；Leu→Ile 改變分支位置",
    "ILE": "Ile→Val 縮短側鏈；Ile→Ala 大幅縮減",
    "VAL": "Val→Ala 縮減疏水接觸；Val→Ile 增加側鏈長度",
    "ALA": "Ala→Gly 增加骨架彈性；Ala→Val 增大體積產生位阻",
    "GLY": "Gly→Ala 引入 Cβ，限制骨架彈性",
    "PRO": "Pro→Ala 移除環狀限制，增加彈性；通常是特殊結構點",
}


# ── PDB 解析 ──────────────────────────────────────────────────────────────────

def parse_pdb(text: str) -> list[dict]:
    """解析 PDB 格式文字，回傳 atom 清單。"""
    atoms = []
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        try:
            atom = {
                "serial":  int(line[6:11]),
                "name":    line[12:16].strip(),
                "altLoc":  line[16].strip(),
                "resName": line[17:20].strip(),
                "chain":   line[21].strip(),
                "resSeq":  int(line[22:26]),
                "iCode":   line[26].strip(),
                "x":       float(line[30:38]),
                "y":       float(line[38:46]),
                "z":       float(line[46:54]),
                "occupancy": float(line[54:60]) if len(line) > 54 else 1.0,
                "bfactor":   float(line[60:66]) if len(line) > 60 else 0.0,
                "element": line[76:78].strip() if len(line) > 76 else "",
                "isHet":   rec == "HETATM",
            }
            # 推斷元素（若欄位為空）
            if not atom["element"]:
                atom["element"] = atom["name"].lstrip("0123456789")[0].upper()
            # 只保留第一個 altLoc
            if atom["altLoc"] and atom["altLoc"] not in ("", "A", " "):
                continue
            atoms.append(atom)
        except (ValueError, IndexError):
            continue
    return atoms


def fetch_pdb(pdb_id: str) -> str:
    """從 RCSB PDB 下載結構文字。"""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"⬇  下載 {url} ...")
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        raise SystemExit(f"HTTP {e.code}：無法取得 {pdb_id}")
    except Exception as e:
        raise SystemExit(f"下載失敗：{e}")


def load_structure(source: str) -> tuple[str, str, list[dict]]:
    """
    source 可為 PDB ID（4 字元）或檔案路徑。
    回傳 (pdb_id, pdb_text, atoms)
    """
    p = Path(source)
    if p.is_file():
        pdb_text = p.read_text(errors="replace")
        pdb_id   = p.stem.upper()
    elif len(source) == 4 and source.isalnum():
        pdb_text = fetch_pdb(source)
        pdb_id   = source.upper()
    else:
        raise SystemExit(f"無法識別輸入：'{source}'（需為 4 字元 PDB ID 或 .pdb 檔案路徑）")
    atoms = parse_pdb(pdb_text)
    if not atoms:
        raise SystemExit("解析失敗：未找到 ATOM/HETATM 記錄")
    print(f"✓  載入 {len(atoms)} 個原子（{pdb_id}）")
    return pdb_id, pdb_text, atoms


# ── Ligand 偵測 ───────────────────────────────────────────────────────────────

def detect_ligands(atoms: list[dict]) -> list[str]:
    """回傳所有非水小分子的 resName 清單。"""
    ligands = set()
    for a in atoms:
        if a["isHet"] and a["resName"] not in WATER_RESNAMES:
            ligands.add(a["resName"])
    return sorted(ligands)


def get_ligand_atoms(atoms: list[dict], lig_resname: str) -> list[dict]:
    return [a for a in atoms if a["resName"] == lig_resname]


def get_protein_atoms(atoms: list[dict]) -> list[dict]:
    return [a for a in atoms if not a["isHet"]]


# ── 距離計算 ──────────────────────────────────────────────────────────────────

def dist3(a, b) -> float:
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)


# ── 氫鍵分析 ──────────────────────────────────────────────────────────────────

def find_hbonds(lig_atoms: list[dict], prot_atoms: list[dict]) -> list[dict]:
    """
    基於重原子 D…A 距離 ≤ 3.5 Å 偵測氫鍵候選對。
    回傳 list of dict，每筆包含兩端原子資訊與距離。
    """
    lig_polar  = [a for a in lig_atoms  if a["element"] in HB_ELEMENTS]
    prot_polar = [a for a in prot_atoms if a["element"] in HB_ELEMENTS]

    hbonds = []
    seen   = set()  # 避免同一殘基重複

    for la in lig_polar:
        for pa in prot_polar:
            d = dist3(la, pa)
            if d <= HB_DA_CUTOFF:
                key = (la["serial"], pa["chain"], pa["resSeq"], pa["resName"])
                if key in seen:
                    # 保留距離最短的
                    for hb in hbonds:
                        if (hb["lig_serial"] == la["serial"] and
                                hb["prot_chain"] == pa["chain"] and
                                hb["prot_resSeq"] == pa["resSeq"]):
                            if d < hb["distance"]:
                                hb.update(_make_hb_entry(la, pa, d))
                    continue
                seen.add(key)
                hbonds.append(_make_hb_entry(la, pa, d))

    hbonds.sort(key=lambda x: x["distance"])
    return hbonds


def _make_hb_entry(la: dict, pa: dict, d: float) -> dict:
    return {
        "lig_atom":    la["name"],
        "lig_element": la["element"],
        "lig_serial":  la["serial"],
        "lig_resName": la["resName"],
        "prot_atom":   pa["name"],
        "prot_element": pa["element"],
        "prot_resName": pa["resName"],
        "prot_chain":  pa["chain"],
        "prot_resSeq": pa["resSeq"],
        "prot_oneLetter": AA_ONE.get(pa["resName"], "?"),
        "distance":    round(d, 3),
        "lig_xyz":     [la["x"], la["y"], la["z"]],
        "prot_xyz":    [pa["x"], pa["y"], pa["z"]],
    }


# ── 近鄰殘基 ──────────────────────────────────────────────────────────────────

def find_nearby_residues(lig_atoms: list[dict], prot_atoms: list[dict],
                         cutoff: float = CONTACT_CUTOFF) -> list[dict]:
    """回傳配體 cutoff Å 內所有蛋白質殘基（去重，附最近距離）。"""
    residues = {}
    for la in lig_atoms:
        for pa in prot_atoms:
            d = dist3(la, pa)
            if d <= cutoff:
                key = (pa["chain"], pa["resSeq"], pa["resName"])
                if key not in residues or d < residues[key]["min_dist"]:
                    residues[key] = {
                        "chain":   pa["chain"],
                        "resSeq":  pa["resSeq"],
                        "resName": pa["resName"],
                        "oneLetter": AA_ONE.get(pa["resName"], "?"),
                        "min_dist": round(d, 3),
                    }
    result = sorted(residues.values(), key=lambda x: x["min_dist"])
    return result


# ── 突變分析 ──────────────────────────────────────────────────────────────────

def mutation_analysis(hbonds: list[dict], nearby: list[dict]) -> list[dict]:
    """針對氫鍵殘基與近鄰殘基提供突變建議。"""
    analyzed = []
    # 氫鍵殘基（高優先）
    hb_keys = {(hb["prot_chain"], hb["prot_resSeq"]) for hb in hbonds}
    for hb in hbonds:
        key = (hb["prot_chain"], hb["prot_resSeq"])
        rn = hb["prot_resName"]
        analyzed.append({
            "residue":     f"{hb['prot_resName']}{hb['prot_resSeq']} (Chain {hb['prot_chain']})",
            "oneLetter":   hb["prot_oneLetter"],
            "role":        "氫鍵交互作用（直接）",
            "priority":    "高",
            "hb_dist":     hb["distance"],
            "suggestion":  MUTATION_HINTS.get(rn, f"{rn}→Ala 移除側鏈功能基"),
        })
    # 近鄰殘基（中優先）
    for res in nearby:
        k = (res["chain"], res["resSeq"])
        if k in hb_keys:
            continue
        rn = res["resName"]
        analyzed.append({
            "residue":     f"{res['resName']}{res['resSeq']} (Chain {res['chain']})",
            "oneLetter":   res["oneLetter"],
            "role":        f"近鄰殘基（{res['min_dist']:.2f} Å）",
            "priority":    "中",
            "hb_dist":     None,
            "suggestion":  MUTATION_HINTS.get(rn, f"{rn}→Ala 探索空間效應"),
        })
    return analyzed


# ── 結構摘要 ──────────────────────────────────────────────────────────────────

def structure_summary(pdb_id: str, atoms: list[dict], lig_resname: str) -> dict:
    chains  = sorted({a["chain"] for a in atoms if not a["isHet"] and a["chain"]})
    n_prot  = len([a for a in atoms if not a["isHet"] and a["name"] == "CA"])
    n_lig   = len([a for a in atoms if a["resName"] == lig_resname])
    ligands = detect_ligands(atoms)
    return {
        "pdb_id":    pdb_id,
        "chains":    chains,
        "n_residues": n_prot,
        "n_lig_atoms": n_lig,
        "lig_resname": lig_resname,
        "all_ligands": ligands,
    }


# ── HTML 報告生成 ─────────────────────────────────────────────────────────────

REPORT_HTML = """\
<!DOCTYPE html>
<html lang="zh-TW">
<head>
<meta charset="UTF-8">
<title>HB Interaction Report – {pdb_id}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  :root {{ --accent:#2563eb; --bg:#f8fafc; --card:#fff;
           --text:#1e293b; --sub:#64748b; --border:#e2e8f0; }}
  *{{ box-sizing:border-box; margin:0; padding:0; }}
  body {{ font-family:'Segoe UI',system-ui,sans-serif; background:var(--bg);
          color:var(--text); line-height:1.6; }}
  header {{ background:linear-gradient(135deg,#1e3a5f,var(--accent));
            color:#fff; padding:2rem; }}
  header h1 {{ font-size:1.8rem; font-weight:700; }}
  header p  {{ opacity:.85; margin-top:.3rem; }}
  .container {{ max-width:1100px; margin:0 auto; padding:2rem 1.5rem; }}
  .grid-2 {{ display:grid; grid-template-columns:1fr 1fr; gap:1.5rem;
             margin-bottom:1.5rem; }}
  @media(max-width:768px){{ .grid-2{{grid-template-columns:1fr;}} }}
  .card {{ background:var(--card); border:1px solid var(--border);
           border-radius:12px; padding:1.25rem; }}
  .card h2 {{ font-size:1rem; font-weight:600; color:var(--accent);
              margin-bottom:.75rem; border-bottom:2px solid var(--border);
              padding-bottom:.5rem; }}
  #viewer {{ width:100%; height:380px; border-radius:8px;
             overflow:hidden; position:relative; }}
  table {{ width:100%; border-collapse:collapse; font-size:.88rem; }}
  th {{ background:#eff6ff; padding:.5rem .75rem; text-align:left;
        font-weight:600; border-bottom:2px solid var(--border); }}
  td {{ padding:.45rem .75rem; border-bottom:1px solid var(--border); }}
  tr:last-child td {{ border-bottom:none; }}
  tr:hover td {{ background:#f0f7ff; }}
  .badge {{ display:inline-block; padding:.15rem .55rem; border-radius:999px;
            font-size:.75rem; font-weight:600; }}
  .badge-high   {{ background:#fee2e2; color:#b91c1c; }}
  .badge-medium {{ background:#fef3c7; color:#92400e; }}
  .stat-grid {{ display:grid; grid-template-columns:repeat(3,1fr); gap:.75rem; }}
  .stat {{ text-align:center; padding:.75rem; background:#eff6ff;
           border-radius:8px; }}
  .stat-value {{ font-size:1.8rem; font-weight:700; color:var(--accent); }}
  .stat-label {{ font-size:.78rem; color:var(--sub); margin-top:.15rem; }}
  .mutation-row {{ margin-bottom:.6rem; padding:.75rem; background:#f8fafc;
                   border-left:3px solid var(--accent); border-radius:0 8px 8px 0; }}
  .mutation-title {{ font-weight:600; font-size:.9rem; }}
  .mutation-role  {{ font-size:.8rem; color:var(--sub); margin:.15rem 0; }}
  .mutation-sugg  {{ font-size:.82rem; color:#374151; }}
  .dist-badge {{ display:inline-block; background:#dcfce7; color:#166534;
                 padding:.1rem .45rem; border-radius:4px; font-size:.78rem;
                 font-weight:600; }}
  footer {{ text-align:center; padding:1.5rem; color:var(--sub);
            font-size:.82rem; border-top:1px solid var(--border); margin-top:2rem; }}
</style>
</head>
<body>
<header>
  <h1>🔬 蛋白質–配體氫鍵交互作用分析報告</h1>
  <p>PDB ID：<strong>{pdb_id}</strong> ｜ 配體：<strong>{lig_resname}</strong> ｜
     鏈：{chains} ｜ 殘基數：{n_residues}</p>
</header>

<div class="container">

  <!-- 統計卡 -->
  <div class="card" style="margin-bottom:1.5rem">
    <h2>📊 分析摘要</h2>
    <div class="stat-grid">
      <div class="stat">
        <div class="stat-value">{n_hbonds}</div>
        <div class="stat-label">氫鍵對數</div>
      </div>
      <div class="stat">
        <div class="stat-value">{n_unique_res}</div>
        <div class="stat-label">氫鍵殘基數</div>
      </div>
      <div class="stat">
        <div class="stat-value">{n_nearby}</div>
        <div class="stat-label">近鄰殘基（5 Å）</div>
      </div>
    </div>
  </div>

  <div class="grid-2">
    <!-- 3D 檢視器 -->
    <div class="card">
      <h2>🧬 3D 結構視覺化</h2>
      <div id="viewer"></div>
      <p style="font-size:.78rem;color:var(--sub);margin-top:.5rem">
        黃色虛線 = 氫鍵 ｜ 棍棒 = 配體 ｜ Cartoon = 蛋白質</p>
    </div>

    <!-- 氫鍵表 -->
    <div class="card">
      <h2>⚡ 氫鍵清單（D–A ≤ 3.5 Å）</h2>
      <table>
        <thead>
          <tr><th>配體原子</th><th>蛋白質殘基</th><th>原子</th>
              <th>距離 (Å)</th><th>類型</th></tr>
        </thead>
        <tbody>{hbond_rows}</tbody>
      </table>
    </div>
  </div>

  <!-- 突變分析 -->
  <div class="card" style="margin-bottom:1.5rem">
    <h2>🧪 突變位點分析</h2>
    {mutation_cards}
  </div>

  <!-- 近鄰殘基 -->
  <div class="card">
    <h2>🔍 5 Å 近鄰殘基</h2>
    <table>
      <thead>
        <tr><th>殘基</th><th>鏈</th><th>最近距離 (Å)</th><th>角色</th></tr>
      </thead>
      <tbody>{nearby_rows}</tbody>
    </table>
  </div>

</div>

<footer>
  PyMOL HB Interaction Analyzer ｜ 分析日期：{date}
</footer>

<script>
const PDB_TEXT = `{pdb_text_js}`;
const HBONDS   = {hbonds_json};

(function(){{
  const viewer = $3Dmol.createViewer(
    document.getElementById('viewer'),
    {{backgroundColor:'#f0f4f8'}}
  );
  viewer.addModel(PDB_TEXT, 'pdb');

  // 蛋白質 cartoon
  viewer.setStyle({{hetflag:false}},
    {{cartoon:{{colorscheme:'spectrum', opacity:0.85}}}});

  // 配體 sticks
  viewer.setStyle({{resn:'{lig_resname}'}},
    {{stick:{{colorscheme:'default', radius:0.25}}}});

  // 水分子隱藏
  viewer.setStyle({{resn:'HOH'}}, {{}});
  viewer.setStyle({{resn:'WAT'}}, {{}});

  // 氫鍵虛線
  for(const hb of HBONDS){{
    viewer.addCylinder({{
      start: {{x:hb.lig_xyz[0], y:hb.lig_xyz[1], z:hb.lig_xyz[2]}},
      end:   {{x:hb.prot_xyz[0], y:hb.prot_xyz[1], z:hb.prot_xyz[2]}},
      radius: 0.06,
      color: '#facc15',
      dashed: true,
      fromCap:1, toCap:1,
    }});
  }}

  // 氫鍵殘基標籤
  const labeled = new Set();
  for(const hb of HBONDS){{
    const key = hb.prot_chain + hb.prot_resSeq;
    if(labeled.has(key)) continue;
    labeled.add(key);
    viewer.addLabel(
      hb.prot_resName + hb.prot_resSeq,
      {{
        position:{{x:hb.prot_xyz[0], y:hb.prot_xyz[1], z:hb.prot_xyz[2]}},
        backgroundColor:'rgba(30,58,95,0.85)',
        fontColor:'#fff',
        fontSize:11,
        borderRadius:3,
        padding:2,
      }}
    );
  }}

  viewer.zoomTo({{resn:'{lig_resname}'}}, 500);
  viewer.render();
}})();
</script>
</body>
</html>
"""


def _hbond_type(la_elem: str, pa_elem: str, dist: float) -> str:
    if dist <= 2.5:
        return "強氫鍵"
    if dist <= 3.2:
        return "中強氫鍵"
    return "弱氫鍵"


def generate_report(pdb_id: str, pdb_text: str, lig_resname: str,
                    hbonds: list[dict], nearby: list[dict],
                    mutations: list[dict], summary: dict) -> str:
    import datetime

    # H-bond table rows
    hbond_rows = ""
    for hb in hbonds:
        htype = _hbond_type(hb["lig_element"], hb["prot_element"], hb["distance"])
        hbond_rows += (
            f"<tr>"
            f"<td>{hb['lig_resName']}:{hb['lig_atom']} ({hb['lig_element']})</td>"
            f"<td>{hb['prot_resName']}{hb['prot_resSeq']} ({hb['prot_chain']})</td>"
            f"<td>{hb['prot_atom']} ({hb['prot_element']})</td>"
            f"<td><span class='dist-badge'>{hb['distance']:.3f}</span></td>"
            f"<td>{htype}</td>"
            f"</tr>\n"
        )

    # Nearby residue rows
    hb_keys = {(hb["prot_chain"], hb["prot_resSeq"]) for hb in hbonds}
    nearby_rows = ""
    for res in nearby:
        role = "✦ 氫鍵殘基" if (res["chain"], res["resSeq"]) in hb_keys else "接觸殘基"
        nearby_rows += (
            f"<tr>"
            f"<td>{res['resName']}{res['resSeq']} ({res['oneLetter']})</td>"
            f"<td>{res['chain']}</td>"
            f"<td>{res['min_dist']:.3f}</td>"
            f"<td>{role}</td>"
            f"</tr>\n"
        )

    # Mutation cards
    mutation_cards = ""
    for m in mutations:
        badge_cls = "badge-high" if m["priority"] == "高" else "badge-medium"
        dist_html = (f" &nbsp;<span class='dist-badge'>{m['hb_dist']:.3f} Å</span>"
                     if m["hb_dist"] else "")
        mutation_cards += (
            f"<div class='mutation-row'>"
            f"<div class='mutation-title'>{m['residue']}"
            f"{dist_html}"
            f" &nbsp;<span class='badge {badge_cls}'>{m['priority']}優先</span></div>"
            f"<div class='mutation-role'>{m['role']}</div>"
            f"<div class='mutation-sugg'>💡 {m['suggestion']}</div>"
            f"</div>\n"
        )

    # Escape PDB text for JS
    pdb_js = pdb_text.replace("`", "\\`").replace("\\", "\\\\")

    n_unique = len({(hb["prot_chain"], hb["prot_resSeq"]) for hb in hbonds})

    html = REPORT_HTML.format(
        pdb_id       = pdb_id,
        lig_resname  = lig_resname,
        chains       = ", ".join(summary["chains"]),
        n_residues   = summary["n_residues"],
        n_hbonds     = len(hbonds),
        n_unique_res = n_unique,
        n_nearby     = len(nearby),
        hbond_rows   = hbond_rows,
        nearby_rows  = nearby_rows,
        mutation_cards = mutation_cards,
        hbonds_json  = json.dumps(hbonds, ensure_ascii=False),
        pdb_text_js  = pdb_js,
        date         = datetime.date.today().isoformat(),
    )
    return html


# ── CLI 進入點 ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="PyMOL HB Interaction Analyzer – 氫鍵分析與報告生成")
    parser.add_argument("source",
        help="PDB ID（如 1HSG）或 .pdb 檔案路徑")
    parser.add_argument("--ligand", "-l", default=None,
        help="指定配體 resname（預設自動偵測第一個非水 HETATM）")
    parser.add_argument("--output", "-o", default=None,
        help="輸出 HTML 報告路徑（預設：<PDB_ID>_hb_report.html）")
    parser.add_argument("--cutoff", "-c", type=float, default=3.5,
        help="D–A 氫鍵距離上限 Å（預設 3.5）")
    args = parser.parse_args()

    global HB_DA_CUTOFF
    HB_DA_CUTOFF = args.cutoff

    pdb_id, pdb_text, atoms = load_structure(args.source)

    ligands = detect_ligands(atoms)
    if not ligands:
        raise SystemExit("未偵測到非水 HETATM 記錄（配體）")

    lig_resname = args.ligand or ligands[0]
    if lig_resname not in ligands:
        raise SystemExit(f"指定的配體 '{lig_resname}' 不在結構中。可用：{ligands}")

    print(f"🔍 分析配體：{lig_resname}（可用配體：{ligands}）")

    lig_atoms  = get_ligand_atoms(atoms, lig_resname)
    prot_atoms = get_protein_atoms(atoms)

    hbonds  = find_hbonds(lig_atoms, prot_atoms)
    nearby  = find_nearby_residues(lig_atoms, prot_atoms)
    summary = structure_summary(pdb_id, atoms, lig_resname)
    muts    = mutation_analysis(hbonds, nearby)

    print(f"✓  氫鍵：{len(hbonds)} 對 ｜ 近鄰殘基：{len(nearby)} 個")

    out_path = args.output or f"{pdb_id}_hb_report.html"
    report   = generate_report(pdb_id, pdb_text, lig_resname,
                                hbonds, nearby, muts, summary)
    Path(out_path).write_text(report, encoding="utf-8")
    print(f"✅ 報告已儲存：{out_path}")

    # 同時輸出 JSON 摘要
    json_path = out_path.replace(".html", "_data.json")
    data = {"summary": summary, "hbonds": hbonds,
            "nearby": nearby, "mutations": muts}
    Path(json_path).write_text(json.dumps(data, indent=2, ensure_ascii=False),
                                encoding="utf-8")
    print(f"✅ JSON 資料：{json_path}")


if __name__ == "__main__":
    main()
