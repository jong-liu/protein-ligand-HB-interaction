"""
pymol_automation.py
====================
PyMOL 自動化分析腳本 – 在 PyMOL 環境內執行

用法（在 PyMOL GUI 或命令列）：
  pymol -c pymol_automation.py -- 1HSG
  pymol -c pymol_automation.py -- 1HSG MK1          # 指定配體
  pymol -c pymol_automation.py -- /path/to/file.pdb # 本地檔案

或在 PyMOL 指令列：
  run pymol_automation.py
  analyze 1HSG
  analyze 1HSG, ligand=MK1
"""

import sys
import os
import json
import math
import datetime
from pathlib import Path

try:
    from pymol import cmd, stored
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    print("警告：未偵測到 PyMOL，某些功能將無法使用")

# ── 全域設定 ──────────────────────────────────────────────────────────────────
OUTPUT_DIR   = Path("pymol_output")
IMG_WIDTH    = 2400
IMG_HEIGHT   = 1800
HB_CUTOFF    = 3.5   # Å

WATER_RES = {"HOH", "WAT", "H2O", "DOD"}

# 色彩方案
COLOR_PROTEIN  = "slate"
COLOR_LIGAND   = "yellow"
COLOR_HB_RESI  = "tv_orange"
COLOR_SURFACE  = "lightblue"

# ── 輔助函數 ──────────────────────────────────────────────────────────────────

def setup_output_dir() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUTPUT_DIR


def detect_ligand(pdb_id: str) -> str | None:
    """自動偵測結構中第一個非水 HETATM 的 resname。"""
    stored.hets = []
    cmd.iterate(f"({pdb_id}) and hetatm and not resn HOH+WAT+H2O+DOD",
                "stored.hets.append(resn)")
    if not stored.hets:
        return None
    # 回傳最多原子數的配體
    counts = {}
    for r in stored.hets:
        counts[r] = counts.get(r, 0) + 1
    return max(counts, key=counts.get)


def get_hbond_residues(pdb_id: str, lig_sel: str) -> list[dict]:
    """
    找出與配體距離 ≤ HB_CUTOFF Å 的極性蛋白質原子對，
    回傳 [{resn, resi, chain, atom, dist}, ...]。
    """
    stored.contacts = []
    cmd.iterate_state(
        1,
        f"polymer and elem N+O+S within {HB_CUTOFF} of ({lig_sel})",
        "stored.contacts.append({'resn':resn,'resi':resi,'chain':chain,"
        "'atom':name,'x':x,'y':y,'z':z})"
    )

    stored.lig_atoms = []
    cmd.iterate_state(1, f"({lig_sel}) and elem N+O+S",
                      "stored.lig_atoms.append({'atom':name,'x':x,'y':y,'z':z})")

    hbonds = []
    for la in stored.lig_atoms:
        for pa in stored.contacts:
            d = math.sqrt((la["x"]-pa["x"])**2 +
                          (la["y"]-pa["y"])**2 +
                          (la["z"]-pa["z"])**2)
            if d <= HB_CUTOFF:
                hbonds.append({**pa, "dist": round(d, 3), "lig_atom": la["atom"]})

    # 去重（保留最短距離）
    best = {}
    for hb in hbonds:
        k = (hb["chain"], hb["resi"], hb["resn"])
        if k not in best or hb["dist"] < best[k]["dist"]:
            best[k] = hb
    return sorted(best.values(), key=lambda x: x["dist"])


# ── PyMOL 視覺化核心 ─────────────────────────────────────────────────────────

def setup_visualization(pdb_id: str, lig_sel: str, hb_residues: list[dict]):
    """設定 PyMOL 視覺化：cartoon + ligand sticks + H-bond annotations。"""

    # ① 初始化：隱藏所有顯示
    cmd.hide("everything", "all")

    # ② 蛋白質 cartoon（spectrum 配色）
    cmd.show("cartoon", f"polymer and {pdb_id}")
    cmd.color("spectrum", f"polymer and {pdb_id}")

    # ③ 配體 sticks（黃色）
    cmd.show("sticks", lig_sel)
    cmd.color(COLOR_LIGAND, lig_sel)
    cmd.show("spheres", f"({lig_sel}) and elem C")
    cmd.set("sphere_scale", 0.15, lig_sel)

    # ④ 氫鍵殘基：高亮橘色 sticks
    if hb_residues:
        hb_sel_list = " or ".join(
            f"(chain {r['chain']} and resi {r['resi']})"
            for r in hb_residues
        )
        hb_sel = f"hb_residues_{pdb_id}"
        cmd.select(hb_sel, f"polymer and ({hb_sel_list})")
        cmd.show("sticks", hb_sel)
        cmd.color(COLOR_HB_RESI, hb_sel)
        cmd.deselect()

    # ⑤ 氫鍵距離線（黃色虛線）
    for i, hb in enumerate(hb_residues):
        dist_name = f"hb_{i+1:02d}_{hb['resn']}{hb['resi']}"
        cmd.distance(dist_name, lig_sel,
                     f"(chain {hb['chain']} and resi {hb['resi']} "
                     f"and name {hb['atom']})")
        cmd.color("yellow", dist_name)

    cmd.hide("labels", "dist*")   # 隱藏距離標籤（後面自訂）

    # ⑥ 殘基標籤（CA 原子）
    if hb_residues:
        cmd.set("label_size", 14)
        cmd.set("label_color", "white")
        cmd.set("label_bg_color", "marine")
        cmd.label(
            f"({hb_sel}) and name CA",
            '"%s%s\\n%.2f Å" % (resn, resi, '
            # 透過 stored 傳入距離
            '0.0)'
        )
        # 使用更好的標籤：自訂 Python 字串
        _add_residue_labels(hb_residues)

    # ⑦ 整體美化
    cmd.bg_color("black")
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.6)
    cmd.set("spec_reflect", 0.3)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("stick_ball", 1)
    cmd.set("stick_ball_ratio", 1.5)

    # ⑧ 聚焦配體
    cmd.orient(lig_sel)
    cmd.zoom(lig_sel, 8)


def _add_residue_labels(hb_residues: list[dict]):
    """替每個氫鍵殘基加上帶距離的 CA 標籤。"""
    for hb in hb_residues:
        sel = f"(chain {hb['chain']} and resi {hb['resi']} and name CA)"
        cmd.label(sel, f'"{hb["resn"]}{hb["resi"]}\\n{hb["dist"]:.2f}Å"')


# ── 圖像匯出 ──────────────────────────────────────────────────────────────────

def export_images(pdb_id: str, lig_sel: str, out_dir: Path) -> list[str]:
    """匯出多角度 PNG 圖像，回傳路徑清單。"""
    paths = []

    views = {
        "front": None,  # 目前視角
        "top":   "cmd.turn('x', 90)",
        "side":  "cmd.turn('y', 90)",
    }

    for view_name, rotate_cmd in views.items():
        if rotate_cmd:
            eval(rotate_cmd)

        img_path = str(out_dir / f"{pdb_id}_{view_name}.png")
        cmd.ray(IMG_WIDTH, IMG_HEIGHT)
        cmd.png(img_path, dpi=300, quiet=1)
        paths.append(img_path)
        print(f"  ✓ 圖像：{img_path}")

    # 恢復原始視角
    cmd.orient(lig_sel)
    cmd.zoom(lig_sel, 8)

    # ── 額外：表面視角 ──
    cmd.show("surface", f"polymer and {pdb_id}")
    cmd.set("transparency", 0.65, f"polymer and {pdb_id}")
    cmd.color(COLOR_SURFACE, f"polymer and {pdb_id}")
    cmd.set("surface_quality", 1)

    surf_path = str(out_dir / f"{pdb_id}_surface.png")
    cmd.ray(IMG_WIDTH, IMG_HEIGHT)
    cmd.png(surf_path, dpi=300, quiet=1)
    paths.append(surf_path)
    print(f"  ✓ 圖像（表面）：{surf_path}")

    cmd.hide("surface", f"polymer and {pdb_id}")
    return paths


def save_session(pdb_id: str, out_dir: Path):
    """儲存 PyMOL .pse 工作階段檔。"""
    pse_path = str(out_dir / f"{pdb_id}_session.pse")
    cmd.save(pse_path)
    print(f"  ✓ PSE 工作階段：{pse_path}")
    return pse_path


# ── 文字輸出 ──────────────────────────────────────────────────────────────────

def print_hbond_report(pdb_id: str, lig_resname: str, hb_residues: list[dict]):
    """在 PyMOL 終端機印出氫鍵摘要。"""
    sep = "─" * 55
    print(f"\n{sep}")
    print(f"  PyMOL H-Bond Analysis  |  {pdb_id}  |  Ligand: {lig_resname}")
    print(sep)
    print(f"  {'殘基':<12} {'鏈':^4} {'原子':<6} {'距離(Å)':>8}")
    print(sep)
    for hb in hb_residues:
        print(f"  {hb['resn']}{hb['resi']:<8} {hb['chain']:^4} "
              f"{hb['atom']:<6} {hb['dist']:>8.3f}")
    print(sep)
    print(f"  共 {len(hb_residues)} 個氫鍵殘基")
    print(f"{sep}\n")


def save_json(pdb_id: str, lig_resname: str,
              hb_residues: list[dict], out_dir: Path) -> str:
    data = {
        "pdb_id": pdb_id,
        "ligand": lig_resname,
        "analysis_date": datetime.date.today().isoformat(),
        "hb_cutoff_A": HB_CUTOFF,
        "hbonds": hb_residues,
    }
    json_path = str(out_dir / f"{pdb_id}_hbonds.json")
    Path(json_path).write_text(json.dumps(data, indent=2, ensure_ascii=False))
    print(f"  ✓ JSON 資料：{json_path}")
    return json_path


# ── 主分析函數（可從 PyMOL 指令列呼叫）──────────────────────────────────────

def analyze(source: str, ligand: str = "", cutoff: float = HB_CUTOFF):
    """
    主要分析入口。
    source  : PDB ID 或 .pdb 檔路徑
    ligand  : 配體 resname（空白則自動偵測）
    cutoff  : H-bond D–A 距離上限（Å）
    """
    global HB_CUTOFF
    HB_CUTOFF = float(cutoff)

    out_dir = setup_output_dir()
    print(f"\n{'='*55}")
    print(f"  PyMOL H-Bond Analyzer")
    print(f"  來源：{source}  配體：{ligand or '自動偵測'}  上限：{HB_CUTOFF} Å")
    print(f"{'='*55}\n")

    # ── 載入結構 ──
    p = Path(source)
    if p.is_file():
        pdb_id = p.stem.upper()
        cmd.load(str(p), pdb_id)
        print(f"✓ 載入本地檔案：{p}")
    else:
        pdb_id = source.upper()
        print(f"⬇ 從 RCSB 下載：{pdb_id}")
        cmd.fetch(pdb_id, async_=0)

    # 移除水分子（保留蛋白質 + 配體）
    cmd.remove(f"({pdb_id}) and (resn HOH+WAT+H2O+DOD)")

    # ── 偵測配體 ──
    lig_resname = ligand.upper() if ligand else detect_ligand(pdb_id)
    if not lig_resname:
        print("⚠ 警告：未找到配體，僅執行蛋白質視覺化")
        cmd.show("cartoon", pdb_id)
        cmd.color("spectrum", pdb_id)
        cmd.bg_color("black")
        return

    lig_sel = f"({pdb_id}) and resn {lig_resname}"
    print(f"✓ 配體：{lig_resname}")

    # ── 計算氫鍵 ──
    hb_residues = get_hbond_residues(pdb_id, lig_sel)
    print_hbond_report(pdb_id, lig_resname, hb_residues)

    # ── 視覺化 ──
    print("⚙ 設定視覺化...")
    setup_visualization(pdb_id, lig_sel, hb_residues)

    # ── 匯出 ──
    print("📸 匯出圖像...")
    img_paths = export_images(pdb_id, lig_sel, out_dir)

    print("💾 儲存工作階段...")
    save_session(pdb_id, out_dir)
    save_json(pdb_id, lig_resname, hb_residues, out_dir)

    print(f"\n✅ 分析完成！輸出目錄：{out_dir.resolve()}")
    return hb_residues


# ── PyMOL 指令註冊 ───────────────────────────────────────────────────────────

if PYMOL_AVAILABLE:
    cmd.extend("analyze", analyze)
    print("✓ 指令 'analyze' 已載入")
    print("  用法：analyze <PDB_ID_或_檔案路徑> [, ligand=<RESNAME>] [, cutoff=3.5]")


# ── 命令列進入點 ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # pymol -c pymol_automation.py -- 1HSG [ligand] [cutoff]
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if not args:
        print("用法：pymol -c pymol_automation.py -- <PDB_ID_或_檔案> [配體] [距離上限]")
        sys.exit(0)

    src   = args[0]
    lig   = args[1] if len(args) > 1 else ""
    cut   = float(args[2]) if len(args) > 2 else 3.5

    if PYMOL_AVAILABLE:
        analyze(src, lig, cut)
    else:
        print("請在 PyMOL 環境中執行此腳本")
        print("  pymol -c pymol_automation.py -- 1HSG")
