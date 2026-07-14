# Protein-Ligand HB Interaction Analyzer

蛋白質–配體 **氫鍵（Hydrogen Bond）交互作用** 分析與突變位點建議工具。
輸入 PDB ID 或上傳結構檔，即可在瀏覽器中完成 3D 視覺化、氫鍵偵測與突變設計參考，全程無需安裝任何軟體。

**線上使用：** https://jong-liu.github.io/protein-ligand-HB-interaction/

---

## 功能特色

- **結構載入** — 輸入 PDB ID（自動從 RCSB 下載）或上傳 `.pdb` / `.ent` / `.cif` 檔案
- **3D 視覺化** — 由 [py3Dmol](https://3dmol.org/) 驅動，支援 Cartoon / Stick / Surface / Sphere 風格切換、標籤、連續旋轉與視角重置
- **氫鍵分析** — 偵測蛋白質與配體間 N/O/S 重原子的 Donor–Acceptor 距離 ≤ cutoff（預設 **3.6 Å**），並逐殘基去重取最短距離
- **突變位點建議** — 依 20 種胺基酸的突變策略字典評估功能影響，僅計入 **side chain** 交互（排除 backbone N/CA/C/O/OXT）
- **多格式匯出** — HTML、JSON、CSV、PDF（html2canvas + jsPDF，支援中文）、PyMOL 指令碼
- **格式支援** — PDB 固定欄位格式與 mmCIF `_atom_site` loop（auth_* 優先）

---

## 快速開始

線上版（推薦）直接開啟 [https://jong-liu.github.io/protein-ligand-HB-interaction/](https://jong-liu.github.io/protein-ligand-HB-interaction/) 即可使用。

本機使用只需下載 `index.html` 後用瀏覽器開啟——它是單一檔案、純瀏覽器端 JavaScript，無伺服器依賴。所有外部函式庫皆透過 CDN 動態載入。

操作流程：

1. 輸入 PDB ID（例如 `1AKE`）或上傳結構檔
2. 填入配體 resname（大小寫不拘）與氫鍵 cutoff（預設 3.6 Å）
3. 點擊分析，於右側檢視 3D 結構、氫鍵清單與突變位點建議
4. 依需求匯出 HTML / JSON / CSV / PDF / PyMOL 指令碼

---

## 突變位點判定原則

| 優先級 | 條件 | 意義 |
|--------|------|------|
| **高優先**（紅） | 蛋白質 **side chain** 原子與配體形成氫鍵 | 突變直接破壞氫鍵，對結合影響最大 |
| **中優先**（橘） | 5 Å 內近鄰、無 side chain 氫鍵（含僅 backbone 者） | 可能有疏水接觸或空間效應，可探索構形影響 |

> **Backbone 排除規則：** backbone 原子（`N`、`CA`、`C`、`O`、`OXT`）為所有胺基酸共有，無論突變成何種殘基都不會消失，因此其形成的氫鍵不會被突變破壞，不列入高優先建議。

---

## 專案結構

```
protein-ligand-HB-interaction/
├── index.html          ← 主交付物（GitHub Pages 唯一追蹤檔）
├── hb_analyzer.py      ← Python CLI（本機離線用，不推送）
├── pymol_automation.py ← PyMOL API 腳本（本機離線用，不推送）
├── CLAUDE.md           ← 開發規範與 Workflow 定義
├── README.md           ← 本文件
└── .gitignore          ← 僅追蹤 index.html / CLAUDE.md / README.md / .gitignore
```

---

## 外部依賴（CDN）

| 函式庫 | 用途 |
|--------|------|
| [py3Dmol](https://3dmol.org/build/3Dmol-min.js) | 3D 分子視覺化 |
| [html2canvas](https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js) | PDF 截圖 |
| [jsPDF](https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js) | PDF 生成 |

---

## 開發

開發規範、Skills、Agents 與標準 Git Workflow 詳見 [`CLAUDE.md`](./CLAUDE.md)。核心原則：

- 只推送 `index.html`、`CLAUDE.md`、`README.md`、`.gitignore`
- 修改 `index.html` 一律使用精準差異編輯，避免整檔改寫造成截斷
- 外部套件只用 CDN 動態載入，不使用 `pip` / `npm`

---

## 授權

本專案供學術與教學研究使用。結構資料來源為 [RCSB Protein Data Bank](https://www.rcsb.org/)。
