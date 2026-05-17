# CLAUDE.md — PyMOL HB Interaction Analyzer

> 本文件定義此專案的 **Skills、Agents 與 Workflow**。  
> Claude 在執行任何修改或新功能任務時，應先讀取此文件，再依據 Workflow 完成並推送。

---

## 1. 專案概述

**GitHub Pages URL：** `https://jong-liu.github.io/pymol-HB-interaction/`  
**本地路徑：** `E:\AI\claude code\projects\pymol\pymol-HB-interaction\`  
**主要交付物：** `index.html`（單一檔案，瀏覽器端純 JavaScript，無伺服器依賴）

### 核心功能

| 功能 | 說明 |
|------|------|
| 蛋白質結構載入 | 輸入 PDB ID（自動從 RCSB 取得）或上傳 `.pdb` / `.ent` / `.cif` 檔案 |
| 3D 視覺化 | py3Dmol（CDN）— cartoon / stick / surface / sphere 風格切換 |
| 氫鍵分析 | D–A 重原子距離 ≤ cutoff（預設 **3.6 Å**），逐殘基去重 |
| 突變位點建議 | 20 種胺基酸的突變策略字典 + 功能影響評估 |
| 報告匯出 | HTML、JSON、CSV、PDF（html2canvas + jsPDF）、PyMOL 指令碼 |
| 格式支援 | PDB 固定欄位格式 + mmCIF `_atom_site` loop 格式（auth_* 優先） |

---

## 2. 檔案結構

```
pymol-HB-interaction/
├── index.html            ← 主交付物（GitHub Pages 唯一追蹤檔案）
├── hb_analyzer.py        ← Python CLI（本機離線用，不推送）
├── pymol_automation.py   ← PyMOL API 腳本（本機離線用，不推送）
├── .gitignore            ← 排除所有 .py、輸出檔、session 等
├── CLAUDE.md             ← 本文件（推送至 git）
└── README.md             ← 簡短說明（推送至 git）
```

### `.gitignore` 策略

只推送 `index.html`、`CLAUDE.md`、`README.md`、`.gitignore`。  
所有 `.py` 檔、`pymol_output/`、`*_hb_report.html`、`*.json` 均排除。

---

## 3. index.html 架構速查

### 全域變數

```javascript
let viewer          // py3Dmol 3D 檢視器物件
let pdbData         // 目前載入的原始結構文字（PDB 或 CIF）
let pdbId           // 當前 PDB ID（大寫）
let ligResname      // 配體 resname（保留原始大小寫，比對時 toLowerCase()）
let allAtoms        // 解析後的原子陣列
let hbonds          // 氫鍵結果陣列
let nearbyRes       // 5Å 內近鄰殘基陣列
let labelsOn        // 標籤開關
let spinning        // 旋轉開關
let currentStyle    // 目前 3D 風格（'cartoon'|'stick'|'surface'|'sphere'）
let surfaceActive   // Boolean：surface layer 是否已 addSurface（防疊加用）
let uploadedPdbText // 上傳檔案的文字內容
let uploadedFileName// 上傳檔案名稱
let fileFormat      // 'pdb' 或 'cif'
```

### 關鍵函數

| 函數 | 說明 |
|------|------|
| `parsePDB(text)` | PDB 固定欄位解析器，回傳原子物件陣列 |
| `parseCIF(text)` | mmCIF `_atom_site` loop 狀態機解析器（auth_* 優先） |
| `cifTokenize(line)` | mmCIF 引號字串 tokenizer |
| `findHBonds(atoms, lig, cutoff)` | D–A 距離 ≤ cutoff，逐殘基去重 |
| `findNearby(atoms, lig, cutoff=5.0)` | 5Å 近鄰殘基 |
| `runAnalysis(id, rawPdb, ligName, cutoff, chainFilter, fmt)` | 分析主流程，配體比對 case-insensitive |
| `renderViewer(pdb, fmt)` | py3Dmol 初始化與格式傳入 |
| `applyStyle(style)` | 套用 3D 風格；執行前先移除舊 surface（防疊加）；以 id 更新按鈕 active 狀態 |
| `setStyle(s)` | 互斥切換：若已 active 同一樣式則回預設 cartoon；否則套用新樣式 |
| `toggleLabels()` | 切換 H-Bond 殘基標籤顯示／隱藏；再按一次 → 關閉（labelsOn 預設 true） |
| `toggleSpin()` | 連續旋轉開關：再按停止；同步更新 `btn-spin` active 狀態 |
| `resetView()` | 停止旋轉 + 視角歸位：有配體以配體為中心，否則 zoom all |
| `fetchAndAnalyze()` | 主按鈕入口：PDB ID 下載 或 使用 uploadedPdbText |
| `readFile(file)` | 上傳處理，偵測 .cif，更新 span 文字，重置 input.value |
| `buildPrintContainer()` | 建立隱藏 794px 報告 div 供 PDF 截圖 |
| `exportPDF()` | 動態載入 html2canvas + jsPDF，分頁截圖輸出 |
| `loadScriptIfNeeded(url, checkFn)` | 動態 script 載入（避免重複） |

### H-Bond 參數

- 預設 cutoff：**3.6 Å**（共 4 處 JS fallback + 1 處 HTML input value）
- 篩選原子：蛋白質 `N`/`O`/`S` 與配體 `N`/`O`/`S` 之間的重原子距離
- 去重邏輯：同一殘基保留最短距離的配對

### 3D 視覺化按鈕對照表（2026-05-17 v2 — toggle + reset-first 版）

| 按鈕 | id | onclick | 行為 | 再按一次 |
|------|----|---------|------|---------|
| Cartoon | `btn-cartoon` | `setStyle('cartoon')` | ribbon 螺旋與 sheet（預設） | 已是預設，維持不變 |
| Stick | `btn-stick` | `setStyle('stick')` | 所有殘基鍵棒顯示 | 回預設 cartoon |
| 表面 | `btn-surface` | `setStyle('surface')` | VDW 表面電位圖（防重複 addSurface） | 回預設 cartoon，移除 surface |
| Sphere | `btn-sphere` | `setStyle('sphere')` | 球體 VDW 半徑顯示 | 回預設 cartoon |
| 標籤 | `btn-labels` | `toggleLabels()` | H-Bond 殘基標籤顯示（預設開） | 關閉標籤 |
| 旋轉 | `btn-spin` | `toggleSpin()` | 連續旋轉 `viewer.spin('y',1)` | 停止旋轉 |
| Reset | `btn-reset` | `resetView()` | 停止旋轉 + zoom 至配體（無配體則 zoom all） | — |

> **按鈕設計原則（2026-05-17 v2）：**
> 1. **不疊加**：`applyStyle()` 執行前先 `viewer.removeAllSurfaces()` + `viewer.setStyle({},{})` 清除舊狀態
> 2. **toggle 解除**：再按同一個 active 按鈕 → 回預設 cartoon（style 類）或關閉（labels/spin 類）
> 3. **ID 精準定位**：所有按鈕改以 `id="btn-*"` 管理 active class，避免舊版 index 錯位問題
> 4. 配體無論在何種模式下，固定以 `stick+sphere` 雙重風格顯示

---

## 4. Skills（技能定義）

### Skill A — `add_feature`（新增功能）

**觸發時機：** 使用者說「新增 X 功能」、「加入 Y」、「支援 Z」

**步驟：**
1. 讀取 `index.html` 相關段落（用 `offset` + `limit` 分段讀取）
2. 規劃修改範圍（影響哪些函數、HTML 區塊、CSS）
3. 使用 `Edit` 工具進行精準差異修改（避免整檔改寫）
4. 執行 [Workflow → Git Push](#6-標準-workflow)

**常見新功能清單：**
- 新增 3D 視覺化風格
- 新增匯出格式（Excel、SVG）
- 新增鏈（Chain）篩選 UI 元件
- 新增水分子顯示開關
- 新增殘基搜尋功能

---

### Skill B — `fix_bug`（修正錯誤）

**觸發時機：** 使用者貼上截圖或說明錯誤行為

**步驟：**
1. 用 `Grep` 找到問題所在行號
2. `Read` 讀取前後 20 行上下文
3. 診斷根因後，用 `Edit` 進行最小範圍修正
4. 執行 [Workflow → Git Push](#6-標準-workflow)

**歷史已知問題（勿重複）：**
- 上傳後再分析：確認 `fetchAndAnalyze()` 有 `uploadedPdbText` fallback 分支
- 重複上傳不生效：確認 `readFile()` 更新的是 `<span>` 而非整個 label 的 innerHTML，且有 `input.value = ""`
- 配體大小寫：`runAnalysis()` 使用 `toLowerCase()` 比對，`ligResname` 儲存原始大小寫
- PDF 截圖文字：使用 html2canvas 而非純文字 jsPDF（支援中文）
- **3D 按鈕錯誤（2026-05-17 修正）**：舊版按鈕缺少 Stick / Sphere，多了「Cartoon+表面」且不在預期清單中；`toggleSpin()` 誤用 `viewer.rotate()`（單次旋轉）而非 `viewer.spin('y',1)`（連續旋轉）。已全數修正，按鈕順序改為 Cartoon → Stick → 表面 → Sphere → 標籤 → 旋轉 → 重置視角。
- **按鈕功能疊加 + toggle + Reset 無效（2026-05-17 v2 修正）**：① Surface 多次點擊會 addSurface 堆疊 → 加入 `surfaceActive` flag 防重複；② style 按鈕再按一次無法回預設 → `setStyle()` 加 toggle 邏輯（`currentStyle===s` 時回 cartoon）；③ 標籤/旋轉按鈕缺乏 active 視覺狀態 → 改以 `id="btn-*"` 精準更新；④ 「重置視角」改名為 Reset，`resetView()` 修正為先停旋轉再 zoom（有配體用 `zoomTo({resn:ligResname})` 否則 `zoomTo({})` zoom all）。

---

### Skill C — `update_param`（修改參數/預設值）

**觸發時機：** 使用者說「把 X 改為 Y」、「預設值改成」

**步驟：**
1. 用 `Grep` 找出所有受影響的出現位置
2. 依 `replace_all` 或逐一 `Edit` 修改（若語境不同則逐一處理）
3. 執行 [Workflow → Git Push](#6-標準-workflow)

**H-Bond cutoff 修改位置（共 5 處）：**

| 位置 | 程式碼片段 |
|------|-----------|
| HTML input | `<input ... value="3.6" ...>` |
| `fetchAndAnalyze()` | `... \|\| 3.6;` |
| `readFile()` | `... \|\| 3.6;` |
| `buildPrintContainer()` | `D–A ≤ 3.6 Å` 標題文字 |
| `copyPyMOL()` | `...)\|\|3.6;` |

---

### Skill D — `style_change`（介面樣式調整）

**觸發時機：** 使用者說「改變顏色」、「調整排版」、「字體太小」

**步驟：**
1. 定位 `<style>` 區段（通常在 index.html 前 350 行）
2. 使用 `Edit` 修改對應 CSS class
3. 執行 [Workflow → Git Push](#6-標準-workflow)

---

### Skill E — `add_export`（新增匯出格式）

**觸發時機：** 使用者說「新增匯出 X 格式」

**步驟：**
1. 在 `<div class="export-bar">` 新增按鈕 HTML
2. 在 JS 區段新增對應 `exportXXX()` 函數
3. 若需外部函式庫，使用 `loadScriptIfNeeded()` 動態載入
4. 執行 [Workflow → Git Push](#6-標準-workflow)

---

## 5. Agents（代理人定義）

### Agent 1 — Feature Developer

**負責：** 新增、修改功能程式碼  
**工具：** `Read`, `Edit`, `Write`, `Grep`, `mcp__workspace__bash`  
**限制：** 只修改 `index.html`（不動 .py 檔，不動 .gitignore）

### Agent 2 — Git Publisher

**負責：** 確認修改正確後，執行 git commit + push  
**工具：** `mcp__workspace__bash`  
**前置條件：** Feature Developer 完成 Edit，且 Claude 確認語法無誤

### Agent 3 — QA Verifier

**負責：** 修改後的自我審查  
**工具：** `Read`, `Grep`, `mcp__workspace__bash`  
**執行內容：**
- 確認所有修改的行號語法正確
- `grep` 驗證不再含有舊值
- 確認 HTML 標籤無未閉合

---

## 6. 標準 Workflow

每次修改任務遵循以下 6 步驟：

```
┌─────────────────────────────────────────────────────┐
│  STEP 1  理解需求                                    │
│  閱讀使用者描述，對應到 Skill A–E，確認修改範圍       │
├─────────────────────────────────────────────────────┤
│  STEP 2  定位程式碼                                  │
│  Grep 關鍵字 → Read 相關行 → 確認上下文              │
├─────────────────────────────────────────────────────┤
│  STEP 3  執行修改（Feature Developer Agent）         │
│  Edit 工具精準差異修改，避免整檔改寫                 │
├─────────────────────────────────────────────────────┤
│  STEP 4  自我驗證（QA Verifier Agent）               │
│  Grep 確認舊值已消失，Read 確認語法，                │
│  bash wc -l 確認行數合理                             │
├─────────────────────────────────────────────────────┤
│  STEP 5  Git Push（Git Publisher Agent）             │
│  執行下方標準 git 指令                               │
├─────────────────────────────────────────────────────┤
│  STEP 6  回報結果                                    │
│  提供 GitHub Pages URL 與 commit message             │
└─────────────────────────────────────────────────────┘
```

### 標準 Git Push 指令

```powershell
cd "E:\AI\claude code\projects\pymol\pymol-HB-interaction"
git add index.html
git commit -m "<類型>: <說明>"
git push
```

**Commit 類型前綴：**

| 類型 | 使用時機 |
|------|---------|
| `feat` | 新增功能 |
| `fix` | 修正錯誤 |
| `refactor` | 重構（不影響功能） |
| `style` | 純介面樣式調整 |
| `docs` | 文件更新 |
| `chore` | 維護性更改（參數調整等） |

### Bash 路徑對照（Linux Sandbox）

| Windows 路徑 | Linux Sandbox 路徑 |
|-------------|------------------|
| `E:\AI\claude code\projects\pymol\pymol-HB-interaction\` | `/sessions/eager-nice-hamilton/mnt/pymol-HB-interaction/` |
| 輸出暫存 | `/sessions/eager-nice-hamilton/mnt/outputs/` |

---

## 7. 功能擴充路線圖（Roadmap）

以下為潛在功能，可用對應 Skill 直接執行：

| 優先級 | 功能 | Skill | 難度 |
|--------|------|-------|------|
| 高 | 多配體同時分析 | A | 中 |
| 高 | 氫鍵網絡圖（D3.js 或 Mermaid） | A/E | 中 |
| 中 | 水分子（HOH）氫鍵顯示開關 | A | 低 |
| 中 | 殘基序列搜尋高亮 | A | 低 |
| 中 | Excel (.xlsx) 匯出 | E | 中 |
| 中 | 疏水接觸分析（≤ 4.0 Å C-C） | A | 中 |
| 低 | 多結構疊合比較 | A | 高 |
| 低 | AlphaFold DB 整合（直接輸入 UniProt ID） | A | 中 |
| 低 | 介面語言切換（中/英） | D | 低 |

---

## 8. 外部依賴（CDN）

| 函式庫 | URL | 用途 |
|--------|-----|------|
| py3Dmol | `https://3dmol.org/build/3Dmol-min.js` | 3D 分子視覺化 |
| html2canvas | `https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js` | PDF 截圖 |
| jsPDF | `https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js` | PDF 生成 |

所有外部函式庫均透過 `loadScriptIfNeeded(url, checkFn)` 動態載入，不預先阻塞頁面。

---

## 9. 開發注意事項

### 修改 index.html 的原則

1. **不整檔改寫**：使用 `Edit` 工具精準差異修改，避免 `Write` 覆蓋全檔（風險：截斷）
2. **分段讀取**：大型讀取用 `offset` + `limit` 分段，避免超出上下文
3. **不破壞函數邊界**：每次 Edit 的 `old_string` 需包含足夠上下文使其唯一
4. **CDN 優先**：禁用 `pip install` / `npm install` 方式引入套件，只用 CDN

### HTML 文件結構

```
index.html
├── <head>        行  1–  30  meta, title, CDN preload
├── <style>       行 30– 370  CSS 變數、RWD、動畫
├── <body>        行 370– 400  layout 骨架
│   ├── .sidebar  行 400– 760  左側控制面板
│   └── #viewer-wrapper  右側 3D 視窗 + 分析面板
├── <script>      行 760–1395  所有 JS 邏輯
└── </html>       行 1395–1407
```

### 測試 Checklist（每次修改後）

- [ ] `grep -c "syntax error" index.html` = 0
- [ ] 開啟 GitHub Pages URL，Console 無 JS 錯誤
- [ ] 輸入 `1HSG`，點擊分析，3D 結構正常顯示
- [ ] 氫鍵表格正確列出 MK1 配體的 HB 殘基
- [ ] PDF 匯出正常（按鈕點擊 → 下載）
- [ ] 上傳 `.cif` 檔案後可正常分析

---

## 10. 快速指令參考

### 查看當前 H-Bond cutoff 值

```bash
grep -n "3\.[0-9]" /sessions/eager-nice-hamilton/mnt/pymol-HB-interaction/index.html | grep -i cutoff
```

### 查看所有匯出函數

```bash
grep -n "^function export" /sessions/eager-nice-hamilton/mnt/pymol-HB-interaction/index.html
```

### 確認最新 git 狀態

```bash
cd /sessions/eager-nice-hamilton/mnt/pymol-HB-interaction && git log --oneline -5 && git status
```

### 行數確認（正常約 1400–1500 行）

```bash
wc -l /sessions/eager-nice-hamilton/mnt/pymol-HB-interaction/index.html
```
