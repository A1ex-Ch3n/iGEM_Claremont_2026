# Factor 1 — GC Content / GC 含量差異

**Owner / 負責人:** Alex Chen  
**Output / 輸出:** `03_feature_weighting/outputs/per_factor/factor1/f01_gc_content.csv`

---

## Setup / 環境安裝

```bash
# Create env (one-time) / 建立環境（一次性）
conda env create -f shared/env/environment.yml
conda activate igem2026

# Register Jupyter kernel / 註冊 Jupyter kernel
python -m ipykernel install --user --name igem2026 --display-name "Python (igem2026)"

# Launch notebook / 開啟筆記本
jupyter lab 03_feature_weighting/processes/factor1_gc_content/factor1_gc_content.ipynb
```

Select kernel **Python (igem2026)**, then **Run All Cells**.  
選擇核心 **Python (igem2026)**，然後 **Run All Cells**。

---

## First-run notes / 首次執行說明

- **Cell 3** downloads 10 host bacterial genomes (~250 MB total) via `ncbi-datasets-cli`. This takes ~5 minutes and is only needed once — subsequent runs use the cached `.fna` files.
- **Cell 3** 會透過 `ncbi-datasets-cli` 下載 10 個宿主細菌基因組（約 250 MB）。首次約需 5 分鐘，後續執行將使用快取檔案，無需重複下載。

---

## Updating the truth table / 更新真值表

Edit **Cell 2 (Config)**, change `TRUTH_TABLE_PATH` to point to the new file, then re-run all cells.  
修改 **Cell 2（Config）** 中的 `TRUTH_TABLE_PATH` 指向新檔案，重新執行所有 Cell 即可。

---

## Expected outputs / 預期輸出

| File | Description / 說明 |
|------|-------------------|
| `per_genome_gc.csv` | GC fraction + length per genome / 每個基因組的 GC 比例與長度 |
| `f01_gc_content.csv` | Per-pair `x_gc` and `x_gc_zscore` / 每對配對的 GC 差值與 Z 分數 |
| `gc_density.png` | Distribution + scatter + box plots / 分布圖、散佈圖與箱型圖 |
