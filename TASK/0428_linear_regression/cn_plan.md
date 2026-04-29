## 1. 專案概述

本研究旨在利用生物資訊計算方法，提取細菌（Host）與噬菌體（Phage）的序列特徵，並透過迴歸分析（Regression）建立預測模型。模型之自變數（$x$）為兩者特徵的關聯值，因變數（$y$）為感染可能性（Likelihood of Infection）。

## 2. 團隊分工與因子定義表

|**編號**|**負責人**|**因子名稱 (Factor)**|**數據來源**|**生物學意義**|**x 的計算方式 (Phage p vs. Host b)**|
|---|---|---|---|---|---|
|**1**|**Alex**|**GC Content**|DNA Genome|基因組組成的匹配度。GC 越接近通常代表長期共演化。|**絕對差值：** $x = \\|GC_p - GC_b$|
|**2**|**Olivia**|**pI & Acidity**|**.faa (Protein)**|蛋白質組電荷分佈與酸性胺基酸比例。影響蛋白在宿主體內的穩定性。|**Acidity 差值：** $x = \\|Acid\%_p - Acid\%_b$<br><br>  <br><br>（計算 D, E 殘基佔比）|
|**3**|**Weitao**|**Size (Length)**|.faa (Protein)|蛋白質組平均長度。反映噬菌體基因組精簡化與複製策略。|**平均值比值：** $x = \text{MeanLength}_p / \text{MeanLength}_b$|
|**4**|**Sarah**|**GRAVY**|.faa (Protein)|蛋白質平均疏水性。反映兩者對環境（如溫度、膜結構）的適應。|**絕對差值：** $x = \\|GRAVY_p - GRAVY_b$|
|**5**|**Angela**|**Similarity %**|.faa (Protein)|序列同源性。檢測水平基因轉移（HGT）或輔助代謝基因（AMGs）。|**歸一化比例：** $x = \frac{\text{Shared Genes}}{\text{Total Phage Genes}}$|
|**6**|**Carol**|**CAI**|DNA + .faa|密碼子適應指數。衡量噬菌體對宿主翻譯系統的徵用效率。|**直接計算值：** $x = \text{Mean}(CAI_{\text{phage based on host}})$|

---

## 3. 數據提取詳細規範

### A. Olivia 的 Acidity 與 pI 計算 (來自 .faa)

- **pI：** 使用 `Biopython` 或 `R (Peptides package)` 計算所有蛋白質的等電點，取中位數。
    
- **Acidity (酸鹼度)：** 統計序列中天門冬胺酸 (D) 與麩胺酸 (E) 的總和百分比。
    
    - $x$ 計算邏輯：若噬菌體的酸性胺基酸佔比與宿主高度一致，代表其蛋白質在宿主細胞質的溶解度與電荷環境適應性較佳。
        

### B. Sarah 的 GRAVY 計算 (來自 .faa)

- 計算公式為：$\frac{\sum \text{殘基疏水性值}}{\text{序列總長度}}$。
    
- 計算後建議將 Phage 與 Host 的所有蛋白數值做成 **Density Plot**，觀察高峰位移情況。
    

### C. Carol 的 CAI 計算 (關鍵因子)

- **Step 1:** 從宿主基因組中找出核糖體蛋白（Ribosomal Proteins）作為 High-expression Set。
    
- **Step 2:** 計算該 Set 的 Codon Usage Table。
    
- **Step 3:** 拿噬菌體的每一條基因與此 Table 比對，算出 CAI。
    

---

## 4. 迴歸分析實作建議

### 模型目標

我們希望找到一個 $f(x)$ 使得 $y = ax + b$，其中 $y$ 是感染力的機率（0 到 1 之間）。

### 建議流程

1. **數據標準化 (Normalization)：** 由於 GC 是百分比（0.1 級別），而 Size 是長度（100 級別），在跑多元回歸前，必須將所有 $x$ 進行 Z-score 標準化。
    
2. **散佈圖檢查：** 在跑 Regression 之前，各組組員應先畫出各自因子與 $y$ 的散佈圖，確認是否為線性關係。
    
