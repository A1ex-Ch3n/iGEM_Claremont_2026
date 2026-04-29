## Step 1: Exporting Data from NCBI | 第一步：从 NCBI 导出数据

- [ ] **Visit the NCBI Virus Database.** 访问 NCBI Virus 数据库。
- [ ] **Search for your target phages.** 搜索你的目标噬菌体（例如 "Xanthomonas phage"）。
- [ ] **Click the "Export" button.** 点击结果表格中的 "Export"（导出）按钮。
- [ ] **Select "CSV" format.** 选择 "CSV" 格式，并确保包含 "Accession"（编号）列。
- [ ] **Rename and Move.** 将下载的文件重命名为 `sequences.csv` 并移动到你的项目文件夹中。

---

## Step 2: Environment Setup | 第二步：环境准备

1. **Open your macOS Terminal.** 打开你的 Mac 终端。
2. **Activate your environment.** 激活你的 iGEM 专属环境。

```bash
conda activate igem_drylab
```

3. **Install libraries.** 安装必要的 Python 库。

```bash
pip install pandas biopython pyrodigal
```

---

## Step 3: Preparing the Master Script | 第三步：准备主脚本

在项目文件夹中创建一个名为 `master_pipeline.py` 的文件，并粘贴以下代码：

```python
import os
import pandas as pd
import pyrodigal
from Bio import Entrez, SeqIO
from pathlib import Path

# Configuration | 配置
Entrez.email = "your_email@example.com"  
INPUT_CSV = "sequences.csv"              
BASE_DATA_DIR = Path("./ncbi_dataset/data")

def run_pipeline():
    # Create directory | 创建目录
    BASE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    # Read CSV | 读取 CSV
    df = pd.read_csv(INPUT_CSV)
    accessions = df['Accession'].unique()
    print(f"🚀 Processing {len(accessions)} sequences...")

    # Initialize Engines | 初始化引擎
    single_finder = pyrodigal.GeneFinder(meta=False)
    meta_finder = pyrodigal.GeneFinder(meta=True)

    for acc in accessions:
        folder = BASE_DATA_DIR / acc
        folder.mkdir(exist_ok=True)
        fna_path = folder / f"{acc}.fna"
        faa_path = folder / "proteins.faa"

        # Download Stage | 下载阶段
        if not fna_path.exists():
            print(f"📥 [DOWNLOADING] {acc}...")
            try:
                handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                with open(fna_path, "w") as f:
                    f.write(handle.read())
                handle.close()
            except Exception as e:
                print(f"❌ {acc} Download failed: {e}")
                continue

        # Annotation Stage | 注释阶段
        if not faa_path.exists():
            print(f"🧬 [ANNOTATING] {acc}...")
            try:
                for record in SeqIO.parse(fna_path, "fasta"):
                    dna_seq = str(record.seq)
                    # Length check to prevent crashes | 长度检查防止崩溃
                    if len(dna_seq) >= 20000:
                        single_finder.train(dna_seq)
                        genes = single_finder.find_genes(dna_seq)
                    else:
                        genes = meta_finder.find_genes(dna_seq)
                    
                    with open(faa_path, "w") as out:
                        genes.write_translations(out, sequence_id=record.id)
                print(f"✨ [SUCCESS] {acc} finished.")
            except Exception as e:
                print(f"❌ {acc} Error: {e}")
        else:
            print(f"✅ [SKIPPED] {acc} already exists.")

if __name__ == "__main__":
    run_pipeline()
```

---

## Step 4: Execution | 第四步：执行脚本

1. 在终端中进入你的文件夹。
2. 运行主脚本：

```bash
python master_pipeline.py
```

> [!NOTE]
> 如果网络断开，只需再次运行脚本即可断点续传。

---

## Step 5: Verification of Outputs | 第五步：结果验证

- **Location:** 打开 `ncbi_dataset/data` 目录。
- **Contents:** 包含 `[Accession].fna` (DNA) 和 `proteins.faa` (蛋白质)。
- **Tools:** 你可以使用 VS Code 或任何文本编辑器打开 `.faa` 文件。