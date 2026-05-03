# 📘 iGEM Drylab: The Complete "Zero-to-Hero" Manual
# 📘 iGEM 干实验：从零到一全能手册

> **Guide Intro**: This document is your master reference. It contains everything we discussed: setup, logic, the script, and how to verify your success.
> **手册简介**: 这是你的核心参考文档。它包含了我们讨论过的所有内容：环境搭建、逻辑说明、代码脚本以及如何验证结果。

---

## 🛠️ Phase 1: Environment Configuration (The Engine)
## 🛠️ 第一阶段：环境配置（安装引擎）

Before running any code, your Mac needs to have the right "tools" installed in a specific room called `igem_drylab`. 
在运行代码前，你的 Mac 需要在一个名为 `igem_drylab` 的独立空间里装好工具。

### 1. Create and Activate / 创建并激活
Open your Terminal and run these lines. This ensures we don't mess up your computer's default settings.
打开终端并运行以下命令。这能确保我们不会搞乱你电脑的默认设置。

```bash
# Create a fresh environment with Python 3.12
# 创建一个带有 Python 3.12 的全新环境
conda create -n igem_drylab python=3.12 -y

# Enter the environment (You should see the name on the left)
# 进入环境（你应该能在左侧看到环境名）
conda activate igem_drylab
```

### 2. Install the Biological "Brain" / 安装生物学“大脑”
We use `pyrodigal` for gene prediction and `biopython` for data reading. Since you are on a Mac (M1/M2/M3), we use the `conda-forge` channel for better compatibility.
我们使用 `pyrodigal` 进行基因预测，用 `biopython` 读取数据。因为你使用的是 Mac (M系列芯片)，我们加入 `conda-forge` 频道以获得更好的兼容性。

```bash
# One-line install with all necessary channels
# 一行命令安装所有频道
conda install -c conda-forge -c bioconda pyrodigal biopython -y
```

---

## 🚀 Phase 2: The Core Python Script (The Magic)
## 🚀 第二阶段：核心 Python 脚本（自动化魔法）

**Instruction**: Save the code below as `batch_prodigal.py`. It is designed to scan folders and extract proteins automatically.
**说明**：将以下代码保存为 `batch_prodigal.py`。它会自动扫描文件夹并提取蛋白质。

```python
import pyrodigal
from Bio import SeqIO
from pathlib import Path

def run_prodigal_pipeline(target_directory):
    """
    Main Logic: Find DNA (.fna) -> Train AI -> Predict Genes -> Save Proteins (.faa)
    核心逻辑：找 DNA -> 训练 AI -> 预测基因 -> 保存蛋白质
    """
    base_path = Path(target_directory)
    # Recursively find all genome files (.fna) in subfolders
    fna_files = list(base_path.rglob("*.fna"))
    
    if not fna_files:
        print(f"❌ Error: No .fna files found in: {target_directory}")
        return

    print(f"🔍 Found {len(fna_files)} genome(s). Starting annotation...")
    
    # Initialize GeneFinder (Specific for pyrodigal v3.0+)
    gene_finder = pyrodigal.GeneFinder(meta=False)

    for fna_path in fna_files:
        # Define output path: it creates 'proteins.faa' in the SAME folder
        output_faa = fna_path.parent / "proteins.faa"
        
        # Smart Check: Skip if already done to save time
        if output_faa.exists():
            print(f"✅ [SKIPPED] {fna_path.name} already processed.")
            continue
            
        print(f"⏳ [WORKING] Extracting proteins from: {fna_path.name} ...")
        
        with open(output_faa, "w") as out_faa:
            for record in SeqIO.parse(fna_path, "fasta"):
                dna_seq = str(record.seq)
                
                # CRITICAL: Self-training is required for single-species mode
                # 关键：单物种模式必须先进行自训练
                gene_finder.train(dna_seq)
                
                # Predict genes and write the amino acid translations
                genes = gene_finder.find_genes(dna_seq)
                genes.write_translations(out_faa, sequence_id=record.id)
                
        print(f"✨ [SUCCESS] Proteins successfully saved to: {output_faa}")

if __name__ == "__main__":
    # >>> ACTION REQUIRED: Change the path below to your DATA FOLDER path <<<
    MY_DATA_PATH = "/Users/oliviayan/Desktop/ASM714v1/ncbi_dataset/data/GCF_000007145.1/"
    
    run_prodigal_pipeline(MY_DATA_PATH)
```

---

## 🔍 Phase 3: How to Run & Witness the Results
## 🔍 第三阶段：如何运行与见证成果

### 1. Execution Steps (In VS Code) / 运行步骤
* **Open the file**: Open `batch_prodigal.py` in your editor. / **打开文件**：在编辑器中打开代码。
* **Select Environment**: Look at the **bottom-right corner** of VS Code. Click the Python version and select `igem_drylab (3.12.x)`. This is like putting on the right glasses to see the code. / **选择环境**：看 VS Code 右下角，点击并选择 `igem_drylab`。这能确保编辑器使用你装好插件的那个 Python。
* **Run**: Click the **Triangle (Play)** icon in the top-right corner. / **点击运行**：点击右上角的三角形图标。

### 2. Witnessing the Results / 见证成果
* **Go to your data folder**: You should see a brand new file named `proteins.faa`. / **去你的数据文件夹**：你会看到一个全新的文件 `proteins.faa`。
* **What's inside?**: Open it with any text editor (even VS Code). You should see sequences starting with `>` followed by letters like `MKRIST...`. These are the protein building blocks of your bacterium! / **里面是什么？**：用文本编辑器打开它。你应该看到以 `>` 开头的序列，这些就是你细菌的蛋白质序列！

---

## 💡 Troubleshooting (If things turn red)
## 💡 故障排除（如果报错了）

* **"No module named pyrodigal"**: It means VS Code is still using the default Python. Double check the bottom-right corner! / **说明环境没选对**：确保右下角显示的是 `igem_drylab`。
* **"No .fna files found"**: Check your path in the script. Ensure it's the folder path. / **找不到文件**：检查代码里的路径，确保它指向的是文件夹。

🎉 **Happy Coding for iGEM! 🚀**