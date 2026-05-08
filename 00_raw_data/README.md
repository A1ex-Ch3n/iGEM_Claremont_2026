# 00_raw_data — Raw Genome Database / 原始基因组数据库

This folder is the **single source of truth** for all raw nucleotide sequences used in the iGEM 2026 active-learning phage–host engineering pipeline. Every downstream step reads from here; nothing in this folder is generated — it is either downloaded from NCBI or validated as read-only.

本目录是 iGEM 2026 主动学习噬菌体-宿主工程 pipeline 中所有原始核苷酸序列的**单一来源（single source of truth）**。所有下游步骤都从这里读取数据；本目录中的文件全部来自 NCBI 下载或验证，绝不由其他步骤生成。

---

## Role in the Active-Learning Pipeline / 在主动学习 Pipeline 中的角色

Module 00 fills three roles simultaneously:

1. **Primary annotation substrate** — `phage/*/genome.fna` feeds PHANOTATE (Module 02); `bacteria/*/genome.fna` feeds Prodigal (Module 02).
2. **ESM-2 embedding substrate** — `phage/*/proteins.faa` and `bacteria/*/proteins.faa` feed the protein embedding step (Module 04) as a reference.
3. **Contrastive prior corpus (Layer 3 of the 6-layer data-scarcity strategy)** — The 443 non-*Xanthomonas* phages in the negative pool provide a structural diversity prior, helping the model distinguish Xanthomonas-specificity from general phage sequence features.

模块 00 同时承担三个角色：

1. **注释底层输入** —— `phage/*/genome.fna` 输入 PHANOTATE（模块 02）；`bacteria/*/genome.fna` 输入 Prodigal（模块 02）。
2. **ESM-2 嵌入底层输入** —— `phage/*/proteins.faa` 和 `bacteria/*/proteins.faa` 作为参考提供给蛋白嵌入步骤（模块 04）。
3. **对比先验语料库（6 层数据稀缺策略中的第 3 层）** —— 443 个非 *Xanthomonas* 宿主噬菌体（negative pool）提供结构多样性先验，帮助模型区分 Xanthomonas 特异性与普通噬菌体序列特征。

---

## Directory Layout / 目录结构

```
00_raw_data/
├── README.md                    ← this file / 本文件
├── MANIFEST.csv                 ← global checksum index (INTERFACE.md §Universal conventions)
├── data_needs.md                ← what still needs to be fetched / 尚需补充的数据
├── bacteria_list.csv            ← all bacteria with accessions, GCF mappings, and notes
├── phage_list.csv               ← all 777 phages with source tags (canonical / negative_pool)
│
├── bacteria/                    ← one subdirectory per bacterium / 每个细菌一个子目录
│   └── <seq_accession>/
│       ├── genome.fna           ← whole-genome nucleotide sequence
│       ├── proteins.faa         ← NCBI-annotated protein sequences
│       ├── genes.gff            ← gene coordinates and annotation
│       └── cds.fna              ← CDS nucleotide sequences
│
├── phage/                       ← one subdirectory per phage / 每个噬菌体一个子目录
│   └── <accession>/             ← versioned accession, e.g. EU717894.1
│       ├── genome.fna           ← whole-genome nucleotide sequence
│       ├── proteins.faa         ← NCBI-annotated protein sequences (reference copy)
│       └── cds.fna              ← CDS nucleotide sequences
│
└── processes/
    ├── 01_verify_dataset.ipynb  ← audit notebook; generates MANIFEST.csv
    ├── fetch_bacteria.py        ← download script for bacteria
    ├── fetch_phages.py          ← download script for phages (777 total)
    ├── accession_to_gcf.json    ← lookup: seq accession → GCF assembly accession
    └── tests/
        ├── test_schema.py       ← MANIFEST.csv column/format validation
        ├── test_smoke.py        ← parse 3 smallest genomes with SeqIO
        └── test_sanity.py       ← phiL7 length check; list vs dir count audit
```

For the full data contract (column schemas, filename conventions), see [`INTERFACE.md`](../INTERFACE.md) §Module 00.

完整的数据契约（列名规范、文件命名约定），见 [`INTERFACE.md`](../INTERFACE.md) §Module 00。

---

## Canonical Reference Accessions / 规范参考 Accession

These three are the highest-priority genomes — all downstream dry-lab modules are anchored to them:

| Role | Organism | Accession |
|------|----------|-----------|
| RBP scaffold phage | *Xanthomonas* phage phiL7 | `EU717894.1` |
| Host receptor source | *X. campestris* pv. *campestris* ATCC 33913 (Xcc) | `GCF_000007145.1` |
| ELISA RBP positive control | T7 phage | `GCF_000840885.1` |

The phiL7 genome (44,080 bp) is the **reference scaffold** for the RBP active-learning loop. Its receptor system (TonB-ExbB-ExbD1D2, Wang et al. 2003) is experimentally validated — Module 03 should confirm that phiL7's tail-spike gp25 scores as the top RBP candidate.

phiL7 基因组（44,080 bp）是 RBP 主动学习回路的**参考 scaffold**。其受体系统（TonB-ExbB-ExbD1D2，Wang et al. 2003）已通过实验验证——模块 03 应确认 phiL7 尾刺蛋白 gp25 被评为最高分 RBP 候选。

---

## Dataset Summary / 数据集概览

| Category | Count | Notes |
|----------|-------|-------|
| Phage directories on disk | 774 | 777 in phage_list.csv; 3 missing (see data_needs.md) |
| Bacteria directories on disk | 34 | 49 rows in bacteria_list.csv; 3 invalid + unresolved remainder |
| Canonical Xanthomonas phages | 334 | Appear in interaction matrix |
| Negative-pool phages | 443 | Non-Xanthomonas; used as structural negatives for ML |
| Total files in MANIFEST.csv | 2398 | Verified SHA-256; see MANIFEST.csv |

---

## Organisms / 物种

### Bacteria — 34 directories / 细菌——34 个目录

| Source tag | Count | Description |
|------------|-------|-------------|
| `confirmed_matrix` | 34 | Have accessions directly in the interaction matrix |
| `newly_resolved` | 4 | Had `[No Complete Genome Found]`; accessions found via literature |
| `accession_invalid` | 2 | `KY000037` (Ti plasmid) and `PY746849` (patent sequence) — not genome assemblies |
| `unresolved_TODO` | 11 groups | Cannot be resolved without team input; see `data_needs.md` |

### Phages — 777 in list, 774 on disk / 噬菌体——清单 777 条，磁盘 774 个

| Source tag | Count | Description |
|------------|-------|-------------|
| `canonical_xanthomonas` | 334 | Core Xanthomonas-infecting phages; in the interaction matrix |
| `negative_pool` | 443 | Non-Xanthomonas phages (*E. coli*, *Pseudomonas*, *Salmonella*, etc.) |

> **Important / 注意:** The 443 negative-pool phages are NOT yet rows in the interaction matrix. Before ML training, explicit `0` labels must be added for all (negative_pool phage × Xanthomonas host) pairs. Do not treat their absence as missing data — it is structurally different.
>
> 443 个 negative-pool 噬菌体在互作矩阵中尚无对应行。在 ML 训练前，必须为所有（negative_pool 噬菌体 × Xanthomonas 宿主）组合添加显式 `0` 标签。不要将其缺失视为缺失数据——这在结构上是不同的。

---

## Known Issues / 已知问题

1. **3 missing phage directories** — `NC_013971.1`, `NZ_CP007800.1`, `NZ_CP008698.1` are in `phage_list.csv` but not on disk. Re-download with `fetch_phages.py`.
2. **Invalid bacteria accessions** — `KY000037` (Ti plasmid) and `PY746849` (patent sequence) cannot be downloaded as genomes. Requires matrix correction by Module 01 / Sarah.
3. **Duplicate accessions in interaction matrix** — `NZ_JBWJFR000000000` and `NZ_CP150073` each appear under two different organism names; this produces identical feature vectors. See `data_needs.md` for resolution options.
4. **Legacy `phage/MANIFEST.csv`** — A pre-existing MANIFEST.csv file in the `phage/` subdirectory uses a different (non-INTERFACE) schema. The canonical MANIFEST.csv is now at `00_raw_data/MANIFEST.csv`.

---

## Running the Verification Notebook / 运行验证 Notebook

```bash
conda activate igem2026
cd 00_raw_data/processes
jupyter lab 01_verify_dataset.ipynb
```

This regenerates `MANIFEST.csv` and prints any audit issues. Run the tests to confirm all checks pass:

重新生成 `MANIFEST.csv` 并打印审计问题。运行测试以确认所有检查通过：

```bash
python3 -m pytest 00_raw_data/processes/tests/ -v
```

---

## Running the Fetch Scripts / 运行下载脚本

All scripts are run from the **repository root / 从仓库根目录运行**:

```bash
# Download all bacteria (34 resolved genomes)
python 00_raw_data/processes/fetch_bacteria.py

# Download all 777 phages (can take 20–40 min)
python 00_raw_data/processes/fetch_phages.py

# Re-download a single failed accession
python 00_raw_data/processes/fetch_phages.py --accession NC_013971.1
```

Scripts are **resumable**: already-complete directories are silently skipped.

脚本支持**断点续传**：已完成的目录会被静默跳过。

---

## Data Flow / 数据流向

```
00_raw_data/phage/<acc>/genome.fna       ──→  02_annotation/  (PHANOTATE input)
00_raw_data/bacteria/<acc>/genome.fna    ──→  02_annotation/  (Prodigal input)
00_raw_data/phage/<acc>/proteins.faa     ──→  04_protein_embedding/  (reference)
00_raw_data/bacteria/<acc>/proteins.faa  ──→  04_protein_embedding/  (receptor ref)
00_raw_data/MANIFEST.csv                 ──→  01_data_ground_truth/  (download index)
```
