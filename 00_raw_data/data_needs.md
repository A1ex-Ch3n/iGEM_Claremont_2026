# Data Needs — 00_raw_data / 数据需求

## Current Status (as of 2026-05-08 audit) / 当前状态（截至 2026-05-08 审计）

| Category | On Disk | In List | Gap |
|----------|---------|---------|-----|
| Phage directories | 774 | 777 | 3 missing (see §Missing phages) |
| Bacteria directories | 34 | 49 rows (49 - 3 invalid - 11 unresolved = ~35 real targets) | Several unresolved (see §Unresolved bacteria) |
| MANIFEST.csv | ✅ present | — | 2398 files indexed |

**Pointer:** Actual genome fetching is orchestrated by Module 01 (`01_data_ground_truth/`). This file documents *what* still needs to be fetched and *why*; see `01_data_ground_truth/processes/` for the download scripts.

**指针**：实际的基因组下载由模块 01（`01_data_ground_truth/`）统筹。本文件记录*还需要什么*以及*原因*；下载脚本见 `01_data_ground_truth/processes/`。

---

## File Types: What, Why, and How / 文件类型：内容、用途与获取方式

### For Bacteria (Hosts) / 细菌（宿主）

| File | Why / 用途 | How Obtained / 获取方式 |
|------|-----------|------------------------|
| `genome.fna` | Whole-genome nucleotide. Required for Prodigal (Module 02), GC content. / 全基因组核苷酸，Prodigal 输入（模块 02）、GC 含量计算所需。 | NCBI datasets CLI: `--include genome` |
| `proteins.faa` | NCBI-annotated proteins. ESM-2 receptor embeddings (Module 04). / NCBI 注释蛋白，用于 ESM-2 受体嵌入（模块 04）。 | `--include protein` |
| `genes.gff` | Gene coordinates for receptor identification. / 基因坐标，用于受体鉴别。 | `--include gff3` |
| `cds.fna` | CDS nucleotides. Ribosomal gene CDS for codon usage analysis. / CDS 核苷酸，用于密码子使用分析。 | `--include cds` |

### For Phages / 噬菌体

| File | Why / 用途 | How Obtained / 获取方式 |
|------|-----------|------------------------|
| `genome.fna` | **Required.** PHANOTATE input (Module 02), GC content. / **必须。** PHANOTATE 输入（模块 02）、GC 含量。 | NCBI datasets CLI: `--include genome` |
| `proteins.faa` | NCBI reference copy (Module 02 re-derives via PHANOTATE). / NCBI 参考副本（模块 02 用 PHANOTATE 重新推导）。 | `--include protein` |
| `cds.fna` | CDS nucleotides for codon frequency analysis. / CDS 核苷酸，用于密码子频率分析。 | `--include cds` |

---

## Missing Phages (3) / 缺失噬菌体（3 个）

These accessions appear in `phage_list.csv` but have no directory on disk. Re-download with:

```bash
python 00_raw_data/processes/fetch_phages.py --accession NC_013971.1
python 00_raw_data/processes/fetch_phages.py --accession NZ_CP007800.1
python 00_raw_data/processes/fetch_phages.py --accession NZ_CP008698.1
```

以下 accession 在 `phage_list.csv` 中有记录但磁盘上无对应目录，需重新下载：

| Accession | Source tag | Action |
|-----------|------------|--------|
| `NC_013971.1` | TBD | Re-download |
| `NZ_CP007800.1` | TBD | Re-download |
| `NZ_CP008698.1` | TBD | Re-download |

---

## Canonical Reference Accessions (Highest Priority) / 规范参考 Accession（最高优先级）

These three are required before any dry-lab module can run its reference scaffold validation:

| Role / 角色 | Organism | Accession | Status |
|------------|----------|-----------|--------|
| RBP scaffold phage | *Xanthomonas* phage phiL7 | `EU717894.1` | ✅ Present; 44,080 bp verified |
| Host receptor source | *X. campestris* pv. *campestris* ATCC 33913 | `GCF_000007145.1` / `NC_007086.1` | ✅ Present |
| ELISA positive-control phage | T7 phage | `GCF_000840885.1` | ✅ Present (fetched by Module 01) |

---

## Unresolved Bacteria / 未解析细菌

The following entries are in `bacteria_list.csv` but cannot yet be downloaded. Resolution is required by Module 01 / Sarah before these organisms can be used in ML training.

以下条目在 `bacteria_list.csv` 中，但尚无法下载，需由模块 01 / Sarah 解析后才能用于 ML 训练。

### Invalid accessions (not genome assemblies) / 无效 accession（非基因组 assembly）

| Accession | Problem | Recommended Fix |
|-----------|---------|-----------------|
| `KY000037` | Ti plasmid sequence, not a genome assembly / Ti 质粒序列，非基因组 assembly | Find proper *Agrobacterium sp.* genome accession and update interaction matrix |
| `PY746849` | Patent sequence (US 12571058 B2), not a genome assembly / 专利序列，非基因组 assembly | Find proper *Pseudomonas sp.* genome accession and update interaction matrix |

### Duplicate accessions (two matrix columns → same genome) / 重复 accession（两个矩阵列指向同一基因组）

| Accession | Organism A | Organism B | Recommended Fix |
|-----------|------------|------------|-----------------|
| `NZ_JBWJFR000000000` | *X. arboricola* | *X. arboricola* pv. *pruni* | Merge columns or find distinct genome for species-level entry |
| `NZ_CP150073` | *X. oryzae* | *X. oryzae* pv. *oryzae* | Same as above |

### Strain too vague or no deposited genome / 菌株描述过于笼统或无公开基因组

See original `data_needs.md` (pre-pivot version archived in `archive/`) for the full list of Category A–D strains. Resolution strategies are documented there.

详见 `archive/` 中归档的原始 `data_needs.md`（主动学习转型前版本）中的 A–D 类菌株完整列表。

---

## Do We Need More Data? / 还需要补充哪些数据？

Based on the active-learning pipeline's current needs:

| Priority | What | Why | Owner |
|----------|------|-----|-------|
| HIGH | Re-download 3 missing phages | Module 02 cannot annotate them | Alex / automated |
| HIGH | Fix KY000037 + PY746849 in interaction matrix | These bacteria have no genome | Sarah (Module 01) |
| MEDIUM | Add *X. hortorum* pv. *vitians* CFBP 8686PT (`GCA_012922135.1`) | Only unresolved high-value *Xanthomonas* host | Sarah |
| MEDIUM | Add explicit `0` labels for 443 negative-pool phages × Xanthomonas hosts | Required before ML training (Module 06) | Sarah / Alex |
| LOW | Resolve duplicate accessions (`NZ_JBWJFR000000000`, `NZ_CP150073`) | Produces identical feature vectors | Sarah (Module 01) |
| LOW | Resolve Category A–C unresolved strains | Expands interaction matrix coverage | Sarah |

---

## Negative-Pool Phages: Labeling Requirement / Negative-Pool 噬菌体的标注要求

The 443 `negative_pool` phages have no rows in the interaction matrix. Before Module 06 (uncertainty model) can train:

- All `(negative_pool phage × Xanthomonas host)` pairs must be added with label `0`.
- These are **not** missing data — they are structural negatives (no known cross-infection).
- Pairs involving non-*Xanthomonas* bacteria can also receive label `0` without wet-lab confirmation.

443 个 `negative_pool` 噬菌体在互作矩阵中无记录。模块 06（不确定性模型）训练前：

- 所有 `（negative_pool 噬菌体 × Xanthomonas 宿主）` 组合必须添加标签 `0`。
- 这些**不是**缺失数据——它们是结构性阴性（无已知交叉感染）。
- 涉及非 *Xanthomonas* 细菌的组合也可以无需湿实验确认地赋予标签 `0`。
