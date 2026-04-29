# -*- coding: utf-8 -*-
"""
2026 iGEM Dry Lab — Negative Data Fetching
==========================================
Purpose:
    Generate Affinity = 0 (negative) training samples for the phage–host
    infection prediction ML model.

Strategy:
    Two complementary approaches are implemented:

    [Module A] Cross-genus negatives
        Fetch phages whose annotated host is a bacterium from a DIFFERENT
        genus than Xanthomonas (e.g. Pseudomonas, Escherichia, Bacillus …).
        These phages are highly unlikely to infect Xanthomonas, so every
        (phage, Xanthomonas_strain) pair becomes a confident Affinity = 0
        row in the training set.

    [Module B] Within-Xanthomonas negatives from host-range tables
        Some Xanthomonas phage papers deposit GenBank records that contain
        structured notes / publications describing which strains were TESTED
        but NOT infected.  We mine those records for explicit "0" entries.
        This is rarer but gives high-quality within-genus negatives.

Output files:
    negative_cross_genus.csv          — Module A results (long format)
    negative_within_xanthomonas.csv   — Module B results (long format)
    negative_data_combined.csv        — Union of both, deduplicated, ready
                                        to merge with your positive matrix

Each output row has the same schema as phage_host_matrix_with_ids.csv:
    Phage, Phage_Accession, Host_Name, Host_Accession, Affinity, Source
"""

from Bio import Entrez, SeqIO
import pandas as pd
import time
import re

# ── 基本設定 ────────────────────────────────────────────────────────────────
Entrez.email = "cchen29@cmc.edu"   # ← 保持和原始 notebook 一致

# ── 快取（避免重複查 NCBI）──────────────────────────────────────────────────
host_acc_cache: dict[str, str] = {}


# ════════════════════════════════════════════════════════════════════════════
# 共用工具函數
# ════════════════════════════════════════════════════════════════════════════

def get_host_accession(host_name: str) -> str:
    """
    給定宿主細菌名稱，回傳其 RefSeq 完整基因組 Accession ID。
    與原始 notebook 中的同名函數完全相容，共享 host_acc_cache。
    """
    if not host_name or host_name in ("Unknown", "ERROR"):
        return "Unknown"
    if host_name in host_acc_cache:
        return host_acc_cache[host_name]

    for query in [
        f'"{host_name}"[Organism] AND "complete genome"[All Fields] AND srcdb_refseq[PROP]',
        f'"{host_name}"[Organism] AND "complete genome"[All Fields]',
    ]:
        try:
            sh = Entrez.esearch(db="nucleotide", term=query, retmax=1)
            sr = Entrez.read(sh); sh.close()
            if sr["IdList"]:
                sumh = Entrez.esummary(db="nucleotide", id=sr["IdList"][0])
                sumr = Entrez.read(sumh); sumh.close()
                acc = sumr[0]["Caption"]
                host_acc_cache[host_name] = acc
                return acc
        except Exception as e:
            print(f"  ⚠️  get_host_accession({host_name}): {e}")
        finally:
            time.sleep(0.4)

    host_acc_cache[host_name] = "No Complete Genome Found"
    return "No Complete Genome Found"


def fetch_gb_records_batch(accession_ids: list[str]) -> list:
    """批次下載 GenBank records，自動重試一次。"""
    query = ",".join(accession_ids)
    for attempt in range(2):
        try:
            h = Entrez.efetch(db="nucleotide", id=query,
                              rettype="gb", retmode="text")
            records = list(SeqIO.parse(h, "genbank"))
            h.close()
            return records
        except Exception as e:
            print(f"  ⚠️  efetch attempt {attempt+1} failed: {e}")
            time.sleep(3)
    return []


def extract_host_from_record(record) -> str:
    """從 GenBank record 的 source feature 抽取 host 欄位。"""
    for feat in record.features:
        if feat.type == "source":
            info = feat.qualifiers.get(
                "host", feat.qualifiers.get("lab_host", ["Unknown"])
            )
            return info[0] if info else "Unknown"
    return "Unknown"


# ════════════════════════════════════════════════════════════════════════════
# Module A：跨屬負樣本
# ════════════════════════════════════════════════════════════════════════════
#
# 抓取感染下列「非 Xanthomonas」屬細菌的完整噬菌體基因組。
# 這些噬菌體對所有 Xanthomonas 菌株的 Affinity 設為 0。
#
# 選菌邏輯：
#   - 優先選 Xanthomonadaceae 同科的近緣屬（Stenotrophomonas、Lysobacter）
#     → 系統發育較近，但仍有宿主特異性，是高品質的困難負樣本
#   - 再選常用 phage model host（Pseudomonas、Escherichia、Bacillus）
#     → 系統發育較遠，提供容易的負樣本，幫助模型學習基本特徵
# ════════════════════════════════════════════════════════════════════════════

# 每個屬最多抓幾筆（平衡各屬樣本數，避免 Escherichia phage 壓倒性多數）
NON_XANTHOMONAS_GENERA = {
    # --- 近緣屬（困難負樣本） ---
    "Stenotrophomonas": 80,   # 同屬 Xanthomonadaceae，高品質困難負樣本
    "Lysobacter":       40,
    "Xylella":          40,   # 同科但不同屬，植物病原，生態位相似
    # --- 常見模式細菌（容易負樣本） ---
    "Pseudomonas":      80,
    "Escherichia":      80,
    "Bacillus":         60,
    "Klebsiella":       40,
    "Salmonella":       40,
    "Staphylococcus":   40,
    "Mycobacterium":    30,
}

# 你在 Module A 最終想要的 Xanthomonas 代表菌株
# （這些是 cross-genus phage 的 negative host column）
XANTHOMONAS_REPRESENTATIVE_STRAINS = [
    "Xanthomonas oryzae pv. oryzae",
    "Xanthomonas campestris pv. campestris",
    "Xanthomonas axonopodis",
    "Xanthomonas citri",
]


def fetch_non_xanthomonas_phage_accessions(genus: str, retmax: int) -> list[str]:
    """
    搜尋某個屬細菌的噬菌體完整基因組，回傳 Accession 清單。
    """
    query = (
        f'("{genus}"[All Fields] AND '
        f'("phage"[All Fields] OR "bacteriophage"[All Fields]) AND '
        f'"complete genome"[All Fields] AND "Viruses"[Organism])'
    )
    try:
        sh = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)
        sr = Entrez.read(sh); sh.close()
        id_list = sr["IdList"]
        time.sleep(0.5)

        if not id_list:
            return []

        # 轉換 UID → Caption (Accession)
        accessions = []
        for start in range(0, len(id_list), 200):
            batch = id_list[start:start+200]
            sumh = Entrez.esummary(db="nucleotide", id=",".join(batch))
            sumr = Entrez.read(sumh); sumh.close()
            accessions += [item["Caption"] for item in sumr
                           if int(item.get("Length", 0)) > 10000]
            time.sleep(0.8)

        return accessions

    except Exception as e:
        print(f"  ⚠️  搜尋 {genus} phage 失敗: {e}")
        return []


def run_module_a(xanthomonas_strains: list[str]) -> pd.DataFrame:
    """
    執行 Module A：抓跨屬噬菌體，對所有 Xanthomonas 代表菌株設 Affinity = 0。

    Returns
    -------
    pd.DataFrame  columns: Phage, Phage_Accession, Host_Name,
                           Host_Accession, Affinity, Source
    """
    print("\n" + "═"*60)
    print("MODULE A：跨屬負樣本抓取")
    print("═"*60)

    all_phage_rows = []   # 先存每個 phage 自身的 host 資訊
    all_accessions = []   # 所有跨屬 phage accession，之後展開成 negative pairs

    for genus, retmax in NON_XANTHOMONAS_GENERA.items():
        print(f"\n▶  正在抓取 {genus} phage（最多 {retmax} 筆）...")
        acc_list = fetch_non_xanthomonas_phage_accessions(genus, retmax)
        print(f"   找到 {len(acc_list)} 筆符合條件的 Accession")

        batch_size = 50
        for i in range(0, len(acc_list), batch_size):
            batch = acc_list[i:i+batch_size]
            records = fetch_gb_records_batch(batch)

            for rec in records:
                host_name = extract_host_from_record(rec)
                # 只保留明確標注非 Xanthomonas 宿主的記錄
                if "Xanthomonas" in host_name:
                    continue  # 萬一混進來，跳過

                all_phage_rows.append({
                    "Phage":           rec.description,
                    "Phage_Accession": rec.id,
                    "True_Host":       host_name,   # 這個 phage 真正感染的宿主
                    "Genus":           genus,
                })

            time.sleep(1.2)

    print(f"\n共收集到 {len(all_phage_rows)} 隻跨屬噬菌體，正在展開 negative pairs...")

    # 預先查詢所有 Xanthomonas 代表菌株的 Accession
    xantho_acc_map = {}
    for strain in xanthomonas_strains:
        print(f"  🔍 查詢 Xanthomonas 代表菌株 Accession: {strain}")
        xantho_acc_map[strain] = get_host_accession(strain)

    # 展開：每隻跨屬 phage × 每個 Xanthomonas 代表菌株 = 一筆 Affinity=0 資料
    negative_rows = []
    for phage in all_phage_rows:
        for strain, strain_acc in xantho_acc_map.items():
            negative_rows.append({
                "Phage":          phage["Phage"],
                "Phage_Accession": phage["Phage_Accession"],
                "Host_Name":      strain,
                "Host_Accession": strain_acc,
                "Affinity":       0,
                "Source":         f"Module_A_cross_genus:{phage['Genus']}|NCBI:{phage['Phage_Accession']}",
            })

    df_a = pd.DataFrame(negative_rows)
    print(f"✅  Module A 完成，共產出 {len(df_a)} 筆 negative pairs")
    return df_a


# ════════════════════════════════════════════════════════════════════════════
# Module B：Xanthomonas 噬菌體的文獻已知負感染記錄
# ════════════════════════════════════════════════════════════════════════════
#
# 部分 Xanthomonas phage 的 GenBank record 在 comment / references 中
# 明確記載了哪些菌株被測試但「不感染」。
# 這裡的策略：
#   1. 從你已有的 xanthomonas_phages_accession_list.csv 讀取所有 phage
#   2. 全文下載 GenBank record
#   3. 用正規表達式在 comment 和 structured_comment 欄位搜尋
#      「resistant / no lysis / not infected / 0 / 不感染」等關鍵詞
#   4. 提取對應菌株名稱，輸出 Affinity = 0
#
# 注意：這個方法覆蓋率較低（大多數論文並不把 negative 結果上傳 NCBI），
#       但精確度高，可作為 Module A 的高品質補充。
# ════════════════════════════════════════════════════════════════════════════

# 正則：偵測「不感染」的上下文
NEGATIVE_KEYWORDS_RE = re.compile(
    r"(no\s+lysis|resistant|not\s+infected|unable\s+to\s+infect|"
    r"no\s+plaques?|insensitive|no\s+infection|did\s+not\s+infect|"
    r"could\s+not\s+infect|negative\s+result)",
    re.IGNORECASE
)

# 正則：從「不感染」語境中抽取 Xanthomonas 菌株名稱
XANTHOMONAS_STRAIN_RE = re.compile(
    r"(Xanthomonas\s+[a-z]+(?:\s+pv\.\s+\w+)?(?:\s+\w+)?)",
    re.IGNORECASE
)


def parse_negative_hosts_from_record(record) -> list[str]:
    """
    從一個 GenBank record 的自由文字欄位（comment、references）
    嘗試提取被明確標注「不感染」的 Xanthomonas 菌株名稱。
    """
    # 合并所有自由文字欄位
    text_sources = [record.description]

    # GenBank comment 欄位
    if "comment" in record.annotations:
        text_sources.append(record.annotations["comment"])

    # structured_comment (有些 WGS 會有)
    if "structured_comment" in record.annotations:
        sc = record.annotations["structured_comment"]
        text_sources.append(str(sc))

    # References abstract / title
    for ref in record.annotations.get("references", []):
        text_sources.append(getattr(ref, "title", ""))
        text_sources.append(getattr(ref, "comment", ""))

    full_text = " ".join(text_sources)

    # 先確認文字中有「不感染」關鍵詞
    if not NEGATIVE_KEYWORDS_RE.search(full_text):
        return []

    # 抽取被提及的 Xanthomonas 菌株
    candidates = XANTHOMONAS_STRAIN_RE.findall(full_text)

    # 去重並清理空白
    return list({s.strip() for s in candidates if s.strip()})


def run_module_b(positive_csv: str = "xanthomonas_phages_accession_list.csv") -> pd.DataFrame:
    """
    執行 Module B：從 Xanthomonas phage 的 GenBank records 中挖掘
    已知負感染記錄。

    Parameters
    ----------
    positive_csv : str
        你在 Step 1 生成的 xanthomonas_phages_accession_list.csv 路徑

    Returns
    -------
    pd.DataFrame  同 Module A 的 schema
    """
    print("\n" + "═"*60)
    print("MODULE B：Xanthomonas 噬菌體文獻負感染記錄挖掘")
    print("═"*60)

    try:
        df = pd.read_csv(positive_csv)
        acc_list = df["Accession"].tolist()
    except FileNotFoundError:
        print(f"  ⚠️  找不到 {positive_csv}，跳過 Module B。")
        print("     請先執行原始 notebook 的 Step 1 生成該檔案。")
        return pd.DataFrame()

    print(f"  載入 {len(acc_list)} 個 Xanthomonas phage accessions")

    negative_rows = []
    batch_size = 50

    for i in range(0, len(acc_list), batch_size):
        batch = acc_list[i:i+batch_size]
        print(f"  處理進度：{i+1} ~ {min(i+batch_size, len(acc_list))} / {len(acc_list)}")

        records = fetch_gb_records_batch(batch)

        for rec in records:
            neg_hosts = parse_negative_hosts_from_record(rec)

            for host_name in neg_hosts:
                host_acc = get_host_accession(host_name)
                negative_rows.append({
                    "Phage":           rec.description,
                    "Phage_Accession": rec.id,
                    "Host_Name":       host_name,
                    "Host_Accession":  host_acc,
                    "Affinity":        0,
                    "Source":          f"Module_B_literature|NCBI:{rec.id}",
                })

        time.sleep(1.2)

    df_b = pd.DataFrame(negative_rows)
    if df_b.empty:
        print("  ℹ️  Module B 未找到任何明確的負感染記錄（這是正常的，大多數 NCBI 記錄不包含此資訊）")
    else:
        print(f"✅  Module B 完成，共產出 {len(df_b)} 筆 negative pairs")
    return df_b


# ════════════════════════════════════════════════════════════════════════════
# Module C：Xanthomonas phage 感染同屬其他 pv. 的負樣本（基於已知陽性推斷）
# ════════════════════════════════════════════════════════════════════════════
#
# 補充思路：
#   某些 phage 已知感染 Xanthomonas oryzae pv. oryzae，
#   但文獻指出它們不感染 X. oryzae pv. oryzicola。
#   這類「同種不同 pathovar」的特異性資訊可以從你的 positive matrix 推斷：
#   已知 Affinity=1 的 (phage, hostA) pair，
#   若 hostA 和 hostB 是同種不同 pathovar，
#   且 (phage, hostB) 沒有任何記錄 → 可視為 Affinity=0。
#
#   這個 Module 不需要 NCBI，直接從你的 positive matrix 計算。
# ════════════════════════════════════════════════════════════════════════════

# 已知 Xanthomonas pv. 同種關係（你可以根據自己的知識擴展這個表）
XANTHOMONAS_PV_GROUPS = [
    # 同種 pathovar 組合，phage 感染其中一個但可能不感染另一個
    ["Xanthomonas oryzae pv. oryzae", "Xanthomonas oryzae pv. oryzicola"],
    ["Xanthomonas campestris pv. campestris", "Xanthomonas campestris pv. vesicatoria"],
    ["Xanthomonas axonopodis pv. citri", "Xanthomonas axonopodis pv. phaseoli"],
    ["Xanthomonas citri subsp. citri", "Xanthomonas citri subsp. malvacearum"],
]


def run_module_c(positive_matrix_csv: str = "phage_host_matrix_with_ids.csv") -> pd.DataFrame:
    """
    執行 Module C：從現有 positive matrix 推斷同屬不同 pathovar 的負樣本。

    Parameters
    ----------
    positive_matrix_csv : str
        你在原始 notebook Step 3 生成的 phage_host_matrix_with_ids.csv 路徑

    Returns
    -------
    pd.DataFrame  同 Module A 的 schema
    """
    print("\n" + "═"*60)
    print("MODULE C：Xanthomonas 同屬 pathovar 特異性負樣本推斷")
    print("═"*60)

    try:
        df = pd.read_csv(positive_matrix_csv)
    except FileNotFoundError:
        print(f"  ⚠️  找不到 {positive_matrix_csv}，跳過 Module C。")
        return pd.DataFrame()

    df_pos = df[df["Affinity"] == 1].copy()
    print(f"  讀取到 {len(df_pos)} 筆 Affinity=1 資料")

    negative_rows = []

    for group in XANTHOMONAS_PV_GROUPS:
        for i, host_a in enumerate(group):
            for host_b in group[i+1:]:
                # 找出感染 host_a 的 phage 且沒有 host_b 的任何記錄
                phages_infect_a = set(
                    df_pos[df_pos["Host_Name"].str.contains(host_a, na=False)]["Phage_Accession"]
                )
                phages_with_b_data = set(
                    df[df["Host_Name"].str.contains(host_b, na=False)]["Phage_Accession"]
                )

                # 感染 A、但對 B 完全無記錄的 phage → negative for B
                candidates = phages_infect_a - phages_with_b_data

                if candidates:
                    print(f"  {host_a} → {host_b}: 找到 {len(candidates)} 個候選負樣本 phage")

                for phage_acc in candidates:
                    phage_row = df_pos[df_pos["Phage_Accession"] == phage_acc].iloc[0]
                    host_b_acc = get_host_accession(host_b)
                    negative_rows.append({
                        "Phage":           phage_row["Phage"],
                        "Phage_Accession": phage_acc,
                        "Host_Name":       host_b,
                        "Host_Accession":  host_b_acc,
                        "Affinity":        0,
                        "Source":          f"Module_C_pv_inference|positive_host:{host_a}",
                    })

                # 反向：感染 B、但對 A 完全無記錄
                phages_infect_b = set(
                    df_pos[df_pos["Host_Name"].str.contains(host_b, na=False)]["Phage_Accession"]
                )
                phages_with_a_data = set(
                    df[df["Host_Name"].str.contains(host_a, na=False)]["Phage_Accession"]
                )
                candidates_rev = phages_infect_b - phages_with_a_data

                if candidates_rev:
                    print(f"  {host_b} → {host_a}: 找到 {len(candidates_rev)} 個候選負樣本 phage")

                for phage_acc in candidates_rev:
                    phage_row = df_pos[df_pos["Phage_Accession"] == phage_acc].iloc[0]
                    host_a_acc = get_host_accession(host_a)
                    negative_rows.append({
                        "Phage":           phage_row["Phage"],
                        "Phage_Accession": phage_acc,
                        "Host_Name":       host_a,
                        "Host_Accession":  host_a_acc,
                        "Affinity":        0,
                        "Source":          f"Module_C_pv_inference|positive_host:{host_b}",
                    })

    df_c = pd.DataFrame(negative_rows)
    if df_c.empty:
        print("  ℹ️  Module C 未產出任何負樣本（可能是 positive matrix 中 pathovar 欄位命名不一致）")
    else:
        print(f"✅  Module C 完成，共產出 {len(df_c)} 筆 negative pairs")
    return df_c


# ════════════════════════════════════════════════════════════════════════════
# 主流程：執行所有 Module 並合併輸出
# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":

    # ── Module A ──────────────────────────────────────────────────────────
    df_a = run_module_a(XANTHOMONAS_REPRESENTATIVE_STRAINS)
    df_a.to_csv("negative_cross_genus.csv", index=False, encoding="utf-8-sig")
    print(f"\n💾  Module A 輸出已儲存：negative_cross_genus.csv  ({len(df_a)} 行)")

    # ── Module B ──────────────────────────────────────────────────────────
    df_b = run_module_b("xanthomonas_phages_accession_list.csv")
    if not df_b.empty:
        df_b.to_csv("negative_within_xanthomonas.csv", index=False, encoding="utf-8-sig")
        print(f"💾  Module B 輸出已儲存：negative_within_xanthomonas.csv  ({len(df_b)} 行)")

    # ── Module C ──────────────────────────────────────────────────────────
    df_c = run_module_c("phage_host_matrix_with_ids.csv")
    if not df_c.empty:
        df_c.to_csv("negative_pv_inference.csv", index=False, encoding="utf-8-sig")
        print(f"💾  Module C 輸出已儲存：negative_pv_inference.csv  ({len(df_c)} 行)")

    # ── 合併三個 Module 的結果 ────────────────────────────────────────────
    frames = [df for df in [df_a, df_b, df_c] if not df.empty]
    if frames:
        df_combined = pd.concat(frames, ignore_index=True)

        # 去重：同一個 (Phage_Accession, Host_Name) pair 只保留一筆
        df_combined = df_combined.drop_duplicates(
            subset=["Phage_Accession", "Host_Name"], keep="first"
        )

        # 最終確保 Affinity 都是數字 0
        df_combined["Affinity"] = 0

        df_combined.to_csv("negative_data_combined.csv", index=False, encoding="utf-8-sig")
        print(f"\n🎉  合併完成！negative_data_combined.csv  共 {len(df_combined)} 筆 Affinity=0 資料")
        print(f"    Schema：{list(df_combined.columns)}")
        print("\n=== 預覽 ===")
        print(df_combined.head(10).to_string(index=False))
    else:
        print("\n⚠️  所有 Module 均未產出資料，請檢查網路連線或輸入檔案。")