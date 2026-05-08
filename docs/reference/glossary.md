# Glossary — iGEM Claremont 2026

A reference vocab list covering biology, molecular biology, lab equipment, computational tools, machine learning, regulatory, and project-specific terminology used across this project's documentation.

> **How to use:** Sections are organized by domain. Each entry gives the expansion (for acronyms), a plain-English definition, and *italicized context* on how it shows up in our project.

---

## Quick A–Z index of most-used acronyms

| Acronym | Stands for | Section |
|---|---|---|
| AL | Active Learning | ML |
| AF3 | AlphaFold 3 | Bioinformatics |
| ALDE | Active Learning-assisted Directed Evolution | ML |
| APHIS | Animal and Plant Health Inspection Service | Regulatory |
| ATCC | American Type Culture Collection | Databases |
| AUPR | Area Under Precision-Recall curve | ML |
| BALD | Bayesian Active Learning by Disagreement | ML |
| BLI | Biolayer Interferometry | Binding assays |
| BOED | Bayesian Optimal Experimental Design | ML |
| BSL | Biosafety Level | Regulatory |
| CDFA | California Department of Food and Agriculture | Regulatory |
| CFU | Colony-Forming Unit | Microbiology |
| DBTL | Design–Build–Test–Learn | Methodology |
| DMS | Deep Mutational Scanning | Methodology |
| DSMZ | Deutsche Sammlung von Mikroorganismen | Databases |
| ELISA | Enzyme-Linked Immunosorbent Assay | Binding assays |
| EIG | Expected Information Gain | ML |
| EPPO | European and Mediterranean Plant Protection Org. | Regulatory |
| ESM-2 | Evolutionary Scale Modeling v2 | Bioinformatics |
| FPLC | Fast Protein Liquid Chromatography | Equipment |
| GP | Gaussian Process | ML |
| HMC | Harvey Mudd College | Institutions |
| HMM | Hidden Markov Model | Bioinformatics |
| HPC | High-Performance Computing | Equipment |
| HRP | Horseradish Peroxidase | Binding assays |
| IMAC | Immobilized Metal Affinity Chromatography | Purification |
| IPTG | Isopropyl β-D-1-thiogalactopyranoside | Cloning |
| LPS | Lipopolysaccharide | Microbiology |
| MALDI-TOF | Matrix-Assisted Laser Desorption Ionization, Time-of-Flight | Equipment |
| MLP | Multilayer Perceptron | ML |
| MLST | Multilocus Sequence Typing | Bioinformatics |
| MOI | Multiplicity of Infection | Microbiology |
| NCBI | National Center for Biotechnology Information | Databases |
| NCPPB | National Collection of Plant Pathogenic Bacteria | Databases |
| OMP | Outer Membrane Protein | Microbiology |
| ORF | Open Reading Frame | Bioinformatics |
| PBS | Phosphate-Buffered Saline | Materials |
| PDB | Protein Data Bank | Databases |
| PFU | Plaque-Forming Unit | Microbiology |
| PHANOTATE | Phage gene-calling tool | Bioinformatics |
| PLM | Protein Language Model | ML |
| PPI | Protein-Protein Interaction | ML / Biology |
| PPQ-526 | Permit to Move Live Plant Pests (USDA form) | Regulatory |
| RBP | Receptor-Binding Protein | Microbiology |
| SDS-PAGE | Sodium Dodecyl Sulfate–Polyacrylamide Gel Electrophoresis | Purification |
| SLA | Service-Level Agreement | Project mgmt |
| SPR | Surface Plasmon Resonance | Binding assays |
| TEM | Transmission Electron Microscopy | Equipment |
| UQ | Uncertainty Quantification | ML |
| USDA | U.S. Department of Agriculture | Regulatory |
| WGS | Whole-Genome Sequencing | Bioinformatics |
| Xcc | *Xanthomonas campestris* pv. *campestris* | Microbiology |
| YPGA | Yeast Peptone Glycerol Agar | Materials |

---

## 1. Bacteriology & phage biology

- **Xanthomonas** — Genus of Gram-negative plant-pathogenic bacteria; many species cause major crop diseases. *Our wet-lab target organism.*
- **Pathovar (pv.)** — A formal subdivision of a bacterial species defined by host range or symptoms (e.g. *X. campestris* pv. *campestris* infects brassicas; pv. *pelargonii* infects pelargoniums). *Different pathovars often need different phages.*
- **Xcc** — Shorthand for *Xanthomonas campestris* pv. *campestris*; the canonical model strain (ATCC 33913) causes black rot of brassicas. *Our reference host strain.*
- **Strain** vs **isolate** — A *strain* is a cultivated lineage with a fixed identifier; an *isolate* is what you just pulled from the environment before identification. *Once we identify our isolate by 16S/MALDI-TOF, it becomes a named strain.*
- **Type strain** — The reference exemplar of a species, deposited in multiple culture collections; everyone refers to it. *Xcc ATCC 33913 is the type strain for X. campestris.*
- **BSL-1 / BSL-2** — Biosafety Level: BSL-1 is for organisms with no known human disease risk (basic precautions); BSL-2 needs containment cabinets etc. *We work BSL-1 only; all Xanthomonas plant pathogens are BSL-1 in California academic labs.*
- **Bacteriophage** (or **phage**) — A virus that infects bacteria. *Our project uses lytic phages as biocontrol agents against Xanthomonas.*
- **Lytic** vs **temperate** — Lytic phages always kill the host on infection; temperate phages can integrate as a prophage and lie dormant (lysogeny). *We need lytic phages — temperate ones don't reliably kill.*
- **Filamentous phage (Inoviridae)** — Long thread-like phages (e.g., M13, CP1) that extrude through the host without killing it. *We exclude these — wrong life cycle for biocontrol.*
- **Podoviridae / Siphoviridae / Myoviridae** — Three classical families of tailed phages by tail morphology: Podo = short tail, Sipho = long non-contractile tail, Myo = long contractile tail. *phiL7 is a Siphoviridae T7-like phage.*
- **Capsid** — The protein shell that packages the phage's genome.
- **Tail fiber / tail spike** — The protein appendages that recognize and bind the bacterial surface receptor. *These are the RBPs — the central focus of our design work.*
- **RBP** — Receptor-Binding Protein. The phage protein that physically binds the bacterial receptor; it determines host range. *The molecule we engineer.*
- **Receptor** — A bacterial surface molecule (protein or sugar) that the phage RBP grabs to initiate infection. *For phiL7, the receptor is the TonB-ExbB-ExbD1D2 system.*
- **TonB-dependent receptor** — A class of outer-membrane bacterial proteins (typical β-barrel structure) that import nutrients and double as phage receptors; they require energy from the inner-membrane TonB complex. *phiL7's identified receptor system.*
- **LPS** — Lipopolysaccharide, the major glycan component of the Gram-negative outer membrane; many phages bind LPS sugars. *A possible alternative receptor type for our isolated phages.*
- **OMP** — Outer Membrane Protein; β-barrel proteins embedded in the bacterial outer membrane. *Protein-class receptors are easier to model than glycan ones.*
- **Glycan** — A carbohydrate chain; phage RBPs binding glycans are one of the major modeling gaps in the field. *One of the open problems our work touches on.*
- **Plaque** — A clear zone of bacterial lysis on a lawn of cells, formed by one phage seeding successive infection rounds. *Plaque assay is our basic infectivity readout.*
- **Plaque assay** — Quantitative method: serially dilute phage, plate on bacterial lawn, count plaques. *Sarah's existing Benchling protocol.*
- **MOI** — Multiplicity of Infection; phage particles per bacterial cell at the start of an infection. *Standard 0.01–0.1 for our liquid amplification.*
- **PFU / mL** — Plaque-Forming Units per mL; the standard unit for phage titer.
- **CFU / mL** — Colony-Forming Units per mL; the standard unit for live bacterial count.
- **Lysate** — The cell-free supernatant after phages have lysed a culture; contains high-titer phage particles. *Stock material for downstream assays.*
- **Lysis** — The bursting of bacterial cells caused by lytic phage replication.
- **Adsorption** — The first physical attachment of phage to host cell surface; necessary but not sufficient for productive infection. *A confounder — see Farquharson 2021.*
- **Productive infection** — Infection that successfully produces progeny phage and lyses the cell, vs binding-only or aborted infection. *What we ultimately care about, beyond binding.*
- **Defense systems** (CRISPR-Cas, R-M, Abi, etc.) — Bacterial systems that abort phage infection after entry; can decouple binding from lysis.
- **R-M** — Restriction-Modification; bacterial defense system that cuts foreign DNA (incl. phage DNA) at unmethylated sites.
- **Abi** — Abortive infection; defense systems that kill the infected cell to prevent phage spread.
- **Knockout** vs **knockdown** — Knockout = gene completely removed/disrupted; knockdown = gene expression reduced (e.g., CRISPRi). *We're doing knockout (markerless deletion) for receptor validation.*
- **Markerless deletion** — Gene deletion that leaves no antibiotic resistance scar; uses two-step selection (e.g., sacB counter-selection on sucrose).
- **pK18mobsacB** — A suicide-vector plasmid (Schäfer 1994) for markerless gene deletion in Gram-negative bacteria; carries kanamycin resistance + sacB (sucrose lethality) for double selection. *Addgene #87097, our planned receptor-knockout vector.*
- **sacB** — Bacillus subtilis levansucrase gene; kills Gram-negatives in the presence of sucrose. *Counter-selection marker for markerless deletion.*
- **Electroporation** — Transferring DNA into bacterial cells by brief electric pulses opening the membrane. *14 kV/cm in Sarah's Benchling Xcc protocol.*
- **Transformation** — General term for introducing exogenous DNA into bacteria; electroporation is one method.
- **Glycerol stock** — Bacterial culture mixed with glycerol (final ~20%) and stored at −80°C for long-term preservation.

## 2. Molecular biology, cloning, expression

- **Gibson assembly** — One-pot isothermal cloning method using exonuclease + polymerase + ligase to seamlessly join overlapping DNA fragments (Gibson 2009). *Our standard cloning method for RBP variant libraries.*
- **pET vector** — Expression vector family driven by T7 promoter; pET-28a is the common N-terminal His6-tag version. *Our RBP-His6 expression backbone.*
- **BL21(DE3)** — *E. coli* expression strain carrying T7 RNA polymerase under IPTG control; the workhorse for pET-based expression. *Our RBP production host.*
- **IPTG** — Isopropyl β-D-1-thiogalactopyranoside; non-metabolizable lactose analog that induces T7 polymerase expression in BL21(DE3). *Standard ~0.5 mM, 18°C overnight for soluble trimer expression.*
- **His6 / His-tag** — Six-histidine tag attached to a protein for nickel-affinity purification. *N-terminal on our RBP constructs.*
- **N-terminal / C-terminal** — The two ends of a protein chain; tail-fiber RBPs typically have the N-term anchor and C-term host-specificity head.
- **Trimerization** — Many phage RBPs assemble as homotrimers; correct trimer formation is required for binding activity. *Risk: misfolded monomers don't bind.*
- **GCN4 leucine zipper** — A short helical peptide that forces trimerization; can be fused to RBPs as a backup trimerization aid. *Our backup tag if natural trimerization fails.*
- **Inclusion body** — Misfolded protein aggregates inside *E. coli*; insoluble, biologically inactive. *Common failure mode for RBPs over-expressed at 37°C.*
- **Codon optimization** — Recoding a gene's DNA sequence to use codons preferred by the expression host (e.g., *E. coli*).
- **Gene synthesis** — Commercial DNA synthesis (IDT, Twist) of designer sequences; ~1–2 weeks turnaround.
- **ORF** — Open Reading Frame; a stretch of DNA encoding a continuous protein.
- **HMM** — Hidden Markov Model; a probabilistic model used for protein domain detection (e.g., Pfam HMMs).
- **Pfam** — A database of protein-family HMMs; the standard way to assign function to unknown proteins.
- **MLST** — Multilocus Sequence Typing; species/strain identification by sequencing 6–8 conserved housekeeping genes (e.g., *gyrB*, *rpoD*).
- **gyrB / rpoD** — Two of the standard Xanthomonas MLST loci; used to confirm pathovar identity.
- **16S rRNA gene** — The universal bacterial gene used for genus/species identification by Sanger sequencing.
- **WGS** — Whole-Genome Sequencing; we plan Illumina MiSeq for our isolates and phages.
- **Illumina MiSeq** — A small benchtop short-read sequencer; the standard for bacterial WGS in academic cores.

## 3. Protein purification & quality control

- **Ni-NTA** — Nickel-Nitrilotriacetic Acid resin; binds His-tagged proteins via Ni²⁺ chelation. *Our first-step purification matrix.*
- **IMAC** — Immobilized Metal Affinity Chromatography; the umbrella term for Ni-NTA-style purifications.
- **Imidazole** — Small molecule that competes with His-tag for nickel; used as gradient elution to release bound protein.
- **FPLC** — Fast Protein Liquid Chromatography; benchtop chromatography system used for column-based purifications.
- **SDS-PAGE** — Sodium Dodecyl Sulfate Polyacrylamide Gel Electrophoresis; separates proteins by size for purity / size verification.
- **Western blot** — Transfer SDS-PAGE proteins to a membrane and detect with antibody; we use anti-His6 to confirm RBP identity.
- **Bradford assay** — Colorimetric protein concentration assay using Coomassie dye binding.
- **Size-exclusion chromatography (SEC)** — Separates proteins by hydrodynamic radius; used to confirm trimer formation by elution volume.

## 4. Binding assays & their abbreviations

- **ELISA** — Enzyme-Linked Immunosorbent Assay; the workhorse plate-based protein-binding assay. *Our primary readout for the AL cycle: serial dilutions of His6-RBP on Xanthomonas-coated plates → HRP-anti-His6 detection → OD450.*
- **HRP** — Horseradish Peroxidase; an enzyme that produces a colored signal from substrates like TMB.
- **TMB** — 3,3',5,5'-Tetramethylbenzidine; the substrate that turns blue (then yellow when stopped) under HRP catalysis.
- **Anti-His6 antibody** — Antibody that recognizes the six-His tag; HRP-conjugated form gives the ELISA signal.
- **OD450** — Optical density at 450 nm; the wavelength read after TMB development with sulfuric-acid stop.
- **4-parameter logistic (4PL) fit** — The standard sigmoidal curve fit for ELISA dose-response data; parameters give EC50, Hill slope, baseline, asymptote.
- **EC50** — The concentration at which a binding signal reaches 50% of maximum; functional surrogate for binding affinity in ELISA.
- **Kd** — Dissociation constant; the equilibrium binding affinity (lower = tighter binding). *We report apparent Kd from ELISA, since true Kd needs SPR/BLI.*
- **Apparent Kd** — Kd measured under conditions that may not be equilibrium (e.g., ELISA wash steps); approximation but reproducible.
- **BSA blocking** — Coating ELISA plates with bovine serum albumin to prevent non-specific protein binding to the plastic.
- **Dot blot** — Pipetting protein samples directly onto a membrane (no gel); cheaper than ELISA, less quantitative.
- **Immunoblot** — Generic term for antibody-based detection on membranes (Western blot, dot blot are special cases).
- **SPR** — Surface Plasmon Resonance; gold-standard label-free binding assay measuring real-time binding kinetics on a sensor chip (instrument: Biacore). *Capital-intensive; we may validate one ELISA point against SPR at HMC core.*
- **BLI** — Biolayer Interferometry; alternative label-free binding tech (instrument: Octet); cheaper and easier than SPR. *Backup option if ELISA dynamic range is insufficient.*

## 5. Computational & bioinformatics tools

- **PHANOTATE** — Phage-specific ORF caller (McNair 2019); handles overlapping genes that defeat general-purpose callers.
- **Prodigal** — Gold-standard bacterial ORF caller (Hyatt 2010); we use it for Xanthomonas genomes.
- **pharokka** — Integrated phage annotation pipeline (Bouras 2023); combines PHANOTATE + functional databases.
- **PhageRBPdetect** — RBP identification tool (Boeckaerts 2022) using Pfam HMMs + custom HMMs + XGBoost; AUPR 93.8%.
- **AlphaFold 3 (AF3)** — Latest DeepMind model (Abramson 2024) predicting protein and protein-complex structures; we use it for both phage RBPs and host receptors.
- **Boltz-2** — 2025 protein-complex predictor with explicit affinity head; we use it as zero-shot binding-affinity prior to seed our model.
- **ESM-2** — Evolutionary Scale Modeling v2 (Lin 2023); a transformer-based protein language model trained on 65M sequences. *Our embedding backbone — 650M-parameter for routine, 3B for benchmarking.*
- **PLM-interact** — ESM-2 fine-tuned on human PPI data (Liu 2025); we may use as a transfer-learning prior to phage-bacteria binding.
- **BLAST** — Basic Local Alignment Search Tool; classic sequence-similarity search against NCBI databases.
- **HMMER** — HMM-search engine used with Pfam domain libraries.
- **MALDI-TOF** — Matrix-Assisted Laser Desorption Ionization, Time-of-Flight mass spectrometry; used to identify bacterial isolates by protein-mass fingerprint within minutes. *HMC core service.*
- **TEM** — Transmission Electron Microscopy; used optionally to image phage particle morphology (capsid + tail).

## 6. Machine learning & statistics

- **Active Learning (AL)** — ML paradigm where the model selects which data points to acquire next, maximizing information gain per experiment. *Our project's core methodology.*
- **BOED** — Bayesian Optimal Experimental Design (Lindley 1956); the formal mathematical framework underlying active learning.
- **BALD** — Bayesian Active Learning by Disagreement (Houlsby 2011); an acquisition function that selects the variant where ensemble members maximally disagree. *Our primary acquisition function.*
- **Acquisition function** — In AL, the function that scores candidate experiments; the next experiment is the one with the highest score (most information).
- **Deep ensemble** — Ensemble of N independently-trained neural networks; variance across members serves as uncertainty estimate (Lakshminarayanan 2017). *Our regression model.*
- **MLP** — Multilayer Perceptron; a feed-forward neural network. *Each ensemble member is a small MLP over ESM-2 embeddings.*
- **Gaussian Process (GP)** — A probabilistic regression model with built-in uncertainty estimates; works well for small datasets. *Alternative to deep ensembles; will compare in later cycles.*
- **MC Dropout** — Monte Carlo Dropout (Gal & Ghahramani 2016); approximates Bayesian uncertainty by running stochastic dropout at inference time.
- **UQ** — Uncertainty Quantification; producing well-calibrated estimates of prediction confidence (a key bottleneck in our pipeline).
- **Calibration** — How well a model's predicted uncertainty matches its actual error rate; a well-calibrated 80% confidence interval contains the truth ~80% of the time.
- **Posterior** — In Bayesian terms, the updated belief distribution after observing data.
- **Prior** — The starting belief distribution before observing data; informative priors leverage existing knowledge (e.g., Boltz-2 affinity predictions).
- **KL divergence** — Kullback–Leibler divergence; a measure of how different two probability distributions are. *Used to quantify information gain per cycle.*
- **EIG** — Expected Information Gain; the expected reduction in posterior entropy from a candidate experiment.
- **Embedding** — A learned vector representation of a sequence (or other input); ESM-2 produces 1280-D (650M) or 2560-D (3B) embeddings per residue.
- **Mean-pooled** — Averaging per-residue embeddings to get one fixed-length vector per protein.
- **Transfer learning** — Reusing a model trained on one task as the starting point for another; saves data on the target task.
- **Fine-tuning** — Training a pre-trained model further on a smaller, task-specific dataset.
- **Pre-training** — The (usually large) initial training of a model on broad data, before any specific task.
- **Zero-shot** — Predicting on a task without any task-specific training data; e.g., Boltz-2 affinity scores on our designed variants.
- **ALDE** — Active Learning-assisted Directed Evolution (Yang 2025); the closest published benchmark of AL methods for protein engineering.
- **AUPR** — Area Under the Precision-Recall Curve; classification metric particularly useful when positives are rare.
- **AUC / AUROC** — Area Under the Receiver Operating Characteristic curve; standard binary classification performance metric.
- **F1-score** — Harmonic mean of precision and recall; a single number summary of classifier performance.
- **R²** — Coefficient of determination; the fraction of variance in the target explained by the model.
- **Cross-validation** — Splitting data into K folds, training on K−1 and testing on the 1 held-out, rotating; gives more robust performance estimates than a single train/test split.
- **Held-out / test set** — Data the model has never seen during training, used to estimate generalization performance.
- **Retrospective replay** — After all cycles complete, simulating what would have happened if a different acquisition strategy (e.g., random) had been used; how we'll quantify AL advantage.
- **DMS** — Deep Mutational Scanning; high-throughput experimental measurement of fitness for thousands of single-mutant variants. *Our project is a low-throughput cousin — small batches, smarter selection.*
- **PPI** — Protein-Protein Interaction; the broader category that phage RBP – receptor binding belongs to.
- **PLM** — Protein Language Model; transformer architectures trained on protein sequences (ESM-2, ProtBERT, etc.).
- **XGBoost** — Extreme Gradient Boosting; tree-based ML algorithm widely used for tabular data; PhageRBPdetect's classifier.

## 7. Genomes, databases & accessions

- **NCBI** — National Center for Biotechnology Information; hosts GenBank, PubMed, BLAST, and most public bioinformatic resources.
- **GenBank accession** — Unique ID assigned to each sequence record at NCBI (e.g., AE008922 for Xcc, EU717894 for phiL7).
- **ATCC** — American Type Culture Collection; the largest U.S. bacterial culture collection. *We can't ship from ATCC due to USDA permits (see Regulatory).*
- **DSMZ** — Deutsche Sammlung von Mikroorganismen und Zellkulturen; the German equivalent of ATCC.
- **NCPPB** — National Collection of Plant Pathogenic Bacteria (UK); plant-pathology specialist culture collection.
- **BCCM/LMG** — Belgian Coordinated Collections of Microorganisms (Laboratory of Microbiology Ghent collection).
- **Phage Directory** — A community-driven phage exchange platform (phage.directory) connecting researchers with phage hunters.
- **PDB** — Protein Data Bank (rcsb.org); the central repository of solved 3D protein structures. *No PDB structure exists for phiL7 RBP — hence our reliance on AF3.*
- **AE008922** — Xcc ATCC 33913 complete genome accession.
- **EU717894** — phiL7 (Xanthomonas campestris phage) complete genome accession.

## 8. Lab facilities & equipment

- **HMC core** — Harvey Mudd College's shared instrumentation core facility; provides Sanger sequencing, MALDI-TOF, and (potentially) BLI/SPR access for the team.
- **Laguna** — Cluster supercomputer accessible to Claremont Colleges researchers; we use it for ESM-2 3B embedding, AlphaFold 3 batch jobs, and Boltz-2 inference.
- **HPC** — High-Performance Computing; the umbrella term for shared compute clusters like Laguna.
- **GPU** — Graphics Processing Unit; required for transformer-model inference (ESM-2, AF3); we need ≥12 GB VRAM for ESM-2 650M, ≥24 GB for 3B.
- **Spectrophotometer** — Instrument that reads OD600 (bacterial growth) or OD450 (ELISA color); standard benchtop kit.
- **Electroporator** — Pulses bacterial cells with a high-voltage spike to permit DNA uptake; we use Eppendorf Multiporator-class.
- **Ni-NTA column** — Disposable cartridge or refillable column packed with Ni-NTA resin for IMAC purification.

## 9. Regulatory, permits & institutional

- **USDA APHIS** — U.S. Department of Agriculture's Animal and Plant Health Inspection Service; regulates plant pest movement.
- **PPQ Form 526** — Permit to Move Live Plant Pests; required by APHIS before any plant-pathogen culture (incl. Xcc) can be shipped between facilities. *Average processing time 127 days — this is why we self-isolate instead.*
- **Plant pest** — APHIS regulatory term covering all plant-pathogenic microbes; Xanthomonas plant pathogens fall here.
- **Quarantine pest** — A subset of plant pests under stricter movement controls in specific regions; *X. citri* (citrus canker) is a CDFA quarantine concern in California.
- **Select Agent** — The strictest tier of USDA/CDC-regulated organisms (bioterrorism / high-consequence pathogens); we work with NONE of these.
- **CDFA** — California Department of Food and Agriculture; state-level plant-health regulator.
- **EPPO** — European and Mediterranean Plant Protection Organization; publishes standardized diagnostic protocols (e.g., PM 7/110 for Xanthomonas isolation) used worldwide.

## 10. Project methodology & management

- **DBTL cycle** — Design–Build–Test–Learn; the standard iGEM engineering cycle our active-learning loop instantiates.
- **Closed-loop** — A system where outputs feed back as inputs, enabling iterative refinement; our AL pipeline is closed-loop between dry and wet lab.
- **Cycle 0 / Cycle 1 / Cycle 2** — Our internal labeling for the first three rounds of the AL loop. Cycle 0 = cold start with seed variants; Cycles 1–2 = closed-loop with acquisition-function recommendations.
- **SLA** — Service-Level Agreement; informal commitment that dry lab will turn around new ELISA data into next-cycle recommendations within 48 hours.
- **Control arm** — A parallel set of experiments performed under a baseline condition (e.g., random selection) to compare against the experimental condition (AL selection).
- **Positive / negative control** — Reference samples with known outcomes used to validate the assay (positive = expected to bind; negative = expected NOT to bind).
- **Replicate** — Repeat measurement: *technical replicate* = same sample re-measured; *biological replicate* = independent biological sample.
- **Blind test** — A held-out set of samples whose true values are kept hidden from the predictor until after predictions are made.
- **Sanity check** — A simple verification that the system behaves as expected before trusting more complex results (e.g., positive-control binding before testing variants).
- **Epistasis** — Non-additive interaction between two mutations or conditions: combined effect ≠ sum of individual effects. *We test for epistasis between motif variation and receptor knockouts.*

## 11. Materials, media & reagents

- **NYG medium** — Nutrient–Yeast–Glycerol medium; standard rich growth medium for Xanthomonas (used in all Sarah's existing protocols).
- **YPGA + cycloheximide** — Yeast Peptone Glycerol Agar + cycloheximide (250 mg/L) to suppress fungal contamination; standard semi-selective medium for Xanthomonas isolation from plant tissue.
- **Cycloheximide** — Antifungal antibiotic used at 250 mg/L to suppress fungal overgrowth during bacterial isolation.
- **LB medium** — Luria-Bertani broth/agar; the standard *E. coli* growth medium (used for BL21 expression).
- **SOC medium** — Super Optimal broth with Catabolite repression; rich recovery medium used immediately after electroporation to maximize transformation yield.
- **PBS** — Phosphate-Buffered Saline; the standard physiological-pH wash buffer.
- **Kanamycin (kan)** — Aminoglycoside antibiotic; selection marker on pET-28a, pK18mobsacB, etc. (typically 20–50 µg/mL for plasmid maintenance).
- **Pluronic F-127** — A thermo-responsive triblock copolymer; forms a gel above ~20°C that can encapsulate phages for stability under field conditions.
- **Hydrogel** — A water-swollen polymer network; here used as a phage-protective delivery matrix for biocontrol applications.
- **Glycerol stock** — 20% glycerol + bacterial culture mix for −80°C storage; the standard long-term archive format.

## 12. iGEM-specific terminology

- **iGEM** — International Genetically Engineered Machine; the synthetic biology competition we're entering.
- **Wiki** — Each iGEM team's documentation website; the primary deliverable judges score.
- **Composite Part** — An iGEM Registry submission combining multiple BioBrick parts into a functional unit (e.g., our planned RBP-His6 expression library).
- **Best Agriculture Project** — One of iGEM's village awards; a key target for our project's framing.
- **Best Model** — Special prize for outstanding modeling work; the active-learning framework is our pitch here.
- **Best Composite Part** — Special prize for the best registered composite part; our RBP-His6 library entry.
- **Jamboree** — The annual iGEM in-person competition where projects are presented and judged (Grand Jamboree 2026-11-13).

---

**Last updated:** 2026-05-07
**Suggested place for additions:** when a new term shows up in any project doc, paste its abbreviation + 1-line definition here.
