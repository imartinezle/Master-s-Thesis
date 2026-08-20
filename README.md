# Master-s-Thesis
═══════════════════════════════════════════════════════════════════
 0. WET LAB — Muestra de tejido HGSOC
═══════════════════════════════════════════════════════════════════
INPUT:   Tejido tumoral HGSOC (4 pacientes: GIN63, GIN65, GIN67, GIN71)
PROCESO: Corte y montaje en Visium slide → hibridación mRNA →
         cDNA → secuenciación (Illumina NovaSeq)
OUTPUT:  Archivos FASTQ (GEX) + imagen H&E de la sección Visium (.tif)


═══════════════════════════════════════════════════════════════════
 0b. ANOTACIÓN PATOLÓGICA (paralelo, sección adyacente)
═══════════════════════════════════════════════════════════════════
INPUT:   Sección H&E adyacente (NO la misma sección de Visium)
PROCESO: Anotación manual en QuPath (tumor / stroma / stroma_linfos)
         → exportación → registro de imagen (Elastix) para
         transferir anotaciones a la sección Visium
OUTPUT:  PIT_17/20/22/28.csv (barcode → anotación patólogo)


═══════════════════════════════════════════════════════════════════
 SPACE RANGER — Procesamiento inicial (cluster HPC, Slurm)
═══════════════════════════════════════════════════════════════════
INPUT:   FASTQs + imagen H&E + referencia GRCh38 + alineación
         manual revisada (Loupe Browser, .json)
PROCESO: spaceranger count (alineamiento, matriz de conteo espacial,
         vinculación con imagen)
OUTPUT:  filtered_feature_bc_matrix.h5, tissue_positions_list.csv,
         scalefactors_json.json, imagen escalada, .cloupe


═══════════════════════════════════════════════════════════════════
 SCRIPT 00 — Data wrangling patólogos (Python)
═══════════════════════════════════════════════════════════════════
INPUT:   Output Space Ranger (barcodes, tissue_positions) +
         PIT_XX.csv (anotación patólogo)
PROCESO: Cruce barcode ↔ anotación, unificación de categorías
         (tumor / stroma / stroma_linfos / unannotated)
OUTPUT:  GIN{XX}_patologist_annotation.h5ad


═══════════════════════════════════════════════════════════════════
 SCRIPT 01 — Filtrado y normalización (Python, notebook por paciente)
═══════════════════════════════════════════════════════════════════
INPUT:   GIN{XX}_patologist_annotation.h5ad
PROCESO: QC (nº genes, UMIs, %MT) → filtrado por spot y gen
         (thresholds visuales por paciente) → normalización
         (10,000 reads/spot, log1p) → HVGs (top 3,000)
OUTPUT:  GIN{XX}_filtered_normalized.h5ad


═══════════════════════════════════════════════════════════════════
 SCRIPT 02a — Clustering Leiden / Louvain / GraphST (Python)
═══════════════════════════════════════════════════════════════════
INPUT:   GIN{XX}_filtered_normalized.h5ad
PROCESO: PCA (30 PCs) → grafo vecinos (k=15) → Leiden/Louvain
         (res=0.6) → GraphST (embedding espacial + GMM, q=7)
OUTPUT:  GIN{XX}_louvain_leiden_graphst.h5ad


═══════════════════════════════════════════════════════════════════
 SCRIPT 02b — Clustering BayesSpace (R/.qmd)
═══════════════════════════════════════════════════════════════════
INPUT:   GIN{XX}_louvain_leiden_graphst.h5ad (convertido a SCE)
PROCESO: Modelo bayesiano con prior espacial (Potts model),
         MCMC (10,000 iter + 1,000 burn-in), q=7
OUTPUT:  GIN{XX}_bayesspace.h5ad


═══════════════════════════════════════════════════════════════════
 SCRIPT 03 — ESTIMATE (R/.qmd)
═══════════════════════════════════════════════════════════════════
INPUT:   GIN{XX}_bayesspace.h5ad
PROCESO: ssGSEA sobre expresión log-normalizada → StromalScore,
         ImmuneScore, ESTIMATEScore, TumorPurity (proxy) →
         3 thresholds de background (A/B/C) + filtro espacial
         RANN (k=6, máx 1 vecino tumoral)
OUTPUT:  GIN{XX}_estimate.h5ad (incluye columnas booleanas de
         background por método)


═══════════════════════════════════════════════════════════════════
 SCRIPT 04 — InferCNV (R/.qmd)
═══════════════════════════════════════════════════════════════════
INPUT:   GIN{XX}_estimate.h5ad + gencode_v19_gene_post.txt +
         counts crudos (layer "counts")
PROCESO: Grupos reference/observation (según threshold ESTIMATE) →
         moving average por ventana genómica → denoise (cutoff=0.1)
         → subclustering Leiden (res 0.01 o 0.03) → HMM (i6)
OUTPUT:  GIN{XX}_infercnv.h5ad + heatmap InferCNV + CNV score/spot


═══════════════════════════════════════════════════════════════════
 SCRIPT 05 — Cell2location (Python, notebook)
═══════════════════════════════════════════════════════════════════
INPUT:   GIN{XX}_infercnv.h5ad (counts crudos) +
         ATLAS_MSK.h5ad (referencia scRNA-seq, 37 subtipos)
PROCESO: Regresión binomial negativa (reference signatures) →
         modelo espacial (GPU, batch_size=2500) → posterior
         q05_cell_abundance_w_sf
OUTPUT:  GIN{XX}_cell2location.h5ad (abundancia por tipo celular
         y spot)


═══════════════════════════════════════════════════════════════════
 SCRIPT 06 — CellChat (R/.qmd) — SOLO LECTURA, sin output h5ad
═══════════════════════════════════════════════════════════════════
INPUT:   GIN{XX}_cell2location.h5ad (o clustering + coordenadas)
PROCESO: Comunicación ligando-receptor entre spots/dominios
         (identifyOverExpressedGenes, smoothData + PPI.human,
         contact.dependent=FALSE)
OUTPUT:  Redes de interacción (netVisual_circle), lista de genes
         por cluster de interés (consulta STRING)


