# Batch 初步审阅优先级

优先级用于决定先核查哪批元数据与文献，不是已经证明的新颖性排名或可立即下载的最终名单。

门槛：扣除 low_reference 成员后至少 100 个有 WGS FASTQ 的去重 sample accession，至少 2 个不同 study accession。E. coli 的原记录保留。

文献指标为样本加权的 log10(1+PubMed 命中数)，越低表示此口径下关注度相对较低；不是“未被挖掘”的概率。

| 优先级 | 科/待核查分组 | 合格 WGS 样本 | GENOMIC WGS 样本 | 项目数 | 最大项目占比 | 文献关注度指标 | FASTQ GB |
|---|---|---:|---:|---:|---:|---:|---:|
| P1 先核查 | Microbacteriaceae | 196 | 144 | 64 | 14.3% | 2.479 | 120.55 |
| P1 先核查 | Actinomycetaceae | 119 | 43 | 24 | 62.2% | 2.679 | 67.71 |
| P1 先核查 | Flavobacteriaceae | 118 | 108 | 54 | 9.3% | 2.754 | 156.57 |
| P1 先核查 | Eggerthellaceae | 120 | 119 | 5 | 95.0% | 2.766 | 73.21 |
| P1 先核查 | Lachnospiraceae | 132 | 132 | 50 | 15.2% | 2.783 | 172.49 |
| P1 先核查 | Comamonadaceae | 247 | 238 | 54 | 24.7% | 2.792 | 230.58 |
| P1 先核查 | Nitrobacteraceae | 155 | 152 | 96 | 11.6% | 3.037 | 190.97 |
| P1 先核查 | Micrococcaceae | 177 | 164 | 95 | 5.6% | 3.097 | 135.98 |
| P1 先核查 | Metamycoplasmataceae | 176 | 175 | 21 | 26.1% | 3.189 | 90.59 |
| P1 先核查 | Erwiniaceae | 207 | 206 | 68 | 16.9% | 3.462 | 214.77 |
| P2 后续核查 | Alcaligenaceae | 286 | 280 | 64 | 23.1% | 3.547 | 247.91 |
| P2 后续核查 | Rhizobiaceae | 392 | 390 | 107 | 9.7% | 3.566 | 318.33 |
| P2 后续核查 | Morganellaceae | 436 | 413 | 135 | 11.0% | 3.593 | 358.78 |
| P2 后续核查 | Pasteurellaceae | 609 | 605 | 88 | 29.6% | 3.643 | 349.30 |
| P2 后续核查 | Thermaceae | 208 | 40 | 15 | 80.8% | 3.653 | 68.29 |
| P2 后续核查 | Paenibacillaceae | 331 | 331 | 62 | 61.9% | 3.676 | 372.04 |
| P2 后续核查 | Burkholderiaceae | 608 | 596 | 124 | 6.2% | 3.714 | 703.46 |
| P2 后续核查 | Lysobacteraceae | 820 | 790 | 178 | 11.5% | 3.721 | 666.24 |
| P2 后续核查 | Legionellaceae | 184 | 179 | 23 | 21.2% | 3.936 | 112.39 |
| P2 后续核查 | Lactobacillaceae | 439 | 398 | 182 | 8.4% | 3.976 | 388.09 |
| P2 后续核查 | Yersiniaceae | 591 | 588 | 132 | 11.3% | 4.038 | 470.38 |
| P2 后续核查 | Brucellaceae | 118 | 116 | 36 | 20.3% | 4.043 | 78.17 |
| P2 后续核查 | Aeromonadaceae | 131 | 130 | 42 | 22.1% | 4.070 | 126.50 |
| P2 后续核查 | Peptostreptococcaceae | 425 | 425 | 52 | 33.2% | 4.227 | 253.66 |
| P2 后续核查 | Neisseriaceae | 412 | 359 | 58 | 21.1% | 4.261 | 684.73 |
| P2 后续核查 | Bacteroidaceae | 776 | 606 | 45 | 45.2% | 4.272 | 825.02 |
| P2 后续核查 | Campylobacteraceae | 972 | 968 | 69 | 24.4% | 4.321 | 388.19 |
| P2 后续核查 | Streptomycetaceae | 1704 | 1703 | 114 | 50.8% | 4.382 | 1450.85 |
| P2 后续核查 | Clostridiaceae | 101 | 97 | 50 | 22.8% | 4.415 | 130.85 |
| P2 后续核查 | Listeriaceae | 829 | 822 | 61 | 13.5% | 4.420 | 650.23 |
| P3 低优先级参考 | Moraxellaceae | 1496 | 1475 | 200 | 24.3% | 4.422 | 897.19 |
| P3 低优先级参考 | Vibrionaceae | 1679 | 1646 | 143 | 27.8% | 4.489 | 1107.58 |
| P3 低优先级参考 | Enterococcaceae | 2112 | 1995 | 186 | 17.9% | 4.573 | 1798.13 |
| P3 低优先级参考 | Bacillaceae | 1911 | 1889 | 366 | 20.4% | 4.597 | 2697.66 |
| P3 低优先级参考 | Enterobacteriaceae | 14471 | 14195 | 860 | 4.1% | 4.666 | 12841.15 |
| P3 低优先级参考 | Helicobacteraceae | 229 | 229 | 21 | 45.0% | 4.744 | 137.16 |
| P3 低优先级参考 | Mycobacteriaceae | 2726 | 2521 | 123 | 15.7% | 4.861 | 1767.14 |
| P3 低优先级参考 | Streptococcaceae | 1753 | 1687 | 210 | 14.9% | 4.895 | 1228.87 |
| P3 低优先级参考 | Pseudomonadaceae | 2539 | 2488 | 449 | 4.8% | 5.091 | 1795.04 |
| P3 低优先级参考 | Staphylococcaceae | 3758 | 3524 | 340 | 22.2% | 5.205 | 2031.21 |
| 暂缓：数据不足 | Corynebacteriaceae | 94 | 79 | 31 | 22.3% |  | 79.82 |
| 暂缓：数据不足 | Nocardiaceae | 86 | 82 | 52 | 29.1% |  | 92.92 |
| 暂缓：数据不足 | Borreliaceae | 84 | 84 | 13 | 56.0% |  | 41.80 |
| 暂缓：数据不足 | Anaplasmataceae | 79 | 77 | 23 | 25.3% |  | 196.68 |
| 暂缓：数据不足 | Weeksellaceae | 78 | 75 | 42 | 9.0% |  | 43.61 |
| 暂缓：数据不足 | Bifidobacteriaceae | 75 | 67 | 29 | 14.7% |  | 68.98 |
| 暂缓：数据不足 | Prevotellaceae | 68 | 62 | 17 | 42.6% |  | 53.80 |
| 暂缓：数据不足 | Pectobacteriaceae | 67 | 67 | 36 | 16.4% |  | 123.66 |
| 暂缓：数据不足 | Pseudonocardiaceae | 64 | 64 | 20 | 40.6% |  | 41.17 |
| 暂缓：数据不足 | Halomonadaceae | 64 | 63 | 33 | 14.1% |  | 74.93 |
| 暂缓：数据不足 | Caryophanaceae | 61 | 61 | 13 | 78.7% |  | 66.81 |
| 暂缓：数据不足 | Prochlorococcaceae | 60 | 6 | 4 | 61.7% |  | 66.94 |
| 暂缓：数据不足 | Carnobacteriaceae | 58 | 58 | 17 | 53.4% |  | 38.96 |
| 暂缓：数据不足 | Roseobacteraceae | 58 | 53 | 21 | 22.4% |  | 69.30 |
| 暂缓：数据不足 | Shewanellaceae | 55 | 55 | 28 | 14.5% |  | 54.65 |
| 暂缓：数据不足 | Oscillospiraceae | 52 | 51 | 19 | 25.0% |  | 80.98 |
| 暂缓：数据不足 | Candidatus Pelagibacteraceae | 50 | 45 | 5 | 38.0% |  | 8.28 |
| 暂缓：数据不足 | Oxalobacteraceae | 50 | 49 | 29 | 20.0% |  | 36.11 |
| 暂缓：数据不足 | Synechococcaceae | 49 | 5 | 4 | 89.8% |  | 74.32 |
| 暂缓：数据不足 | Micromonosporaceae | 49 | 49 | 11 | 75.5% |  | 40.30 |
| 暂缓：数据不足 | Phyllobacteriaceae | 47 | 45 | 23 | 17.0% |  | 62.84 |
| 暂缓：数据不足 | Treponemataceae | 46 | 45 | 7 | 87.0% |  | 103.22 |
| 暂缓：数据不足 | Colwelliaceae | 45 | 45 | 5 | 80.0% |  | 168.74 |
| 暂缓：数据不足 | Acetobacteraceae | 45 | 45 | 20 | 22.2% |  | 37.39 |
| 暂缓：数据不足 | Caulobacteraceae | 43 | 41 | 25 | 25.6% |  | 28.23 |
| 暂缓：数据不足 | Rhodanobacteraceae | 42 | 42 | 18 | 23.8% |  | 28.15 |
| 暂缓：数据不足 | Mycoplasmoidaceae | 42 | 39 | 8 | 38.1% |  | 27.71 |
| 暂缓：数据不足 | Bartonellaceae | 42 | 42 | 12 | 47.6% |  | 29.65 |
| 暂缓：数据不足 | Tannerellaceae | 41 | 41 | 7 | 73.2% |  | 68.81 |
| 暂缓：数据不足 | Francisellaceae | 41 | 41 | 22 | 22.0% |  | 34.74 |
| 暂缓：数据不足 | Pseudoalteromonadaceae | 40 | 39 | 26 | 22.5% |  | 68.07 |
| 暂缓：数据不足 | Alteromonadaceae | 40 | 40 | 14 | 27.5% |  | 49.94 |
| 暂缓：数据不足 | Hafniaceae | 39 | 39 | 19 | 28.2% |  | 27.02 |
| 暂缓：数据不足 | Paracoccaceae | 39 | 30 | 22 | 23.1% |  | 36.42 |
| 暂缓：数据不足 | Sphingobacteriaceae | 36 | 35 | 22 | 25.0% |  | 21.59 |
| 暂缓：数据不足 | Propionibacteriaceae | 35 | 29 | 16 | 20.0% |  | 48.44 |
| 暂缓：数据不足 | Anoxybacillaceae | 34 | 33 | 12 | 29.4% |  | 41.79 |
| 暂缓：数据不足 | Methylobacteriaceae | 33 | 33 | 20 | 24.2% |  | 29.07 |
| 暂缓：数据不足 | Sphingomonadaceae | 33 | 31 | 24 | 12.1% |  | 18.34 |
| 暂缓：数据不足 | Azospirillaceae | 31 | 31 | 9 | 38.7% |  | 25.95 |
| 暂缓：数据不足 | Acholeplasmataceae | 28 | 19 | 20 | 17.9% |  | 123.65 |
| 暂缓：数据不足 | Gordoniaceae | 28 | 22 | 15 | 21.4% |  | 28.18 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2518612 | 26 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | Aerococcaceae | 25 | 20 | 10 | 48.0% |  | 15.72 |
| 暂缓：数据不足 | Akkermansiaceae | 24 | 23 | 5 | 83.3% |  | 19.64 |
| 暂缓：数据不足 | Fusobacteriaceae | 23 | 20 | 11 | 34.8% |  | 7.21 |
| 暂缓：数据不足 | Chlamydiaceae | 21 | 21 | 7 | 38.1% |  | 37.03 |
| 暂缓：数据不足 | Erysipelotrichaceae | 20 | 18 | 11 | 20.0% |  | 14.13 |
| 暂缓：数据不足 | Deinococcaceae | 19 | 19 | 7 | 68.4% |  | 6.56 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2212473 | 19 | 0 | 1 | 100.0% |  | 0.02 |
| 暂缓：数据不足 | Thermoanaerobacteraceae | 18 | 18 | 8 | 50.0% |  | 20.16 |
| 暂缓：数据不足 | Desulfovibrionaceae | 18 | 18 | 13 | 16.7% |  | 24.22 |
| 暂缓：数据不足 | Streptosporangiaceae | 18 | 18 | 2 | 94.4% |  | 6.36 |
| 暂缓：数据不足 | Marinobacteraceae | 18 | 16 | 13 | 16.7% |  | 15.75 |
| 暂缓：数据不足 | Sphaerotilaceae | 18 | 16 | 15 | 16.7% |  | 46.89 |
| 暂缓：数据不足 | unplaced_genus_33986 | 18 | 18 | 13 | 22.2% |  | 20.64 |
| 暂缓：数据不足 | Methylococcaceae | 17 | 17 | 14 | 17.6% |  | 22.08 |
| 暂缓：数据不足 | Peptoniphilaceae | 16 | 13 | 9 | 31.2% |  | 14.47 |
| 暂缓：数据不足 | Orbaceae | 15 | 15 | 3 | 60.0% |  | 13.89 |
| 暂缓：数据不足 | Cellvibrionaceae | 15 | 15 | 8 | 33.3% |  | 18.98 |
| 暂缓：数据不足 | Chroococcidiopsidaceae | 15 | 11 | 4 | 40.0% |  | 12.29 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_77133 | 15 | 2 | 5 | 46.7% |  | 197.68 |
| 暂缓：数据不足 | Leptotrichiaceae | 14 | 13 | 5 | 57.1% |  | 11.30 |
| 暂缓：数据不足 | Eubacteriaceae | 14 | 13 | 14 | 7.1% |  | 17.63 |
| 暂缓：数据不足 | Kribbellaceae | 14 | 14 | 1 | 100.0% |  | 6.75 |
| 暂缓：数据不足 | Alloboseaceae | 14 | 14 | 7 | 57.1% |  | 7.86 |
| 暂缓：数据不足 | Nocardioidaceae | 14 | 14 | 11 | 21.4% |  | 11.70 |
| 暂缓：数据不足 | Rikenellaceae | 13 | 10 | 5 | 46.2% |  | 30.06 |
| 暂缓：数据不足 | Chitinophagaceae | 13 | 12 | 5 | 30.8% |  | 11.07 |
| 暂缓：数据不足 | Brevibacteriaceae | 13 | 10 | 12 | 15.4% |  | 22.76 |
| 暂缓：数据不足 | Nostocaceae | 12 | 11 | 9 | 33.3% |  | 40.44 |
| 暂缓：数据不足 | Leptospiraceae | 12 | 12 | 6 | 50.0% |  | 14.29 |
| 暂缓：数据不足 | Microcoleaceae | 12 | 1 | 3 | 66.7% |  | 12.50 |
| 暂缓：数据不足 | Thermomonosporaceae | 12 | 12 | 3 | 58.3% |  | 15.48 |
| 暂缓：数据不足 | Mycoplasmataceae | 12 | 12 | 4 | 33.3% |  | 6.05 |
| 暂缓：数据不足 | Geminocystaceae | 12 | 12 | 10 | 16.7% |  | 11.60 |
| 暂缓：数据不足 | Stappiaceae | 12 | 12 | 6 | 41.7% |  | 15.38 |
| 暂缓：数据不足 | Caldicellulosiruptoraceae | 12 | 12 | 4 | 58.3% |  | 13.79 |
| 暂缓：数据不足 | Veillonellaceae | 12 | 11 | 7 | 41.7% |  | 9.03 |
| 暂缓：数据不足 | Promicromonosporaceae | 12 | 9 | 9 | 25.0% |  | 10.50 |
| 暂缓：数据不足 | Sporomusaceae | 11 | 11 | 11 | 9.1% |  | 11.24 |
| 暂缓：数据不足 | Nitrospiraceae | 11 | 2 | 4 | 72.7% |  | 27.46 |
| 暂缓：数据不足 | Sphingopyxidaceae | 11 | 11 | 3 | 54.5% |  | 5.69 |
| 暂缓：数据不足 | Dermabacteraceae | 11 | 11 | 7 | 45.5% |  | 16.61 |
| 暂缓：数据不足 | Chromatiaceae | 10 | 7 | 7 | 30.0% |  | 15.18 |
| 暂缓：数据不足 | Oceanospirillaceae | 10 | 10 | 7 | 20.0% |  | 17.80 |
| 暂缓：数据不足 | Microbulbiferaceae | 10 | 10 | 3 | 60.0% |  | 15.84 |
| 暂缓：数据不足 | Merismopediaceae | 10 | 10 | 3 | 40.0% |  | 23.86 |
| 暂缓：数据不足 | Pseudothioglobaceae (SeqCode) | 10 | 10 | 5 | 60.0% |  | 5.00 |
| 暂缓：数据不足 | Cellulomonadaceae | 10 | 10 | 10 | 10.0% |  | 5.90 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2715253 | 10 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | Piscirickettsiaceae | 9 | 6 | 6 | 33.3% |  | 30.36 |
| 暂缓：数据不足 | Zoogloeaceae | 9 | 4 | 5 | 33.3% |  | 16.30 |
| 暂缓：数据不足 | Spiroplasmataceae | 9 | 7 | 9 | 11.1% |  | 26.33 |
| 暂缓：数据不足 | Ilumatobacteraceae | 9 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | Pirellulaceae | 9 | 0 | 1 | 100.0% |  | 0.02 |
| 暂缓：数据不足 | Neomoorellaceae | 9 | 9 | 9 | 11.1% |  | 8.37 |
| 暂缓：数据不足 | Methylocystaceae | 9 | 9 | 9 | 11.1% |  | 11.22 |
| 暂缓：数据不足 | Cyclobacteriaceae | 9 | 3 | 4 | 66.7% |  | 5.31 |
| 暂缓：数据不足 | Rickettsiaceae | 9 | 6 | 8 | 22.2% |  | 63.69 |
| 暂缓：数据不足 | Dermatophilaceae | 9 | 9 | 2 | 77.8% |  | 6.64 |
| 暂缓：数据不足 | Chromobacteriaceae | 8 | 8 | 8 | 12.5% |  | 16.48 |
| 暂缓：数据不足 | Idiomarinaceae | 8 | 3 | 5 | 62.5% |  | 1.20 |
| 暂缓：数据不足 | Turicibacteraceae | 8 | 7 | 5 | 50.0% |  | 2.21 |
| 暂缓：数据不足 | Thalassospiraceae | 8 | 8 | 7 | 25.0% |  | 5.90 |
| 暂缓：数据不足 | Thermosynechococcaceae | 8 | 8 | 8 | 12.5% |  | 7.28 |
| 暂缓：数据不足 | Methylophilaceae | 8 | 3 | 3 | 62.5% |  | 1.69 |
| 暂缓：数据不足 | Xanthobacteraceae | 8 | 8 | 2 | 87.5% |  | 5.24 |
| 暂缓：数据不足 | Erythrobacteraceae | 8 | 8 | 6 | 25.0% |  | 5.53 |
| 暂缓：数据不足 | Nitriliruptoraceae | 8 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Hyphomonadaceae | 8 | 3 | 3 | 62.5% |  | 3.68 |
| 暂缓：数据不足 | Intrasporangiaceae | 8 | 8 | 6 | 25.0% |  | 4.25 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2052180 | 8 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | Desulfobacteraceae | 7 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | Pseudobdellovibrionaceae | 7 | 6 | 4 | 42.9% |  | 8.14 |
| 暂缓：数据不足 | Candidatus Nanopelagicaceae | 7 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Devosiaceae | 7 | 5 | 4 | 28.6% |  | 3.60 |
| 暂缓：数据不足 | Sphingobiaceae | 7 | 7 | 6 | 28.6% |  | 7.28 |
| 暂缓：数据不足 | Saprospiraceae | 7 | 3 | 3 | 57.1% |  | 3.74 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3120121 | 7 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Tolypothrichaceae | 6 | 6 | 6 | 16.7% |  | 4.31 |
| 暂缓：数据不足 | Dermacoccaceae | 6 | 6 | 5 | 33.3% |  | 10.04 |
| 暂缓：数据不足 | Selenomonadaceae | 6 | 6 | 6 | 16.7% |  | 12.09 |
| 暂缓：数据不足 | Desulfobulbaceae | 6 | 1 | 5 | 33.3% |  | 53.66 |
| 暂缓：数据不足 | Coprobacillaceae | 6 | 6 | 5 | 33.3% |  | 2.39 |
| 暂缓：数据不足 | Nodosilineaceae | 6 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | Acetivibrionaceae | 6 | 6 | 2 | 83.3% |  | 5.57 |
| 暂缓：数据不足 | Fictibacillaceae | 6 | 6 | 3 | 66.7% |  | 3.02 |
| 暂缓：数据不足 | Rhodospirillaceae | 6 | 6 | 5 | 33.3% |  | 6.19 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2201156 | 6 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2497613 | 6 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Brachyspiraceae | 5 | 5 | 2 | 80.0% |  | 6.67 |
| 暂缓：数据不足 | Allobogoriellaceae | 5 | 5 | 3 | 60.0% |  | 3.54 |
| 暂缓：数据不足 | Halieaceae | 5 | 3 | 3 | 40.0% |  | 1.20 |
| 暂缓：数据不足 | Aphanizomenonaceae | 5 | 3 | 4 | 40.0% |  | 3.36 |
| 暂缓：数据不足 | Nitrosomonadaceae | 5 | 5 | 3 | 60.0% |  | 11.54 |
| 暂缓：数据不足 | Geobacteraceae | 5 | 3 | 4 | 40.0% |  | 2.11 |
| 暂缓：数据不足 | Alcanivoracaceae | 5 | 4 | 5 | 20.0% |  | 3.84 |
| 暂缓：数据不足 | Cryomorphaceae | 5 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | Sulfurovaceae | 5 | 1 | 2 | 80.0% |  | 0.09 |
| 暂缓：数据不足 | Spirosomataceae | 5 | 5 | 5 | 20.0% |  | 2.55 |
| 暂缓：数据不足 | Tichowtungiaceae | 5 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Candidatus Izemoplasmataceae | 5 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Entomoplasmataceae | 5 | 5 | 4 | 40.0% |  | 3.05 |
| 暂缓：数据不足 | Nocardiopsidaceae | 5 | 5 | 5 | 20.0% |  | 4.16 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1898207 | 5 | 5 | 2 | 60.0% |  | 8.87 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2201154 | 5 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_91347 | 5 | 5 | 1 | 100.0% |  | 1.50 |
| 暂缓：数据不足 | Coxiellaceae | 4 | 2 | 4 | 25.0% |  | 22.58 |
| 暂缓：数据不足 | Spirochaetaceae | 4 | 4 | 1 | 100.0% |  | 1.89 |
| 暂缓：数据不足 | Thermoactinomycetaceae | 4 | 4 | 3 | 50.0% |  | 2.44 |
| 暂缓：数据不足 | Leptolyngbyaceae | 4 | 3 | 3 | 50.0% |  | 5.97 |
| 暂缓：数据不足 | Verrucomicrobiaceae | 4 | 4 | 3 | 50.0% |  | 22.96 |
| 暂缓：数据不足 | Acidithiobacillaceae | 4 | 4 | 4 | 25.0% |  | 4.71 |
| 暂缓：数据不足 | Arcobacteraceae | 4 | 4 | 1 | 100.0% |  | 1.55 |
| 暂缓：数据不足 | Thalassobaculaceae | 4 | 1 | 2 | 75.0% |  | 0.88 |
| 暂缓：数据不足 | Desulfocapsaceae | 4 | 1 | 2 | 75.0% |  | 0.97 |
| 暂缓：数据不足 | Acutalibacteraceae | 4 | 2 | 4 | 25.0% |  | 3.21 |
| 暂缓：数据不足 | Roseomonadaceae | 4 | 4 | 4 | 25.0% |  | 1.93 |
| 暂缓：数据不足 | Candidatus Njordibacteraceae | 4 | 4 | 1 | 100.0% |  | 1.56 |
| 暂缓：数据不足 | Aquificaceae | 4 | 4 | 4 | 25.0% |  | 3.34 |
| 暂缓：数据不足 | Candidatus Rhabdochlamydiaceae | 4 | 4 | 1 | 100.0% |  | 0.92 |
| 暂缓：数据不足 | Ectothiorhodospiraceae | 4 | 4 | 3 | 50.0% |  | 2.85 |
| 暂缓：数据不足 | Coriobacteriaceae | 4 | 3 | 4 | 25.0% |  | 11.01 |
| 暂缓：数据不足 | Gallionellaceae | 4 | 4 | 4 | 25.0% |  | 14.13 |
| 暂缓：数据不足 | Acidaminococcaceae | 4 | 4 | 3 | 50.0% |  | 7.05 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2026780 | 4 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2026801 | 4 | 0 | 1 | 100.0% |  | 0.01 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2030927 | 4 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Scytonemataceae | 3 | 2 | 3 | 33.3% |  | 5.30 |
| 暂缓：数据不足 | Porticoccaceae | 3 | 3 | 2 | 66.7% |  | 0.70 |
| 暂缓：数据不足 | Odoribacteraceae | 3 | 3 | 3 | 33.3% |  | 5.35 |
| 暂缓：数据不足 | Oscillatoriaceae | 3 | 3 | 3 | 33.3% |  | 2.39 |
| 暂缓：数据不足 | Kangiellaceae | 3 | 3 | 3 | 33.3% |  | 4.17 |
| 暂缓：数据不足 | Dysgonomonadaceae | 3 | 3 | 2 | 66.7% |  | 6.99 |
| 暂缓：数据不足 | Endozoicomonadaceae | 3 | 3 | 3 | 33.3% |  | 3.28 |
| 暂缓：数据不足 | Salinarimonadaceae | 3 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Saccharospirillaceae | 3 | 1 | 2 | 66.7% |  | 2.55 |
| 暂缓：数据不足 | Moritellaceae | 3 | 3 | 3 | 33.3% |  | 11.85 |
| 暂缓：数据不足 | Sulfurimonadaceae | 3 | 0 | 2 | 66.7% |  | 3.09 |
| 暂缓：数据不足 | Vicingaceae | 3 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Maricaulaceae | 3 | 1 | 2 | 66.7% |  | 0.76 |
| 暂缓：数据不足 | Allostellaceae | 3 | 3 | 1 | 100.0% |  | 1.83 |
| 暂缓：数据不足 | Zymomonadaceae | 3 | 3 | 2 | 66.7% |  | 0.17 |
| 暂缓：数据不足 | Chitinibacteraceae | 3 | 3 | 2 | 66.7% |  | 1.41 |
| 暂缓：数据不足 | Salinicoccaceae | 3 | 3 | 2 | 66.7% |  | 3.23 |
| 暂缓：数据不足 | Myxococcaceae | 3 | 3 | 2 | 66.7% |  | 7.64 |
| 暂缓：数据不足 | Anaerotignaceae | 3 | 3 | 3 | 33.3% |  | 1.42 |
| 暂缓：数据不足 | Puniceispirillaceae (SeqCode) | 3 | 3 | 2 | 66.7% |  | 0.58 |
| 暂缓：数据不足 | Litorivicinaceae | 3 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Anaerovoracaceae | 3 | 3 | 2 | 66.7% |  | 3.23 |
| 暂缓：数据不足 | Kineosporiaceae | 3 | 3 | 3 | 33.3% |  | 0.69 |
| 暂缓：数据不足 | Geodermatophilaceae | 3 | 3 | 2 | 66.7% |  | 3.82 |
| 暂缓：数据不足 | unplaced_genus_1392389 | 3 | 3 | 3 | 33.3% |  | 7.36 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2026788 | 3 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Dehalococcoidaceae | 2 | 2 | 2 | 50.0% |  | 5.22 |
| 暂缓：数据不足 | Beutenbergiaceae | 2 | 0 | 1 | 100.0% |  | 0.10 |
| 暂缓：数据不足 | Sanguibacteraceae | 2 | 2 | 2 | 50.0% |  | 0.91 |
| 暂缓：数据不足 | Petrotogaceae | 2 | 2 | 2 | 50.0% |  | 0.69 |
| 暂缓：数据不足 | Pseudanabaenaceae | 2 | 0 | 2 | 50.0% |  | 4.58 |
| 暂缓：数据不足 | Microcystaceae | 2 | 1 | 2 | 50.0% |  | 2.56 |
| 暂缓：数据不足 | Chlorobiaceae | 2 | 2 | 1 | 100.0% |  | 3.46 |
| 暂缓：数据不足 | Candidatus Hepatincolaceae | 2 | 2 | 2 | 50.0% |  | 8.04 |
| 暂缓：数据不足 | Desulfuromonadaceae | 2 | 2 | 1 | 100.0% |  | 1.75 |
| 暂缓：数据不足 | Candidatus Nanogingivalaceae | 2 | 0 | 1 | 100.0% |  | 1.47 |
| 暂缓：数据不足 | Desulfonatronaceae | 2 | 1 | 2 | 50.0% |  | 2.37 |
| 暂缓：数据不足 | Aurantimonadaceae | 2 | 2 | 2 | 50.0% |  | 12.88 |
| 暂缓：数据不足 | Steroidobacteraceae | 2 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Roseivirgaceae | 2 | 2 | 2 | 50.0% |  | 3.04 |
| 暂缓：数据不足 | Robiginitomaculaceae | 2 | 2 | 1 | 100.0% |  | 0.58 |
| 暂缓：数据不足 | Bruguierivoracaceae | 2 | 2 | 2 | 50.0% |  | 20.71 |
| 暂缓：数据不足 | Syntrophotaleaceae | 2 | 2 | 2 | 50.0% |  | 1.20 |
| 暂缓：数据不足 | Desulfobaccaceae | 2 | 2 | 2 | 50.0% |  | 6.24 |
| 暂缓：数据不足 | Aristaeellaceae | 2 | 2 | 1 | 100.0% |  | 3.35 |
| 暂缓：数据不足 | Desulfothermaceae | 2 | 2 | 2 | 50.0% |  | 0.50 |
| 暂缓：数据不足 | Sedimenticolaceae | 2 | 0 | 1 | 100.0% |  | 0.42 |
| 暂缓：数据不足 | Desulfatiglandaceae | 2 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Aerosakkonemataceae | 2 | 2 | 1 | 100.0% |  | 10.40 |
| 暂缓：数据不足 | Nodulariaceae | 2 | 2 | 2 | 50.0% |  | 2.54 |
| 暂缓：数据不足 | Pseudohongiellaceae | 2 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Butyricicoccaceae | 2 | 2 | 2 | 50.0% |  | 0.37 |
| 暂缓：数据不足 | Nitrospinaceae | 2 | 0 | 1 | 100.0% |  | 1.05 |
| 暂缓：数据不足 | Hyphomicrobiaceae | 2 | 2 | 2 | 50.0% |  | 2.03 |
| 暂缓：数据不足 | Gemellaceae | 2 | 1 | 2 | 50.0% |  | 3.40 |
| 暂缓：数据不足 | Rhodothermaceae | 2 | 2 | 1 | 100.0% |  | 5.06 |
| 暂缓：数据不足 | Synergistaceae | 2 | 2 | 2 | 50.0% |  | 6.12 |
| 暂缓：数据不足 | Rubrobacteraceae | 2 | 2 | 2 | 50.0% |  | 0.75 |
| 暂缓：数据不足 | Sphaerobacteraceae | 2 | 2 | 2 | 50.0% |  | 5.86 |
| 暂缓：数据不足 | Tsukamurellaceae | 2 | 2 | 2 | 50.0% |  | 0.31 |
| 暂缓：数据不足 | Christensenellaceae | 2 | 2 | 1 | 100.0% |  | 0.36 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1909294 | 2 | 1 | 2 | 50.0% |  | 1.62 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2026772 | 2 | 1 | 2 | 50.0% |  | 5.53 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2044939 | 2 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2053541 | 2 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2053569 | 2 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3060593 | 2 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3122005 | 2 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Chloroflexaceae | 1 | 1 | 1 | 100.0% |  | 4.76 |
| 暂缓：数据不足 | Planctomycetaceae | 1 | 1 | 1 | 100.0% |  | 0.58 |
| 暂缓：数据不足 | Oligosphaeraceae | 1 | 1 | 1 | 100.0% |  | 2.41 |
| 暂缓：数据不足 | Schleiferiaceae | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Opitutaceae | 1 | 1 | 1 | 100.0% |  | 1.66 |
| 暂缓：数据不足 | Thiotrichaceae | 1 | 1 | 1 | 100.0% |  | 0.80 |
| 暂缓：数据不足 | Marinifilaceae | 1 | 1 | 1 | 100.0% |  | 0.08 |
| 暂缓：数据不足 | Atopobiaceae | 1 | 1 | 1 | 100.0% |  | 4.38 |
| 暂缓：数据不足 | Wenzhouxiangellaceae | 1 | 1 | 1 | 100.0% |  | 0.24 |
| 暂缓：数据不足 | Porphyromonadaceae | 1 | 1 | 1 | 100.0% |  | 0.45 |
| 暂缓：数据不足 | Tissierellaceae | 1 | 1 | 1 | 100.0% |  | 1.04 |
| 暂缓：数据不足 | Mesoaciditogaceae | 1 | 1 | 1 | 100.0% |  | 6.86 |
| 暂缓：数据不足 | Tepidisphaeraceae | 1 | 1 | 1 | 100.0% |  | 0.39 |
| 暂缓：数据不足 | Candidatus Paracaedibacteraceae | 1 | 1 | 1 | 100.0% |  | 0.19 |
| 暂缓：数据不足 | Balneolaceae | 1 | 1 | 1 | 100.0% |  | 0.08 |
| 暂缓：数据不足 | Peptococcaceae | 1 | 1 | 1 | 100.0% |  | 0.34 |
| 暂缓：数据不足 | Alicyclobacillaceae | 1 | 1 | 1 | 100.0% |  | 0.17 |
| 暂缓：数据不足 | Gloeobacteraceae | 1 | 1 | 1 | 100.0% |  | 0.48 |
| 暂缓：数据不足 | Acaryochloridaceae | 1 | 1 | 1 | 100.0% |  | 1.17 |
| 暂缓：数据不足 | Pleurocapsaceae | 1 | 1 | 1 | 100.0% |  | 0.83 |
| 暂缓：数据不足 | Chlorogloeopsidaceae | 1 | 1 | 1 | 100.0% |  | 2.45 |
| 暂缓：数据不足 | Hapalosiphonaceae | 1 | 1 | 1 | 100.0% |  | 0.02 |
| 暂缓：数据不足 | Thermomicrobiaceae | 1 | 1 | 1 | 100.0% |  | 0.44 |
| 暂缓：数据不足 | Gemmataceae | 1 | 1 | 1 | 100.0% |  | 0.09 |
| 暂缓：数据不足 | Bryobacteraceae | 1 | 1 | 1 | 100.0% |  | 0.64 |
| 暂缓：数据不足 | Muribaculaceae | 1 | 1 | 1 | 100.0% |  | 0.59 |
| 暂缓：数据不足 | Barnesiellaceae | 1 | 1 | 1 | 100.0% |  | 2.15 |
| 暂缓：数据不足 | Flammeovirgaceae | 1 | 1 | 1 | 100.0% |  | 1.28 |
| 暂缓：数据不足 | Azonexaceae | 1 | 1 | 1 | 100.0% |  | 18.18 |
| 暂缓：数据不足 | Silvanigrellaceae | 1 | 1 | 1 | 100.0% |  | 0.80 |
| 暂缓：数据不足 | Thioalkalibacteraceae | 1 | 1 | 1 | 100.0% |  | 0.03 |
| 暂缓：数据不足 | Fastidiosibacteraceae | 1 | 1 | 1 | 100.0% |  | 1.88 |
| 暂缓：数据不足 | Hydrogenophilaceae | 1 | 1 | 1 | 100.0% |  | 3.53 |
| 暂缓：数据不足 | Abditibacteriaceae | 1 | 1 | 1 | 100.0% |  | 2.25 |
| 暂缓：数据不足 | Desulfomicrobiaceae | 1 | 1 | 1 | 100.0% |  | 1.75 |
| 暂缓：数据不足 | Syntrophobacteraceae | 1 | 1 | 1 | 100.0% |  | 2.98 |
| 暂缓：数据不足 | Saccharimonadaceae (SeqCode) | 1 | 0 | 1 | 100.0% |  | 0.03 |
| 暂缓：数据不足 | Vicinamibacteraceae | 1 | 1 | 1 | 100.0% |  | 0.17 |
| 暂缓：数据不足 | Hydrogenothermaceae | 1 | 1 | 1 | 100.0% |  | 0.22 |
| 暂缓：数据不足 | Hahellaceae | 1 | 1 | 1 | 100.0% |  | 1.75 |
| 暂缓：数据不足 | Nautiliaceae | 1 | 1 | 1 | 100.0% |  | 0.34 |
| 暂缓：数据不足 | Capsulimonadaceae | 1 | 1 | 1 | 100.0% |  | 0.64 |
| 暂缓：数据不足 | Baekduiaceae | 1 | 1 | 1 | 100.0% |  | 3.00 |
| 暂缓：数据不足 | Thermoguttaceae | 1 | 1 | 1 | 100.0% |  | 1.03 |
| 暂缓：数据不足 | Zhaonellaceae | 1 | 1 | 1 | 100.0% |  | 1.42 |
| 暂缓：数据不足 | Fulvivirgaceae | 1 | 1 | 1 | 100.0% |  | 3.35 |
| 暂缓：数据不足 | Marivirgaceae | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Reichenbachiellaceae | 1 | 1 | 1 | 100.0% |  | 5.84 |
| 暂缓：数据不足 | Salibacteraceae | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Kytococcaceae | 1 | 1 | 1 | 100.0% |  | 1.44 |
| 暂缓：数据不足 | Ornithinimicrobiaceae | 1 | 0 | 1 | 100.0% |  | 0.12 |
| 暂缓：数据不足 | Thermodesulfovibrionaceae | 1 | 1 | 1 | 100.0% |  | 0.21 |
| 暂缓：数据不足 | Terasakiellaceae | 1 | 1 | 1 | 100.0% |  | 0.45 |
| 暂缓：数据不足 | Sphingosinicellaceae | 1 | 1 | 1 | 100.0% |  | 0.99 |
| 暂缓：数据不足 | Segnochrobactraceae | 1 | 1 | 1 | 100.0% |  | 0.62 |
| 暂缓：数据不足 | Kaistiaceae | 1 | 1 | 1 | 100.0% |  | 0.96 |
| 暂缓：数据不足 | Phreatobacteraceae | 1 | 1 | 1 | 100.0% |  | 0.49 |
| 暂缓：数据不足 | Anaerolineaceae | 1 | 1 | 1 | 100.0% |  | 1.45 |
| 暂缓：数据不足 | Desulfococcaceae | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Sirenicapillariaceae | 1 | 1 | 1 | 100.0% |  | 3.44 |
| 暂缓：数据不足 | Desulfotomaculaceae | 1 | 1 | 1 | 100.0% |  | 0.73 |
| 暂缓：数据不足 | Desulfonauticaceae | 1 | 1 | 1 | 100.0% |  | 2.41 |
| 暂缓：数据不足 | Sporanaerobacteraceae | 1 | 1 | 1 | 100.0% |  | 0.08 |
| 暂缓：数据不足 | Tepidimicrobiaceae | 1 | 1 | 1 | 100.0% |  | 1.21 |
| 暂缓：数据不足 | Ignatzschineriaceae | 1 | 1 | 1 | 100.0% |  | 0.53 |
| 暂缓：数据不足 | Aminobacteriaceae | 1 | 1 | 1 | 100.0% |  | 0.80 |
| 暂缓：数据不足 | Geoalkalibacteraceae | 1 | 1 | 1 | 100.0% |  | 4.05 |
| 暂缓：数据不足 | Desulfitibacteraceae | 1 | 1 | 1 | 100.0% |  | 0.78 |
| 暂缓：数据不足 | Cloacimonadaceae (SeqCode) | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Thioreductoraceae | 1 | 1 | 1 | 100.0% |  | 0.04 |
| 暂缓：数据不足 | Candidatus Moduliflexaceae | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Thiovibrionaceae | 1 | 1 | 1 | 100.0% |  | 0.54 |
| 暂缓：数据不足 | Anthocerotibacteraceae | 1 | 1 | 1 | 100.0% |  | 19.97 |
| 暂缓：数据不足 | Thiohalobacteraceae | 1 | 1 | 1 | 100.0% |  | 0.51 |
| 暂缓：数据不足 | Filifactoraceae | 1 | 1 | 1 | 100.0% |  | 0.62 |
| 暂缓：数据不足 | Gallaecimonadaceae | 1 | 1 | 1 | 100.0% |  | 0.50 |
| 暂缓：数据不足 | Halobacillaceae | 1 | 1 | 1 | 100.0% |  | 1.18 |
| 暂缓：数据不足 | Candidatus Chloroheliaceae | 1 | 1 | 1 | 100.0% |  | 0.98 |
| 暂缓：数据不足 | Fidelibacteraceae | 1 | 1 | 1 | 100.0% |  | 5.29 |
| 暂缓：数据不足 | Acidocellaceae | 1 | 1 | 1 | 100.0% |  | 0.48 |
| 暂缓：数据不足 | Guptibacillaceae | 1 | 1 | 1 | 100.0% |  | 1.80 |
| 暂缓：数据不足 | Sphingomicrobiaceae | 1 | 1 | 1 | 100.0% |  | 0.65 |
| 暂缓：数据不足 | Sphingorhabdaceae | 1 | 1 | 1 | 100.0% |  | 1.52 |
| 暂缓：数据不足 | Candidatus Levibacteraceae | 1 | 1 | 1 | 100.0% |  | 0.06 |
| 暂缓：数据不足 | Fervidibacteraceae | 1 | 1 | 1 | 100.0% |  | 1.66 |
| 暂缓：数据不足 | Fusibacteraceae | 1 | 1 | 1 | 100.0% |  | 1.94 |
| 暂缓：数据不足 | Desulfarculaceae | 1 | 1 | 1 | 100.0% |  | 0.24 |
| 暂缓：数据不足 | Caldilineaceae | 1 | 1 | 1 | 100.0% |  | 0.83 |
| 暂缓：数据不足 | Polyangiaceae | 1 | 1 | 1 | 100.0% |  | 0.38 |
| 暂缓：数据不足 | Acanthopleuribacteraceae | 1 | 1 | 1 | 100.0% |  | 0.64 |
| 暂缓：数据不足 | Halobacteroidaceae | 1 | 1 | 1 | 100.0% |  | 0.03 |
| 暂缓：数据不足 | Symbiobacteriaceae | 1 | 1 | 1 | 100.0% |  | 2.19 |
| 暂缓：数据不足 | Desulfurobacteriaceae | 1 | 1 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Marinilabiliaceae | 1 | 1 | 1 | 100.0% |  | 0.14 |
| 暂缓：数据不足 | Holophagaceae | 1 | 1 | 1 | 100.0% |  | 0.99 |
| 暂缓：数据不足 | Mariprofundaceae | 1 | 1 | 1 | 100.0% |  | 0.82 |
| 暂缓：数据不足 | Iamiaceae | 1 | 1 | 1 | 100.0% |  | 0.78 |
| 暂缓：数据不足 | Phycisphaeraceae | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | Frankiaceae | 1 | 1 | 1 | 100.0% |  | 3.23 |
| 暂缓：数据不足 | Dietziaceae | 1 | 1 | 1 | 100.0% |  | 0.07 |
| 暂缓：数据不足 | Nakamurellaceae | 1 | 0 | 1 | 100.0% |  | 0.03 |
| 暂缓：数据不足 | Cytophagaceae | 1 | 1 | 1 | 100.0% |  | 1.60 |
| 暂缓：数据不足 | unplaced_genus_114248 | 1 | 1 | 1 | 100.0% |  | 7.82 |
| 暂缓：数据不足 | unplaced_genus_1381133 | 1 | 1 | 1 | 100.0% |  | 0.28 |
| 暂缓：数据不足 | unplaced_genus_1470577 | 1 | 1 | 1 | 100.0% |  | 0.66 |
| 暂缓：数据不足 | unplaced_genus_190972 | 1 | 1 | 1 | 100.0% |  | 1.88 |
| 暂缓：数据不足 | unplaced_genus_204619 | 1 | 1 | 1 | 100.0% |  | 1.62 |
| 暂缓：数据不足 | unplaced_genus_2066914 | 1 | 1 | 1 | 100.0% |  | 0.67 |
| 暂缓：数据不足 | unplaced_genus_3459270 | 1 | 0 | 1 | 100.0% |  | 1.09 |
| 暂缓：数据不足 | unplaced_genus_415014 | 1 | 1 | 1 | 100.0% |  | 2.26 |
| 暂缓：数据不足 | unplaced_genus_495819 | 1 | 1 | 1 | 100.0% |  | 4.12 |
| 暂缓：数据不足 | unplaced_genus_54066 | 1 | 1 | 1 | 100.0% |  | 0.91 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1232452 | 1 | 1 | 1 | 100.0% |  | 2.38 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1232457 | 1 | 1 | 1 | 100.0% |  | 1.45 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1563636 | 1 | 1 | 1 | 100.0% |  | 0.71 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1760 | 1 | 1 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1778434 | 1 | 1 | 1 | 100.0% |  | 0.80 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1879010 | 1 | 1 | 1 | 100.0% |  | 0.30 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1891238 | 1 | 1 | 1 | 100.0% |  | 0.23 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1909295 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1913989 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2026720 | 1 | 0 | 1 | 100.0% |  | 3.61 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2044936 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2044940 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2044943 | 1 | 1 | 1 | 100.0% |  | 0.77 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2052159 | 1 | 0 | 1 | 100.0% |  | 96.85 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2053615 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_205975 | 1 | 0 | 1 | 100.0% |  | 1.78 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2496632 | 1 | 0 | 1 | 100.0% |  | 8.18 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2562987 | 1 | 1 | 1 | 100.0% |  | 0.78 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2572087 | 1 | 1 | 1 | 100.0% |  | 0.06 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2589994 | 1 | 1 | 1 | 100.0% |  | 0.33 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2663014 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2746493 | 1 | 0 | 1 | 100.0% |  | 0.46 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2758719 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2774202 | 1 | 1 | 1 | 100.0% |  | 0.25 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2802972 | 1 | 1 | 1 | 100.0% |  | 3.44 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2829083 | 1 | 1 | 1 | 100.0% |  | 0.03 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2910168 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2918525 | 1 | 1 | 1 | 100.0% |  | 0.98 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3034175 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_305753 | 1 | 1 | 1 | 100.0% |  | 0.40 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3104632 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3108503 | 1 | 1 | 1 | 100.0% |  | 0.92 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3109346 | 1 | 1 | 1 | 100.0% |  | 3.59 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3233533 | 1 | 1 | 1 | 100.0% |  | 3.39 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3235145 | 1 | 1 | 1 | 100.0% |  | 0.09 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3235146 | 1 | 1 | 1 | 100.0% |  | 0.16 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3235148 | 1 | 1 | 1 | 100.0% |  | 0.04 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3243068 | 1 | 1 | 1 | 100.0% |  | 2.09 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3391387 | 1 | 1 | 1 | 100.0% |  | 2.40 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3395322 | 1 | 1 | 1 | 100.0% |  | 0.11 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3395323 | 1 | 1 | 1 | 100.0% |  | 0.12 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3409803 | 1 | 1 | 1 | 100.0% |  | 0.02 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3410811 | 1 | 1 | 1 | 100.0% |  | 2.97 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3410812 | 1 | 1 | 1 | 100.0% |  | 4.02 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3411324 | 1 | 1 | 1 | 100.0% |  | 0.10 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3411325 | 1 | 1 | 1 | 100.0% |  | 0.33 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3411326 | 1 | 1 | 1 | 100.0% |  | 0.73 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3411327 | 1 | 1 | 1 | 100.0% |  | 0.02 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3411328 | 1 | 1 | 1 | 100.0% |  | 0.77 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3411329 | 1 | 1 | 1 | 100.0% |  | 0.51 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3411330 | 1 | 1 | 1 | 100.0% |  | 0.08 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3411331 | 1 | 1 | 1 | 100.0% |  | 0.06 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3417971 | 1 | 1 | 1 | 100.0% |  | 0.08 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3420184 | 1 | 1 | 1 | 100.0% |  | 0.44 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3436344 | 1 | 1 | 1 | 100.0% |  | 1.07 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3447365 | 1 | 1 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3450247 | 1 | 1 | 1 | 100.0% |  | 0.37 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3450248 | 1 | 1 | 1 | 100.0% |  | 0.35 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3450249 | 1 | 1 | 1 | 100.0% |  | 1.16 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3450250 | 1 | 1 | 1 | 100.0% |  | 0.52 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3450251 | 1 | 1 | 1 | 100.0% |  | 0.84 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3458707 | 1 | 1 | 1 | 100.0% |  | 2.27 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3465769 | 1 | 1 | 1 | 100.0% |  | 0.27 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3465770 | 1 | 1 | 1 | 100.0% |  | 0.74 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3472973 | 1 | 0 | 1 | 100.0% |  | 1.12 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3472974 | 1 | 0 | 1 | 100.0% |  | 1.02 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3477787 | 1 | 1 | 1 | 100.0% |  | 2.88 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3681801 | 1 | 0 | 1 | 100.0% |  | 1.24 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3693643 | 1 | 1 | 1 | 100.0% |  | 0.14 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3701362 | 1 | 0 | 1 | 100.0% |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_384638 | 1 | 1 | 1 | 100.0% |  | 2.13 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_713059 | 1 | 1 | 1 | 100.0% |  | 0.08 |
| 暂缓：数据不足 | Ferrovaceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Salinibacteraceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Thermotogaceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Immundisolibacteraceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Sterolibacteriaceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Candidatus Amarolineaceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Sphaerochaetaceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Nitratiruptoraceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Candidatus Scalinduaceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Candidatus Magnetomoraceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Candidatus Phycocordibacteraceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Beijerinckiaceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | Succinivibrionaceae | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_genus_1918454 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1648056 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_1805100 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_207831 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_2608721 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230096 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230097 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230098 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230099 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230100 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230101 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230102 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230103 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230104 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230105 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230106 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230107 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230108 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230109 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230110 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230111 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3230112 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386001 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386002 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386003 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386004 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386005 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386006 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386007 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386008 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386009 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386010 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386011 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386012 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386013 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386014 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386015 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386016 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3386017 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3447660 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3447661 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3447662 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：数据不足 | unplaced_unplaced_taxon_3447663 | 0 | 0 | 0 |  |  | 0.00 |
| 暂缓：分类待核查 | unplaced_unplaced_taxon_1869227 | 861 | 343 | 10 | 57.8% |  | 255.98 |
