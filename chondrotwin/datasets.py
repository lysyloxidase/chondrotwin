"""Registry and downloader for the open datasets used by the knee pipeline (:mod:`chondrotwin.knee`).

Nothing here is redistributed with the package: files are downloaded from the
original providers into a local data directory (default ``data/raw`` in the
repository, override with ``CHONDROTWIN_DATA``), checked against the sizes and
SHA-256 sums below, and must be cited / licensed as stated.

    python scripts/fetch_data.py            # all datasets (~93 MB)
    python scripts/fetch_data.py orthoload  # one source
"""
from __future__ import annotations

import hashlib
import os
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional


@dataclass(frozen=True)
class Source:
    name: str
    license: str
    citation: str
    landing_page: str


@dataclass(frozen=True)
class RemoteFile:
    source: str
    filename: str
    url: str
    size: int
    sha256: Optional[str] = None


SOURCES: Dict[str, Source] = {
    "orthoload": Source(
        "orthoload",
        "OrthoLoad terms: non-commercial use; cite the article AND the OrthoLoad database",
        "Bergmann G, Bender A, Graichen F, et al. Standardized loads acting in knee implants. "
        "PLoS ONE 2014;9(1):e86035. doi:10.1371/journal.pone.0086035; OrthoLoad database, orthoload.com",
        "https://orthoload.com/test-loads/standard-loads-knee-joint/"),
    "gse104782": Source(
        "gse104782",
        "NCBI GEO public data",
        "Ji Q, Zheng Y, Zhang G, et al. Single-cell RNA-seq analysis reveals the progression of human "
        "osteoarthritis. Ann Rheum Dis 2019;78(1):100-110. GEO GSE104782",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104782"),
    "gse152805": Source(
        "gse152805",
        "NCBI GEO public data",
        "Chou CH, Jain V, Gibson J, et al. Synovial cell cross-talk with cartilage plays a major role in the "
        "pathogenesis of osteoarthritis. Sci Rep 2020;10:10868. GEO GSE152805",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152805"),
    "gse57218": Source(
        "gse57218",
        "NCBI GEO public data",
        "Ramos YFM, den Hollander W, Bovee JVMG, et al. Genes involved in the osteoarthritis process identified "
        "through genome wide expression analysis in articular cartilage; the RAAK study. PLoS ONE 2014;9:e103056. "
        "GEO GSE57218 (platform annotation GPL6947)",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57218"),
    "gse114007": Source(
        "gse114007",
        "NCBI GEO public data",
        "Fisch KM, Gamini R, Alvarez-Garcia O, et al. Identification of transcription factors responsible for "
        "dysregulated networks in human osteoarthritis cartilage by global gene expression analysis. "
        "Osteoarthritis Cartilage 2018;26:1531-1538. GEO GSE114007",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114007"),
    "senmayo": Source(
        "senmayo",
        "MSigDB gene set SAUL_SEN_MAYO (CC BY 4.0, MSigDB terms)",
        "Saul D, Kosinsky RL, Atkinson EJ, et al. A new gene set identifies senescent cells and predicts "
        "senescence-associated pathways across tissues. Nat Commun 2022;13:4827",
        "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SAUL_SEN_MAYO.html"),
    "gse117999": Source(
        "gse117999", "NCBI GEO public data",
        "Rai MF, Tycksen ED, Cai L, et al. Distinct degenerative phenotype of articular cartilage from knees with meniscus "
        "tear compared to knees with osteoarthritis. Osteoarthritis Cartilage 2019;27:945-955. GEO GSE117999",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117999"),
    "gse169077": Source(
        "gse169077", "NCBI GEO public data",
        "GEO GSE169077, Dysregulated gene expression in human osteoarthritis cartilage (Attur M, Abramson SB; no primary "
        "publication linked)",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169077"),
    "hallmark_inflammation": Source(
        "hallmark_inflammation", "MSigDB gene set HALLMARK_INFLAMMATORY_RESPONSE (CC BY 4.0, MSigDB terms)",
        "Liberzon A, Birger C, Thorvaldsdottir H, et al. The Molecular Signatures Database hallmark gene set collection. "
        "Cell Syst 2015;1:417-425",
        "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INFLAMMATORY_RESPONSE.html"),
    "gse178557": Source(
        "gse178557", "NCBI GEO public data",
        "GEO GSE178557: Microarray identifies differentially expressed mRNAs of cartilages between osteoarthritis (OA) "
        "patients and control patients (no primary publication linked in GEO)",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178557"),
    "gse55235": Source(
        "gse55235", "NCBI GEO public data",
        "Woetzel D, Huber R, Kupfer P, et al. Identification of rheumatoid arthritis and osteoarthritis patients by "
        "transcriptome-based rule set generation. Arthritis Res Ther 2014;16:R84. GEO GSE55235, GSE55457",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55235"),
    "oaizib_cm": Source(
        "oaizib_cm",
        "CC BY-NC 4.0; citing the CartiMorph paper is compulsory",
        "Yao Y, et al. CartiMorph: a framework for automated knee articular cartilage morphometrics. "
        "Med Image Anal 2024. Dataset OAIZIB-CM, doi:10.5281/zenodo.14934086 (labels derived from OAI-ZIB, "
        "Ambellan F, et al. Med Image Anal 2019; images from the Osteoarthritis Initiative)",
        "https://zenodo.org/records/14934086"),
}

_OL = "https://orthoload.com/wp-content/uploads/"
_GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104782/"
_ZEN = "https://zenodo.org/records/14934086/files/"
_GEOS = "https://ftp.ncbi.nlm.nih.gov/geo/samples/"
_GEO57 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57218/"
_GEO114 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114007/"

FILES: List[RemoteFile] = [
    RemoteFile("hallmark_inflammation", "HALLMARK_INFLAMMATORY_RESPONSE.json", "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_INFLAMMATORY_RESPONSE&fileType=json", 1894,
               "35203e49f5dbc3d0ada06ef005bc555a07666b860a05ac0de12bc6a6c19773ae"),
    RemoteFile("gse178557", "GSE178557_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178557/matrix/GSE178557_series_matrix.txt.gz", 1107674,
               "23762ddcb9828c43e71bad7d5af06b82f72a731777e6ca222012c788d55acefd"),
    RemoteFile("gse178557", "GSM5393765_OA-1.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393765/suppl/GSM5393765_OA-1.txt.gz", 2190554,
               "f3996395c6befc465bf5b60e11964bc32670d23413e8b7518658465fe7b1d7b7"),
    RemoteFile("gse178557", "GSM5393766_OA-2.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393766/suppl/GSM5393766_OA-2.txt.gz", 2200544,
               "753d8bf86e3624790b6ff12fef31d73c94ca90b5a7d65b3b73d30ed534f2d6b0"),
    RemoteFile("gse178557", "GSM5393767_OA-3.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393767/suppl/GSM5393767_OA-3.txt.gz", 2171310,
               "d9a4787542235851dfa7e95bbbbacd8e23c4e1be7a25cc4289e58d96e46dbf3c"),
    RemoteFile("gse178557", "GSM5393768_OA-4.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393768/suppl/GSM5393768_OA-4.txt.gz", 2180803,
               "14d3e66fbf76fed502cc2964bbba051c29764d27afeeed697948881f9a5270fe"),
    RemoteFile("gse178557", "GSM5393769_Con-1.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393769/suppl/GSM5393769_Con-1.txt.gz", 2209449,
               "dc38b75a5692810dd5f762125e4a3f23e151ed345478fcf6290f7ddd2eef24c1"),
    RemoteFile("gse178557", "GSM5393770_Con-2.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393770/suppl/GSM5393770_Con-2.txt.gz", 2228479,
               "4b11f914520f5580a7589d7aca9d78b3a3e4d5123b27f9a67d9e4a7c648c93d9"),
    RemoteFile("gse178557", "GSM5393771_Con-3.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393771/suppl/GSM5393771_Con-3.txt.gz", 2216206,
               "1482a397e4d65c832b6a02632ad019e4a580cbd797038362161c51942af84930"),
    RemoteFile("gse178557", "GSM5393772_Con-4.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393772/suppl/GSM5393772_Con-4.txt.gz", 2183707,
               "d9bf90db7a54bc912037076515f6cf8f532c21452148c5769d8f12762a7f354e"),
    RemoteFile("gse55235", "GSE55235_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55235/matrix/GSE55235_series_matrix.txt.gz", 1374133,
               "a9d2ad2cc673a5a06aa7ea95c4eb830488505e2ec99f78422af5eb1a1047b57d"),
    RemoteFile("gse55235", "GSE55457_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55457/matrix/GSE55457_series_matrix.txt.gz", 1687146,
               "9817a63dc041dd684c1039258bd2f6f0dfa7e1f5d0fb26f09577c294ca508596"),
    RemoteFile("senmayo", "SAUL_SEN_MAYO.json", "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=SAUL_SEN_MAYO&fileType=json", 1271,
               "4716e3238ea4a82e65f09ee05dfe073dcfa428e95c307a7e1c3044431a8559e9"),
    RemoteFile("gse117999", "GSE117999_normalized_matrix.xlsx", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117999/suppl/GSE117999_normalized_matrix.xlsx", 19772578,
               "d2564611995ae035c3b6d81974666572a8d302477533d42f19c1347fb40cabb0"),
    RemoteFile("gse117999", "GSE117999_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117999/matrix/GSE117999_series_matrix.txt.gz", 7078568,
               "fe56852f74992b08e5bfbbe056ce69191501bd73eabaf74b83c1916a91215e0a"),
    RemoteFile("gse169077", "GSE169077_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169077/matrix/GSE169077_series_matrix.txt.gz", 987919,
               "b15fed9991e4d132b7dfeda0b34524a74f644974e1e10cab80a3d1919b4224b7"),
    RemoteFile("gse169077", "GPL96.annot.gz", "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz", 4522748,
               "88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4"),
    RemoteFile("gse114007", "GSE114007_raw_counts.xlsx", _GEO114 + "suppl/GSE114007_raw_counts.xlsx", 4049817,
               "e904f0d75fa6b3c7944589238be003311fbd9b0ab1b0def887f823280ec8c530"),
    RemoteFile("gse57218", "GSE57218_series_matrix.txt.gz", _GEO57 + "matrix/GSE57218_series_matrix.txt.gz", 18676256,
               "d2c631f8b85fa6a65cf134992ed13181f1e63702a602a1c696596bcd38b0ab9a"),
    RemoteFile("gse57218", "GPL6947.annot.gz", "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6947/annot/GPL6947.annot.gz", 6956785,
               "244c4ef6394abd8c87309318f3f27e9e63c45d2d1c3e89303818f1f63a290824"),
    RemoteFile("gse114007", "GSE114007_normal_normalized.counts.txt.gz", _GEO114 + "suppl/GSE114007_normal_normalized.counts.txt.gz", 1477804,
               "12197be538a85359ca548956176af7e609ce8463240a02c086bf266c2c30701a"),
    RemoteFile("gse114007", "GSE114007_OA_normalized.counts.txt.gz", _GEO114 + "suppl/GSE114007_OA_normalized.counts.txt.gz", 1662255,
               "2036d4bc9b526b286c40def57d51380bfbe93a18505ded820db3725ee96f6788"),
    RemoteFile("gse114007", "GSE114007-GPL11154_series_matrix.txt.gz", _GEO114 + "matrix/GSE114007-GPL11154_series_matrix.txt.gz", 3826,
               "01f498484d1060fd58e27b88696956d45815fa1cc8ad393b7832fd339650eb9f"),
    RemoteFile("gse114007", "GSE114007-GPL18573_series_matrix.txt.gz", _GEO114 + "matrix/GSE114007-GPL18573_series_matrix.txt.gz", 3865,
               "d0c66b65010b9f7a710e9033b74550d00eed5b905d76f22d45723dd16a954ca4"),
    RemoteFile("gse152805", "GSM4626766_OA_oLT_113.barcodes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626766/suppl/GSM4626766_OA_oLT_113.barcodes.tsv.gz", 25188,
               "0918d0b850a25c1cb44cb499f84af10c6743e050ded71b2052379e3fe2541a64"),
    RemoteFile("gse152805", "GSM4626766_OA_oLT_113.genes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626766/suppl/GSM4626766_OA_oLT_113.genes.tsv.gz", 257283,
               "7e35722524a00e68eccf7cb27541d4b77c69f562aa58c489226ce353ac5f4642"),
    RemoteFile("gse152805", "GSM4626766_OA_oLT_113.matrix.mtx.gz",
               _GEOS + "GSM4626nnn/GSM4626766/suppl/GSM4626766_OA_oLT_113.matrix.mtx.gz", 44541103,
               "101b752b309ff76eb4d06bcf29a0b9d0c456d4facaf4bc8d07a21ab0e6155db4"),
    RemoteFile("gse152805", "GSM4626767_OA_oLT_116.barcodes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626767/suppl/GSM4626767_OA_oLT_116.barcodes.tsv.gz", 16138,
               "8d88712115e2dee308b3ae606798653a3ab182f415c6f71ee10c48be214d9202"),
    RemoteFile("gse152805", "GSM4626767_OA_oLT_116.genes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626767/suppl/GSM4626767_OA_oLT_116.genes.tsv.gz", 257283,
               "4dde563d7f55abd68ad4d0fb41cd6fc8ece726d4f040f3df4da29f088b321b65"),
    RemoteFile("gse152805", "GSM4626767_OA_oLT_116.matrix.mtx.gz",
               _GEOS + "GSM4626nnn/GSM4626767/suppl/GSM4626767_OA_oLT_116.matrix.mtx.gz", 31459196,
               "1cd1d1237809e5bbec2b98a52a6ec509a55b503890776eb2f09ef3eff6467526"),
    RemoteFile("gse152805", "GSM4626768_OA_oLT_118.barcodes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626768/suppl/GSM4626768_OA_oLT_118.barcodes.tsv.gz", 19379,
               "56ab7b98bd7112b33bdd2bc590d1ff0a929147089e6877e1b3f6840f56049bc0"),
    RemoteFile("gse152805", "GSM4626768_OA_oLT_118.genes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626768/suppl/GSM4626768_OA_oLT_118.genes.tsv.gz", 257283,
               "f3aedd85fe7757a0b0acf9ea49b3cb8a0708a7d5f7cb40d65a8008cc5c80d39d"),
    RemoteFile("gse152805", "GSM4626768_OA_oLT_118.matrix.mtx.gz",
               _GEOS + "GSM4626nnn/GSM4626768/suppl/GSM4626768_OA_oLT_118.matrix.mtx.gz", 45518932,
               "b523e3063ca02bc523140844395c9d3fb4e4757da832bf9416162b347fb59a13"),
    RemoteFile("gse152805", "GSM4626769_OA_MT_113.barcodes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626769/suppl/GSM4626769_OA_MT_113.barcodes.tsv.gz", 13735,
               "03c1c8887815d2fbce1acc304c88de2ecc0330a677f618e05e268010d985036a"),
    RemoteFile("gse152805", "GSM4626769_OA_MT_113.genes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626769/suppl/GSM4626769_OA_MT_113.genes.tsv.gz", 257282,
               "655289925924a0bebba331853259bd9c481610edc48462233412233683ffe10b"),
    RemoteFile("gse152805", "GSM4626769_OA_MT_113.matrix.mtx.gz",
               _GEOS + "GSM4626nnn/GSM4626769/suppl/GSM4626769_OA_MT_113.matrix.mtx.gz", 31733896,
               "d40bace26266c65fa7a59c09d71d64e7c4a8989bcb8e3cfc6d8ea7fe64d9d7c3"),
    RemoteFile("gse152805", "GSM4626770_OA_MT_116.barcodes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626770/suppl/GSM4626770_OA_MT_116.barcodes.tsv.gz", 15192,
               "eca53e5466fa4fc930cb2b47444b2d148177bb11d45cfee06bf3bf256f3a3327"),
    RemoteFile("gse152805", "GSM4626770_OA_MT_116.genes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626770/suppl/GSM4626770_OA_MT_116.genes.tsv.gz", 257282,
               "f103ffb21cf1881c8ccc183e7c8ace1a523cab4bd9ff4e20d12faee490035d5f"),
    RemoteFile("gse152805", "GSM4626770_OA_MT_116.matrix.mtx.gz",
               _GEOS + "GSM4626nnn/GSM4626770/suppl/GSM4626770_OA_MT_116.matrix.mtx.gz", 37945778,
               "5258d4a7885a01aed8816c4691933308228cfa4eceb3abef31c262fe2d9b3745"),
    RemoteFile("gse152805", "GSM4626771_OA_MT_118.barcodes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626771/suppl/GSM4626771_OA_MT_118.barcodes.tsv.gz", 20296,
               "b802e940093395258dd86dff5f38d7c1e11deb781a6dfc0e2039b4cfcfc18b50"),
    RemoteFile("gse152805", "GSM4626771_OA_MT_118.genes.tsv.gz",
               _GEOS + "GSM4626nnn/GSM4626771/suppl/GSM4626771_OA_MT_118.genes.tsv.gz", 257282,
               "e208aa1e4c982cf7deef7c69179bdd7c6c0a58f69242712ab208e257f37f7e4c"),
    RemoteFile("gse152805", "GSM4626771_OA_MT_118.matrix.mtx.gz",
               _GEOS + "GSM4626nnn/GSM4626771/suppl/GSM4626771_OA_MT_118.matrix.mtx.gz", 48154426,
               "4c68ebe923b93cd5294278897c77aab7f5c069d588489b9d4c835b45c12ff02d"),
    RemoteFile("orthoload", "AverageSubject_AllLoadLevels_Data_Diagrams.zip",
               _OL + "AverageSubject_AllLoadLevels_Data_Diagrams.zip", 4688908,
               "f02c615b5b978a455575207d4c7f324bb148d90a9d5ccc916ec8d04e8b4aaed9"),
    RemoteFile("orthoload", "K1L_Walking.xlsx",
               _OL + "K1L_Walking.xlsx", 821573,
               "d00e4b20cd94567ddcd2d58034accc5bff005b43377ba7101614e6aec4d65387"),
    RemoteFile("orthoload", "K2L_Walking.xlsx",
               _OL + "K2L_Walking.xlsx", 807486,
               "f8f50ebf78823b29646527b55c17476666fd4611e0fcbade196db5928eb5c56f"),
    RemoteFile("orthoload", "K3R_Walking.xlsx",
               _OL + "K3R_Walking.xlsx", 863176,
               "8d25b2edb3c632a2afab8bc2e5e0fd7f6f65c4dec6bf4a6a86338f697a8f633b"),
    RemoteFile("orthoload", "K5R_Walking.xlsx",
               _OL + "K5R_Walking.xlsx", 817409,
               "16ebffb796dd883b5c49f152ecb68101835e5cb395f8cf6be5fa319834b22523"),
    RemoteFile("orthoload", "K6L_Walking.xlsx",
               _OL + "K6L_Walking.xlsx", 866083,
               "9b98a9a1b1e8f358b81f0360d7927cacd4c129ae77ac58bf201f00f0ba335861"),
    RemoteFile("orthoload", "K7L_Walking.xlsx",
               _OL + "K7L_Walking.xlsx", 836789,
               "54836fc92a9b9f367a4d02972205193420465b8dfaa298c077e8fc5441839aec"),
    RemoteFile("orthoload", "K8L_Walking.xlsx",
               _OL + "K8L_Walking.xlsx", 856370,
               "28367e56ce6cee716bdecad73ac9fd430ea559a2ca73118d318a030dd95b9a52"),
    RemoteFile("orthoload", "K9L_Walking.xlsx",
               _OL + "K9L_Walking.xlsx", 768721,
               "fe56e27a7c9b2580f32d091e2bd1c32fcbb235e740a214f09ab81cbc0ad88e08"),
    RemoteFile("gse104782", "GSE104782_allcells_UMI_count.txt.gz",
               _GEO + "suppl/GSE104782_allcells_UMI_count.txt.gz", 10227738,
               "544d03f5165cb1c50a413a120110494c2076dd775ae237c916827118023b014f"),
    RemoteFile("gse104782", "GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx",
               _GEO + "suppl/GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx", 112348,
               "2e0b31c40637e1dab36553a1f3bb9abc93d97acf416e49e379a9be623483098c"),
    RemoteFile("gse104782", "GSE104782_series_matrix.txt.gz",
               _GEO + "matrix/GSE104782_series_matrix.txt.gz", 3349,
               "4adeca7d77aff97a3c97ff80604d3333674899592067a821b0ac63ba607b2c61"),
    RemoteFile("oaizib_cm", "info.zip",
               _ZEN + "info.zip?download=1", 51680,
               "bb2007ba6ea67f49c35d14ffd0ed3dfe003a756ce953d89ba658166afb66b05e"),
    RemoteFile("oaizib_cm", "labelsTr.zip",
               _ZEN + "labelsTr.zip?download=1", 57538129,
               "6d202f6cd2905723d59e33510f2f49716a83b11aabc71fa71aa84b17aa978ddd"),
    RemoteFile("oaizib_cm", "labelsTs.zip",
               _ZEN + "labelsTs.zip?download=1", 14122969,
               "72a81672e6af56a1ecff4dbf3ef66c127497288ff8eefe7ac0124994d9623c6d"),
]


def data_dir(root: Optional[os.PathLike] = None) -> Path:
    if root is not None:
        return Path(root)
    env = os.environ.get("CHONDROTWIN_DATA")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[1] / "data" / "raw"


def path(source: str, filename: str, root: Optional[os.PathLike] = None) -> Path:
    return data_dir(root) / source / filename


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(f: RemoteFile, root: Optional[os.PathLike] = None) -> bool:
    p = path(f.source, f.filename, root)
    if not p.exists():
        return False
    if f.size and p.stat().st_size != f.size:
        return False
    return f.sha256 is None or _sha256(p) == f.sha256


def fetch(sources: Optional[Iterable[str]] = None, root: Optional[os.PathLike] = None, quiet: bool = False) -> List[Path]:
    """Download the files of the given sources (default: all) unless already present and valid."""
    wanted = set(sources or SOURCES)
    unknown = wanted - set(SOURCES)
    if unknown:
        raise ValueError(f"unknown sources: {sorted(unknown)}")
    out = []
    for f in FILES:
        if f.source not in wanted:
            continue
        p = path(f.source, f.filename, root)
        if not verify(f, root):
            p.parent.mkdir(parents=True, exist_ok=True)
            tmp = p.with_suffix(p.suffix + ".part")
            req = urllib.request.Request(f.url, headers={"User-Agent": "chondrotwin-fetch/0.4"})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            tmp.replace(p)
            if not verify(f, root):
                raise IOError(f"{f.filename}: size or checksum mismatch (got {p.stat().st_size} bytes, sha256 {_sha256(p)})")
        if not quiet:
            print(f"{f.source:10s} {f.filename:75s} {p.stat().st_size / 1e6:8.2f} MB  sha256 {_sha256(p)[:16]}")
        out.append(p)
    return out


def require(source: str, filename: str, root: Optional[os.PathLike] = None) -> Path:
    p = path(source, filename, root)
    if not p.exists():
        raise FileNotFoundError(f"{p} not found - run `python scripts/fetch_data.py {source}` first "
                                f"(license: {SOURCES[source].license})")
    return p
