# Data Directory

The processed data files are provided separately from this code repository through the [Figshare private link](https://figshare.com/s/4ce6b177cd5f04ea7cba). To keep the repository lightweight, the large processed data files are not committed to GitHub.

After downloading the data archive, place it in the root directory of this repository and run:

```bash
unzip scvolt_data.zip -d .
```

The extracted files should follow the directory structure described in the data documentation.

## Expected Layout

```text
data/
  01_single_cell_prediction/
    451Lu_PLX-4720_anndata.h5ad
    Docetaxel_data.h5ad
    Paclitaxel_data.h5ad
    PLX-4720_data.h5ad
    Sorafenib_data.h5ad
  02_time_series_analysis/
    GSE114461_Bortezomib.csv.gz
  03_spatial_analysis/
    GDSC2_fitted_dose_response_27Oct23.xlsx
    V19L29-097.h5ad
```

## Usage

After extracting the data archive, run the notebooks from either the repository root or the `analyses/` directory:

- `analyses/single_cell_prediction.ipynb` uses `data/01_single_cell_prediction/`.
- `analyses/time_series_analysis.ipynb` uses `data/02_time_series_analysis/GSE114461_Bortezomib.csv.gz`.
- `analyses/spatial_analysis.ipynb` uses `data/03_spatial_analysis/V19L29-097.h5ad` and `data/03_spatial_analysis/GDSC2_fitted_dose_response_27Oct23.xlsx`.

Generated outputs should be written to `results/`, not to this directory.

## Public Sources and Study Roles

The processed files in this directory are derived from the public transcriptomic and pharmacogenomic records listed below. Single-cell sensitive/resistant labels are used as benchmark annotations for evaluation only; they are not treated as directly measured per-cell dose-response ground truth.

| Source | Modality and biological context | Drug/condition and groups | Scale and study role |
| --- | --- | --- | --- |
| GEO GSE140440 | Single-cell total RNA-seq from human prostate cancer docetaxel-resistance models based on DU145 and PC3 cell lines. | Docetaxel; each cell line includes sensitive and resistant variants. | 324 single cells: DU145, 72 sensitive and 72 resistant; PC3, 90 sensitive and 90 resistant. Used as the Docetaxel benchmark task. |
| GEO GSE114461 | 10x Genomics scRNA-seq from an MCF7 breast cancer cell-line drug-treatment time course. | Bortezomib, 500 nM; MCF7 strain AA measured at t0, t12, t48, and t96. The t96 condition corresponds to 72 h exposure followed by drug washout and 24 h recovery. | 7,555 cells: t0, 1,726; t12, 2,743; t48, 1,851; t96, 1,235. Used for time-course response analysis. |
| GEO GSE163836 | scRNA-seq from human inflammatory breast cancer chemotherapy-resistance models. | Paclitaxel; FCIBC02 parental cells and paclitaxel-resistant FCIBC02PR and SUM149PR derivatives. | 17,202 cells: FCIBC02, 1,978; FCIBC02PR, 5,483; SUM149PR, 9,741. Used as the Paclitaxel benchmark task. |
| GEO GSE175716 | scRNA-seq from a hepatocellular carcinoma sorafenib response/sensitivity comparison. | Sorafenib; primary response group (PR), sensitive group A (SA), and sensitive group B (SB). | 16,068 cells: PR, 4,989; SA, 4,774; SB, 6,305. Used as the Sorafenib benchmark task. |
| GEO GSE108383 (451Lu) | SMART-seq/Fluidigm scRNA-seq from a human melanoma BRAF-inhibitor resistance model. | PLX-4720; 451Lu parental/sensitive and resistant states after inhibitor exposure. | 197 single cells: 84 parental/sensitive and 113 resistant. Used as the PLX-4720--451Lu benchmark task. |
| GEO GSE108383 (A375) | SMART-seq/Fluidigm scRNA-seq from a human melanoma BRAF-inhibitor resistance model. | PLX-4720; A375 parental/sensitive and resistant states after inhibitor exposure. | 130 single cells: 78 parental/sensitive and 52 resistant. Used as the PLX-4720--A375 benchmark task. |
| 10x Genomics V19L29-097 | Visium spatial transcriptomics from fresh-frozen invasive ductal carcinoma breast tissue, AJCC/UICC Stage Group IIA, ER positive, PR negative, HER2 positive. | Computational Lapatinib sensitivity mapping and Dasatinib complementary-candidate analysis; histological regions include Healthy, Invasive, Surrounding tumor, and Tumor. | 3,798 tissue-covered spots. The analyzed grouping includes 2,490 tumor-like spots and 1,308 non-tumor spots. Used for spatial concordance and drug-coverage analyses. |
| GDSC2 release 8.5 | Bulk cell-line pharmacogenomic dose-response screen across cancer models. | Multi-drug screen; the Dasatinib subset is used as external pharmacogenomic evidence for the spatial candidate analysis. | 242,036 dose-response records covering 969 cell lines and 286 drugs. The BRCA subset has 13,106 records; the Dasatinib subset has 962 records, including 50 BRCA records. |

Counts refer to the public records analyzed in this study. The GDSC2 row reports cell-line dose-response records, not cells. For V19L29-097, tumor-like spots combine Tumor and Invasive regions, and non-tumor spots combine Healthy and Surrounding tumor regions.

## Manifest

The files below are processed analysis inputs corresponding to the manuscript experiments. Public source records, biological context, drug or treatment groups, and study roles are summarized above.

| File | Used in | Description | Size | SHA256 |
| --- | --- | --- | ---: | --- |
| `data/01_single_cell_prediction/451Lu_PLX-4720_anndata.h5ad` | single-cell benchmark | Processed AnnData benchmark input for PLX-4720 response prediction. | 7.49 MB | `75211DBC47C286A3CBF5AC86035C9CBF290E8AC7967C32E8E1446C61B412A2AE` |
| `data/01_single_cell_prediction/Docetaxel_data.h5ad` | single-cell benchmark | Processed AnnData benchmark input for Docetaxel response prediction. | 32.81 MB | `4C07AB1E1C2E75AD5DA7D1D66A6F81A9DD4D47F32A19B575B1E0CD5C8983C3A4` |
| `data/01_single_cell_prediction/Paclitaxel_data.h5ad` | single-cell benchmark | Processed AnnData benchmark input for Paclitaxel response prediction. | 562.14 MB | `5B2C61E6F7FE917F6AD52DFD99F7E324C03D4D30487993A2F831918E360F1D6F` |
| `data/01_single_cell_prediction/PLX-4720_data.h5ad` | single-cell benchmark | Processed AnnData benchmark input for PLX-4720 response prediction. | 13.69 MB | `C81D71B3660E56CA973850AC1DBBC3826F92C432EEED22629E181FE566509054` |
| `data/01_single_cell_prediction/Sorafenib_data.h5ad` | single-cell benchmark | Processed AnnData benchmark input for Sorafenib response prediction. | 11.33 MB | `75AFB323A3D763F8B142B470E5FB965B9ED50407B8E561210D142B377DDA917A` |
| `data/02_time_series_analysis/GSE114461_Bortezomib.csv.gz` | time-course analysis | Processed expression matrix for the Bortezomib time-course analysis. | 24.94 MB | `1CC65D617214DE9C73E8B4EA4ED889862AC950D94066099FB58175A58B324B73` |
| `data/03_spatial_analysis/GDSC2_fitted_dose_response_27Oct23.xlsx` | spatial analysis | GDSC2 fitted dose-response table used for the external Dasatinib breast-cancer cell-line comparison. | 20.34 MB | `F950A7027BE265F8A7A74220A27FD18CBD368485349BD8C2048E88BB1CD07560` |
| `data/03_spatial_analysis/V19L29-097.h5ad` | spatial analysis | Processed spatial transcriptomics object used for spatial sensitivity and coverage analyses. | 216.57 MB | `D8D53A8937FD9DA1D022C67A0C32A5FB5547A60EE8EF30BDDB64DA383FC06463` |



