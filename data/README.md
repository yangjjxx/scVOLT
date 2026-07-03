# Data Directory

This directory stores the processed input datasets used by the scVOLT paper analyses. Large data files are intentionally not tracked in Git. For review or reproduction, download the processed data archive from the [Figshare private link](https://figshare.com/s/4ce6b177cd5f04ea7cba) and unpack it at the repository root so that the paths below are restored exactly.

Do not rename these files unless the corresponding notebook paths are updated.

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

## Manifest

The files below are processed analysis inputs corresponding to the manuscript experiments. Primary public sources and accessions are described in the manuscript and supplementary dataset table.

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

## Archiving Recommendation

The processed data files are provided separately from the code repository through the [Figshare private link](https://figshare.com/s/4ce6b177cd5f04ea7cba). The code repository should keep this README but should not commit the large processed files.

To restore the data after download:

```bash
unzip scvolt_data.zip -d .
```

If the downloaded archive has a different file name, use that file name in the command. The extracted files should match the layout and checksums listed above.
