# Paper analysis files

This folder stores the scripts and small result artifacts used for the paper
analyses. These files are bundled with the package as documentation/supporting
materials; they are not sourced when the R package is loaded.

After installing the package, the folder can be located with:

```r
system.file("docs", package = "conticox")
```

## Folder layout

- `simulation/`: simulation driver scripts, scenario search scripts, small
  simulation result objects, and simulation figures.
- `real_data/`: TCGA and other real-data analysis scripts, selected summary
  tables, AUC curve figures, and median-risk KM plots.

Most scripts retain the original project-level path assumptions. For exact
reproduction, run them from the project root
`/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls` or adjust the paths at the top of
the script before running from an installed package location.

## Not bundled

Large binary data/cache files are intentionally not copied here, so the package
remains lightweight. These remain in the project root or subfolders:

- `data_cache/`
- large root-level `TCGA_*.RData` validation result files
- large root-level `*_auc_curves_list_*.RData` objects
- `figures/*_fit_results.rds`
- large `candidate_results/*.RData` objects
- downloaded raw expression archives or compressed matrices

The bundled scripts document how those files are produced or consumed.
