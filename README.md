# GeoGenoPlot
A python package for plotting correlations between variants and environmental factors


## Installation
```
pip install geogenoplot
```

## Usage

```
geogenoplot [-h] [-p] [-n ENV_FACTOR_NAME] [-s {pearson,spearman,kruskal,all}] [-t] [--lat LAT] [--lon LON] vcf_file raster_file georef_table var_id

positional arguments:
  vcf_file              input vmo directory
  raster_file           output file
  georef_table          output file
  var_id                output file

optional arguments:
  -h, --help            show this help message and exit
  -p, --phased          phased or not, default is not phased
  -n ENV_FACTOR_NAME, --env_factor_name ENV_FACTOR_NAME
                        environmental factor name
  -s {pearson,spearman,kruskal,all}, --statistic {pearson,spearman,kruskal,all}
                        correlation statistic, choices are pearson, spearman, kruskal, all. Default is all.
  -t, --no_trendline    no trendline, default is to show trendline
  --lat LAT             latitude range, default is -37,51
  --lon LON             longitude range, default is -20,95
```