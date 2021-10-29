# WKREF1

Workshop on ICES reference points

## SOFTWARE installation

Packages are available in both source and binary (Windows) format for R 4.1

```r
R.version
```

To INSTALL necessary packages and all dependencies, both CRAN and FLR, run the following

```
install.packages(c("mse", "mseviz", "FLSRTMB"),
  repos=structure(c(CRAN="https://cloud.r-project.org/",
  FLR="https://flr-project.org/R")))
```

Code in this repository is organized following the ICES TAF guidelines. To ensure you get the precise versions used in the runs, please do

```
install.packages("icesTAF", repos="https://cloud.r-project.org/")
library(icesTAF)

taf.bootstrap()
```

This will install packages from source, so [Rtools](https://cran.r-project.org/bin/windows/Rtools/) needs to be installed if run from Windows.


