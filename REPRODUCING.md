# Running the simulation

## Environment verified on 2026-10-20

- R 4.5.3
- nloptr 2.2.1
- tidyverse 2.0.0

Install the dependencies in an R library for your project environment:

```r
install.packages(c("nloptr", "tidyverse"), repos = "https://cloud.r-project.org")
```

From the repository root, run:

```sh
Rscript evalue.R
```

The script sets seed 123, simulates 10,000 genes, prints comparison tables, and
draws diagnostic plots. In the tested noninteractive Rscript run, the default
graphics device wrote `Rplots.pdf` in the working directory.

The command completed with exit code 0 in the environment above. This confirms
execution, not an independent check of the FDR claims or every statistical method.
The versions above record one tested environment; they are not a complete package
lockfile and numerical results may depend on dependencies and platform.
