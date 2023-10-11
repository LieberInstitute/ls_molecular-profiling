#Goal of this script: Format OUD summart statistics into a snploc file for MAGMA. 
library(here)
library(sessioninfo)

#Load in the OUD statistics
OUD <- read.delim(file = here("MAGMA",
                              "SNP_Data",
                              "GWAS_Tables",
                              "OUD",
                              "OUD.phs001672.pha004954.txt"),
                  comment.char = "#",
                  header = TRUE)

#Snploc files need to be rsID, chromosome, location. 
write.table(x = data.frame(SNP.ID = OUD$SNP.ID,
                           Chr.ID = OUD$Chr.ID,
                           Pos    = OUD$Chr.Position),
            file = here("MAGMA",
                        "SNP_Data",
                        "GWAS_Tables",
                        "OUD",
                        "OUD.phs001672.pha004954.snploc"),
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

Sys.time()
#[1] "2023-10-03 16:48:13 EDT"
options(width = 120)
# session_info()
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-10-03
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version date (UTC) lib source
# cli           3.6.1   2023-03-23 [2] CRAN (R 4.3.1)
# colorout    * 1.2-2   2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# here        * 1.0.1   2020-12-13 [2] CRAN (R 4.3.1)
# rprojroot     2.0.3   2022-04-02 [2] CRAN (R 4.3.1)
# sessioninfo * 1.2.2   2021-12-06 [2] CRAN (R 4.3.1)
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
