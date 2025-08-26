## ----echo = FALSE-------------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
# library(surveysd)
# 
# set.seed(1234)
# eusilc <- demo.eusilc(n = 2, prettyNames = TRUE)
# 
# eusilc[1:5, .(year, povertyRisk, gender, pWeight)]

## ----eval = FALSE-------------------------------------------------------------
# dat_boot_rw <- draw.bootstrap(eusilc,
#                               method = "Rao-Wu",
#                               REP = 10,
#                               hid = "hid",
#                               weights = "pWeight",
#                               strata = "region",
#                               period = "year")

## -----------------------------------------------------------------------------
# dat_boot_calib <- recalib(dat_boot_rw,
#                           conP.var = "gender",
#                           conH.var = "region",
#                           epsP = 1e-2,
#                           epsH = 2.5e-2,
#                           verbose = FALSE)
# dat_boot_calib[1:5, .(year, povertyRisk, gender, pWeight, w1, w2, w3, w4)]

