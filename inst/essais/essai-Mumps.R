load("./inst/essais/mdat01.rda")
dat <- droplevels(mdat01[, -attr(mdat01, "toRemove")])
dat <- dat[, -c(1,2)] # 1=date fout la merde, 2 n'a qu'un level
dat <- dat[, -c(81,82,83)] # potency columns
nas <- unique(which(is.na(dat), arr.ind = TRUE)[,"row"])
dat <- droplevels(dat[-nas,])
gfi <- gfiUltra(
  Potency.by.titration..CVP...Log.CCID50.mL...Two.Decimal.Places. ~ .,
  data = dat
)
head(gfi$models)
