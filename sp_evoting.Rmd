---
title: "Spatial panel e-voting"
author: "Märten Veskimäe"
date: "`r format(Sys.time(), '%B %d, %Y')`"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = T,
                      warning = F,
                      message = F,
                      comment = NA,
                      fig.align="center")
options(width=120)
Sys.setlocale("LC_ALL", "UTF-8")
```

## Data setup
```{r data}
library(readxl)
library(dplyr)
library(stargazer)
library(ggplot2)
library(ggthemes)

library(sp)
library(maptools)
library(spdep)
library(rgdal)

library(splm)
library(plm)
library(SpatialTools)

library(lmtest)
library(car)
library(sandwich)

library(reshape2)
library(tmap)

# Spatial data
epsg3301 = "+proj=lcc +lat_1=59.33333333333334 +lat_2=58 +lat_0=57.51755393055556 +lon_0=24 +x_0=500000 +y_0=6375000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
coord = spTransform(readShapePoly("omavalitsus_shp/omavalitsus.shp", proj4string=CRS(epsg3301)), CRS("+init=epsg:4267"))
coord@data$ONIMI = iconv(coord@data$ONIMI, "latin1", "UTF-8")
coord@data$MNIMI = iconv(coord@data$MNIMI, "latin1", "UTF-8")

# library(rgeos)
# coord = SpatialPolygonsDataFrame(gSimplify(coord, tol=.01, topologyPreserve=T), data=coord@data)




## Panel data
y = 1999
evotes = read_excel("evote.xlsx") %>% filter(valimine=="kov", aasta>=y) %>%
  mutate(ONIMI = ringkond,
         aasta = relevel(as.factor(aasta), as.character(y)),
         turnout = h22letanud/valijad_nimekirjas,
         evote = e_valis/h22letanud,
         evote = ifelse(is.na(evote),0,evote))

## Balancing
tmp = left_join(coord@data["ONIMI"], evotes[c("ONIMI", "aasta", "turnout", "evote", "competition", "candidates", "valijad_nimekirjas")], by="ONIMI")
df = pdata.frame(tmp, index=c("ONIMI","aasta"), drop.index=F, row.names=T)
df = make.pbalanced(df)

full.set = df %>%
  group_by(ONIMI) %>%
  summarize(s = sum(turnout))

df = subset(df, !(df$ONIMI %in% as.character(full.set$ONIMI[is.na(full.set$s)]))) %>%
  mutate(ONIMI = droplevels(ONIMI))

coord = subset(coord, coord$ONIMI %in% levels(df$ONIMI))

# spts = ggplot2::fortify(coord)
# write.csv(spts[c("long", "lat", "order", "id")], "coords.csv")
# write.csv(df,"evote.csv")
```
Iot obtain a balanced panel dataset, districts with missing cells were excluded. Using `balance.type="shared"` in `make.pbalanced(...)` resulted in an empty set.

## Descriptive image
```{r dplot}
df[order(df$aasta),] %>% ggplot() +
  geom_point(aes(evote, turnout, color=ONIMI), show.legend = F, size=.5, alpha=.5) +
  geom_line(aes(evote, turnout, color=ONIMI, group=ONIMI), show.legend = F, size=.5, alpha=.3) +
  theme_bw()
```

## Model specification
```{r formula}
f = as.formula(turnout ~ evote + competition + candidates + log(valijad_nimekirjas) + aasta)
f
```

Year dummies were included iot avoid bias due to low evote diffusion in 2005.

## Panel data approach
### Panel models
```{r panel, results="asis"}
pool = plm(f, data=df, model="pooling")
bw = plm(f, data=df, model="between")
fe = plm(f, data=df, model="within")
re = plm(f, data=df, model="random")

stargazer(pool,bw,fe,re,
          column.labels=c("Pooled OLS",
                          "Between effects",
                          "Fixed effects",
                          "Random effects"),
          star.cutoffs=c(0.05, 0.01, 0.001),
          omit = c("Constant"),
          type="html")
```

### F-test
```{r ftest}
pFtest(fe, pool)
```
FEs are different from zero. FE should preferred.

### Hausman
```{r hausman}
phtest(fe, re)
```
FE should preferred.

###  Groupwise heteroscedasticity and contemporaneous correlation
```{r hcc}
bptest(f, data=df, studentize=F)
pcdtest(fe, test="cd")
```
Heteroscedasticity is present.

### Corrected standard errors
```{r}
coeftest(fe, vcovHC(fe, method="arellano", cluster="group"))
```

### Residuals
```{r feresid}
resid.fe = resid(fe)
names(resid.fe) = paste0(df$ONIMI,"-",df$aasta)
fe.res = data.frame(ONIMI = gsub("-*[[:digit:]]", "", names(resid.fe)),
                    year = gsub("^.*\\-", "", names(resid.fe)),
                    residuals = resid(fe),
                    row.names = NULL) %>%
  dcast(ONIMI ~ year)

coord.res = coord
coord.res@data = left_join(coord.res@data,fe.res, by="ONIMI")
qtm(shp = coord.res, fill = c("1999", "2002", "2005","2009","2013"), fill.palette ="Blues", ncol = 3)
```
FE residuals do seem to have spatial correlation.

## Spatial panel
### Spatial weights
```{r weightmatrix}
dists = as.matrix(dist(cbind(as.vector(coordinates(coord)[,1]), as.vector(coordinates(coord)[,2]))))

W = mat2listw(dists, style="W") # row standardised
B = mat2listw(dists, style="B") # basic binary coding
C = mat2listw(dists, style="C") # globally standardised
U = mat2listw(dists, style="U") # equal to C divided by the number of neighbours
minmax = mat2listw(dists, style="minmax")
S = mat2listw(dists, style="S") # variance-stabilizing coding scheme

wm = W
```
Neighbour based weights did not work for some reason.

### Tests
```{r random}
bsktest(f, df, listw = wm, test = "LM1")
```
Random regional effects are present.

```{r}
bsktest.mod = function (formula, data, index = NULL, listw, standardize=F, ...) 
{
    if (!is.null(index)) {
        data <- plm.data(data, index)
    }
    index <- data[, 1]
    tindex <- data[, 2]
    x <- model.matrix(formula, data = data)
    y <- model.response(model.frame(formula, data = data))
    cl <- match.call()
    names(index) <- row.names(data)
    ind <- index[which(names(index) %in% row.names(x))]
    tind <- tindex[which(names(index) %in% row.names(x))]
    oo <- order(tind, ind)
    x <- x[oo, ]
    y <- y[oo]
    ind <- ind[oo]
    tind <- tind[oo]
    N <- length(unique(ind))
    k <- dim(x)[[2]]
    time <- max(tapply(x[, 1], ind, length))
    NT <- length(ind)
    ols <- lm(y ~ x - 1)
    XpXi <- solve(crossprod(x))
    n <- dim(ols$model)[1]
    indic <- seq(1, time)
    inde <- as.numeric(rep(indic, each = N))
    ind1 <- seq(1, N)
    inde1 <- as.numeric(rep(ind1, time))
    bOLS <- coefficients(ols)
    e <- as.matrix(residuals(ols))
    ee <- crossprod(e)
    Ws <- listw2dgCMatrix(listw)
    Wst <- t(Ws)
    WWp <- (Ws + Wst)/2
    yy <- function(q) {
        wq <- WWp %*% q
        wq <- as.matrix(wq)
    }
    IWWpe <- unlist(tapply(e, inde, yy))
    H <- crossprod(e, IWWpe)/crossprod(e)
    W2 <- Ws %*% Ws
    WW <- crossprod(Ws)
    tr <- function(R) sum(diag(R))
    b <- tr(W2 + WW)
    LM2 <- sqrt((N^2 * time)/b) * as.numeric(H)
    s <- NT - k
    lag <- function(QQ) lag.listw(listw, QQ)
    fun2 <- function(Q) unlist(tapply(Q, inde, lag))
    Wx <- apply(x, 2, fun2)
    WX <- matrix(Wx, NT, k)
    XpWx <- crossprod(x, WX)
    D2M <- XpWx %*% XpXi
    Ed2 <- (time * sum(diag(Ws)) - tr(D2M))/s
    WWx <- apply(WX, 2, fun2)
    WWX <- matrix(WWx, NT, k)
    XpWWX <- crossprod(x, WWX)
    spb <- XpWWX %*% XpXi
    spbb <- tr(spb)
    tpb <- XpWx %*% XpXi %*% XpWx %*% XpXi
    fintr2 <- time * tr(W2) - 2 * spbb + tr(tpb)
    Vd2 <- 2 * (s * fintr2 - (sum(diag(D2M))^2))/(s^2 * (s + 
        2))
    We <- unlist(tapply(e, inde, function(W) lag.listw(listw, 
        W)))
    d2 <- crossprod(e, We)/ee
    SLM2 <- (d2 - Ed2)/sqrt(Vd2)
    statistics <- if (standardize) 
        abs(SLM2)
    else abs(LM2)
    pval <- 2 * pnorm(statistics, lower.tail = FALSE)
    names(statistics) <- if (standardize) 
        "SLM2"
    else "LM2"
    method <- "Baltagi, Song and Koh LM2 marginal test"
    dname <- deparse(formula)
    RVAL <- list(statistic = statistics, method = method, p.value = pval, 
        data.name = deparse(formula), alternative = "Spatial autocorrelation")
    class(RVAL) <- "htest"
    return(RVAL)
}
```

```{r}
bsktest.mod(f, df, listw = wm, test="LM2")
```
Some spatial autocorrelation is detected.

```{r}
sphtest(f, data = df, listw = wm, spatial.model = "error", method = "ML")
```
RE can be used.

### ML Random effects model
```{r ml.models}
RE = spml(f, df, listw=wm, model="random", lag=T, spatial.error="b")
summary(RE)
```

### GM models
```{r gm.models}
GM.RE = spgm(f, df,
          listw = wm,
          model = "random",
          moments = "fullweights",
          lag=T, spatial.error = T)

GM.FE = spgm(f, df,
          listw = wm,
          model = "within",
          moments = "fullweights",
          lag=T, spatial.error = T)

sphtest(GM.RE, GM.FE,mehod="GM")
```
Hausman suggests FE.

```{r}
summary(GM.FE)
```

