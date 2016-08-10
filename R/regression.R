
.waldTest <- function(Est, cov) {
    if (length(Est) == 1) {
        W <- (Est^2)/cov
    } else {
        W <- as.numeric(t(Est) %*% solve(cov) %*% Est)
    }
    pval <- pchisq(W, df=length(Est), lower.tail=FALSE)
    c(W, pval)
}

.runRegression <- function(model.string, model.data, model.type) {
    model.formula <- as.formula(model.string)
    tryCatch({
        if (model.type == "linear") {
            mod <- lm(model.formula, data=model.data)
        } else if (model.type == "logistic") {
            mod <- glm(model.formula, data=model.data, family=binomial())
        }
    
        Est <- unname(coef(mod)["genotype"])
        cov <- vcov(mod)["genotype","genotype"]
        
        c(Est, sqrt(cov), .waldTest(Est, cov))
    }, warning=function(w) rep(NA, 4), error=function(e) rep(NA, 4))
}

.runFirth <- function(model.string, model.data, geno.index=NULL) {
    model.formula <- as.formula(model.string)
    tryCatch({
        mod <- tryCatch({
            logistf(model.formula, data=model.data, plconf=geno.index, dataout=FALSE)
        }, error=function(e) {
            ## test will fail if geno.index is too large
            geno.index <- which(colnames(model.matrix(model.formula, model.data)) == "genotype")
            logistf(model.formula, data=model.data, plconf=geno.index, dataout=FALSE)
        })
        ## test will be wrong if geno.index is too small
        ind <- which(mod$terms == "genotype")
        if (ind != geno.index) {
            mod <- logistf(model.formula, data=model.data, plconf=ind, dataout=FALSE)
        }
    
        Est <- unname(coef(mod)[ind])
        cov <- vcov(mod)[ind,ind]
        pval <- unname(mod$prob[ind])
        Stat <- qchisq(pval, df=1, lower.tail=FALSE)
        
        c(Est, sqrt(cov), Stat, pval)
    }, warning=function(w) rep(NA, 4), error=function(e) rep(NA, 4))
}

.outputNames <- function(model.type) {
    switch(model.type,
           linear=c("Est", "SE", "Wald.Stat", "Wald.Pval"),
           logistic=c("Est", "SE", "Wald.Stat", "Wald.Pval"),
           firth=c("Est", "SE", "PPL.Stat", "PPL.Pval"))
}
    
.freqByOutcome <- function(model.data, model.type, outcome) {
    if (model.type %in% c("logistic", "firth")) {
        c0 <- model.data[[outcome]] == 0
        c1 <- model.data[[outcome]] == 1
        c(n0=sum(!is.na(model.data[["genotype"]][c0])),
          n1=sum(!is.na(model.data[["genotype"]][c1])),
          freq0=0.5*mean(model.data[["genotype"]][c0], na.rm=TRUE),
          freq1=0.5*mean(model.data[["genotype"]][c1], na.rm=TRUE))
    } else {
        c(n=sum(!is.na(model.data[["genotype"]])),
          freq=0.5*mean(model.data[["genotype"]], na.rm=TRUE))
    }
}

.freqFromBinary <- function(freq) {
    freq[is.na(freq)] <- 0
    unname((freq["n0"]*freq["freq0"] + freq["n1"]*freq["freq1"]) /
           (freq["n0"] + freq["n1"]))
}

setMethod("regression",
          "SeqVarData",
          function(gdsobj, outcome, covar=NULL,
                   model.type=c("linear", "logistic", "firth")) {
              model.type <- match.arg(model.type)
              
              ## get covariates
              dat <- pData(sampleData(gdsobj))[,c(outcome, covar),drop=FALSE]
              
              ## create model formula
              model.string <- paste(outcome, "~", paste(c(covar, "genotype"), collapse=" + "))
              
              ## for firth test - determine index of genotype in model matrix
              if (model.type == "firth") {
                  tmp <- cbind(dat, "genotype"=0)
                  geno.index <- which(colnames(model.matrix(as.formula(model.string), tmp)) == "genotype")
                  rm(tmp)
              }

              ## apply function over variants
              res <- seqApply(gdsobj, "genotype", function(x) {
              #res <- seqApply(gdsobj, c(geno="genotype", id="variant.id"), function(x) {
                  #print(x$id)
                  ## assume we want effect of reference allele
                  model.data <- cbind(dat, genotype=colSums(x == 0))
                  model.data <- droplevels(model.data[complete.cases(model.data),])
                  
                  ## don't bother with monomorphic variants
                  freq <- .freqByOutcome(model.data, model.type, outcome)
                  freq.all <- if (model.type == "linear") freq["freq"] else .freqFromBinary(freq)
                  if (freq.all %in% c(0,1)) {
                      reg <- rep(NA, 4)
                  } else {
                      if (model.type %in% c("linear", "logistic")) {
                          reg <- .runRegression(model.string, model.data, model.type)
                      } else if (model.type == "firth") {
                          reg <- .runFirth(model.string, model.data, geno.index)
                      }
                  }
                  c(freq, setNames(reg, .outputNames(model.type)))
              }, margin="by.variant", as.is="list")
              res <- do.call(rbind, res)
              data.frame(variant.id=seqGetData(gdsobj, "variant.id"), res,
                         stringsAsFactors=FALSE)
          })
