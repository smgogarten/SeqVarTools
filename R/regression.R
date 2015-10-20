
.waldTest <- function(Est, cov) {
    if (length(Est) == 1) {
        W <- (Est^2)/cov
    } else {
        W <- as.numeric(t(Est) %*% solve(cov) %*% Est)
    }
    pval <- pchisq(W, df=length(Est), lower.tail=FALSE)
    c(Wald.Stat=W, Wald.pval=pval)
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
        ret <- c(Est=Est, SE=sqrt(cov), .waldTest(Est, cov))

    }, warning=function(w) NA, error=function(e) NA)
}

setMethod("regression",
          "SeqVarData",
          function(gdsobj, outcome, covar=NULL,
                   model.type=c("linear", "logistic")) {
              model.type <- match.arg(model.type)
              
              ## get covariates
              dat <- pData(sampleData(gdsobj))[,c(outcome, covar)]
              
              ## create model formula
              model.string <- paste(outcome, "~", paste(c(covar, "genotype"), collapse=" + "))
              
              ## apply function over variants
              res <- seqApply(gdsobj, "genotype", function(x) {
                  ## assume we want effect of reference allele
                  model.data <- cbind(dat, genotype=colSums(x == 0))
                  .runRegression(model.string, model.data, model.type)
              }, margin="by.variant", as.is="list")
              res <- do.call(rbind, res)
              data.frame(variant.id=seqGetData(gdsobj, "variant.id"), res,
                         stringsAsFactors=FALSE)
          })
