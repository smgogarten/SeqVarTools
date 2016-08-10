library(logistf)

test_lm <- function() {
    n <- 100
    dat <- data.frame(outcome=rnorm(n),
                      covarA=rnorm(n),
                      covarB=sample(letters[1:3], n, replace=TRUE),
                      genotype=sample(0:2, n, replace=TRUE))
    model.string <- "outcome ~ covarA + covarB + genotype"

    res <- SeqVarTools:::.runRegression(model.string, dat, "linear")
    exp <- lm(model.string, dat)
    checkEquals(summary(exp)$coef["genotype",1:2], res[1:2], check.names=FALSE)
}
    
test_glm <- function() {
    n <- 100
    dat <- data.frame(outcome=rbinom(n,1,0.3),
                      covarA=rnorm(n),
                      covarB=sample(letters[1:3], n, replace=TRUE),
                      genotype=sample(0:2, n, replace=TRUE))
    model.string <- "outcome ~ covarA + covarB + genotype"

    res <- SeqVarTools:::.runRegression(model.string, dat, "logistic")
    exp <- glm(model.string, dat, family="binomial")
    checkEquals(summary(exp)$coef["genotype",1:2], res[1:2], check.names=FALSE)
}

test_firth <- function() {
    n <- 100
    dat <- data.frame(outcome=rbinom(n,1,0.3),
                      covarA=rnorm(n),
                      covarB=sample(letters[1:3], n, replace=TRUE),
                      genotype=sample(0:2, n, replace=TRUE))
    model.string <- "outcome ~ covarA + covarB + genotype"

    res <- SeqVarTools:::.runFirth(model.string, dat, geno.index=5)
    exp <- logistf(as.formula(model.string), dat)
    checkEquals(coef(exp)["genotype"], res[1], check.names=FALSE)
}

.testData <- function(binary=FALSE, nv=100) {
    gds <- seqOpen(seqExampleFileName("gds"))
    seqSetFilter(gds, variant.id=1:nv, verbose=FALSE)
    
    require(Biobase)
    sample.id <- seqGetData(gds, "sample.id")
    n <- length(sample.id)
    df <- data.frame(sample.id=seqGetData(gds, "sample.id"),
                     outcome=(if (binary) rbinom(n,1,0.3) else rnorm(n)),
                     covarA=rnorm(n),
                     covarB=factor(sample(letters[1:3], n, replace=TRUE)),
                     stringsAsFactors=FALSE)
    adf <- AnnotatedDataFrame(df)

    SeqVarData(gds, adf)
}

test_regression <- function() {
    gds <- .testData()
    res <- regression(gds, outcome="outcome", covar=c("covarA", "covarB"))

    dat <- cbind(pData(sampleData(gds)), genotype=refDosage(gds)[,1])
    exp <- lm("outcome ~ covarA + covarB + genotype", dat)
    seqClose(gds)
    checkEquals(summary(exp)$coef["genotype",1], res[1,"Est"], check.names=FALSE)
    checkEquals(summary(exp)$coef["genotype",2], res[1,"SE"], check.names=FALSE)
}

test_regression_firth <- function() {
    gds <- .testData(binary=TRUE)
    res <- regression(gds, outcome="outcome", covar=c("covarA", "covarB"),
                      model.type="firth")

    dat <- cbind(pData(sampleData(gds)), genotype=refDosage(gds)[,1])
    seqClose(gds)
    exp <- logistf(as.formula("outcome ~ covarA + covarB + genotype"), dat)
    checkEquals(coef(exp)["genotype"], res[1,"Est"], check.names=FALSE)
}

test_freq <- function() {
    gds <- .testData()
    geno <- refDosage(gds)
    freq <- alleleFrequency(gds)
    for (i in 1:ncol(geno)) {
        model.data <- cbind(pData(sampleData(gds)), genotype=geno[,i])
        fo <- SeqVarTools:::.freqByOutcome(model.data, "linear", "outcome")
        checkEquals(freq[i], fo["freq"], check.names=FALSE)
        checkEquals(sum(!is.na(geno[,i])), fo["n"], check.names=FALSE)
    }
    seqClose(gds)
}

test_freq_binary <- function() {
    gds <- .testData(binary=TRUE)
    status <- sampleData(gds)$outcome
    geno <- refDosage(gds)
    n <- lapply(c(0,1), function(c) colSums(!is.na(geno[(status == c),])))
    samp.id <- sampleData(gds)$sample.id
    freq <- lapply(c(0,1), function(c)
                   applyMethod(gds, alleleFrequency, sample=samp.id[status == c]))
    for (i in 1:ncol(geno)) {
        model.data <- cbind(pData(sampleData(gds)), genotype=geno[,i])
        fo <- SeqVarTools:::.freqByOutcome(model.data, "logistic", "outcome")
        for (c in c(0,1)) {
            checkEquals(freq[[c+1]][i], fo[paste0("freq", c)], check.names=FALSE)
            checkEquals(n[[c+1]][i], fo[paste0("n", c)], check.names=FALSE)
        }
    }
    seqClose(gds)
}

mono <- function(x) {
    (x["freq0"] == 0 & x["freq1"] == 0) | (x["freq0"] == 1 & x["freq1"] == 1)
}

test_freqFromBinary <- function() {
    gds <- .testData(binary=TRUE)
    geno <- refDosage(gds)
    res <- regression(gds, outcome="outcome", covar=c("covarA", "covarB"),
                      model.type="logistic")
    for (i in 1:ncol(geno)) {
        model.data <- cbind(pData(sampleData(gds)), genotype=geno[,i])
        fo <- SeqVarTools:::.freqByOutcome(model.data, "logistic", "outcome")
        fc <- SeqVarTools:::.freqByOutcome(model.data, "linear", "outcome")
        checkEquals(fc["freq"], SeqVarTools:::.freqFromBinary(fo), checkNames=FALSE)
        checkEquals(mono(fo), fc["freq"] %in% c(0,1), checkNames=FALSE)
    }
    seqClose(gds)
}

test_freq_mono <- function() {
    gds <- .testData(binary=FALSE)
    res <- regression(gds, outcome="outcome", covar=c("covarA", "covarB"),
                      model.type="linear")
    for (i in 1:nrow(res)) {
        if (res[i,"freq"] %in% c(0,1)) checkTrue(is.na(res[i,"Est"]))
    }
    seqClose(gds)
    
    gds <- .testData(binary=TRUE)
    res <- regression(gds, outcome="outcome", covar=c("covarA", "covarB"),
                      model.type="logistic")
    for (i in 1:nrow(res)) {
        if (mono(res[i,])) checkTrue(is.na(res[i,"Est"]))
    }
    seqClose(gds)
}

test_firth_badindex <- function() {
    n <- 100
    dat <- data.frame(outcome=rbinom(n,1,0.3),
                      covarA=rnorm(n),
                      covarB=sample(letters[1:3], n, replace=TRUE),
                      genotype=sample(0:2, n, replace=TRUE))
    model.string <- "outcome ~ covarA + covarB + genotype"

    exp <- logistf(as.formula(model.string), dat)
    res <- SeqVarTools:::.runFirth(model.string, dat, geno.index=4) # too small
    checkEquals(coef(exp)["genotype"], res[1], check.names=FALSE)
    res <- SeqVarTools:::.runFirth(model.string, dat, geno.index=6) # too big
    checkEquals(coef(exp)["genotype"], res[1], check.names=FALSE)
}

test_droplevels <- function() {
    gds <- .testData(binary=TRUE)
    seqSetFilter(gds, sample.sel=sampleData(gds)$covarB != "c")
    res <- regression(gds, outcome="outcome", covar=c("covarA", "covarB"),
                      model.type="firth")

    dat <- droplevels(cbind(pData(sampleData(gds)), genotype=refDosage(gds)[,1]))
    seqClose(gds)
    exp <- logistf(as.formula("outcome ~ covarA + covarB + genotype"), dat)
    checkEquals(coef(exp)["genotype"], res[1,"Est"], check.names=FALSE)
}

test_nocovar <- function() {
    gds <- .testData(binary=FALSE)
    res <- regression(gds, outcome="outcome", covar=NULL, model.type="linear")
    
    dat <- cbind(pData(sampleData(gds)), genotype=refDosage(gds)[,1])
    exp <- lm("outcome ~ genotype", dat)
    seqClose(gds)
    checkEquals(summary(exp)$coef["genotype",1], res[1,"Est"], check.names=FALSE)
    checkEquals(summary(exp)$coef["genotype",2], res[1,"SE"], check.names=FALSE)
}
