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

.testData <- function() {
    gds <- seqOpen(seqExampleFileName("gds"))
    seqSetFilter(gds, variant.id=1:100, verbose=FALSE)
    
    require(Biobase)
    sample.id <- seqGetData(gds, "sample.id")
    n <- length(sample.id)
    df <- data.frame(sample.id=seqGetData(gds, "sample.id"),
                     outcome=rnorm(n),
                     covarA=rnorm(n),
                     covarB=sample(letters[1:3], n, replace=TRUE),
                     stringsAsFactors=FALSE)
    adf <- AnnotatedDataFrame(df)

    SeqVarData(gds, adf)
}

test_regression <- function() {
    gds <- .testData()
    res <- regression(gds, outcome="outcome", covar=c("covarA", "covarB"))

    dat <- cbind(pData(sampleData(gds)), genotype=refDosage(gds)[,1])
    exp <- lm("outcome ~ covarA + covarB + genotype", dat)
    checkEquals(summary(exp)$coef["genotype",1], res[1,"Est"], check.names=FALSE)
    checkEquals(summary(exp)$coef["genotype",2], res[1,"SE"], check.names=FALSE)
    seqClose(gds)
}
