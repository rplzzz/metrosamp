context('analysis')

mk_test_ms <- function(vals, samplp, paccept=0, scale=numeric(0)) {
    if(is.matrix(vals)) {
        samples <- vals
    }
    else {
        samples <- matrix(1, nrow=length(vals), ncol=3) * vals
    }
    if(length(scale)==0) {
        scale <- rep(1, ncol(samples))
    }
    structure(
        list(samples=samples, samplp=samplp, accept=paccept, plast=samples[nrow(samples),], scale=scale),
        class=c('metrosamp','list'))
}

ms1 <- mk_test_ms(seq(1,25), seq(-1,-25))
ms2 <- mk_test_ms(seq(26,50), seq(0,-24))
ms3 <- mk_test_ms(seq(51,75), seq(-1,-25))
ms4 <- mk_test_ms(seq(76,100), seq(-1,-25))
pn <- c('a','b','c')
colnames(ms1$samples) <- colnames(ms2$samples) <- colnames(ms3$samples) <- colnames(ms4$samples) <- pn
msruns <- list(ms1, ms2, ms3, ms4)

test_that('MAP returns correct value',{
    mx <- MAP(msruns)
    expect_equal(mx, structure(c(26, 26, 26), names=pn))
})

test_that('EV returns correct value', {
    mx <- EV(msruns)
    expect_equal(mx, structure(c(50.5, 50.5, 50.5), names=pn))
})

test_that('CI returns correct values', {
    mx <- CI(msruns)
    mx_expected <- matrix(c(rep(3.475, 3), rep(97.525, 3)), ncol=2)
    rownames(mx_expected) <- pn
    colnames(mx_expected) <- c('2.5%', '97.5%')
    expect_equal(mx, mx_expected)
})

test_that('rsample returns an appropriately constructed matrix', {
    samps <- rsample(msruns, 5)
    expect_is(samps, 'matrix')
    expect_equal(dim(samps), c(5,3))
    expect_equal(colnames(samps), pn)
})
