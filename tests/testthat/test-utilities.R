varnames.grp <- c('SL', 'SW', 'PL', 'PW')
varnames <- unique(varnames.grp)
ints <- data.frame(c(3, 8), c(7, 4), c(3, 4))

test_that('groupVars works', {
  expect_equal(groupVars(varnames.grp, iris[, 1:4]),
               varnames.grp)
})

test_that('pasteInt works', {
  int <- c('Sepal', 'Petal')
  expect_equal(pasteInt(int),
               'Sepal_Petal')
})

test_that('nameInts works', {
  expected <- c(c.3..8. = 'PL-_PW+',
                c.7..4. = 'PL+_PW-',
                c.3..4. = 'PL-_PW-')
  expect_equal(nameInts(ints, varnames, TRUE),
               expected)
})

test_that('nameInt works', {
  idx <- c(3, 4)
  sgn <- c('+', '-')
  expect_equal(nameInt(varnames, idx, sgn),
               'PL+_PW-')
})

test_that('int2Id works', {
  expect_equal(int2Id('PL+_PW-', varnames.grp,
                      signed=TRUE, split=FALSE),
               c(PL=7, PW=4))
})

test_that('unsign works', {
  expect_equal(unsign('PL+_PW-'),
               'PL_PW')
})

test_that('intSign works', {
  expect_equal(intSign('PL+_PW-', split=FALSE),
               c(1, -1))
})

test_that('intSubsets works', {
  actual <- intSubsets(c(2L, 3L, 5L, 7L))
  actual <- lapply(actual, sort)
  expect <- list(2L, 3L, 5L, 7L,
                 c(2L, 3L, 5L),
                 c(2L, 3L, 7L),
                 c(2L, 5L, 7L),
                 c(3L, 5L, 7L),
                 c(2L, 3L, 5L, 7L))
  expect <- lapply(expect, sort)
  expect_true(all(actual %in% expect) &&
              all(expect %in% actual))
})

test_that('lreplicate works', {
  expect_equal(lreplicate(5, exp(0)),
               list(1, 1, 1, 1, 1))
})
