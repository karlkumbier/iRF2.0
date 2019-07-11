context('test intSubsetss')

test_that('intsSubsets works', {
  skip('This function is not called anywhere in the current code base.')
})

test_that('intSubsets works', {
  actual <- intSubsets('A+_B-_C+_D-')
  actual <- lapply(actual, sort)
  expect <- list("A+", "B-", "C+", "D-",
                 "A+_B-", "A+_C+", "A+_D-",
                 "B-_C+", "B-_D-", "C+_D-",
                 "A+_B-_C+", "A+_B-_D-",
                 "A+_C+_D-", "B-_C+_D-",
                 "A+_B-_C+_D-")
  expect <- lapply(expect, sort)
  expect_true(all(actual %in% expect) &&
              all(expect %in% actual))
})
