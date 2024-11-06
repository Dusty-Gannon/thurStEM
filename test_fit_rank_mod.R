set.seed(5663)
mu <- c(-3, -2, -1, 0)
devtools::load_all()

surveys <- 1000

ranks <- sapply(
  1:surveys,
  function(x){
    s <- rnorm(length(mu), mean = mu)
    return(order(s, decreasing = T))
  }
) |> t()

# give names
colnames(ranks) <- paste0("i", 1:4)

pairs_dat <- ranks_as_pairs(ranks)

test <- fit_rank_mod(pairs_dat)

summary(test)

# Intercept only model looks good

#### Model with two groups ####

mu1 <- c(-3, -2, -1, 0)
mu2 <- c(-3, -2.5, -0.5, 0)
surveys <- 10000

ranks2 <- matrix(nrow = surveys, ncol = length(mu1))

for(n in 1:surveys){
  if(n <= surveys/2){
    ranks2[n, ] <- rnorm(length(mu1), mu1) |> order(decreasing = T)
  } else{
    ranks2[n, ] <- rnorm(length(mu2), mu2) |> order(decreasing = T)
  }
}

colnames(ranks2) <- paste0("i", 1:4)

pairs_dat2 <- ranks_as_pairs(ranks2)

# add grouping factor
pairs_dat2$grp <- factor(rep(1:2, each = surveys/2))

test2 <- fit_rank_mod(pairs_dat2, pairs_cols = 1:6, form = ~ grp)
summary(test2)

pairs_dat2$grp <- ifelse(pairs_dat2$grp == 1, 0, 1)
pairs_dat3 <- pairs_dat2
pairs_dat3$grp <- factor(pairs_dat3$grp)

reg2 <- '
  i4 ~ 0 * 1 + mu4 * 1
  i1i2 ~ d12 * 1 + g12 * grp
  i1i3 ~ d13 * 1 + g13 * grp
  i1i4 ~ d14 * 1 + g14 * grp
  i2i3 ~ d23 * 1 + g23 * grp
  i2i4 ~ d24 * 1 + g24 * grp
  i3i4 ~ d34 * 1 + g34 * grp
'

# fix the factor loadings as pairwise differences
meas2 <- '
  i1 =~ 1 * i1i2 + 1 * i1i3 + 1 * i1i4
  i2 =~ -1 * i1i2 + 1 * i2i3 + 1 * i2i4
  i3 =~ -1 * i1i3 + -1 * i2i3 + 1 * i3i4
  i4 =~ -1 * i1i4 + -1 * i2i4 + -1 * i3i4
'

covars2 <- '
  i1 ~~ 1 * i1
  i2 ~~ 1 * i2
  i3 ~~ 1 * i3
  i4 ~~ 1 * i4
  i1i2 ~~ 2 * i1i2
  i1i3 ~~ 2 * i1i3
  i1i4 ~~ 2 * i1i4
  i2i3 ~~ 2 * i2i3
  i2i4 ~~ 2 * i2i4
  i3i4 ~~ 2 * i3i4
'

derived2 <- '
  mu11 := sqrt(2) * d14
  mu21 := sqrt(2) * d24
  mu31 := sqrt(2) * d34
  b1 := sqrt(2) * g14
  b2 := sqrt(2) * g24
  b3 := sqrt(2) * g34
  mu12 := sqrt(2) * (d14 + g14)
  mu22 := sqrt(2) * (d24 + g24)
  mu32 := sqrt(2) * (d34 + g34)
'

# fit the model
mfit2 <- lavaan(
  model = c(reg2, meas2, covars2, derived2),
  data = pairs_dat3,
  ordered = names(pairs_dat2)[1:6],
  parameterization = "theta",
  meanstructure = TRUE,
  orthogonal = T,
  std.lv = T
)
lavaan::summary(mfit2)



