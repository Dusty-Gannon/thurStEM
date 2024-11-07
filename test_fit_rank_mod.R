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

pairs_dat2$grp2 <- ifelse(pairs_dat2$grp == 1, 0, 1)
pairs_dat2$grp3 <- rep(1:2, each = surveys/2)
pairs_dat3 <- pairs_dat2
pairs_dat3$grp <- factor(pairs_dat3$grp)


# i1i2 ~ d12 * 1 + g12 * grp2
# i1i3 ~ d13 * 1 + g13 * grp2
# i1i4 ~ d14 * 1 + g14 * grp2
# i2i3 ~ d23 * 1 + g23 * grp2
# i2i4 ~ d24 * 1 + g24 * grp2
# i3i4 ~ d34 * 1 + g34 * grp2

reg2 <- '
  i1 ~ mu1 * 1 + g11 * grp2
  i2 ~ mu2 * 1 + g21 * grp2
  i3 ~ mu3 * 1 + g31 * grp2
  i4 ~ 0 * 1 + mu4 * 1
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
  mu11 := sqrt(2) * mu1
  mu21 := sqrt(2) * mu2
  mu31 := sqrt(2) * mu3
  mu12 := sqrt(2) * mu1 + g11
  mu22 := sqrt(2) * mu2 + g21
  mu32 := sqrt(2) * mu3 + g31
'

# fit the model
mfit2 <- lavaan(
  model = c(reg2, meas2, covars2, derived2),
  data = pairs_dat2,
  ordered = names(pairs_dat2)[1:6],
  parameterization = "theta",
  meanstructure = TRUE,
  orthogonal = T,
  std.lv = T
)
lavaan::summary(mfit2, standardized = F)

colMeans(pairs_dat2[1:surveys/2, 1:6])

pnorm(coef(mfit2))


colMeans(pairs_dat2[(surveys/2 + 1):surveys, 1:6])
pnorm(coef(mfit2)[which(1:12 %% 2 == 0)] +
        coef(mfit2)[which(1:12 %% 2 == 1)])


coef(mfit2)


