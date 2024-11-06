set.seed(5663)
mu <- c(-3, -2, -1, 0)

surveys <- 200

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

pairs_dat <- as.data.frame(pairs_dat)

test <- fit_rank_mod(pairs_dat)

summary(test)









