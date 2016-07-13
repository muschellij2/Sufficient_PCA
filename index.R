## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
cacher = FALSE

## ----data_prcomp, cache=cacher-------------------------------------------
set.seed(20160706)
p = 81
R = 4e5
Y = matrix(rnorm(p * R, mean = 6, sd = 4), nrow = R, ncol = p)
Z = scale(Y, center = TRUE, scale = TRUE)

## ----prcomp, cache=cacher, dependson="data_prcomp"-----------------------
system.time({
    pc_comp = prcomp(
      x = Z, 
      center = FALSE,
      scale. = FALSE)
    pcs = pc_comp$x
})

## ----svd_z, cache=cacher, dependson="data_prcomp"------------------------
system.time({
    svd_z = svd(Z, nu = 0)
    pc_svd_z = Z %*% svd_z$v
})

## ----equiv, cache=cacher, dependson="data_prcomp"------------------------
identical(dim(pc_svd_z), dim(pcs))
corr_pcs = diag(cor(pcs, pc_svd_z))
corr_pcs = abs(corr_pcs) # sign variance
all(abs(corr_pcs - 1) < 1e-10)

## ----pc_cov_z, cache=cacher, dependson="data_prcomp"---------------------
system.time({
    cov = t(Z) %*% Z
    svd_cov = svd(cov, nv = 0)
    pc_cov_z = Z %*% svd_cov$u 
})

## ----equiv_cov, cache=cacher---------------------------------------------
identical(dim(pc_svd_z), dim(pc_cov_z))
corr_pcs = diag(cor(pc_cov_z, pc_svd_z))
corr_pcs = abs(corr_pcs) # sign variance
all(abs(corr_pcs - 1) < 1e-10)

## ----left_is_right, cache=cacher, dependson="data_prcomp"----------------
svd_cov_all = svd(cov)
max(abs(svd_cov_all$u - svd_cov_all$v))

## ----funcs---------------------------------------------------------------
pca_stats = function(L){
    N = length(L)
    tmp_l = vector(mode = "list",
        length = N)
    suff_L = list(
        cross_i = tmp_l,
        csums_i = tmp_l,
        csums_sq_i = tmp_l,
        n_i = tmp_l
        )
    iid = 1
    for (iid in seq(N)) {
        mat = L[[iid]]
        mat = as.matrix(mat)
        cross_i = crossprod(mat)
        csums_i = colSums(mat)
        csums_sq_i = colSums(mat ^ 2)
        n_i = nrow(mat)
        suff_L$cross_i[[iid]] = cross_i 
        suff_L$csums_i[[iid]] = csums_i 
        suff_L$csums_sq_i[[iid]] = csums_sq_i 
        suff_L$n_i[[iid]] = n_i 
    }
    return(suff_L)
}

## ----split, cache=cacher, dependson="data_prcomp"------------------------
n_ids = 5
ids = sample(1:n_ids, size = R, replace = TRUE)
X_i = lapply(1:n_ids, function(id) {
  Y[ids == id,]
})
stats = pca_stats(X_i)

## ----combiner------------------------------------------------------------
make_cov = function(suff_L){
    # stopifnot(fslr::same_dims(suff_L$cross_i))
    # mu = mapply(function(cs, n){
    #     cs/n
    # }, suff_L$csums_i, suff_L$n_i,
    #     SIMPLIFY = FALSE)
    s_i = Reduce("+", suff_L$csums_i)
    R = Reduce("+", suff_L$n_i)
    mu = s_i / R
    s2_i = Reduce("+", suff_L$csums_sq_i)
    sigma = sqrt((s2_i - R* mu ^ 2)/(R - 1))
    sigma_mat = outer(sigma, sigma)

    mu_mat = outer(mu, mu)

    C = Reduce("+", suff_L$cross_i)
    C = (C - R *mu_mat)/sigma_mat
    L = list(sigma = sigma,
             mu = mu,
             C =  C)
    return(L)
}

## ----combine, cache=cacher, dependson = "split"--------------------------
res = make_cov(stats)
split_cov = res$C
max(abs(cov - split_cov)) < 1e-7

## ----pc_cov_split, cache=cacher, dependson="data_prcomp"-----------------
system.time({
    svd_split_cov = svd(split_cov, nv = 0)
})
max(abs(svd_split_cov$u - svd_cov$u )) < 1e-7

## ----equiv_cov_split, cache=cacher---------------------------------------
pc_split = Z %*% svd_split_cov$u
identical(dim(pc_split), dim(pc_cov_z))
corr_pcs = diag(cor(pc_cov_z, pc_split))
corr_pcs = abs(corr_pcs) # sign variance
all(abs(corr_pcs - 1) < 1e-10)

