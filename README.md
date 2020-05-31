# The uniform general signed rank test

This R package implements the uniform general signed rank test described in

Howard, S.R. and Pimentel, S. D. (2019), [The uniform general signed rank test
and its design sensitivity](https://arxiv.org/abs/1904.08895), preprint,
arXiv:1904.08895.

You can install this package like so:

```R
install.packages('devtools')
devtools::install_github('gostevehoward/uniformrank')
```

Alternatively, you can find a package tarball under Github releases.

The main entry points are `uniform_max_Gamma()` and `uniform_pval()`. See the R
documentation for details.
