# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

lhoodcpp <- function(eta, beta, doc_ct, mu, siginv) {
  .Call("stm_lhoodcpp", eta, beta, doc_ct, mu, siginv)
}

gradcpp <- function(eta, beta, doc_ct, mu, siginv) {
  .Call("stm_gradcpp", eta, beta, doc_ct, mu, siginv)
}

hpbcpp <- function(eta, beta, doc_ct, mu, siginv, sigmaentropy) {
  .Call("stm_hpbcpp", eta, beta, doc_ct, mu, siginv, sigmaentropy)
}

