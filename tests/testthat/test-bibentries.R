test_that("bibentries are valid", {
  expect_s3_class(bibentries, "bibentry")
  expect_gte(length(bibentries), 5)

  # Minor formatting change around doi between 4.5 and 4.6, not particularly relevant
  if (R.version$major == "4" & as.numeric(R.version$minor) >= 6) {
    reference_print = "Ewald F, Bothmann L, Wright M, Bischl B, Casalicchio G, König G (2024).\n\\dQuote{A Guide to Feature Importance Methods for Scientific Inference.}\nIn Longo L, Lapuschkin S, Seifert C (eds.), \\emph{Explainable Artificial Intelligence}, 440--464.\nISBN 978-3-031-63797-1.\n\\doi{10.1007/978-3-031-63797-1_22}.\n\nKönig G, Molnar C, Bischl B, Grosse-Wentrup M (2021).\n\\dQuote{Relative Feature Importance.}\nIn \\emph{2020 25th International Conference on Pattern Recognition (ICPR)}, 9318--9325.\n\\doi{10.1109/ICPR48806.2021.9413090}."
  } else {
    reference_print = "Ewald F, Bothmann L, Wright M, Bischl B, Casalicchio G, König G (2024).\n\\dQuote{A Guide to Feature Importance Methods for Scientific Inference.}\nIn Longo L, Lapuschkin S, Seifert C (eds.), \\emph{Explainable Artificial Intelligence}, 440--464.\nISBN 978-3-031-63797-1, \\doi{10.1007/978-3-031-63797-1_22}.\n\nKönig G, Molnar C, Bischl B, Grosse-Wentrup M (2021).\n\\dQuote{Relative Feature Importance.}\nIn \\emph{2020 25th International Conference on Pattern Recognition (ICPR)}, 9318--9325.\n\\doi{10.1109/ICPR48806.2021.9413090}."
  }

  printed_bib = print_bib("ewald_2024", "konig_2021")
  expect_identical(printed_bib, reference_print)
})
