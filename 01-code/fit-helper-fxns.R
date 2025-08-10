
# Extract lm coefficients cleanly 
get_lm_coefs <- function(x) {
  mdl <- x %>% 
    # get the lm model object
    extract_fit_engine()
  
  prms <- tidy(mdl)
  
  pp <- x %>% 
    # get the lm model object
    extract_preprocessor() |> 
    summary() |> filter(role == "predictor") |> 
    unnest_wider(type,names_sep = "_") 
  
  emmeans(mdl,'genotype')
  estimate_expectation(mdl)
}

# A helper function that closes any parallel workers
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# Load commonly used and additional libraries
load_libs <- function(...){
  libs <-
    c(
      c(
        "tidyverse",
        "gt",
        "cowplot",
        "patchwork",
        "ggplot2",
        "ggh4x",
        "latex2exp",
        "readxl",
        "tidymodels",
        "furrr",
        "doParallel",
        "doFuture",
        "broom.mixed"
      ),
      ...
    )
  lapply(libs,require, character = TRUE)
  select <- dplyr::select
  map <- purrr::map
  tidymodels_prefer()
  
}

# Fix eggs tag facet
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}
