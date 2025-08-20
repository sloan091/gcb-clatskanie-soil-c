# Harmonizing Clatskanie Datasets for GCB Publication
Brandon Sloan
19 August, 2025

# Introduction

We harmonized several different data sets collected at the Clatskanie
Common Garden Site over the past few years for this *Global Change
Biology* manuscript. We will briefly describe soil, root and aboveground
biomass data below and discuss the pre-processing that occurred to
arrive at the final data set used for our heritability and regression
analysis in the manuscript. The final data set is located in this repo
at */02-data/02-processed/clatskanie-c-fit-data.csv.* The manuscript
analysis using this data can be found in this folder in the
*02-gcb-pub-analysis.qmd* notebook.

# SEED Soil Core Data

Soil cores were collected in 2022 as part of and ORNL SEED project led
by Melanie Mayes. The cores were analyzed for mineral associated,
particulate, and total C and N, organic acids, pH, and soil textural
data. The cores were taken for 23 genotypes with 3 replicates for most
genotypes (66 cores total), with an additional genotype taken as part of
a later project (BESC-13). The cores were taken at 15 and 30 cm, with a
few cores having deeper measurements. As part of the harmonization, I
have identified each tree based on its genotype label and replicate
block, which are represented by the *tree, genotype,* and *block*
columns in all of the following tables.

<details class="code-fold">
<summary>Code</summary>

``` r
# Clatskanie latitude and longitude data from Stan Martins
df_ll <-
  read_csv(
    file.path(d_path,"03-cbi/clatskanie-lat-lon.csv")) |> 
  mutate(genotype = as.factor(CBI_genotype),
         block = as.factor(CBI_block),
         tree = fct_cross(genotype,block,sep = "-")
         ) |> 
    select(tree,latitude,longitude) |> 
    group_by(tree) |> 
    summarize(across(.cols = c(latitude,longitude),first))

# SEED soil core data
df_core <-
  read_xlsx(file.path(d_path, "01-seed/seed-soil-core-update-besc-13.xlsx"),
            sheet = "Soils") |>
  mutate(
    genotype = as.factor(Genotype),
    block = as.factor(Block),
    tree = fct_cross(genotype, block, sep = "-"),
    across(c(TotpN, TotpC), as.numeric)
  ) |>
  left_join(df_ll, by = "tree") |>
  select(tree,
         genotype,
         block,
         Depth,
         latitude,
         longitude,
         Row,
         Column,
         Trait:Citrate) |>
  mutate(
    TotC = TotpC * 10,
    # convert to units tonnes mg C/g soil
    C_chk = abs((MAOMC + POMC - TotC) / TotC * 100),
    C_chk_rel = (MAOMC + POMC - TotC) / TotC * 100,
    C_chk_abs = (MAOMC + POMC - TotC)
  )

display_table <- function(df){
 head(df) |> 
    knitr::kable(digits = 2)
}

display_table(df_core)
```

</details>
<div id="tbl-seed-core">

<div class="cell-output-display">

| tree | genotype | block | Depth | latitude | longitude | Row | Column | Trait | BD | MCwetBD | MCdryBD | MCwetCHEM | MCdryCHEM | pH | Sand | Clay | Silt | POMpN | POMpC | MAOMpN | MAOMpC | POMN | POMC | MAOMN | MAOMC | TotpN | TotpC | OrgAcidMass | Lactate_conc | Acetate_conc | Propionate_conc | Formate_conc | Butyrate_conc | Pyruvate_conc | Succinate_conc | Oxalate_conc | Furmarate_conc | Citrate_conc | Lactate | Acetate | Propionate | Formate | Butyrate | Pyruvate | Succinate | Oxalate | Furmarate | Citrate | TotC | C_chk | C_chk_rel | C_chk_abs |
|:---|:---|:---|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|
| BESC-102-1 | BESC-102 | 1 | 15 | 46.12 | -123.27 | 18 | 8 | HF | 0.88 | 0.37 | 0.60 | 0.42 | 0.72 | 5.11 | 0.00 | 45.90 | 55.75 | 0.89 | 19.61 | 0.31 | 4.05 | 0.53 | 11.56 | 3.05 | 39.52 | 0.41 | 6.04 | 1.21 | 19.89 | 9.44 | 14.05 | 9.24 | 0 | 0 | 0 | 0 | 4.33 | 0 | 26.241100053454279 | 12.447701307407799 | 18.542919117712412 | 12.190435620609227 | 0 | 0 | 0 | 0 | 5.7165754917856901 | 0 | 60.4 | 15.43 | -15.43 | -9.32 |
| BESC-102-1 | BESC-102 | 1 | 30 | 46.12 | -123.27 | 18 | 8 | HF | 0.61 | 0.32 | 0.47 | 0.31 | 0.44 | 5.41 | 0.00 | 50.45 | 50.45 | NA | NA | NA | NA | NA | NA | NA | NA | 0.16 | 2.09 | 1.24 | 26.03 | 0.00 | 0.00 | 0.00 | 0 | 0 | 0 | 0 | 4.24 | 0 | 30.761400137078475 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 5.0045918167003709 | 0 | 20.9 | NA | NA | NA |
| BESC-102-2 | BESC-102 | 2 | 15 | 46.12 | -123.27 | 61 | 21 | HF | 0.62 | 0.41 | 0.69 | 0.43 | 0.77 | 5.09 | 7.85 | 45.15 | 47.00 | 0.81 | 16.53 | 0.30 | 3.40 | 0.40 | 8.11 | 2.83 | 32.30 | 0.53 | 8.15 | 1.29 | 11.48 | 12.30 | 14.13 | 8.95 | 0 | 0 | 0 | 0 | 4.27 | 0 | 15.06244404695275 | 16.14911335695335 | 18.540310799973991 | 11.747315209921938 | 0 | 0 | 0 | 0 | 5.6000216736383539 | 0 | 81.5 | 50.42 | -50.42 | -41.09 |
| BESC-102-2 | BESC-102 | 2 | 30 | 46.12 | -123.27 | 61 | 21 | HF | 0.62 | 0.45 | 0.82 | 0.41 | 0.70 | 6.15 | 0.35 | 52.15 | 47.50 | NA | NA | NA | NA | NA | NA | NA | NA | 0.44 | 6.54 | 1.28 | 8.91 | 0.00 | 0.00 | 0.00 | 0 | 0 | 0 | 0 | 4.30 | 0 | 12.678931857096265 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6.1161551475784091 | 0 | 65.4 | NA | NA | NA |
| BESC-102-3 | BESC-102 | 3 | 15 | 46.12 | -123.27 | 92 | 22 | HF | 0.79 | 0.32 | 0.48 | 0.35 | 0.54 | 5.70 | 8.39 | 38.21 | 51.20 | 0.64 | 13.30 | 0.40 | 5.64 | 0.37 | 7.59 | 3.68 | 51.99 | 0.38 | 5.35 | 1.08 | 9.20 | 0.00 | 0.00 | 0.00 | 0 | 0 | 0 | 0 | 0.00 | 0 | 12.609040922118737 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 53.5 | 11.36 | 11.36 | 6.08 |
| BESC-102-3 | BESC-102 | 3 | 30 | 46.12 | -123.27 | 92 | 22 | HF | 0.83 | 0.33 | 0.50 | 0.35 | 0.54 | 5.34 | 0.00 | 45.90 | 55.00 | NA | NA | NA | NA | NA | NA | NA | NA | 0.29 | 3.27 | 1.01 | 0.00 | 0.00 | 0.00 | 0.00 | 0 | 0 | 0 | 0 | 0.00 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 32.7 | NA | NA | NA |

</div>

Table 1: SEED soil core data collected in 2022.

</div>
<div id="tbl-seed-core-key">

| Variable | Description | Units |
|----|----|----|
| tree | Unique ID for all poplar trees across these data sets (69 trees) |  |
| genotype | ID for genotypes contained across these data sets (24 genotypes) |  |
| block | ID for replicate blocks (3 replicates) |  |
| Depth | Depth that soil core was taken to at each tree (max of 4 depths) | cm |
| Traits | The aboveground phenotypic traits used to select the soil cores |  |
| BD | Soil bulk density | g/cm$^3$ |
| MCwetBD, MCdryBD, MCwetCHEM, MCwetBD | Wet and dry soil moisture content from bulk density and chemistry samples | cm$^3$ water/ cm$^3$ soil |
| pH | Soil acidity |  |
| Sand | Percent Sand | % |
| Clay | Percent Clay | % |
| Silt | Percent Silt | % |
| POMpN, POMpC | Particulate organic matter percent Carbon and Nitrogen | % |
| POMN, POMC | Particulate organic matter percent Carbon and Nitrogen | mg C or N/g soil |
| MAOMpN, MAOMpC | Mineral associated organic matter percent C and N | % |
| MAOMN, MAOMC | Mineral associated organic matter C and N | mg C or N/g soil |
| TotpN, TotpC | Total percent N or C | % |
| OrgAcidMass | Fresh soil mass for organic acid analysis | g |
| Lactate_conc, Acetate_conc, … | Organic acid concentrations | $\mu$moles |
| Lactate, Acetate, … | Specific organic acid concentrations | $\mu$moles/g soil |

Table 2: <a href="#tbl-seed-core" class="quarto-xref">Table 1</a> key

</div>

## C Mass Balance Discrepancy

<details class="code-fold">
<summary>Code</summary>

``` r
# Combined data set for prelim analysis and plotting
df_check <-
  df_core |> filter(Depth == 15) |>
  select(tree:TotpC) |>
  mutate(
    TotC = TotpC * 10, # convert to units tonnes mg C/g soil
    C_chk = abs((MAOMC + POMC - TotC) / TotC * 100),
    C_chk_rel = (MAOMC + POMC - TotC) / TotC * 100,
    C_chk_abs = (MAOMC + POMC - TotC),
    C_chk_id = fct(ifelse(C_chk > 30, "Outlier", "Good")),
    C_chk_plt = fct(ifelse(C_chk > 30, as.character(tree), ""))
  )

outliers <-
  df_check |>
  filter(C_chk_id %in% "Outlier") |>
  pluck("tree") |>
  as.character()
```

</details>

I initially checked the soil C measurements to ensure that the MAOM and
POM soil C fractions matched the overall total C, since they were
measured from separate sub-samples. We considered mass balance errors
greater than 30% as unacceptable.
<a href="#fig-c-chk" class="quarto-xref">Figure 1</a> shows the seven
trees that have unrealistic mass balance (BESC-102-2, BESC-876-2,
BESC-144-3, BESC-192-1, BESC-192-2, BESC-192-3, BESC-145-1). Therefore,
the elemental analysis (and not fractionation) were re-run for these
outliers by Tommy Mead with replicates to characterize the experimental
error.

<details class="code-fold">
<summary>Code</summary>

``` r
# Sanity check on MAOM and POM calcs
g1 <- df_check |>
  ggplot(aes(
    x = MAOMC + POMC,
    y = TotC,
    fill = Clay,
    shape = C_chk_id,
    label = C_chk_plt
  )) +
  geom_point(size = 3) +
  geom_text_repel() +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = 2) +
  scale_shape_manual(name = NULL,
                     values = c("Outlier" = 24, "Good" = 21)) +
  scale_fill_distiller(
    name = "% Clay",
    palette = "RdYlBu",
    type = "div",
    limits = c(30, 50),
    oob = scales::squish
  ) +
  xlab("0-15 cm MAOM C + POM C \n [mg C/g soil]") +
  ylab("0-15 cm Total C \n [mg C/g soil]") +
  theme_cowplot(fs)
plot(g1)
```

</details>
<div id="fig-c-chk">

![](01-harmonize-clatskanie-data-pub_files/figure-commonmark/fig-c-chk-1.png)


Figure 1: Sanity check on soil C mass balance showing the sums of
Mineral Associated Soil C and Particulate Soil C observations versus the
Total Soil C..

</div>

The re-run samples are compared to the original values and checked
against the mass balance

<details class="code-fold">
<summary>Code</summary>

``` r
# Read Original Fractionation data
df_frac <-
  file.path(d_path, "01-seed/seed-maom-pom-fractionation.xlsx") |>
  read_xlsx(sheet = "POM_MAOM_raw_and_calcs",
            skip = 1) |>
  select(`Full ID`, `Soil dry weight (g)`, `POM mass`, `MAOM mass`) |>
  mutate(
    genotype = ifelse(
      str_detect(`Full ID`, "ROAD"),
      "ROAD",
      str_extract(`Full ID`, "[:upper:]+(\\-|\\s|\\S)\\d+\\-?\\d*")
    ),
    block = ifelse(
      str_detect(`Full ID`, "ROAD"),
      1,
      str_extract(`Full ID`, "(?<=(Cl|CL|U))\\d+")
    ),
    tree = paste(genotype, block, sep = "-"),
    across(genotype:tree, fct),
    pom_frac = `POM mass` / `Soil dry weight (g)`,
    maom_frac = `MAOM mass` / `Soil dry weight (g)`
  ) |>
  drop_na()
```

</details>
<details class="code-fold">
<summary>Code</summary>

``` r
# Read in updated samples from Tommy Mead by 01/02/24
df_core_update_1 <-
  list.files(file.path(d_path,
                       "02-ldrd/"), "*Tray*", full.names = TRUE) |>
  map(
    \(x)
    read_xlsx(x,
              sheet = "Output") |>
      filter(Type %in% "Test") |>
      select(`Sample ID`:`Final %C`)
  ) |>
  bind_rows() |>
  mutate(
    genotype = ifelse(
      str_detect(`Sample ID`, "ROAD"),
      "ROAD",
      str_extract(`Sample ID`, "[:alpha:]+(\\-|\\s|\\S)\\d+\\-?\\d*")
    ) |> str_to_upper(),
    block = ifelse(
      str_detect(`Sample ID`, "ROAD"),
      1,
      str_extract(`Sample ID`, "(?<=(Cl|CL|U))\\d+")
    ),
    rep = str_extract(`Sample ID`, "(?<=r)\\d+"),
    tree = paste(genotype, block, sep = "-"),
    fraction = case_when(
      str_detect(`Sample ID`, "MOAM") ~ "MAOM",
      str_detect(`Sample ID`, "POM") ~ "POM",
      .default = "Tot"
    ),
    depth = ifelse(fraction %in% "Tot",
                   ifelse(str_detect(`Sample ID`, "_30"),
                          "30",
                          "15"),
                   "15"),
    across(genotype:depth, fct)
  ) |>
  select(-`Sample ID`) |>
  filter(depth  %in% "15") |>
  rename(pN = `Final %N`,
         pC = `Final %C`,
         weight = `Weight [mg]`) |>
  pivot_wider(
    names_from = fraction,
    values_from = c(weight, pN, pC),
    names_glue = "{fraction}{.value}"
  ) |>
  left_join(df_frac |> select(tree, pom_frac, maom_frac), by = "tree") |>
  drop_na() |>
  mutate(
    MAOMN = maom_frac * MAOMpN * 10,
    POMN = pom_frac * POMpN * 10,
    TotN = TotpN * 10,
    MAOMC = maom_frac * MAOMpC * 10,
    POMC = pom_frac * POMpC * 10,
    TotC = TotpC * 10,
    C_chk = abs((MAOMC + POMC - TotC) / TotC * 100),
    C_chk_abs = (MAOMC + POMC - TotC),
    C_chk_rel = (MAOMC + POMC - TotC) / TotC * 100,
    C_chk_id = fct(ifelse(C_chk > 30, "Outlier", "Good"))
  )

# Read in updated samples from Tommy Mead on 03/07/24
df_core_update_2 <-
  list.files(file.path(d_path,
                       "02-ldrd/"), "*Rerun*", full.names = TRUE) |>
  map(
    \(x)
    read_xlsx(x,
              sheet = "Output") |>
      filter(Type %in% "Test") |>
      select(`Sample ID`:`Final %C`)
  ) |>
  bind_rows() |>
  mutate(
    genotype = 
      str_extract(`Sample ID`, "[:alpha:]+(\\-|\\s|\\S)\\d+") |> 
      str_to_upper() |> 
      str_replace(" ","-"),
    block = str_extract(`Sample ID`, "(?<=(\\-))\\d+"),
    rep = str_extract(`Sample ID`, "\\d+$"),
    tree = paste(genotype, block, sep = "-"),
    fraction = case_when(
      str_detect(`Sample ID`, "MOAM") ~ "MAOM",
      str_detect(`Sample ID`, "Pom") ~ "POM",
      .default = "Tot"
    ),
    depth = ifelse(fraction %in% "Tot",
                   ifelse(str_detect(`Sample ID`, "_30"),
                          "30",
                          "15"),
                   "15"),
    across(genotype:depth, fct)
  ) |>
  select(-`Sample ID`) |>
  filter(depth  %in% "15") |>
  rename(pN = `Final %N`,
         pC = `Final %C`,
         weight = `Weight [mg]`) |>
  pivot_wider(
    names_from = fraction,
    values_from = c(weight, pN, pC),
    names_glue = "{fraction}{.value}"
  ) |>
  left_join(df_frac |> select(tree, pom_frac, maom_frac), by = "tree") |>
  drop_na() |>
  mutate(
    MAOMN = maom_frac * MAOMpN * 10,
    POMN = pom_frac * POMpN * 10,
    TotN = TotpN * 10,
    MAOMC = maom_frac * MAOMpC * 10,
    POMC = pom_frac * POMpC * 10,
    TotC = TotpC * 10,
    C_chk = abs((MAOMC + POMC - TotC) / TotC * 100),
    C_chk_abs = (MAOMC + POMC - TotC),
    C_chk_rel = (MAOMC + POMC - TotC) / TotC * 100,
    C_chk_id = fct(ifelse(C_chk > 30, "Outlier", "Good"))
  )
```

</details>
<details class="code-fold">
<summary>Code</summary>

``` r
# Combine new samples with fractionation data
df_update_plt <-
  bind_rows(list(`01/24 EA` = df_core_update_1,
                 `03/24 EA` = df_core_update_2,
                 "Original EA" = df_check |>
              filter(C_chk_id %in% "Outlier")), .id = "run") |>
  select(run, tree, genotype, MAOMC:C_chk_id)

g1a <-  df_update_plt |>
  ggplot(aes(
    x = tree,
    y = C_chk_rel,
    shape = C_chk_id,
    color = run,
    fill = run
  )) +
  geom_hline(yintercept = 30,
             linetype = 2,
             color = "red") +
  geom_hline(yintercept = -30,
             linetype = 2,
             color = "red") +
  geom_jitter(size = 3,
              position = position_dodge2(width = 0.4)) +
  scale_shape_manual(name = NULL,
                     values = c("Outlier" = 24, "Good" = 21)) +
  scale_fill_manual(
    name = "Run",
    values = c(
      `01/24 EA` = "red",
      `03/24 EA` = "green",
      `Original EA` = "white"
    )
  ) +
  scale_color_manual(
    name = "Run",
    values = c(
      `01/24 EA` = "red",
      `03/24 EA` = "green",
      `Original EA` = "black"
    )
  ) +
  xlab(NULL) +
  ylab("Percent Error,\n (MAOM C + POM C - Total C)/Total C") +
  theme_cowplot(fs) +
  theme(axis.text.x = element_blank()) 

g1b <-  df_update_plt |>
  ggplot(aes(
    x = tree,
    y = C_chk_abs,
    shape = C_chk_id,
    color = run,
    fill = run
  )) +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "red") +
  geom_jitter(size = 3,
              position = position_dodge2(width = 0.4)) +
  scale_shape_manual(name = NULL,
                     values = c("Outlier" = 24, "Good" = 21)) +
  scale_fill_manual(
    name = "Run",
    values = c(
      `01/24 EA` = "red",
      `03/24 EA` = "green",
      `Original EA` = "white"
    )
  ) +
  scale_color_manual(
    name = "Run",
    values = c(
      `01/24 EA` = "red",
      `03/24 EA` = "green",
      `Original EA` = "black"
    )
  ) +
  xlab(NULL) +
  ylab("Absolute Error [mg C/g soil],\n MAOM C + POM C - Total C") +
  theme_cowplot(fs) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))

g1a / g1b + plot_layout(guides = "collect")
```

</details>
<div id="fig-c-re-run">

![](01-harmonize-clatskanie-data-pub_files/figure-commonmark/fig-c-re-run-1.png)


Figure 2: Sanity check on soil C mass balance showing the sums of
Mineral Associated Soil C and Particulate Soil C observations versus the
Total Soil C..

</div>
<details class="code-fold">
<summary>Code</summary>

``` r
df_update_plt |>
  pivot_longer(cols = c(TotC, MAOMC, POMC)) |>
  mutate(name = fct_recode(
    name,
    "Total C" = "TotC",
    "MAOM C" = "MAOMC",
    "POM C" = "POMC"
  )) |>
  ggplot(aes(
    x = tree,
    y = value,
    shape = C_chk_id,
    color = run,
    fill = run
  )) +
  geom_jitter(size = 3,
              position = position_dodge2(width = 0.4)) +
  scale_shape_manual(name = NULL,
                     values = c("Outlier" = 24, "Good" = 21)) +
  scale_fill_manual(
    name = "Run",
    values = c(
      `01/24 EA` = "red",
      `03/24 EA` = "green",
      `Original EA` = "white"
    )
  ) +
  scale_color_manual(
    name = "Run",
    values = c(
      `01/24 EA` = "red",
      `03/24 EA` = "green",
      `Original EA` = "black"
    )
  ) +
  ylab("Soil C [mg C/g soil]") +
  xlab(NULL) +
  facet_wrap( ~ name, ncol = 1, scales = "free_y") +
  theme_cowplot(fs) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))
```

</details>
<div id="fig-c-re-run-frac">

![](01-harmonize-clatskanie-data-pub_files/figure-commonmark/fig-c-re-run-frac-1.png)


Figure 3: Sanity check on soil C mass balance showing the sums of
Mineral Associated Soil C and Particulate Soil C observations versus the
Total Soil C..

</div>

# SEED Root Chemistry Data

Root chemistry was analyzed at each of the SEED soil cores. The
differing nutrients and minerals have units of percent (*p*) and parts
per million (*ppm*). These results can be joined to the soil core data
using the *tree* column.

<details class="code-fold">
<summary>Code</summary>

``` r
# SEED root chemistry data
df_root_ch <-   read_xlsx(
    file.path(d_path,"01-seed/seed-roots-10-10-23.xlsx"),
              sheet = "Sheet1") |> 
  mutate(genotype = as.factor(str_replace(Genotype,"_","-")),
         block = as.factor(Block),
         tree = fct_cross(genotype,block,sep = "-")
         ) |> 
  select(tree,genotype,block,pCa:`Zn ppm`) 


display_table(df_root_ch)
```

</details>
<div id="tbl-seed-rc">

<div class="cell-output-display">

| tree | genotype | block | pCa | pK | pMg | pP | pC | pN | pS | Lig | Alppm | Bppm | Cdppm | Crppm | Cuppm | Feppm | Mn ppm | Mo ppm | Na ppm | Ni ppm | Pb ppm | Zn ppm |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|:---|---:|---:|:---|---:|
| BESC-102-1 | BESC-102 | 1 | 0.72 | 0.44 | 0.16 | 0.08 | 47.49 | 0.96 | 0.08 | 26.58 | 7378.17 | 16.35 | \<1.76 | 73.56 | 15.17 | 6213.76 | 152.44 | 2.2869999999999999 | 285.62 | 7.07 | 8.8350000000000009 | 55.16 |
| BESC-102-2 | BESC-102 | 2 | 0.82 | 0.36 | 0.13 | 0.11 | 51.28 | 1.02 | 0.11 | 32.71 | 5044.21 | 11.64 | \<0.80 | 54.67 | 21.75 | 3243.34 | 283.39 | 2.1030000000000002 | 307.61 | 9.60 | 4.9130000000000003 | 73.66 |
| BESC-102-3 | BESC-102 | 3 | 0.76 | 0.61 | 0.20 | 0.14 | 44.07 | 1.11 | 0.11 | 24.03 | 8403.85 | 16.98 | \<0.80 | 256.51 | 23.01 | 8770.00 | 165.94 | 5.6139999999999999 | 335.89 | 15.06 | 7.0030000000000001 | 51.25 |
| BESC-876-1 | BESC-876 | 1 | 0.88 | 0.37 | 0.14 | 0.11 | 46.33 | 0.90 | 0.12 | 25.67 | 5301.36 | 16.35 | 0.91200000000000003 | 93.37 | 25.61 | 7784.00 | 457.08 | 2.9540000000000002 | 199.36 | 9.89 | 6.1509999999999998 | 56.24 |
| BESC-876-2 | BESC-876 | 2 | 0.82 | 0.47 | 0.15 | 0.10 | 50.21 | 0.74 | 0.11 | 31.78 | 4881.96 | 13.77 | \<0.80 | 66.27 | 22.16 | 3890.71 | 490.52 | 2.036 | 260.95 | 9.88 | 6.226 | 78.03 |
| GW-9776-1 | GW-9776 | 1 | 0.99 | 0.47 | 0.15 | 0.10 | 48.77 | 0.81 | 0.11 | 31.41 | 3859.18 | 27.43 | 1.0920000000000001 | 59.25 | 27.30 | 6618.00 | 703.30 | 1.6080000000000001 | 255.57 | 10.76 | 3.262 | 75.29 |

</div>

Table 3: SEED soil core data collected in 2022.

</div>
<div id="tbl-seed-rc-key">

| Variable | Description | Units |
|----|----|----|
| pCa, pK, pMg, pP, pC, pN, pS | Root Calcium, Potassium, Magnesium, Phosphorus, Sulfur percentages | % |
| Lig | Root Lignin | % |
| Alppm, Bppm, Cdppm, Crppm, Cuppm, Mnppm, Mo ppm, Na ppm, Ni ppm, Pb ppm, Zn ppm | Root Aluminum, Boron, Cadmium, Chromium, Copper, Iron, Manganese, Molibdinum, Sodium, Nickel, Lead, and Zinc concentrations | parts per million (ppm) |

Table 4: <a href="#tbl-seed-rc" class="quarto-xref">Table 3</a> key

</div>

# CBI Aboveground Biomass

A time series of abovegound biomass measurements from 2009-2019 was
collected from *Poplar Innovations* as part of Center for Bioenergy
Innovation (CBI) research. The primary measurements are tree diameters
that were converted to aboveground woody biomass (trunk and branches)
based on the allometric equations from Truax et al. (2014).

<details class="code-fold">
<summary>Code</summary>

``` r
# Calculate growth rate
df_agb_rate <- read_csv(file.path(d_path, "03-cbi/cbi-agb.csv"))|>
  mutate(across(c(tree, genotype, block), as.factor),
         agb = 0.071 * dbh ^ (2.4055)) |> 
  group_by(tree, coppice_number) |> 
  summarize(agb = agb[which.max(year)]/10,
            year = year[which.max(year)],
            cop_id = coppice_number[which.max(year)],
            height = max(height, na.rm = TRUE),
            d20 = max(d20, na.rm = TRUE),
            d50 = max(d50, na.rm = TRUE),
            dbh = max(dbh, na.rm = TRUE),
            ) |> 
  group_by(tree) |> 
  summarize(across(c(agb,cop_id), sum),
            across(c(height,d20,d50,dbh), \(x) max(x, na.rm = TRUE))) |> 
  ungroup() |> 
  relocate(agb,.after=last_col())

display_table(df_agb_rate)
```

</details>
<div id="tbl-cbi-agb">

<div class="cell-output-display">

| tree       | cop_id | height |  d20 |  d50 |  dbh |   agb |
|:-----------|-------:|-------:|-----:|-----:|-----:|------:|
| BESC-102-1 |      1 |    930 |  7.9 |  5.0 | 20.7 | 10.75 |
| BESC-102-2 |      0 |   1000 |  8.9 |  8.0 | 18.3 |  7.73 |
| BESC-102-3 |      0 |    710 | -Inf | -Inf | 24.8 | 16.06 |
| BESC-119-1 |      1 |    720 |  5.5 |  4.9 | 16.8 |  6.69 |
| BESC-119-2 |      0 |   1260 |  8.4 |  7.4 | 19.1 |  8.57 |
| BESC-119-3 |      0 |    790 | -Inf | -Inf |  8.9 |  1.36 |

</div>

Table 5: LDRD aboveground biomass data collected from 2009 -2019

</div>
<div id="tbl-2">

| Variable | Description | Units |
|----|----|----|
| cop_id | Flag to indicate if a given year is before or after coppicing. |  |
| height | Tree height | cm |
| d20 | Diameter at 20 cm height | cm |
| d50 | Diameter at 50 cm height | cm |
| dbh | Diameter at breast height | cm |
| agb | Estimated aboveground biomass growth rate from Truax et al. (2014) | kg/yr |

Table 6: <a href="#tbl-cbi-agb" class="quarto-xref">Table 5</a> key

</div>

# LDRD Soil Chemistry Data

Soil chemistry data for the 69 SEED soil cores were measured at UGA in
May 2024 in order to check whether the root chemistry influence on soil
C was an active or passive effect.

<div id="tbl-ldrd-sc">

<details class="code-fold">
<summary>Code</summary>

``` r
df_soil_chem <-
  read_xlsx(file.path(d_path, "02-ldrd/soil-chem-uga-5-25-24.xlsx"),
            skip = 9) |>
  drop_na(Sample) |>
  mutate(
    genotype =
      str_extract(`Sample`, "(SLMC-28-2|[:alpha:]+(\\-|\\s|\\S)(\\d+))") |>
      str_to_upper() |>
      str_replace(" ", "-"),
    block = str_extract(`Sample`, "(?<=(CL\\s))\\d+"),
    tree = paste(genotype, block, sep = "-"),
    across(.cols = c("genotype", "block","tree"), fct)
    ) |> 
  rename(
    pH = "pH 2",
    LBC = "LBC 1",
    base_sat = "Base\r\nSatur-\r\nation"
  ) |> 
  mutate(across(c(Ca,K,Mg,P), \(x) x/1e4, .names = "p{.col}")) |> 
  select(-c(Ca,K,Mg,P)) |> 
  rename_with(\(x) paste0(x, "_ppm"), c(Cd:Zn,Al)) |> 
  relocate(tree, genotype, block) |> 
  select(-Lab, -Sample)
```

</details>


Table 7: LDRD soil chemistry data from 2022 SEED soil cores measured in
2024.

</div>
<div id="tbl-ldrd-sc-key">

| Variable | Description | Units |
|----|----|----|
| LBC, LBCeq | Lime buffer capacity | ppm CaCO3/pH |
| pH | Soil pH |  |
| CEC | Cation exchange capacity | milliequivalents per 100 grams |
| base_sat | Soil base saturation of CEC | % |
| pCa, pK, pMg, pP | Root Calcium, Potassium, Magnesium, Phosphoruspercentages | % |
| Al_ppm, Cd_ppm, Cr_ppm, Cu_ppm, Fe_ppm, Mn_ppm, Mo_ppm, Na_ppm, Ni_ppm, Pb_ppm, Zn_ppm | Aluminum, Cadmium, Chromium, Copper, Iron, Manganese, Molibdinum, Sodium, Nickel, Lead, and Zinc concentrations. | Parts per million |

Table 8: <a href="#tbl-ldrd-sc" class="quarto-xref">Table 7</a> key

</div>

# Create GCB Manuscript Analysis Dataset

We stitched together the complete tree, soil and root data set covering
69 trees (24 genotypes), focusing on primarily the 0-15 cm depth for the
manuscript analysis. The coding was tricky given we had to update
multiple original soil core samples with two differing updates (1/24/24
and 3/24/24) to the SOC fractionation as well as add soil chemistry data
taken on 5/24/24. We have excluded Cadmium, Chromium, Molibdinum, and
Lead from the root and/or soil chemistry data since many of these were
either not measured for both root and soil or they had only trace
amounts for many of the samples. The final data set is stored at
*./02-data/02-processed/clatskanie-c-fit-data.csv* and described below
in <a href="#tbl-gcb-data" class="quarto-xref">Table 9</a> -
<a href="#tbl-gcb-data-key" class="quarto-xref">Table 10</a>.

<details class="code-fold">
<summary>Code</summary>

``` r
# Initially define the fit data set.
df_fit_orig <-
  df_core |> 
  filter(Depth == 15) |>
  select(tree:TotpC, TotC) |>
  left_join(
    df_root_ch |>
      select(
        tree:block,
        pCa:Lig,
        Alppm:Bppm,
        Crppm:`Mn ppm`,
        `Na ppm`,
        `Ni ppm`,
        `Zn ppm`
      ),
    by = c("tree","genotype","block")
  )

# Create update table to select BESC-13 samples from the 1/24/24 EA analysis, and the remaining 7 outlier samples from the 3/24/24 EA runs.
to_update <- bind_rows(df_core_update_1,df_core_update_2, .id = "run") |>
  select(run,tree,genotype,TotpN:TotC) |> 
  group_by(run,tree,genotype) |> 
  summarize(across(.cols = everything(),mean)) |> 
  ungroup() |> 
  filter(genotype %in% "BESC-13" | tree %in% "BESC-102-2" | run == 2) |> 
  select(tree, TotpN:TotC) |> 
  select(-pom_frac,-maom_frac,-TotN) |> 
  mutate(Depth = 15,
         TotpC = TotC/10)

# Update the rows
df_fit_update <- df_fit_orig |>
  rows_update(
    to_update,
    by  = "tree",
    unmatched = "ignore"
  )


# Create a 30 cm concentration and stock variable to add. Note: BESC-145-1 is missing BD at 30 cm, while BESC-351-3 is missing TotpC at 30 cm. I have filled these based on the 15 cm values. 
df_core_update <- df_core |>
  rows_update(
    to_update,
    by  = c("tree", "Depth"),
    unmatched = "ignore"
  )

df_core_30 <- df_core_update |> 
  filter(Depth <= 30) |> 
  group_by(tree) |> 
  fill(c(BD, TotpC, TotC), .direction = "down") |> 
  select(tree, Depth, BD, TotpC, TotC) |> 
  summarize(TotC_st_30 = sum(TotpC * BD * 15),
            TotC_30 = sum(TotC * BD)/sum(BD))



# Create final data set for Clatskanie soil C manuscript.
df_fit <- df_fit_update |>
  rename_with(.cols = contains("ppm"), \(x) gsub("(\\sppm)|(ppm)","_ppm",x)) |> 
  mutate(root_CN = pC / pN) |>
  left_join(df_agb_rate, by = "tree") |>
  arrange(genotype, block) |>
  rowid_to_column("ord") |>
  mutate(
    tree = fct_reorder(tree, ord, .desc = TRUE),
    SiCl = Silt + Clay,
    # convert to stock units tonnes C/ha
    TotC_st = TotpC * BD * 15,
    MAOMC_st = MAOMC * 15 * BD / 10,
    POMC_st = POMC * 15 * BD / 10,
    agb = agb / 9 * 10 * 0.45,    # convert to units tonnes C/ha/yr (see above)
    pMAOM = MAOMC / TotC * 100,
    pPOM = POMC / TotC * 100,
    MAOM_POM = MAOMC / POMC,
    TotC_agb = TotC / agb,
    MAOM_CN = MAOMC/MAOMN,
    POM_CN = POMC/POMN,
    soil_CN = TotpC / TotpN,
    cop_id = as.factor(ifelse(cop_id == 1, "Coppiced", "Non-coppiced")),
    Trait = as.factor(Trait),
    # Check for outliers
    C_chk = abs((MAOMC + POMC - TotC) / TotC * 100),
    C_chk_abs = (MAOMC + POMC - TotC),
    C_chk_rel = (MAOMC + POMC - TotC) / TotC * 100,
    C_chk_id = fct(ifelse(C_chk > 30, "Outlier", "Good"))
  ) |> 
  # Add 0-30 cm SOC
  left_join(df_core_30, by = "tree") |>
  # Add in new soil chemistry data from 5/24/24
  left_join(df_soil_chem |> 
                        select(-genotype,-block, -where(is.character)),
                      by = "tree",
                      suffix = c("_root", "_soil"))

# Create a variable that calculates the differences between soil and root
root_soil_diff <- df_fit |> select(tree, contains(c("_root", "_soil"))) |>
  pivot_longer(
    pH_root:pP_soil,
    names_to = c("mineral", "organ"),
    names_pattern = "(.*)_(root|soil)"
  ) |>
  group_by(tree, mineral) |>
  summarize(diff = value[organ %in% "root"] - value[organ %in% "soil"]) |>
  ungroup() |>
  pivot_wider(names_from = "mineral",
              values_from = "diff",
              names_glue = "{mineral}_diff")

# Add differences and save
df_fit <- df_fit |> left_join(root_soil_diff, by = "tree") |> 
  rename(B_ppm_root = B_ppm,
         pS_root = pS) |> 
        select(ord:longitude,cop_id, height, dbh, agb, root_CN, Lig, BD,Sand:Silt, SiCl, pH_soil,CEC,base_sat, POMC, MAOMC, TotC, TotC_30, POMC_st, MAOMC_st, TotC_st, TotC_st_30, POM_CN, MAOM_CN, soil_CN, C_chk_rel, C_chk_id, contains(c("Al_ppm", "B_ppm", "Cu_ppm", "Fe_ppm","Mn_ppm","Na_ppm","Ni_ppm","Zn_ppm")), contains(c("pCa","pK","pMg","pP","pS")))

write_csv(df_fit, "./02-data/02-processed/clatskanie-c-fit-data.csv") 
```

</details>
<details class="code-fold">
<summary>Code</summary>

``` r
display_table(df_fit)
```

</details>
<div id="tbl-gcb-data">

<div class="cell-output-display">

| ord | tree | genotype | block | Depth | latitude | longitude | cop_id | height | dbh | agb | root_CN | Lig | BD | Sand | Clay | Silt | SiCl | pH_soil | CEC | base_sat | POMC | MAOMC | TotC | TotC_30 | POMC_st | MAOMC_st | TotC_st | TotC_st_30 | POM_CN | MAOM_CN | soil_CN | C_chk_rel | C_chk_id | Al_ppm_root | Al_ppm_soil | Al_ppm_diff | B_ppm_root | Cu_ppm_root | Cu_ppm_soil | Cu_ppm_diff | Fe_ppm_root | Fe_ppm_soil | Fe_ppm_diff | Mn_ppm_root | Mn_ppm_soil | Mn_ppm_diff | Na_ppm_root | Na_ppm_soil | Na_ppm_diff | Ni_ppm_root | Ni_ppm_soil | Ni_ppm_diff | Zn_ppm_root | Zn_ppm_soil | Zn_ppm_diff | pCa_root | pCa_soil | pCa_diff | pK_root | pK_soil | pK_diff | pMg_root | pMg_soil | pMg_diff | pP_root | Cr_ppm | pPOM | pP_soil | pP_diff | pS_root |
|---:|:---|:---|:---|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | BESC-102-1 | BESC-102 | 1 | 15 | 46.12 | -123.27 | Coppiced | 930 | 20.7 | 5.37 | 49.73 | 26.58 | 0.88 | 0.00 | 45.90 | 55.75 | 101.65 | 5.26 | 32.23 | 42.16 | 11.56 | 39.52 | 60.40 | 44.22 | 15.33 | 52.40 | 80.09 | 99.33 | 21.81 | 12.96 | 14.88 | -15.43 | Good | 7378.17 | 4.47 | 7373.70 | 16.35 | 15.17 | 4.10 | 11.07 | 6213.76 | 267.51 | 5946.25 | 152.44 | 67.24 | 85.20 | 285.62 | 37.81 | 247.82 | 7.07 | 2.11 | 4.96 | 55.16 | 3.30 | 51.86 | 0.72 | 0.21 | 0.51 | 0.44 | 0.02 | 0.42 | 0.16 | 0.03 | 0.13 | 0.08 | 73.56 | 19.14 | 0.01 | 0.07 | 0.08 |
| 2 | BESC-102-2 | BESC-102 | 2 | 15 | 46.12 | -123.27 | Non-coppiced | 1000 | 18.3 | 3.86 | 50.42 | 32.71 | 0.62 | 7.85 | 45.15 | 47.00 | 92.15 | 5.27 | 32.70 | 47.78 | 9.37 | 52.97 | 76.37 | 70.86 | 8.70 | 49.21 | 70.95 | 132.22 | 20.34 | 13.83 | 14.60 | -18.37 | Good | 5044.21 | 2.92 | 5041.29 | 11.64 | 21.75 | 6.08 | 15.67 | 3243.34 | 363.77 | 2879.57 | 283.39 | 59.00 | 224.39 | 307.61 | 53.09 | 254.52 | 9.60 | 3.03 | 6.56 | 73.66 | 7.04 | 66.62 | 0.82 | 0.24 | 0.58 | 0.36 | 0.01 | 0.35 | 0.13 | 0.04 | 0.09 | 0.11 | 54.67 | 12.26 | 0.01 | 0.09 | 0.11 |
| 3 | BESC-102-3 | BESC-102 | 3 | 15 | 46.12 | -123.27 | Non-coppiced | 710 | 24.8 | 8.03 | 39.81 | 24.03 | 0.79 | 8.39 | 38.21 | 51.20 | 89.41 | 5.48 | 28.92 | 73.51 | 7.59 | 51.99 | 53.50 | 42.85 | 9.01 | 61.69 | 63.48 | 104.23 | 20.51 | 14.13 | 13.90 | 11.36 | Good | 8403.85 | 0.14 | 8403.71 | 16.98 | 23.01 | 5.05 | 17.95 | 8770.00 | 215.49 | 8554.51 | 165.94 | 93.89 | 72.05 | 335.89 | 43.53 | 292.36 | 15.06 | 2.82 | 12.24 | 51.25 | 9.72 | 41.53 | 0.76 | 0.33 | 0.43 | 0.61 | 0.02 | 0.59 | 0.20 | 0.05 | 0.15 | 0.14 | 256.51 | 14.19 | 0.02 | 0.13 | 0.11 |
| 4 | BESC-119-1 | BESC-119 | 1 | 15 | 46.12 | -123.27 | Coppiced | 720 | 16.8 | 3.34 | 51.94 | 35.71 | 0.82 | 0.00 | 37.90 | 62.87 | 100.77 | 4.78 | 30.19 | 43.10 | 7.07 | 32.04 | 42.43 | 29.88 | 8.66 | 39.26 | 52.00 | 72.04 | 20.79 | 12.09 | 12.78 | -7.83 | Good | 2109.76 | 4.64 | 2105.12 | 11.84 | 24.69 | 4.75 | 19.95 | 1833.12 | 411.20 | 1421.92 | 149.32 | 57.42 | 91.90 | 107.56 | 34.50 | 73.06 | 4.88 | 1.88 | 3.00 | 57.83 | 4.72 | 53.12 | 0.90 | 0.20 | 0.71 | 0.72 | 0.03 | 0.68 | 0.10 | 0.03 | 0.07 | 0.12 | 14.03 | 16.66 | 0.01 | 0.10 | 0.12 |
| 5 | BESC-13-1 | BESC-13 | 1 | 15 | 46.12 | -123.27 | Non-coppiced | 670 | 12.6 | 1.57 | 40.69 | 27.05 | 0.99 | 1.89 | 34.90 | 63.21 | 98.11 | 5.25 | 25.40 | 67.60 | 2.70 | 21.61 | 24.70 | 24.70 | 4.01 | 32.12 | 36.72 | 69.09 | 11.93 | 11.03 | 9.32 | -1.59 | Good | 1047.54 | 0.41 | 1047.13 | 15.29 | 10.38 | 5.76 | 4.62 | 952.77 | 350.07 | 602.70 | 73.41 | 74.21 | -0.80 | 372.08 | 33.46 | 338.62 | 3.91 | 2.85 | 1.06 | 65.39 | 7.45 | 57.94 | 0.54 | 0.26 | 0.28 | 0.83 | 0.03 | 0.80 | 0.10 | 0.04 | 0.06 | 0.12 | NA | 10.92 | 0.02 | 0.10 | 0.11 |
| 6 | BESC-13-2 | BESC-13 | 2 | 15 | 46.12 | -123.27 | Non-coppiced | 1150 | 23.0 | 6.70 | 52.09 | 31.56 | 0.84 | 8.27 | 36.83 | 54.90 | 91.73 | 5.27 | 26.33 | 71.23 | 14.70 | 28.90 | 49.10 | 49.10 | 18.55 | 36.50 | 61.99 | 114.76 | 23.83 | 10.84 | 13.18 | -11.20 | Good | 2326.70 | 0.53 | 2326.17 | 14.38 | 17.18 | 6.37 | 10.81 | 2202.72 | 447.65 | 1755.07 | 243.11 | 84.15 | 158.96 | 102.42 | 38.25 | 64.18 | 8.19 | 3.13 | 5.07 | 93.19 | 9.29 | 83.90 | 0.70 | 0.30 | 0.41 | 0.48 | 0.03 | 0.45 | 0.12 | 0.04 | 0.08 | 0.14 | NA | 29.93 | 0.01 | 0.12 | 0.09 |

</div>

Table 9: The final processed data set used for the analysis in the GCB
paper.

</div>
<div id="tbl-gcb-data-key">

| Variable | Description | Units |
|----|----|----|
| ord | Order of gentoypes used for manuscript plotting |  |
| tree | Unique ID for all poplar trees across these data sets (max of 66 trees) |  |
| genotype | ID for genotypes contained across these data sets (max of 23 genotypes) |  |
| block | ID for replicate blocks (max of 3 replicates) |  |
| Depth | Soil core depth at each tree. Primarily use 15 cm for manuscript, see raw data for deeper samples (max 60 cm for select trees) | cm |
| latitude | Latitude of tree | degrees |
| longitude | Longitude of tree | degrees |
| cop_id | Flag to indicate if a given year is before or after coppicing. |  |
| height | Tree height | cm |
| dbh | Tree diameter at breast height | cm |
| agb | Estimated aboveground biomass growth rate from Truax et al. (2014) based on dbh | kg/yr |
| root_CN | Root carbon to nitrogen ratio |  |
| Lig | Root lignin percentage | % |
| BD | Soil bulk density | g/cm$^3$ |
| pH_soil | Soil acidity taken from LDRD soil chemistry data. |  |
| CEC | Cation exchange capacity | milliequivalents per 100 grams |
| base_sat | Soil base saturation of CEC | % |
| Sand, Silt, Clay | Percent sand, silt or clay | % |
| POMC, MAOMC, TotC | Particulate (POM), mineral-associated (MAOM), and Total soil organic matter C concentrations | mg C/g soil |
| POMC_st, MAOMC_st, TotC_st | Particulate (POM), mineral-associated (MAOM), and Total soil organic matter C stocks | tonnes C/ha |
| POMC_CN, MAOMC_CN, soil_CN | Particulate (POM), mineral-associated (MAOM), and Total soil organic matter carbon to nitrogen ratios | % |
| TotC_30, TotC_st_30 | Total SOC concentration and stock from 0-30 cm rather than 0-15 cm as the other observations. | mg C/g soil or tonnes C/ha |
| C_chk_rel, C_chk_id | Percentage mass balance error in POMC + MAOMC versus TotC, with less than 30% error being acceptable. | % |
| Al_ppm, B_ppm, Cu_ppm, Fe_ppm, Mn_ppm, Mo_ppm, Na_ppm, Ni_ppm, Zn_ppm | Aluminum, Boron, Copper, Iron, Manganese, Molibdinum, Sodium, Nickel, and Zinc concentrations. | parts per million (ppm) |
| pCa, pK, pMg, pP, pS | Calcium, Potassium, Magnesium, Phosphorus, Sulfur percentages | % |
| \_root, \_soil, \_diff | Suffix indicating either the root or soil elemental or nutrient quantity or the difference in root and soil quantity. |  |

Table 10: <a href="#tbl-gcb-data" class="quarto-xref">Table 9</a> key

</div>

# Supplemental Figures

Additionally, we created supplemental Figures S4 and S6 used in the
manuscript below. These files were created here as they relied on soil
observations that were not included in the final output data.

<details class="code-fold">
<summary>Code</summary>

``` r
g1 <- df_core_update |>
  filter(Depth <= 30) |> 
  mutate(Depth = fct(as.character(Depth))) |> 
  ggplot(aes(y=TotC,x=BD, color=Depth)) +
  geom_point() +
  geom_smooth(method=MASS::rlm,se=FALSE)  +
  xlab("Soil Bulk Density (g/cm3)") +
  ylab("Total Soil C \n [mg C/g soil]") +
  ylim(c(0,100)) +
  theme_cowplot(9) +
  theme(legend.position="inside",legend.position.inside = c(0.8,0.8))

g2 <- df_core_update |>
  select(tree,Depth, BD, Clay, TotC) |> 
  filter(Depth <= 30) |> 
  mutate(Depth = fct(as.character(Depth))) |> 
  pivot_longer(cols = c(BD, Clay, TotC)) |>
  pivot_wider(names_from=Depth, values_from=value) |> 
  ggplot(aes(y=`15`,x=`30`)) +
  geom_abline(linetype="dashed", color="red") +
  geom_point() +
  #geom_smooth(method=MASS::rlm,se=FALSE)  +
  xlab("15-30 cm") +
  ylab("0-15 cm") +
  facet_wrap(~name, scales="free", ncol = 3,
      labeller = as_labeller(
      c(
        BD = "Bulk Density [g/cm3]",
        Clay = " % Clay",
        TotC = "Total Soil C [mg C/g soil]"
      )
    )) +
  theme_cowplot(9)


(g1 / g2) + plot_annotation(tag_levels = 'a', tag_suffix = ')')

ggsave("03-figs/fig-s4.jpg",
        width = 6,
        height = 4,
        units = "in")
```

</details>
<div id="fig-s4">

![](01-harmonize-clatskanie-data-pub_files/figure-commonmark/fig-s4-1.png)


Figure 4: Figure S4 in the manuscript looking at depth dependence soil
core observations.

</div>
<details class="code-fold">
<summary>Code</summary>

``` r
df_fit |> 
  select(-contains(c("B_ppm","pS_"))) |> 
  select(tree, contains(c("_root","_soil"))) |> 
  pivot_longer(c(-tree,-pH_soil),
               names_to = c("mineral", "organ"),
               names_pattern = "(.*)_(root|soil)") |> 
  pivot_wider(
    names_from = organ,
    values_from = value,
    values_fn = mean
  ) |> 
  ggplot(aes(x = soil, y = root, color = mineral)) +
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~mineral, scales = "free") + 
  theme_cowplot(8) + 
  theme(legend.position = "bottom",
        legend.justification = "center")
fig_path <- "./03-figs/"
ggsave2(file.path(fig_path, "fig-s6-soilvroot.jpg"),
        width = 6, height = 6, dpi = 600)
```

</details>
<div id="fig-s6">

![](01-harmonize-clatskanie-data-pub_files/figure-commonmark/fig-s6-1.png)


Figure 5: Figure S6 in the manuscript of soil versus root chemistry
measurements for 69 trees at Clatskanie

</div>

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-truax2014" class="csl-entry">

Truax, Benoit, Daniel Gagnon, Julien Fortier, and France Lambert. 2014.
“Biomass and Volume Yield in Mature Hybrid Poplar Plantations on
Temperate Abandoned Farmland.” *Forests* 5 (12): 3107–30.
<https://doi.org/10.3390/f5123107>.

</div>

</div>
