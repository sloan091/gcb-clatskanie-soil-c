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

<div id="lmkysutyao" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#lmkysutyao table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#lmkysutyao thead, #lmkysutyao tbody, #lmkysutyao tfoot, #lmkysutyao tr, #lmkysutyao td, #lmkysutyao th {
  border-style: none;
}
&#10;#lmkysutyao p {
  margin: 0;
  padding: 0;
}
&#10;#lmkysutyao .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#lmkysutyao .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#lmkysutyao .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#lmkysutyao .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#lmkysutyao .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#lmkysutyao .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#lmkysutyao .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#lmkysutyao .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#lmkysutyao .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#lmkysutyao .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#lmkysutyao .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#lmkysutyao .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#lmkysutyao .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#lmkysutyao .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#lmkysutyao .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#lmkysutyao .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#lmkysutyao .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#lmkysutyao .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#lmkysutyao .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#lmkysutyao .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#lmkysutyao .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#lmkysutyao .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#lmkysutyao .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#lmkysutyao .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#lmkysutyao .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#lmkysutyao .gt_left {
  text-align: left;
}
&#10;#lmkysutyao .gt_center {
  text-align: center;
}
&#10;#lmkysutyao .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#lmkysutyao .gt_font_normal {
  font-weight: normal;
}
&#10;#lmkysutyao .gt_font_bold {
  font-weight: bold;
}
&#10;#lmkysutyao .gt_font_italic {
  font-style: italic;
}
&#10;#lmkysutyao .gt_super {
  font-size: 65%;
}
&#10;#lmkysutyao .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#lmkysutyao .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#lmkysutyao .gt_indent_1 {
  text-indent: 5px;
}
&#10;#lmkysutyao .gt_indent_2 {
  text-indent: 10px;
}
&#10;#lmkysutyao .gt_indent_3 {
  text-indent: 15px;
}
&#10;#lmkysutyao .gt_indent_4 {
  text-indent: 20px;
}
&#10;#lmkysutyao .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="tree">tree</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="genotype">genotype</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="block">block</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pCa">pCa</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pK">pK</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pMg">pMg</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pP">pP</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pC">pC</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pN">pN</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pS">pS</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Lig">Lig</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Alppm">Alppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Bppm">Bppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Cdppm">Cdppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Crppm">Crppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Cuppm">Cuppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Feppm">Feppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Mn ppm">Mn ppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Mo ppm">Mo ppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Na ppm">Na ppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Ni ppm">Ni ppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Pb ppm">Pb ppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Zn ppm">Zn ppm</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="tree" class="gt_row gt_center">BESC-102-1</td>
<td headers="genotype" class="gt_row gt_center">BESC-102</td>
<td headers="block" class="gt_row gt_center">1</td>
<td headers="pCa" class="gt_row gt_right">0.72</td>
<td headers="pK" class="gt_row gt_right">0.44</td>
<td headers="pMg" class="gt_row gt_right">0.16</td>
<td headers="pP" class="gt_row gt_right">0.08</td>
<td headers="pC" class="gt_row gt_right">47.49</td>
<td headers="pN" class="gt_row gt_right">0.95</td>
<td headers="pS" class="gt_row gt_right">0.08</td>
<td headers="Lig" class="gt_row gt_right">26.58</td>
<td headers="Alppm" class="gt_row gt_right">7,378.17</td>
<td headers="Bppm" class="gt_row gt_right">16.35</td>
<td headers="Cdppm" class="gt_row gt_left">&lt;1.76</td>
<td headers="Crppm" class="gt_row gt_right">73.56</td>
<td headers="Cuppm" class="gt_row gt_right">15.17</td>
<td headers="Feppm" class="gt_row gt_right">6,213.76</td>
<td headers="Mn ppm" class="gt_row gt_right">152.44</td>
<td headers="Mo ppm" class="gt_row gt_right">2.2869999999999999</td>
<td headers="Na ppm" class="gt_row gt_right">285.62</td>
<td headers="Ni ppm" class="gt_row gt_right">7.07</td>
<td headers="Pb ppm" class="gt_row gt_right">8.8350000000000009</td>
<td headers="Zn ppm" class="gt_row gt_right">55.16</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-102-2</td>
<td headers="genotype" class="gt_row gt_center">BESC-102</td>
<td headers="block" class="gt_row gt_center">2</td>
<td headers="pCa" class="gt_row gt_right">0.82</td>
<td headers="pK" class="gt_row gt_right">0.36</td>
<td headers="pMg" class="gt_row gt_right">0.13</td>
<td headers="pP" class="gt_row gt_right">0.11</td>
<td headers="pC" class="gt_row gt_right">51.28</td>
<td headers="pN" class="gt_row gt_right">1.02</td>
<td headers="pS" class="gt_row gt_right">0.11</td>
<td headers="Lig" class="gt_row gt_right">32.71</td>
<td headers="Alppm" class="gt_row gt_right">5,044.21</td>
<td headers="Bppm" class="gt_row gt_right">11.64</td>
<td headers="Cdppm" class="gt_row gt_left">&lt;0.80</td>
<td headers="Crppm" class="gt_row gt_right">54.67</td>
<td headers="Cuppm" class="gt_row gt_right">21.75</td>
<td headers="Feppm" class="gt_row gt_right">3,243.34</td>
<td headers="Mn ppm" class="gt_row gt_right">283.39</td>
<td headers="Mo ppm" class="gt_row gt_right">2.1030000000000002</td>
<td headers="Na ppm" class="gt_row gt_right">307.61</td>
<td headers="Ni ppm" class="gt_row gt_right">9.60</td>
<td headers="Pb ppm" class="gt_row gt_right">4.9130000000000003</td>
<td headers="Zn ppm" class="gt_row gt_right">73.66</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-102-3</td>
<td headers="genotype" class="gt_row gt_center">BESC-102</td>
<td headers="block" class="gt_row gt_center">3</td>
<td headers="pCa" class="gt_row gt_right">0.76</td>
<td headers="pK" class="gt_row gt_right">0.61</td>
<td headers="pMg" class="gt_row gt_right">0.20</td>
<td headers="pP" class="gt_row gt_right">0.14</td>
<td headers="pC" class="gt_row gt_right">44.07</td>
<td headers="pN" class="gt_row gt_right">1.11</td>
<td headers="pS" class="gt_row gt_right">0.11</td>
<td headers="Lig" class="gt_row gt_right">24.03</td>
<td headers="Alppm" class="gt_row gt_right">8,403.85</td>
<td headers="Bppm" class="gt_row gt_right">16.98</td>
<td headers="Cdppm" class="gt_row gt_left">&lt;0.80</td>
<td headers="Crppm" class="gt_row gt_right">256.51</td>
<td headers="Cuppm" class="gt_row gt_right">23.01</td>
<td headers="Feppm" class="gt_row gt_right">8,770.00</td>
<td headers="Mn ppm" class="gt_row gt_right">165.94</td>
<td headers="Mo ppm" class="gt_row gt_right">5.6139999999999999</td>
<td headers="Na ppm" class="gt_row gt_right">335.89</td>
<td headers="Ni ppm" class="gt_row gt_right">15.06</td>
<td headers="Pb ppm" class="gt_row gt_right">7.0030000000000001</td>
<td headers="Zn ppm" class="gt_row gt_right">51.25</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-876-1</td>
<td headers="genotype" class="gt_row gt_center">BESC-876</td>
<td headers="block" class="gt_row gt_center">1</td>
<td headers="pCa" class="gt_row gt_right">0.88</td>
<td headers="pK" class="gt_row gt_right">0.37</td>
<td headers="pMg" class="gt_row gt_right">0.14</td>
<td headers="pP" class="gt_row gt_right">0.11</td>
<td headers="pC" class="gt_row gt_right">46.33</td>
<td headers="pN" class="gt_row gt_right">0.90</td>
<td headers="pS" class="gt_row gt_right">0.12</td>
<td headers="Lig" class="gt_row gt_right">25.67</td>
<td headers="Alppm" class="gt_row gt_right">5,301.36</td>
<td headers="Bppm" class="gt_row gt_right">16.35</td>
<td headers="Cdppm" class="gt_row gt_left">0.91200000000000003</td>
<td headers="Crppm" class="gt_row gt_right">93.37</td>
<td headers="Cuppm" class="gt_row gt_right">25.61</td>
<td headers="Feppm" class="gt_row gt_right">7,784.00</td>
<td headers="Mn ppm" class="gt_row gt_right">457.08</td>
<td headers="Mo ppm" class="gt_row gt_right">2.9540000000000002</td>
<td headers="Na ppm" class="gt_row gt_right">199.36</td>
<td headers="Ni ppm" class="gt_row gt_right">9.89</td>
<td headers="Pb ppm" class="gt_row gt_right">6.1509999999999998</td>
<td headers="Zn ppm" class="gt_row gt_right">56.24</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-876-2</td>
<td headers="genotype" class="gt_row gt_center">BESC-876</td>
<td headers="block" class="gt_row gt_center">2</td>
<td headers="pCa" class="gt_row gt_right">0.82</td>
<td headers="pK" class="gt_row gt_right">0.47</td>
<td headers="pMg" class="gt_row gt_right">0.15</td>
<td headers="pP" class="gt_row gt_right">0.10</td>
<td headers="pC" class="gt_row gt_right">50.21</td>
<td headers="pN" class="gt_row gt_right">0.74</td>
<td headers="pS" class="gt_row gt_right">0.11</td>
<td headers="Lig" class="gt_row gt_right">31.78</td>
<td headers="Alppm" class="gt_row gt_right">4,881.96</td>
<td headers="Bppm" class="gt_row gt_right">13.77</td>
<td headers="Cdppm" class="gt_row gt_left">&lt;0.80</td>
<td headers="Crppm" class="gt_row gt_right">66.27</td>
<td headers="Cuppm" class="gt_row gt_right">22.16</td>
<td headers="Feppm" class="gt_row gt_right">3,890.71</td>
<td headers="Mn ppm" class="gt_row gt_right">490.52</td>
<td headers="Mo ppm" class="gt_row gt_right">2.036</td>
<td headers="Na ppm" class="gt_row gt_right">260.95</td>
<td headers="Ni ppm" class="gt_row gt_right">9.88</td>
<td headers="Pb ppm" class="gt_row gt_right">6.226</td>
<td headers="Zn ppm" class="gt_row gt_right">78.03</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">GW-9776-1</td>
<td headers="genotype" class="gt_row gt_center">GW-9776</td>
<td headers="block" class="gt_row gt_center">1</td>
<td headers="pCa" class="gt_row gt_right">0.99</td>
<td headers="pK" class="gt_row gt_right">0.47</td>
<td headers="pMg" class="gt_row gt_right">0.15</td>
<td headers="pP" class="gt_row gt_right">0.10</td>
<td headers="pC" class="gt_row gt_right">48.77</td>
<td headers="pN" class="gt_row gt_right">0.81</td>
<td headers="pS" class="gt_row gt_right">0.11</td>
<td headers="Lig" class="gt_row gt_right">31.41</td>
<td headers="Alppm" class="gt_row gt_right">3,859.18</td>
<td headers="Bppm" class="gt_row gt_right">27.43</td>
<td headers="Cdppm" class="gt_row gt_left">1.0920000000000001</td>
<td headers="Crppm" class="gt_row gt_right">59.25</td>
<td headers="Cuppm" class="gt_row gt_right">27.30</td>
<td headers="Feppm" class="gt_row gt_right">6,618.00</td>
<td headers="Mn ppm" class="gt_row gt_right">703.30</td>
<td headers="Mo ppm" class="gt_row gt_right">1.6080000000000001</td>
<td headers="Na ppm" class="gt_row gt_right">255.57</td>
<td headers="Ni ppm" class="gt_row gt_right">10.76</td>
<td headers="Pb ppm" class="gt_row gt_right">3.262</td>
<td headers="Zn ppm" class="gt_row gt_right">75.29</td></tr>
  </tbody>
  &#10;  
</table>
</div>

</div>

Table 3: SEED soil core data collected in 2022.

</div>
<div id="tbl-seed-rc-key">

| Variable | Description | Units |
|----|----|----|
| pCa, pK, pMg, pP, pC, pN, pS | Root nutrient percentages | % |
| Lig | Root Lignin | % |
| Alppm, Bppm, Cdppm, Crppm, Cuppm, Mnppm, Mo ppm, Na ppm, Ni ppm, Pb ppm, Zn ppm | Root metals and nutrient concentrations | Parts per million |

Table 4: <a href="#tbl-seed-rc" class="quarto-xref">Table 3</a> key

</div>

# CBI Aboveground Biomass

A time series of abovegound biomass measurements from 2009-2019 was
collected from *Poplar Innovations* as part of Center for Bioenergy
Innovation (CBI) research. The primary measurements are tree diameters
that were converted to aboveground woody biomass (trunk and branches)
based on the allometric equations from @truax2014.

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

<div id="ovcjdeqdvj" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ovcjdeqdvj table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#ovcjdeqdvj thead, #ovcjdeqdvj tbody, #ovcjdeqdvj tfoot, #ovcjdeqdvj tr, #ovcjdeqdvj td, #ovcjdeqdvj th {
  border-style: none;
}
&#10;#ovcjdeqdvj p {
  margin: 0;
  padding: 0;
}
&#10;#ovcjdeqdvj .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#ovcjdeqdvj .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#ovcjdeqdvj .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#ovcjdeqdvj .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#ovcjdeqdvj .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#ovcjdeqdvj .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#ovcjdeqdvj .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#ovcjdeqdvj .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#ovcjdeqdvj .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#ovcjdeqdvj .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#ovcjdeqdvj .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#ovcjdeqdvj .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#ovcjdeqdvj .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#ovcjdeqdvj .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#ovcjdeqdvj .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ovcjdeqdvj .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#ovcjdeqdvj .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#ovcjdeqdvj .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#ovcjdeqdvj .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ovcjdeqdvj .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#ovcjdeqdvj .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ovcjdeqdvj .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#ovcjdeqdvj .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ovcjdeqdvj .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#ovcjdeqdvj .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ovcjdeqdvj .gt_left {
  text-align: left;
}
&#10;#ovcjdeqdvj .gt_center {
  text-align: center;
}
&#10;#ovcjdeqdvj .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#ovcjdeqdvj .gt_font_normal {
  font-weight: normal;
}
&#10;#ovcjdeqdvj .gt_font_bold {
  font-weight: bold;
}
&#10;#ovcjdeqdvj .gt_font_italic {
  font-style: italic;
}
&#10;#ovcjdeqdvj .gt_super {
  font-size: 65%;
}
&#10;#ovcjdeqdvj .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#ovcjdeqdvj .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#ovcjdeqdvj .gt_indent_1 {
  text-indent: 5px;
}
&#10;#ovcjdeqdvj .gt_indent_2 {
  text-indent: 10px;
}
&#10;#ovcjdeqdvj .gt_indent_3 {
  text-indent: 15px;
}
&#10;#ovcjdeqdvj .gt_indent_4 {
  text-indent: 20px;
}
&#10;#ovcjdeqdvj .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="tree">tree</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="cop_id">cop_id</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="height">height</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="d20">d20</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="d50">d50</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="dbh">dbh</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="agb">agb</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="tree" class="gt_row gt_center">BESC-102-1</td>
<td headers="cop_id" class="gt_row gt_right">1.00</td>
<td headers="height" class="gt_row gt_right">930.00</td>
<td headers="d20" class="gt_row gt_right">7.90</td>
<td headers="d50" class="gt_row gt_right">5.00</td>
<td headers="dbh" class="gt_row gt_right">20.70</td>
<td headers="agb" class="gt_row gt_right">10.75</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-102-2</td>
<td headers="cop_id" class="gt_row gt_right">0.00</td>
<td headers="height" class="gt_row gt_right">1,000.00</td>
<td headers="d20" class="gt_row gt_right">8.90</td>
<td headers="d50" class="gt_row gt_right">8.00</td>
<td headers="dbh" class="gt_row gt_right">18.30</td>
<td headers="agb" class="gt_row gt_right">7.73</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-102-3</td>
<td headers="cop_id" class="gt_row gt_right">0.00</td>
<td headers="height" class="gt_row gt_right">710.00</td>
<td headers="d20" class="gt_row gt_right">−Inf</td>
<td headers="d50" class="gt_row gt_right">−Inf</td>
<td headers="dbh" class="gt_row gt_right">24.80</td>
<td headers="agb" class="gt_row gt_right">16.06</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-119-1</td>
<td headers="cop_id" class="gt_row gt_right">1.00</td>
<td headers="height" class="gt_row gt_right">720.00</td>
<td headers="d20" class="gt_row gt_right">5.50</td>
<td headers="d50" class="gt_row gt_right">4.90</td>
<td headers="dbh" class="gt_row gt_right">16.80</td>
<td headers="agb" class="gt_row gt_right">6.69</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-119-2</td>
<td headers="cop_id" class="gt_row gt_right">0.00</td>
<td headers="height" class="gt_row gt_right">1,260.00</td>
<td headers="d20" class="gt_row gt_right">8.40</td>
<td headers="d50" class="gt_row gt_right">7.40</td>
<td headers="dbh" class="gt_row gt_right">19.10</td>
<td headers="agb" class="gt_row gt_right">8.57</td></tr>
    <tr><td headers="tree" class="gt_row gt_center">BESC-119-3</td>
<td headers="cop_id" class="gt_row gt_right">0.00</td>
<td headers="height" class="gt_row gt_right">790.00</td>
<td headers="d20" class="gt_row gt_right">−Inf</td>
<td headers="d50" class="gt_row gt_right">−Inf</td>
<td headers="dbh" class="gt_row gt_right">8.90</td>
<td headers="agb" class="gt_row gt_right">1.36</td></tr>
  </tbody>
  &#10;  
</table>
</div>

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
| agb | Estimated aboveground biomass growth rate from @truax2014 | kg/yr |

Table 6: <a href="#tbl-cbi-agb" class="quarto-xref">Table 5</a> key

</div>

# LDRD Soil Chemistry Data

Soil chemistry data for the 69 SEED soil cores were measured at UGA in
May 2024 in order to check whether the root chemistry influence on soil
C was an active or passive effect.

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

# Create Manuscript Analysis Dataset

Creating the overall data set for the Clatskanie soil C manuscript is
bit messy given we have to update multiple original soil core samples
with two differing updates (1/24/24 and 3/24/24) to the SOC
fractionation as well as add soil chemistry data taken on 5/24/24.

Furthermore, I have excluded Cadmium, Chromium, Molibdinum, and Lead
from the root and/or soil chemistry data since many of these were either
not measured for both root and soil or they had only trace amounts for
many of the samples.

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
  
  
  
         # select(-pH_root, -pC, -pN, -Cr_ppm, -Row, -Column, -MCwetCHEM,-MCdryCHEM,-MCwetBD,-MCdryBD)

write_csv(df_fit, "./02-data/02-processed/clatskanie-c-fit-data.csv")
display_table(df_fit)
```

</details>
<div id="vquauxglkj" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#vquauxglkj table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#vquauxglkj thead, #vquauxglkj tbody, #vquauxglkj tfoot, #vquauxglkj tr, #vquauxglkj td, #vquauxglkj th {
  border-style: none;
}
&#10;#vquauxglkj p {
  margin: 0;
  padding: 0;
}
&#10;#vquauxglkj .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#vquauxglkj .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#vquauxglkj .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#vquauxglkj .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#vquauxglkj .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#vquauxglkj .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#vquauxglkj .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#vquauxglkj .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#vquauxglkj .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#vquauxglkj .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#vquauxglkj .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#vquauxglkj .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#vquauxglkj .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#vquauxglkj .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#vquauxglkj .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#vquauxglkj .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#vquauxglkj .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#vquauxglkj .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#vquauxglkj .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#vquauxglkj .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#vquauxglkj .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#vquauxglkj .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#vquauxglkj .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#vquauxglkj .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#vquauxglkj .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#vquauxglkj .gt_left {
  text-align: left;
}
&#10;#vquauxglkj .gt_center {
  text-align: center;
}
&#10;#vquauxglkj .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#vquauxglkj .gt_font_normal {
  font-weight: normal;
}
&#10;#vquauxglkj .gt_font_bold {
  font-weight: bold;
}
&#10;#vquauxglkj .gt_font_italic {
  font-style: italic;
}
&#10;#vquauxglkj .gt_super {
  font-size: 65%;
}
&#10;#vquauxglkj .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#vquauxglkj .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#vquauxglkj .gt_indent_1 {
  text-indent: 5px;
}
&#10;#vquauxglkj .gt_indent_2 {
  text-indent: 10px;
}
&#10;#vquauxglkj .gt_indent_3 {
  text-indent: 15px;
}
&#10;#vquauxglkj .gt_indent_4 {
  text-indent: 20px;
}
&#10;#vquauxglkj .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ord">ord</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="tree">tree</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="genotype">genotype</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="block">block</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Depth">Depth</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="latitude">latitude</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="longitude">longitude</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="cop_id">cop_id</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="height">height</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="dbh">dbh</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="agb">agb</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="root_CN">root_CN</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Lig">Lig</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="BD">BD</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Sand">Sand</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Clay">Clay</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Silt">Silt</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="SiCl">SiCl</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pH_soil">pH_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="CEC">CEC</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="base_sat">base_sat</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="POMC">POMC</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="MAOMC">MAOMC</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="TotC">TotC</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="TotC_30">TotC_30</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="POMC_st">POMC_st</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="MAOMC_st">MAOMC_st</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="TotC_st">TotC_st</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="TotC_st_30">TotC_st_30</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="POM_CN">POM_CN</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="MAOM_CN">MAOM_CN</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="soil_CN">soil_CN</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="C_chk_rel">C_chk_rel</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="C_chk_id">C_chk_id</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Al_ppm_root">Al_ppm_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Al_ppm_soil">Al_ppm_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Al_ppm_diff">Al_ppm_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="B_ppm_root">B_ppm_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Cu_ppm_root">Cu_ppm_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Cu_ppm_soil">Cu_ppm_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Cu_ppm_diff">Cu_ppm_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Fe_ppm_root">Fe_ppm_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Fe_ppm_soil">Fe_ppm_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Fe_ppm_diff">Fe_ppm_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Mn_ppm_root">Mn_ppm_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Mn_ppm_soil">Mn_ppm_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Mn_ppm_diff">Mn_ppm_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Na_ppm_root">Na_ppm_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Na_ppm_soil">Na_ppm_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Na_ppm_diff">Na_ppm_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Ni_ppm_root">Ni_ppm_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Ni_ppm_soil">Ni_ppm_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Ni_ppm_diff">Ni_ppm_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Zn_ppm_root">Zn_ppm_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Zn_ppm_soil">Zn_ppm_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Zn_ppm_diff">Zn_ppm_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pCa_root">pCa_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pCa_soil">pCa_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pCa_diff">pCa_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pK_root">pK_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pK_soil">pK_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pK_diff">pK_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pMg_root">pMg_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pMg_soil">pMg_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pMg_diff">pMg_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pP_root">pP_root</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Cr_ppm">Cr_ppm</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pPOM">pPOM</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pP_soil">pP_soil</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pP_diff">pP_diff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pS_root">pS_root</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="ord" class="gt_row gt_right">1.00</td>
<td headers="tree" class="gt_row gt_center">BESC-102-1</td>
<td headers="genotype" class="gt_row gt_center">BESC-102</td>
<td headers="block" class="gt_row gt_center">1</td>
<td headers="Depth" class="gt_row gt_right">15.00</td>
<td headers="latitude" class="gt_row gt_right">46.12</td>
<td headers="longitude" class="gt_row gt_right">−123.27</td>
<td headers="cop_id" class="gt_row gt_center">Coppiced</td>
<td headers="height" class="gt_row gt_right">930.00</td>
<td headers="dbh" class="gt_row gt_right">20.70</td>
<td headers="agb" class="gt_row gt_right">5.37</td>
<td headers="root_CN" class="gt_row gt_right">49.73</td>
<td headers="Lig" class="gt_row gt_right">26.58</td>
<td headers="BD" class="gt_row gt_right">0.88</td>
<td headers="Sand" class="gt_row gt_right">0.00</td>
<td headers="Clay" class="gt_row gt_right">45.90</td>
<td headers="Silt" class="gt_row gt_right">55.75</td>
<td headers="SiCl" class="gt_row gt_right">101.65</td>
<td headers="pH_soil" class="gt_row gt_right">5.26</td>
<td headers="CEC" class="gt_row gt_right">32.23</td>
<td headers="base_sat" class="gt_row gt_right">42.16</td>
<td headers="POMC" class="gt_row gt_right">11.56</td>
<td headers="MAOMC" class="gt_row gt_right">39.52</td>
<td headers="TotC" class="gt_row gt_right">60.40</td>
<td headers="TotC_30" class="gt_row gt_right">44.22</td>
<td headers="POMC_st" class="gt_row gt_right">15.33</td>
<td headers="MAOMC_st" class="gt_row gt_right">52.40</td>
<td headers="TotC_st" class="gt_row gt_right">80.09</td>
<td headers="TotC_st_30" class="gt_row gt_right">99.33</td>
<td headers="POM_CN" class="gt_row gt_right">21.81</td>
<td headers="MAOM_CN" class="gt_row gt_right">12.96</td>
<td headers="soil_CN" class="gt_row gt_right">14.88</td>
<td headers="C_chk_rel" class="gt_row gt_right">−15.43</td>
<td headers="C_chk_id" class="gt_row gt_center">Good</td>
<td headers="Al_ppm_root" class="gt_row gt_right">7,378.17</td>
<td headers="Al_ppm_soil" class="gt_row gt_right">4.47</td>
<td headers="Al_ppm_diff" class="gt_row gt_right">7,373.70</td>
<td headers="B_ppm_root" class="gt_row gt_right">16.35</td>
<td headers="Cu_ppm_root" class="gt_row gt_right">15.17</td>
<td headers="Cu_ppm_soil" class="gt_row gt_right">4.10</td>
<td headers="Cu_ppm_diff" class="gt_row gt_right">11.07</td>
<td headers="Fe_ppm_root" class="gt_row gt_right">6,213.76</td>
<td headers="Fe_ppm_soil" class="gt_row gt_right">267.51</td>
<td headers="Fe_ppm_diff" class="gt_row gt_right">5,946.25</td>
<td headers="Mn_ppm_root" class="gt_row gt_right">152.44</td>
<td headers="Mn_ppm_soil" class="gt_row gt_right">67.24</td>
<td headers="Mn_ppm_diff" class="gt_row gt_right">85.20</td>
<td headers="Na_ppm_root" class="gt_row gt_right">285.62</td>
<td headers="Na_ppm_soil" class="gt_row gt_right">37.81</td>
<td headers="Na_ppm_diff" class="gt_row gt_right">247.82</td>
<td headers="Ni_ppm_root" class="gt_row gt_right">7.07</td>
<td headers="Ni_ppm_soil" class="gt_row gt_right">2.11</td>
<td headers="Ni_ppm_diff" class="gt_row gt_right">4.96</td>
<td headers="Zn_ppm_root" class="gt_row gt_right">55.16</td>
<td headers="Zn_ppm_soil" class="gt_row gt_right">3.30</td>
<td headers="Zn_ppm_diff" class="gt_row gt_right">51.86</td>
<td headers="pCa_root" class="gt_row gt_right">0.72</td>
<td headers="pCa_soil" class="gt_row gt_right">0.21</td>
<td headers="pCa_diff" class="gt_row gt_right">0.51</td>
<td headers="pK_root" class="gt_row gt_right">0.44</td>
<td headers="pK_soil" class="gt_row gt_right">0.02</td>
<td headers="pK_diff" class="gt_row gt_right">0.42</td>
<td headers="pMg_root" class="gt_row gt_right">0.16</td>
<td headers="pMg_soil" class="gt_row gt_right">0.03</td>
<td headers="pMg_diff" class="gt_row gt_right">0.13</td>
<td headers="pP_root" class="gt_row gt_right">0.08</td>
<td headers="Cr_ppm" class="gt_row gt_right">73.56</td>
<td headers="pPOM" class="gt_row gt_right">19.14</td>
<td headers="pP_soil" class="gt_row gt_right">0.01</td>
<td headers="pP_diff" class="gt_row gt_right">0.07</td>
<td headers="pS_root" class="gt_row gt_right">0.08</td></tr>
    <tr><td headers="ord" class="gt_row gt_right">2.00</td>
<td headers="tree" class="gt_row gt_center">BESC-102-2</td>
<td headers="genotype" class="gt_row gt_center">BESC-102</td>
<td headers="block" class="gt_row gt_center">2</td>
<td headers="Depth" class="gt_row gt_right">15.00</td>
<td headers="latitude" class="gt_row gt_right">46.12</td>
<td headers="longitude" class="gt_row gt_right">−123.27</td>
<td headers="cop_id" class="gt_row gt_center">Non-coppiced</td>
<td headers="height" class="gt_row gt_right">1,000.00</td>
<td headers="dbh" class="gt_row gt_right">18.30</td>
<td headers="agb" class="gt_row gt_right">3.86</td>
<td headers="root_CN" class="gt_row gt_right">50.42</td>
<td headers="Lig" class="gt_row gt_right">32.71</td>
<td headers="BD" class="gt_row gt_right">0.62</td>
<td headers="Sand" class="gt_row gt_right">7.85</td>
<td headers="Clay" class="gt_row gt_right">45.15</td>
<td headers="Silt" class="gt_row gt_right">47.00</td>
<td headers="SiCl" class="gt_row gt_right">92.15</td>
<td headers="pH_soil" class="gt_row gt_right">5.27</td>
<td headers="CEC" class="gt_row gt_right">32.70</td>
<td headers="base_sat" class="gt_row gt_right">47.78</td>
<td headers="POMC" class="gt_row gt_right">9.37</td>
<td headers="MAOMC" class="gt_row gt_right">52.97</td>
<td headers="TotC" class="gt_row gt_right">76.37</td>
<td headers="TotC_30" class="gt_row gt_right">70.86</td>
<td headers="POMC_st" class="gt_row gt_right">8.70</td>
<td headers="MAOMC_st" class="gt_row gt_right">49.21</td>
<td headers="TotC_st" class="gt_row gt_right">70.95</td>
<td headers="TotC_st_30" class="gt_row gt_right">132.22</td>
<td headers="POM_CN" class="gt_row gt_right">20.34</td>
<td headers="MAOM_CN" class="gt_row gt_right">13.83</td>
<td headers="soil_CN" class="gt_row gt_right">14.60</td>
<td headers="C_chk_rel" class="gt_row gt_right">−18.37</td>
<td headers="C_chk_id" class="gt_row gt_center">Good</td>
<td headers="Al_ppm_root" class="gt_row gt_right">5,044.21</td>
<td headers="Al_ppm_soil" class="gt_row gt_right">2.92</td>
<td headers="Al_ppm_diff" class="gt_row gt_right">5,041.29</td>
<td headers="B_ppm_root" class="gt_row gt_right">11.64</td>
<td headers="Cu_ppm_root" class="gt_row gt_right">21.75</td>
<td headers="Cu_ppm_soil" class="gt_row gt_right">6.08</td>
<td headers="Cu_ppm_diff" class="gt_row gt_right">15.67</td>
<td headers="Fe_ppm_root" class="gt_row gt_right">3,243.34</td>
<td headers="Fe_ppm_soil" class="gt_row gt_right">363.77</td>
<td headers="Fe_ppm_diff" class="gt_row gt_right">2,879.57</td>
<td headers="Mn_ppm_root" class="gt_row gt_right">283.39</td>
<td headers="Mn_ppm_soil" class="gt_row gt_right">59.00</td>
<td headers="Mn_ppm_diff" class="gt_row gt_right">224.39</td>
<td headers="Na_ppm_root" class="gt_row gt_right">307.61</td>
<td headers="Na_ppm_soil" class="gt_row gt_right">53.09</td>
<td headers="Na_ppm_diff" class="gt_row gt_right">254.52</td>
<td headers="Ni_ppm_root" class="gt_row gt_right">9.60</td>
<td headers="Ni_ppm_soil" class="gt_row gt_right">3.03</td>
<td headers="Ni_ppm_diff" class="gt_row gt_right">6.56</td>
<td headers="Zn_ppm_root" class="gt_row gt_right">73.66</td>
<td headers="Zn_ppm_soil" class="gt_row gt_right">7.04</td>
<td headers="Zn_ppm_diff" class="gt_row gt_right">66.62</td>
<td headers="pCa_root" class="gt_row gt_right">0.82</td>
<td headers="pCa_soil" class="gt_row gt_right">0.24</td>
<td headers="pCa_diff" class="gt_row gt_right">0.58</td>
<td headers="pK_root" class="gt_row gt_right">0.36</td>
<td headers="pK_soil" class="gt_row gt_right">0.01</td>
<td headers="pK_diff" class="gt_row gt_right">0.35</td>
<td headers="pMg_root" class="gt_row gt_right">0.13</td>
<td headers="pMg_soil" class="gt_row gt_right">0.04</td>
<td headers="pMg_diff" class="gt_row gt_right">0.09</td>
<td headers="pP_root" class="gt_row gt_right">0.11</td>
<td headers="Cr_ppm" class="gt_row gt_right">54.67</td>
<td headers="pPOM" class="gt_row gt_right">12.26</td>
<td headers="pP_soil" class="gt_row gt_right">0.01</td>
<td headers="pP_diff" class="gt_row gt_right">0.09</td>
<td headers="pS_root" class="gt_row gt_right">0.11</td></tr>
    <tr><td headers="ord" class="gt_row gt_right">3.00</td>
<td headers="tree" class="gt_row gt_center">BESC-102-3</td>
<td headers="genotype" class="gt_row gt_center">BESC-102</td>
<td headers="block" class="gt_row gt_center">3</td>
<td headers="Depth" class="gt_row gt_right">15.00</td>
<td headers="latitude" class="gt_row gt_right">46.12</td>
<td headers="longitude" class="gt_row gt_right">−123.27</td>
<td headers="cop_id" class="gt_row gt_center">Non-coppiced</td>
<td headers="height" class="gt_row gt_right">710.00</td>
<td headers="dbh" class="gt_row gt_right">24.80</td>
<td headers="agb" class="gt_row gt_right">8.03</td>
<td headers="root_CN" class="gt_row gt_right">39.81</td>
<td headers="Lig" class="gt_row gt_right">24.03</td>
<td headers="BD" class="gt_row gt_right">0.79</td>
<td headers="Sand" class="gt_row gt_right">8.39</td>
<td headers="Clay" class="gt_row gt_right">38.21</td>
<td headers="Silt" class="gt_row gt_right">51.20</td>
<td headers="SiCl" class="gt_row gt_right">89.41</td>
<td headers="pH_soil" class="gt_row gt_right">5.48</td>
<td headers="CEC" class="gt_row gt_right">28.92</td>
<td headers="base_sat" class="gt_row gt_right">73.51</td>
<td headers="POMC" class="gt_row gt_right">7.59</td>
<td headers="MAOMC" class="gt_row gt_right">51.99</td>
<td headers="TotC" class="gt_row gt_right">53.50</td>
<td headers="TotC_30" class="gt_row gt_right">42.85</td>
<td headers="POMC_st" class="gt_row gt_right">9.01</td>
<td headers="MAOMC_st" class="gt_row gt_right">61.69</td>
<td headers="TotC_st" class="gt_row gt_right">63.48</td>
<td headers="TotC_st_30" class="gt_row gt_right">104.23</td>
<td headers="POM_CN" class="gt_row gt_right">20.51</td>
<td headers="MAOM_CN" class="gt_row gt_right">14.13</td>
<td headers="soil_CN" class="gt_row gt_right">13.90</td>
<td headers="C_chk_rel" class="gt_row gt_right">11.36</td>
<td headers="C_chk_id" class="gt_row gt_center">Good</td>
<td headers="Al_ppm_root" class="gt_row gt_right">8,403.85</td>
<td headers="Al_ppm_soil" class="gt_row gt_right">0.14</td>
<td headers="Al_ppm_diff" class="gt_row gt_right">8,403.71</td>
<td headers="B_ppm_root" class="gt_row gt_right">16.98</td>
<td headers="Cu_ppm_root" class="gt_row gt_right">23.01</td>
<td headers="Cu_ppm_soil" class="gt_row gt_right">5.05</td>
<td headers="Cu_ppm_diff" class="gt_row gt_right">17.95</td>
<td headers="Fe_ppm_root" class="gt_row gt_right">8,770.00</td>
<td headers="Fe_ppm_soil" class="gt_row gt_right">215.49</td>
<td headers="Fe_ppm_diff" class="gt_row gt_right">8,554.51</td>
<td headers="Mn_ppm_root" class="gt_row gt_right">165.94</td>
<td headers="Mn_ppm_soil" class="gt_row gt_right">93.89</td>
<td headers="Mn_ppm_diff" class="gt_row gt_right">72.05</td>
<td headers="Na_ppm_root" class="gt_row gt_right">335.89</td>
<td headers="Na_ppm_soil" class="gt_row gt_right">43.53</td>
<td headers="Na_ppm_diff" class="gt_row gt_right">292.36</td>
<td headers="Ni_ppm_root" class="gt_row gt_right">15.06</td>
<td headers="Ni_ppm_soil" class="gt_row gt_right">2.82</td>
<td headers="Ni_ppm_diff" class="gt_row gt_right">12.24</td>
<td headers="Zn_ppm_root" class="gt_row gt_right">51.25</td>
<td headers="Zn_ppm_soil" class="gt_row gt_right">9.72</td>
<td headers="Zn_ppm_diff" class="gt_row gt_right">41.53</td>
<td headers="pCa_root" class="gt_row gt_right">0.76</td>
<td headers="pCa_soil" class="gt_row gt_right">0.33</td>
<td headers="pCa_diff" class="gt_row gt_right">0.43</td>
<td headers="pK_root" class="gt_row gt_right">0.61</td>
<td headers="pK_soil" class="gt_row gt_right">0.02</td>
<td headers="pK_diff" class="gt_row gt_right">0.59</td>
<td headers="pMg_root" class="gt_row gt_right">0.20</td>
<td headers="pMg_soil" class="gt_row gt_right">0.05</td>
<td headers="pMg_diff" class="gt_row gt_right">0.15</td>
<td headers="pP_root" class="gt_row gt_right">0.14</td>
<td headers="Cr_ppm" class="gt_row gt_right">256.51</td>
<td headers="pPOM" class="gt_row gt_right">14.19</td>
<td headers="pP_soil" class="gt_row gt_right">0.02</td>
<td headers="pP_diff" class="gt_row gt_right">0.13</td>
<td headers="pS_root" class="gt_row gt_right">0.11</td></tr>
    <tr><td headers="ord" class="gt_row gt_right">4.00</td>
<td headers="tree" class="gt_row gt_center">BESC-119-1</td>
<td headers="genotype" class="gt_row gt_center">BESC-119</td>
<td headers="block" class="gt_row gt_center">1</td>
<td headers="Depth" class="gt_row gt_right">15.00</td>
<td headers="latitude" class="gt_row gt_right">46.12</td>
<td headers="longitude" class="gt_row gt_right">−123.27</td>
<td headers="cop_id" class="gt_row gt_center">Coppiced</td>
<td headers="height" class="gt_row gt_right">720.00</td>
<td headers="dbh" class="gt_row gt_right">16.80</td>
<td headers="agb" class="gt_row gt_right">3.34</td>
<td headers="root_CN" class="gt_row gt_right">51.94</td>
<td headers="Lig" class="gt_row gt_right">35.71</td>
<td headers="BD" class="gt_row gt_right">0.82</td>
<td headers="Sand" class="gt_row gt_right">0.00</td>
<td headers="Clay" class="gt_row gt_right">37.90</td>
<td headers="Silt" class="gt_row gt_right">62.87</td>
<td headers="SiCl" class="gt_row gt_right">100.77</td>
<td headers="pH_soil" class="gt_row gt_right">4.78</td>
<td headers="CEC" class="gt_row gt_right">30.19</td>
<td headers="base_sat" class="gt_row gt_right">43.10</td>
<td headers="POMC" class="gt_row gt_right">7.07</td>
<td headers="MAOMC" class="gt_row gt_right">32.04</td>
<td headers="TotC" class="gt_row gt_right">42.43</td>
<td headers="TotC_30" class="gt_row gt_right">29.88</td>
<td headers="POMC_st" class="gt_row gt_right">8.66</td>
<td headers="MAOMC_st" class="gt_row gt_right">39.26</td>
<td headers="TotC_st" class="gt_row gt_right">52.00</td>
<td headers="TotC_st_30" class="gt_row gt_right">72.04</td>
<td headers="POM_CN" class="gt_row gt_right">20.79</td>
<td headers="MAOM_CN" class="gt_row gt_right">12.09</td>
<td headers="soil_CN" class="gt_row gt_right">12.78</td>
<td headers="C_chk_rel" class="gt_row gt_right">−7.83</td>
<td headers="C_chk_id" class="gt_row gt_center">Good</td>
<td headers="Al_ppm_root" class="gt_row gt_right">2,109.76</td>
<td headers="Al_ppm_soil" class="gt_row gt_right">4.64</td>
<td headers="Al_ppm_diff" class="gt_row gt_right">2,105.12</td>
<td headers="B_ppm_root" class="gt_row gt_right">11.84</td>
<td headers="Cu_ppm_root" class="gt_row gt_right">24.69</td>
<td headers="Cu_ppm_soil" class="gt_row gt_right">4.75</td>
<td headers="Cu_ppm_diff" class="gt_row gt_right">19.95</td>
<td headers="Fe_ppm_root" class="gt_row gt_right">1,833.12</td>
<td headers="Fe_ppm_soil" class="gt_row gt_right">411.20</td>
<td headers="Fe_ppm_diff" class="gt_row gt_right">1,421.92</td>
<td headers="Mn_ppm_root" class="gt_row gt_right">149.32</td>
<td headers="Mn_ppm_soil" class="gt_row gt_right">57.42</td>
<td headers="Mn_ppm_diff" class="gt_row gt_right">91.90</td>
<td headers="Na_ppm_root" class="gt_row gt_right">107.56</td>
<td headers="Na_ppm_soil" class="gt_row gt_right">34.50</td>
<td headers="Na_ppm_diff" class="gt_row gt_right">73.06</td>
<td headers="Ni_ppm_root" class="gt_row gt_right">4.88</td>
<td headers="Ni_ppm_soil" class="gt_row gt_right">1.88</td>
<td headers="Ni_ppm_diff" class="gt_row gt_right">3.00</td>
<td headers="Zn_ppm_root" class="gt_row gt_right">57.83</td>
<td headers="Zn_ppm_soil" class="gt_row gt_right">4.72</td>
<td headers="Zn_ppm_diff" class="gt_row gt_right">53.12</td>
<td headers="pCa_root" class="gt_row gt_right">0.90</td>
<td headers="pCa_soil" class="gt_row gt_right">0.20</td>
<td headers="pCa_diff" class="gt_row gt_right">0.71</td>
<td headers="pK_root" class="gt_row gt_right">0.72</td>
<td headers="pK_soil" class="gt_row gt_right">0.03</td>
<td headers="pK_diff" class="gt_row gt_right">0.68</td>
<td headers="pMg_root" class="gt_row gt_right">0.10</td>
<td headers="pMg_soil" class="gt_row gt_right">0.03</td>
<td headers="pMg_diff" class="gt_row gt_right">0.07</td>
<td headers="pP_root" class="gt_row gt_right">0.12</td>
<td headers="Cr_ppm" class="gt_row gt_right">14.03</td>
<td headers="pPOM" class="gt_row gt_right">16.66</td>
<td headers="pP_soil" class="gt_row gt_right">0.01</td>
<td headers="pP_diff" class="gt_row gt_right">0.10</td>
<td headers="pS_root" class="gt_row gt_right">0.12</td></tr>
    <tr><td headers="ord" class="gt_row gt_right">5.00</td>
<td headers="tree" class="gt_row gt_center">BESC-13-1</td>
<td headers="genotype" class="gt_row gt_center">BESC-13</td>
<td headers="block" class="gt_row gt_center">1</td>
<td headers="Depth" class="gt_row gt_right">15.00</td>
<td headers="latitude" class="gt_row gt_right">46.12</td>
<td headers="longitude" class="gt_row gt_right">−123.27</td>
<td headers="cop_id" class="gt_row gt_center">Non-coppiced</td>
<td headers="height" class="gt_row gt_right">670.00</td>
<td headers="dbh" class="gt_row gt_right">12.60</td>
<td headers="agb" class="gt_row gt_right">1.57</td>
<td headers="root_CN" class="gt_row gt_right">40.69</td>
<td headers="Lig" class="gt_row gt_right">27.05</td>
<td headers="BD" class="gt_row gt_right">0.99</td>
<td headers="Sand" class="gt_row gt_right">1.89</td>
<td headers="Clay" class="gt_row gt_right">34.90</td>
<td headers="Silt" class="gt_row gt_right">63.21</td>
<td headers="SiCl" class="gt_row gt_right">98.11</td>
<td headers="pH_soil" class="gt_row gt_right">5.25</td>
<td headers="CEC" class="gt_row gt_right">25.40</td>
<td headers="base_sat" class="gt_row gt_right">67.60</td>
<td headers="POMC" class="gt_row gt_right">2.70</td>
<td headers="MAOMC" class="gt_row gt_right">21.61</td>
<td headers="TotC" class="gt_row gt_right">24.70</td>
<td headers="TotC_30" class="gt_row gt_right">24.70</td>
<td headers="POMC_st" class="gt_row gt_right">4.01</td>
<td headers="MAOMC_st" class="gt_row gt_right">32.12</td>
<td headers="TotC_st" class="gt_row gt_right">36.72</td>
<td headers="TotC_st_30" class="gt_row gt_right">69.09</td>
<td headers="POM_CN" class="gt_row gt_right">11.93</td>
<td headers="MAOM_CN" class="gt_row gt_right">11.03</td>
<td headers="soil_CN" class="gt_row gt_right">9.32</td>
<td headers="C_chk_rel" class="gt_row gt_right">−1.59</td>
<td headers="C_chk_id" class="gt_row gt_center">Good</td>
<td headers="Al_ppm_root" class="gt_row gt_right">1,047.54</td>
<td headers="Al_ppm_soil" class="gt_row gt_right">0.41</td>
<td headers="Al_ppm_diff" class="gt_row gt_right">1,047.13</td>
<td headers="B_ppm_root" class="gt_row gt_right">15.29</td>
<td headers="Cu_ppm_root" class="gt_row gt_right">10.38</td>
<td headers="Cu_ppm_soil" class="gt_row gt_right">5.76</td>
<td headers="Cu_ppm_diff" class="gt_row gt_right">4.62</td>
<td headers="Fe_ppm_root" class="gt_row gt_right">952.77</td>
<td headers="Fe_ppm_soil" class="gt_row gt_right">350.07</td>
<td headers="Fe_ppm_diff" class="gt_row gt_right">602.70</td>
<td headers="Mn_ppm_root" class="gt_row gt_right">73.41</td>
<td headers="Mn_ppm_soil" class="gt_row gt_right">74.21</td>
<td headers="Mn_ppm_diff" class="gt_row gt_right">−0.80</td>
<td headers="Na_ppm_root" class="gt_row gt_right">372.08</td>
<td headers="Na_ppm_soil" class="gt_row gt_right">33.46</td>
<td headers="Na_ppm_diff" class="gt_row gt_right">338.62</td>
<td headers="Ni_ppm_root" class="gt_row gt_right">3.91</td>
<td headers="Ni_ppm_soil" class="gt_row gt_right">2.85</td>
<td headers="Ni_ppm_diff" class="gt_row gt_right">1.06</td>
<td headers="Zn_ppm_root" class="gt_row gt_right">65.39</td>
<td headers="Zn_ppm_soil" class="gt_row gt_right">7.45</td>
<td headers="Zn_ppm_diff" class="gt_row gt_right">57.94</td>
<td headers="pCa_root" class="gt_row gt_right">0.54</td>
<td headers="pCa_soil" class="gt_row gt_right">0.26</td>
<td headers="pCa_diff" class="gt_row gt_right">0.28</td>
<td headers="pK_root" class="gt_row gt_right">0.83</td>
<td headers="pK_soil" class="gt_row gt_right">0.03</td>
<td headers="pK_diff" class="gt_row gt_right">0.80</td>
<td headers="pMg_root" class="gt_row gt_right">0.10</td>
<td headers="pMg_soil" class="gt_row gt_right">0.04</td>
<td headers="pMg_diff" class="gt_row gt_right">0.06</td>
<td headers="pP_root" class="gt_row gt_right">0.12</td>
<td headers="Cr_ppm" class="gt_row gt_right">NA</td>
<td headers="pPOM" class="gt_row gt_right">10.92</td>
<td headers="pP_soil" class="gt_row gt_right">0.02</td>
<td headers="pP_diff" class="gt_row gt_right">0.10</td>
<td headers="pS_root" class="gt_row gt_right">0.11</td></tr>
    <tr><td headers="ord" class="gt_row gt_right">6.00</td>
<td headers="tree" class="gt_row gt_center">BESC-13-2</td>
<td headers="genotype" class="gt_row gt_center">BESC-13</td>
<td headers="block" class="gt_row gt_center">2</td>
<td headers="Depth" class="gt_row gt_right">15.00</td>
<td headers="latitude" class="gt_row gt_right">46.12</td>
<td headers="longitude" class="gt_row gt_right">−123.27</td>
<td headers="cop_id" class="gt_row gt_center">Non-coppiced</td>
<td headers="height" class="gt_row gt_right">1,150.00</td>
<td headers="dbh" class="gt_row gt_right">23.00</td>
<td headers="agb" class="gt_row gt_right">6.70</td>
<td headers="root_CN" class="gt_row gt_right">52.09</td>
<td headers="Lig" class="gt_row gt_right">31.56</td>
<td headers="BD" class="gt_row gt_right">0.84</td>
<td headers="Sand" class="gt_row gt_right">8.27</td>
<td headers="Clay" class="gt_row gt_right">36.83</td>
<td headers="Silt" class="gt_row gt_right">54.90</td>
<td headers="SiCl" class="gt_row gt_right">91.73</td>
<td headers="pH_soil" class="gt_row gt_right">5.27</td>
<td headers="CEC" class="gt_row gt_right">26.33</td>
<td headers="base_sat" class="gt_row gt_right">71.23</td>
<td headers="POMC" class="gt_row gt_right">14.70</td>
<td headers="MAOMC" class="gt_row gt_right">28.90</td>
<td headers="TotC" class="gt_row gt_right">49.10</td>
<td headers="TotC_30" class="gt_row gt_right">49.10</td>
<td headers="POMC_st" class="gt_row gt_right">18.55</td>
<td headers="MAOMC_st" class="gt_row gt_right">36.50</td>
<td headers="TotC_st" class="gt_row gt_right">61.99</td>
<td headers="TotC_st_30" class="gt_row gt_right">114.76</td>
<td headers="POM_CN" class="gt_row gt_right">23.83</td>
<td headers="MAOM_CN" class="gt_row gt_right">10.84</td>
<td headers="soil_CN" class="gt_row gt_right">13.18</td>
<td headers="C_chk_rel" class="gt_row gt_right">−11.20</td>
<td headers="C_chk_id" class="gt_row gt_center">Good</td>
<td headers="Al_ppm_root" class="gt_row gt_right">2,326.70</td>
<td headers="Al_ppm_soil" class="gt_row gt_right">0.53</td>
<td headers="Al_ppm_diff" class="gt_row gt_right">2,326.17</td>
<td headers="B_ppm_root" class="gt_row gt_right">14.38</td>
<td headers="Cu_ppm_root" class="gt_row gt_right">17.18</td>
<td headers="Cu_ppm_soil" class="gt_row gt_right">6.37</td>
<td headers="Cu_ppm_diff" class="gt_row gt_right">10.81</td>
<td headers="Fe_ppm_root" class="gt_row gt_right">2,202.72</td>
<td headers="Fe_ppm_soil" class="gt_row gt_right">447.65</td>
<td headers="Fe_ppm_diff" class="gt_row gt_right">1,755.07</td>
<td headers="Mn_ppm_root" class="gt_row gt_right">243.11</td>
<td headers="Mn_ppm_soil" class="gt_row gt_right">84.15</td>
<td headers="Mn_ppm_diff" class="gt_row gt_right">158.96</td>
<td headers="Na_ppm_root" class="gt_row gt_right">102.42</td>
<td headers="Na_ppm_soil" class="gt_row gt_right">38.25</td>
<td headers="Na_ppm_diff" class="gt_row gt_right">64.18</td>
<td headers="Ni_ppm_root" class="gt_row gt_right">8.19</td>
<td headers="Ni_ppm_soil" class="gt_row gt_right">3.13</td>
<td headers="Ni_ppm_diff" class="gt_row gt_right">5.07</td>
<td headers="Zn_ppm_root" class="gt_row gt_right">93.19</td>
<td headers="Zn_ppm_soil" class="gt_row gt_right">9.29</td>
<td headers="Zn_ppm_diff" class="gt_row gt_right">83.90</td>
<td headers="pCa_root" class="gt_row gt_right">0.70</td>
<td headers="pCa_soil" class="gt_row gt_right">0.30</td>
<td headers="pCa_diff" class="gt_row gt_right">0.41</td>
<td headers="pK_root" class="gt_row gt_right">0.48</td>
<td headers="pK_soil" class="gt_row gt_right">0.03</td>
<td headers="pK_diff" class="gt_row gt_right">0.45</td>
<td headers="pMg_root" class="gt_row gt_right">0.12</td>
<td headers="pMg_soil" class="gt_row gt_right">0.04</td>
<td headers="pMg_diff" class="gt_row gt_right">0.08</td>
<td headers="pP_root" class="gt_row gt_right">0.14</td>
<td headers="Cr_ppm" class="gt_row gt_right">NA</td>
<td headers="pPOM" class="gt_row gt_right">29.93</td>
<td headers="pP_soil" class="gt_row gt_right">0.01</td>
<td headers="pP_diff" class="gt_row gt_right">0.12</td>
<td headers="pS_root" class="gt_row gt_right">0.09</td></tr>
  </tbody>
  &#10;  
</table>
</div>

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
| agb | Estimated aboveground biomass growth rate from @truax2014 based on dbh | kg/yr |
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

Additionally, we create supplemental Figures S4 and S6 used in the
manuscript below.

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
