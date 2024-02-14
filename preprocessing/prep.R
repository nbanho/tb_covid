#### Libraries ####

library(tidyverse)


#### Raw data ####

raw_df <- read.csv("data-raw/TB_COVID19_survey_FUP_2023_cleaned.csv") %>%
  filter(country_ap != "Japan")

#### Site characteristics ####

site <- raw_df %>%
  mutate(
    site_region = factor(
      region,
      levels = c(
        "Asia-Pacific",
        "Central Africa",
        "East Africa",
        "Southern Africa",
        "West Africa"
      )
    )
  ) %>%
  mutate(
    site_country = trimws(paste0(
      country_ap,
      country_ca,
      country_ea,
      country_sa,
      country_wa
    )),
    site_country = factor(site_country, levels = sort(unique(site_country)))
  ) %>%
  mutate(
    site_setting = ifelse(grepl("Peri", clinic_setting),
      "peri-urban", clinic_setting
    ),
    site_setting = tolower(site_setting),
    site_setting = factor(
      site_setting,
      levels = c("urban", "peri-urban", "rural")
    )
  ) %>%
  mutate(
    site_owner = ifelse(grepl("Public", clinic_own), "public", "other"),
    site_owner = factor(site_owner, levels = c("public", "other"))
  ) %>%
  mutate(
    site_carelevel = stringi::stri_extract(level_of_care, regex = "\\w*"),
    site_carelevel = ifelse(site_carelevel == "Teritary",
      "Tertiary", site_carelevel
    ),
    site_carelevel = tolower(site_carelevel),
    site_carelevel = factor(
      site_carelevel,
      levels = c("primary", "secondary", "tertiary", "other")
    )
  ) %>%
  mutate(
    site_patients = ifelse(grepl("Both", patients), "in and out", "only out"),
    site_patients = factor(site_patients, levels = c("in and out", "only out"))
  ) %>%
  mutate(
    site_treated = ifelse(grepl("Both", treated),
      "adults and child", "only adults"
    ),
    site_treated = factor(
      site_treated,
      levels = c("adults and child", "only adults")
    )
  ) %>%
  mutate(
    site_integration =
      ifelse(grepl("Other", integration), "partial",
        ifelse(grepl("Full-integration", integration), "full",
          ifelse(grepl("Partial", integration), "partial", "not integrated")
        )
      ),
    site_integration = factor(
      site_integration,
      levels = c("full", "partial", "not integrated")
    )
  ) %>%
  dplyr::select(matches("site_"))

#### IPC measures ####

ipc <- raw_df %>%
  mutate(
    across(
      c(
        matches("distance_"),
        matches("envmt_"),
        matches("msr_")
      ),
      ~ ifelse(grepl("Already in place", .x), "before",
        ifelse(grepl("Only in place", .x), "during",
          ifelse(grepl("maintained", .x), "since",
            ifelse(grepl("Never", .x), "never",
              ifelse(grepl("Information", .x), "unknown", "not applicable")
            )
          )
        )
      )
    ),
    across(
      c(
        matches("distance_"),
        matches("envmt_"),
        matches("msr_"),
        matches("changes_")
      ),
      ~ factor(.x, levels = c(
        "never", "before",
        "during", "since",
        "unknown", "not applicable"
      ))
    )
  ) %>%
  rename(
    ipc_fixed_app = distance_fixed_app,
    ipc_sep_wr = distance_wr_sep,
    ipc_open_wr = distance_wr_open,
    ipc_triage = distance_triage,
    ipc_vent = envmt_vent,
    ipc_phys_dist = distance_phys,
    ipc_masks_staff_comm = msr_staff_comm_masks,
    ipc_masks_pat_comm = msr_patients_comm_masks,
    ipc_masks_staff_surg = msr_staff_surg_masks,
    ipc_masks_pat_surg = msr_patients_surg_masks,
    ipc_masks_staff_ffp2 = msr_staff_ffp2,
    ipc_masks_pat_ffp2 = msr_patients_ffp2
  ) %>%
  dplyr::select(matches("ipc_"))


#### Shortages ####

short <- raw_df %>%
  mutate(
    across(
      c(changes_staff, changes_hrs),
      ~ ifelse(grepl("Only during", .x), "during",
        ifelse(grepl("Started during", .x), "since",
          ifelse(grepl("Started after", .x), "after",
            ifelse(grepl("Not observed", .x), "never",
              ifelse(grepl("Information", .x), "unknown", "not applicable")
            )
          )
        )
      )
    ),
    across(
      c(restr_tb_tr, restr_tb_ct),
      ~ ifelse(.x == "Yes", "yes", ifelse(.x == "No", "no", "not applicable"))
    ),
    across(
      c(
        short_drugs_fla, short_drugs_sla,
        short_drugs_pra, short_drugs_pa,
        short_drugs_fltb, short_drugs_sltb
      ),
      ~ ifelse(.x == "Yes", "yes", ifelse(.x == "No", "no", "not applicable"))
    ),
    across(
      c(
        restr_tb_tr, restr_tb_ct,
        short_drugs_fla, short_drugs_sla,
        short_drugs_pra, short_drugs_pa,
        short_drugs_fltb, short_drugs_sltb
      ),
      ~ factor(.x, levels = c("yes", "no", "not applicable"))
    )
  ) %>%
  rename(
    short_changes_staff = changes_staff,
    short_changes_hrs = changes_hrs,
    short_restrict_tb_tr = restr_tb_tr,
    short_restrict_tb_ct = restr_tb_ct
  ) %>%
  dplyr::select(
    short_changes_hrs, short_changes_staff,
    short_restrict_tb_ct, short_restrict_tb_tr,
    short_drugs_fla, short_drugs_sla,
    short_drugs_pra, short_drugs_pa,
    short_drugs_fltb, short_drugs_sltb
  )


#### Care provision ####

care <- raw_df %>%
  mutate(
    across(
      c(
        changes_dot_virt_cl,
        tb_restr_addr_drugs,
        tb_restr_addr_tele,
        tb_restr_addr_diff
      ),
      ~ ifelse(grepl("Only in place", .x), "during",
        ifelse(grepl("maintained", .x), "since",
          ifelse(grepl("Never", .x), "never",
            ifelse(grepl("Information", .x), "unknown", "not applicable")
          )
        )
      )
    ),
    across(
      c(
        changes_dot_virt_cl,
        tb_restr_addr_drugs,
        tb_restr_addr_tele,
        tb_restr_addr_diff
      ),
      ~ factor(.x, levels = c(
        "never",
        "during", "since",
        "unknown", "not applicable"
      ))
    )
  ) %>%
  rename(
    care_dot_virt_cl = changes_dot_virt_cl,
    care_restr_addr_drugs = tb_restr_addr_drugs,
    care_restr_addr_tele = tb_restr_addr_tele,
    care_restr_addr_diff = tb_restr_addr_diff
  ) %>%
  dplyr::select(matches("care_"))


#### Combine characteristics ####

# combine
site_descr <- cbind(record_id = 1:nrow(raw_df), site, ipc, short, care)

# save
saveRDS(site_descr, "data-clean/site-characteristics.rds")

#### Case outcomes ####

# long format
df_long <- raw_df %>%
  dplyr::select(
    record_id,
    matches("_new_"),
    matches("_art_"),
    matches("_tpt_"),
    matches("_genex_")
  ) %>%
  mutate(record_id = 1:n()) %>%
  reshape2::melt("record_id") %>%
  mutate(
    quarter = as.numeric(
      gsub(
        "q", "",
        stringi::stri_extract(variable, regex = "q\\d")
      )
    ),
    year = as.numeric(
      paste0("20", gsub(
        "_", "",
        stringi::stri_extract(variable, regex = "_\\d{2}")
      ))
    ),
    outcome = gsub("_", "", stringi::stri_extract(variable, regex = "_\\w*_"))
  )

# quarterly data
df_quarterly <- df_long %>%
  dplyr::select(-variable) %>%
  mutate(year_quart = year + quarter / 4) %>%
  reshape2::dcast(record_id + year + quarter + year_quart ~ outcome) %>%
  mutate(across(c(art, genex, new, tpt), as.numeric))

# only annual data
df_annual <- df_quarterly %>%
  filter(is.na(quarter)) %>%
  dplyr::select(-quarter)

# filter years in quarterly data
df_quarterly <- df_quarterly %>%
  filter(!is.na(quarter))

# save data
saveRDS(df_quarterly, "data-clean/quarterly-data.rds")
saveRDS(df_annual, "data-clean/annual-data.rds")
