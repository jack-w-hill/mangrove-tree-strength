# Mangrove tree strength estimated with field experiments
# Jack W Hill, Vicki Bennion, Catherine E Lovelock

##### 0. set up workspace #####

# packages may need to be installed using 
# install.packages("package-name")
library(patchwork)
library(ragg)
library(readxl)
library(rjson)
library(tidyverse)

set.seed(170823)

# negated inclusive operator
"%!in%" <- Negate("%in%")

# calculate basal area from raw circumference data
calc_ba <- function (dat) {
  max_circs <- max(str_count(dat$circ, ","), 
                   na.rm = TRUE) + 1
  
  circ_cols <- paste0("circ_", seq(1:max_circs))
  
  ba_dat <- dat %>% 
    separate(col = circ,
             into = circ_cols,
             sep = ", ",
             fill = "right") %>% 
    pivot_longer(cols = starts_with("circ")) %>% 
    group_by(tree_id) %>% 
    mutate(value = as.numeric(value),
           value = pi * ((value / pi / 2)^2)) %>% 
    summarise(ba = sum(value, na.rm = TRUE)) %>% 
    ungroup()
  
  return(ba_dat)
}

# custom ggplot theme
theme_strength <- function() {
  theme_classic(base_size = 12, base_family = "sans") %+replace%
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(size = rel(0.75)),
      axis.text.x = element_text(margin = margin(2, 0, 3, 0)),
      axis.text.y = element_text(margin = margin(0, 2, 0, 1)),
      axis.ticks = element_line(colour = "black"),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      plot.tag = element_text(face = "bold")
    )
}

# feature colours for plots
triL_col <- "#7D31ED" # purple
triM_col <- "#ED7D31" # orange
triR_col <- "#13d863" # green

#
##### 1. import data  #####

raw_tree_dat <- read_xlsx("data/static-pulling-trees.xlsx",
                          sheet = "trees")

ba_dat <- calc_ba(raw_tree_dat)

tree_dat <- left_join(raw_tree_dat, ba_dat,
                      by = "tree_id") %>% 
  mutate(lean_from_vertical = case_when(
    trunk_lean_angle == "vertical" ~ 0,
    trunk_lean_angle != 0 ~ abs(90 - as.numeric(trunk_lean_angle)),
    .default = NA),
    tree_species = fct_relevel(as_factor(tree_species),
                               c("am",
                                 "rs",
                                 "ca")),
    dbh = case_when(str_detect(circ, ",") ~ mean(12.1, 11.9),
                    .default = as.numeric(circ) / pi))

read_pul_file <- function (.x) {
  pul_file_location <- paste0("data/pull_tests/", .x, ".pul")
  
  pulling_json_raw <- jsonlite::fromJSON(pul_file_location)
  
  clino1_tib <- as_tibble(pulling_json_raw$data$inclino) %>% 
    mutate(instrument = "clino1",
           tare = pulling_json_raw$inclino_tara)
  
  clino2_tib <- as_tibble(pulling_json_raw$data$inclino2) %>% 
    mutate(instrument = "clino2",
           tare = pulling_json_raw$inclino2_tara)
  
  elasto1_tib <- as_tibble(pulling_json_raw$data$dlen) %>% 
    mutate(instrument = "elasto1",
           tare = pulling_json_raw$dlen_tara) %>% 
    rename(val = val1) %>% 
    select(-val2)
  
  if (str_detect(.x, "gladstone")) {
    # a second elastometer was never used on Gladstone trees
    elasto2_tib <- NULL
  } else {
    elasto2_tib <- as_tibble(pulling_json_raw$data$dlen) %>% 
      mutate(instrument = "elasto2",
             tare = pulling_json_raw$dlen_tara) %>% 
      rename(val = val2) %>% 
      select(-val1)
  }
  
  force_tib <- as_tibble(pulling_json_raw$data$force) %>% 
    mutate(instrument = "force_resp",
           tare = pulling_json_raw$force_tara) %>% 
    mutate(filename = pulling_json_raw$name)
  
  pulling_tib_raw <- bind_rows(clino1_tib,
                               clino2_tib,
                               elasto1_tib,
                               elasto2_tib) %>% 
    mutate(filename = pulling_json_raw$name)
  
  return(list(instruments = pulling_tib_raw, force = force_tib))
}

# skip gladstone_10 (file not relevant because laptop died mid-test)
# and halloran_5 (file not relevant because cable got caught on a prop root
# during test)
relevant_pull_tests <- c(paste0("gladstone_", seq(1, 9)),
                         paste0("gladstone_", seq(11, 33)),
                         paste0("halloran_", seq(1, 4)),
                         paste0("halloran_", seq(6, 10)))

pull_dat  <- map(.x = relevant_pull_tests,
                 .f = read_pul_file) %>% 
  set_names(relevant_pull_tests)

# import terrestrial tree overturning data
dat_terrestrial_om <- read_excel("data/Detter et al 2023_terrestrial-pull-tests.xlsx",
                             sheet = "data") %>% 
  select(-"species_code") %>% 
  rename(scaled_moment = scaled_om)

#
##### 2. calculate overturning moments #####

# given a tree_id, calculate the overturning moment
calc_tree_responses <- function (.x, .y,
                                 special_rs_mode = FALSE,
                                 target_angle = 0.25) {
  focal_tree_id <- .x
  focal_tree_dat <- tree_dat %>% 
    filter(tree_id == focal_tree_id)
  
  pull_id <- focal_tree_dat %>% 
    pull(pul_file)
  
  # instrument data from pull test (clinos, elastos)
  test_inst <- pull_dat[[pull_id]]$instruments %>% 
    select(-c(raw, tare, timestamp)) %>% 
    # identify each row for pivoting
    group_by(instrument) %>% 
    mutate(id = row_number()) %>% 
    ungroup() %>% 
    pivot_wider(names_from = "instrument",
                values_from = "val") %>% 
    select(-id)
  
  # force data from pull test
  test_force <- pull_dat[[pull_id]]$force %>% 
    rename_with(~ paste0(., "_force")) %>% 
    select(force = val_force)
  
  test_dat <- bind_cols(test_inst, test_force)
  
  # heights and distances in m
  sling_height <- focal_tree_dat %>% 
    pull(sling_height_test_tree) 
  
  anchor_sling_height <- focal_tree_dat %>% 
    pull(sling_height_anchor_tree)
  
  tree_dist <- focal_tree_dat %>% 
    pull(dist_test_to_anchor)
  
  # extract which clino and which elasto were used in which locations on tree -- 
  # must be as an rlang "symbol" for use in a pipe, later in this function
  root_clino <- focal_tree_dat %>% 
    pull(root_clino_id) %>% 
    paste0("clino", .) %>% 
    sym(.)
  
  trunk_elasto <- focal_tree_dat %>% 
    pull(trunk_elasto_id) %>% 
    paste0("elasto", .) %>% 
    sym(.)
  
  if (is.na(focal_tree_dat$root_elasto_position)) {
    root_elasto = NULL
  } else if (trunk_elasto == "elasto1") {
    root_elasto = sym("elasto2")
  } else {
    root_elasto = sym("elasto1")
  }
  
  # calculate horizontal force applied, using trigonometry given the 
  # angular winch force
  dat <- test_dat %>%
    # also convert kgf to kN
    mutate(force_h = force * 0.00980665 * 
             cos(atan((sling_height - anchor_sling_height) / tree_dist)))
  
  # get row of interest given angle of interest (normally 0.25 deg)
  interest_dat <- dat %>% 
    filter(!is.na(!!root_clino)) %>% 
    filter(abs(!!root_clino - target_angle) == 
             min(abs(!!root_clino - target_angle))) %>% 
    distinct(!!root_clino, .keep_all = TRUE) 
  
  angle_extracted <- interest_dat %>% 
    pull(!!root_clino)
  print(paste0(angle_extracted, " angle for tree ", .x))
  
  force_h_interest <- interest_dat %>% 
    pull(force_h)
  
  # estimate overturning moment using multiplicative factor from literature
  overturn_moment <- 2.5 * force_h_interest * sling_height
  
  elasto_trunk <- interest_dat %>% 
    pull(!!trunk_elasto)
  
  if (!is.null(root_elasto)) {
    elasto_root <- interest_dat %>% 
      pull(!!root_elasto)
  } else {
    elasto_root <- NULL
  }
  
  return(tibble(tree_id = .x,
                overturn_moment,
                elasto_trunk,
                elasto_root,
                pulling_dat = list(dat)))
}

raw_full_resp_dat <- map(.x = seq(1, max(tree_dat$tree_id)),
                         .f = calc_tree_responses) %>% 
  list_rbind() 

ided_pull_dat <- raw_full_resp_dat %>% 
  select(tree_id, pulling_dat) %>% 
  left_join(., tree_dat,
            by = "tree_id")

raw_resp_dat <- raw_full_resp_dat %>% 
  select(-pulling_dat)

resp_dat <- tree_dat %>%  
  left_join(., raw_resp_dat,
            by = "tree_id") %>% 
  mutate(# convert basal area to m2, from cm2
    tree_size = ba / 10000 * sling_height_test_tree,
    scaled_moment = overturn_moment / tree_size) %>% 
  # remove tree 16 because cable became caught on prop roots during pull test
  filter(tree_id != "16")

#
##### 3. analyse overturning moments #####

# test if scaled moment differs across species
summary(aov(scaled_moment ~ tree_species,
            data = resp_dat))

# test if un-scaled moment differs across species or tree sizes
summary(aov(overturn_moment ~ tree_size + tree_species,
            data = resp_dat))
# (interaction previously tested and was not significant)

# test if fit is better with a log-log transformation
mod_raw_moment <- lm(overturn_moment ~ tree_size,
                     data = resp_dat)
mod_log_raw_moment <- lm(log(overturn_moment) ~ log(tree_size),
                         data = resp_dat)
summary(mod_raw_moment)$adj.r.squared
summary(mod_log_raw_moment)$adj.r.squared

# combine terrestrial and mangrove moments and
dat_plot_grand_om <- resp_dat %>% 
  select(scaled_moment) %>% 
  mutate(study = "mangroves") %>% 
  bind_rows(., dat_terrestrial_om)

# test if moments differ between the groups
t.test(x = dat_plot_grand_om %>% 
         filter(study == "mangroves") %>% 
         pull(scaled_moment),
       y = dat_plot_grand_om %>% 
         filter(study == "detter_2023") %>% 
         pull(scaled_moment),
       alternative = "greater")

#
##### 4. plot overturning moments #####

dat_sum_terrestrial_om <- dat_terrestrial_om %>% 
  summarise(mean_scaled_moment = mean(scaled_moment),
            se_scaled_moment = sd(scaled_moment)/sqrt(n()))

dat_sum_mang_om <- resp_dat %>% 
  summarise(mean_scaled_moment = mean(scaled_moment),
            se_scaled_moment = sd(scaled_moment)/sqrt(n()))

dat_sum_plot_grand_om <- bind_rows(dat_sum_terrestrial_om,
                                   dat_sum_mang_om) %>% 
  mutate(study = c("detter_2023", "mangroves"))

dat_plot_om <- resp_dat %>% 
  group_by(tree_species) %>%
  summarise(mean_scaled_moment = mean(scaled_moment),
            se_scaled_moment = sd(scaled_moment)/sqrt(n()))

plot_grand_scaled_moment <- ggplot(
  mapping = aes(x = fct_relevel(as_factor(study), 
                                c("mangroves", "detter_2023")))
) +
  geom_point(data = dat_plot_grand_om,
             aes(y = scaled_moment),
             colour = "grey70",
             position = position_jitter(width = 0.07)) +
  geom_pointrange(data = dat_sum_plot_grand_om,
                  aes(y = mean_scaled_moment,
                      ymin = mean_scaled_moment - se_scaled_moment,
                      ymax = mean_scaled_moment + se_scaled_moment),
                  size = 0.4,
                  shape = "square") + 
  scale_y_continuous(name = expression("          Scaled predicted\noverturning moment (kNm m"^-3*")"),
                     limits = c(0, NA),
                     expand = expansion(add = c(0, 3))) +
  scale_x_discrete(name = NULL,
                   labels = c("Mangroves\n",
                              "Terrestrial\ntrees")) +
  labs(tag = "a") +
  theme_strength() %+replace%
  theme(plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.3, 0.99))

plot_scaled_moment <- ggplot(mapping = aes(x = tree_species)) +
  geom_point(data = resp_dat,
             aes(y = scaled_moment),
             colour = "grey70") +
  geom_pointrange(data = dat_plot_om,
                  aes(y = mean_scaled_moment,
                      ymin = mean_scaled_moment - se_scaled_moment,
                      ymax = mean_scaled_moment + se_scaled_moment),
                  size = 0.4,
                  shape = "square") +
  scale_y_continuous(name = NULL,
                     limits = c(0, NA),
                     expand = expansion(add = c(0, 3))) +
  scale_x_discrete(name = NULL,
                   labels = c("Avicennia\nmarina",
                              "Rhizophora\nstylosa",
                              "Ceriops\naustralis")) +  
  labs(tag = "b") +
  theme_strength() %+replace%
  theme(axis.text.x = element_text(face = "italic",
                                   margin = margin(2, 0, 3, 0)),
        plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.11, 0.99))

plot_raw_moment <- ggplot(data = resp_dat,
                          mapping = aes(x = log(tree_size),
                                        y = log(overturn_moment),
                                        colour = tree_species,
                                        shape = tree_species)) +
  geom_abline(intercept = mod_log_raw_moment$coefficients[["(Intercept)"]],
              slope = mod_log_raw_moment$coefficients[["log(tree_size)"]],
              size = 1, colour = "grey80") +
  geom_point(size = 1.7) +
  scale_y_continuous(name = "Predicted overturning moment (log kNm)",
                     expand = expansion(add = c(0.07, 0.05))) +
  scale_x_continuous(name = expression("Tree size (log m"^3*")"),
                     expand = expansion(add = c(0.05, 0.05))) +
  scale_colour_manual(name = NULL,
                      labels = c("A. marina",
                                 "R. stylosa",
                                 "C. australis"),
                      values = c(triL_col,
                                 triM_col,
                                 triR_col)) +
  scale_shape_manual(name = NULL,
                     labels = c("A. marina",
                                "R. stylosa",
                                "C. australis"),
                     values = c("square",
                                "circle",
                                "triangle")) +
  labs(tag = "c") +
  theme_strength() +
  theme(legend.position = c(1.02, -0.01),
        legend.justification = c("right", "bottom"),
        legend.text = element_text(size = rel(0.75),
                                   face = "italic",
                                   margin = margin(0,0,0,-7)),
        legend.key.height = unit(11, "points"),
        legend.background = element_blank(),
        plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.112, 0.99))

layout <- c(
  patchwork::area(1, 1, 1, 1),
  patchwork::area(1, 2, 1, 3),
  patchwork::area(2, 1, 2, 3)
)

agg_png("fig-output/figure-2.png",
        res = 1080, scaling = 1, units = "cm",
        width = 15, height = 18)
plot_grand_scaled_moment + plot_scaled_moment + plot_raw_moment +
  plot_layout(design = layout) &
  theme(plot.background = element_blank())
dev.off()

#
##### 5. tidy, analyse, and plot elasticity data #####

dat_long_elasto <- resp_dat %>% 
  mutate(elasto_trunk = case_when(trunk_elasto_position == "extension" & 
                                    elasto_trunk < 0 ~ elasto_trunk * -1,
                                  trunk_elasto_position == "compression" & 
                                    elasto_trunk > 0 ~ elasto_trunk * -1,
                                  TRUE ~ elasto_trunk),
         elasto_root = case_when(root_elasto_position == "extension" & 
                                   elasto_root < 0 ~ elasto_root * -1,
                                 root_elasto_position == "compression" & 
                                   elasto_root > 0 ~ elasto_root * -1,
                                 TRUE ~ elasto_root),
         tree_id = as.character(tree_id)) %>% 
  filter(# remove two trees where the elastometer jaws were not released -- 
    # i.e. no measurement as elastometer could not extend/compress
    tree_id %!in% c(25, 26), 
    trunk_elasto_position != "compression") %>% 
  pivot_longer(cols = starts_with("elasto"),
               names_to = "elasto",
               values_to = "delta_length") %>% 
  # convert elastometer units -> per hundred (%)
  mutate(delta_length = delta_length / 1000) 

dat_plot_trunk_elasto <- dat_long_elasto %>% 
  filter(elasto == "elasto_trunk") 

# test if trunk elasticity differs across species
summary(aov(delta_length ~ tree_species, 
            dat_plot_trunk_elasto))

dat_sum_plot_trunk_elasto <- dat_plot_trunk_elasto %>% 
  group_by(tree_species) %>% 
  summarise(mean_dlen = mean(delta_length),
            se_dlen = sd(delta_length)/sqrt(n())) %>% 
  ungroup()

plot_trunk_elasto <- ggplot(mapping = aes(x = tree_species)) +
  geom_point(data = dat_plot_trunk_elasto,
             aes(y = delta_length),
             colour = "grey70") +
  geom_pointrange(data = dat_sum_plot_trunk_elasto,
                  aes(y = mean_dlen,
                      ymin = mean_dlen - se_dlen,
                      ymax = mean_dlen + se_dlen),
                  size = 0.4,
                  shape = "square") +
  scale_y_continuous(name = "Trunk Δlength (%)",
                     limits = c(0, NA),
                     expand = expansion(add = c(0, 0.1))) +
  scale_x_discrete(name = NULL,
                   labels = c("Avicennia\nmarina",
                              "Rhizophora\nstylosa",
                              "Ceriops\naustralis")) +
  labs(tags = "a") +
  theme_strength() %+replace%
  theme(axis.text.x = element_text(face = "italic",
                                   margin = margin(2, 0, 3, 0)),
        plot.tag.position = c(0.16, 0.99))

# only trees at the Halloran site had two elastometers
dat_plot_both_elasto <- dat_long_elasto %>% 
  filter(site == "halloran")

dat_plot_elasto_compare <- dat_plot_both_elasto %>% 
  pivot_wider(names_from = elasto,
              values_from = delta_length) %>% 
  mutate(dlen_ratio = abs(elasto_root / elasto_trunk)) 

dat_sum_plot_elasto_compare <- dat_plot_elasto_compare %>% 
  group_by(root_elasto_position) %>% 
  summarise(mean_dlen_ratio = mean(dlen_ratio),
            se_dlen_ratio = sd(dlen_ratio) / sqrt(n())) %>% 
  ungroup()

plot_elasto_compare <- ggplot() +
  # geom_rect() requires a tibble
  geom_rect(data = tibble(xmin = 0, xmax = 3,
                          ymin = 0, ymax = 1),
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            fill = "grey90") +
  geom_point(data = dat_plot_elasto_compare,
             aes(x = root_elasto_position,
                 y = dlen_ratio),
             colour = "grey70") +
  geom_pointrange(data = dat_sum_plot_elasto_compare,
                  aes(x = root_elasto_position, 
                      y = mean_dlen_ratio,
                      ymin = mean_dlen_ratio - se_dlen_ratio, 
                      ymax = mean_dlen_ratio + se_dlen_ratio),
                  size = 0.4,
                  shape = "square",
                  colour = "black") +
  scale_y_continuous(name = "Ratio of root Δlength : trunk Δlength",
                     limits = c(0, NA),
                     expand = expansion(add = c(0, 0.1))) +
  scale_x_discrete(name = "Force experienced\nby prop root",
                   labels = c("Compression ",
                              "  Extension")) +
  labs(tags = "b") +
  theme_strength() %+replace%
  theme(plot.tag.position = c(0.2, 0.99),
        axis.title.x = element_text(margin = margin(-7, 0, 0, 0)))

agg_png("fig-output/figure-3.png",
        res = 1080, scaling = 1, units = "cm",
        width = 15, height = 12)
plot_trunk_elasto + plot_elasto_compare +
  plot_layout(widths = c(3, 2)) &
  theme(plot.background = element_blank())
dev.off()

#
##### 6. estimate whole-tree modulus of elasticity #####

stress_dat <- resp_dat %>% 
  mutate(overturn_moment = overturn_moment * 1000000, # kNm -> Nmm
         dbh = dbh * 10, # cm -> mm
         sect_mod = dbh**3 * pi/32, # calculate section modulus
         stress = overturn_moment / sect_mod) %>% # in Nmm2
  filter(!is.na(dbh))

moe_dat <- dat_plot_trunk_elasto %>% 
  select(tree_id, delta_length) %>%
  mutate(tree_id = as.numeric(tree_id)) %>% 
  left_join(stress_dat, .,
            by = "tree_id") %>%
  mutate(strain = delta_length, # % (functionally unitless)
         moe = stress / strain) # Nmm2 (== MPa)

moe_dat %>% 
  filter(!is.na(moe)) %>% 
  group_by(tree_species) %>% 
  summarise(mean(moe),
            sd(moe)/sqrt(n()))

stress_dat %>% 
  group_by(tree_species) %>% 
  summarise(mean(stress),
            sd(stress)/sqrt(n()))
# gives stress at predicted overturning moment, i.e. critical stress

# fragility function from Yanagisawa et al. (2010):
  # prob_damage = pnorm((log(stress) - 2.54) / 0.90)
# rearranged for stress:
  # stress = exp(qnorm(prob_damage) * 0.9 + 2.54)
exp(qnorm(.75) * 0.9 + 2.54)
# gives stress in Nmm2

#
##### 7. make wood density table #####

dat_raw_wood_dens <- read_excel("data/Chave et al 2009_wood-density-database.xls",
                                sheet = "Data") 

mang_genera <- c("Avicennia", "Rhizophora", "Ceriops")

terrestrial_genera <- c("Betula", "Picea", "Fraxinus", "Pinus", "Fagus", "Populus",
                        "Platanus", "Acer", "Prunus", "Quercus", "Tillia", "Robina",
                        "Alnus", "Pseudotsuga")

dat_wood_dens <- dat_raw_wood_dens %>% 
  rename(wood_dens = `Wood density (g/cm^3), oven dry mass/fresh volume`) %>% 
  mutate(genus = str_split_i(Binomial, " ", 1),
         study = case_when(genus %in% mang_genera ~ "mang",
                           genus %in% terrestrial_genera ~ "terrestrial",
                           .default = "neither")) %>% 
  filter(study != "neither",
         !is.na(wood_dens))

dat_wood_dens %>% 
  group_by(study) %>% 
  summarise(min(wood_dens),
            max(wood_dens),
            median(wood_dens),
            sd(wood_dens)/sqrt(n()),
            n()) %>% 
  view()

#
##### 8. plot pulling test trajectories #####

# function modified from calc_tree_reponses()
plot_pull_test <- function (focal_tree,
                            pull_filename = NULL) {
  focal_tree_id <- focal_tree
  focal_tree_dat <- tree_dat %>% 
    filter(tree_id == focal_tree_id)
  
  if (!is.null(pull_filename)) {
    pull_id <- pull_filename
  } else {
    pull_id <- focal_tree_dat %>% 
      pull(pul_file)
  }
  
  test_inst <- pull_dat[[pull_id]]$instruments %>% 
    select(-c(raw, tare, timestamp)) %>% 
    group_by(instrument) %>% 
    mutate(id = row_number()) %>% 
    ungroup() %>% 
    pivot_wider(names_from = "instrument",
                values_from = "val") %>% 
    select(-id)
  
  test_force <- pull_dat[[pull_id]]$force %>% 
    rename_with(~ paste0(., "_force")) %>% 
    select(force = val_force)
  
  test_dat <- bind_cols(test_inst, test_force)
  
  sling_height <- focal_tree_dat %>% 
    pull(sling_height_test_tree) 
  
  anchor_sling_height <- focal_tree_dat %>% 
    pull(sling_height_anchor_tree)
  
  tree_dist <- focal_tree_dat %>% 
    pull(dist_test_to_anchor)
  
  root_clino <- focal_tree_dat %>% 
    pull(root_clino_id) %>% 
    paste0("clino", .) %>% 
    sym(.)
  
  trunk_elasto <- focal_tree_dat %>% 
    pull(trunk_elasto_id) %>% 
    paste0("elasto", .) %>% 
    sym(.)
  
  dbh_raw <- focal_tree_dat %>% 
    pull(dbh)
  dbh <- dbh_raw * 10 # cm -> mm
  
  if (is.na(focal_tree_dat$root_elasto_position)) {
    root_elasto = NA
  } else if (trunk_elasto == "elasto1") {
    root_elasto = sym("elasto2")
  } else {
    root_elasto = sym("elasto1")
  }
  
  if (is.na(focal_tree_dat$height_trunk_clino)) {
    trunk_clino = NA
  } else if (focal_tree_dat$root_clino_id == "1") {
    trunk_clino = sym("clino2")
  } else {
    trunk_clino = sym("clino1")
  }
  
  title <- paste0("Tree ID ", focal_tree_id, " —")
  
  dat <- test_dat %>%
    mutate(force_h = force * 0.00980665 * 
             cos(atan((sling_height - anchor_sling_height) / tree_dist)),
           across(starts_with("elasto"),
                  ~ . / 1000))%>% 
    rownames_to_column(var = "time") %>% 
    mutate(time = as.numeric(time),
           # normalise time to between 0 and 1 so plots line up in Patchwork
           time = (time - min(time)) / (max(time) - min(time))) 
  
  pull_plots <- map2(.x = c("force", 
                            root_clino, trunk_clino,
                            trunk_elasto, root_elasto),
                     .y = c(NA,
                            "Tilt at root plate (°)",
                            "Tilt at pull height (°)",
                            "Trunk Δlength (%)",
                            "Root Δlength (%)"),
                     .f = function(.x, .y) {
                       
                       if (is.na(.x)) {
                         return(plot_spacer())
                       } else if (.x == "force") {
                         force_plot <- ggplot(data = dat) +
                           geom_point(aes(x = time,
                                          y = force_h),
                                      size = 1) +
                           scale_y_continuous(limits = c(0, NA),
                                              name = "Horizontal force (kN)",
                                              expand = expansion(add = c(0,
                                                                         0.01))) +
                           scale_x_discrete(name = "Test time",
                                            labels = NULL) + 
                           ggtitle(title) +
                           theme_strength() %+replace%
                           theme(axis.ticks.x = element_blank(),
                                 plot.title = element_text(size = 10,
                                                           hjust = 0))
                         
                         return(force_plot)
                       }
                       
                       current_plot <- ggplot(data = dat) +
                         geom_point(aes(x = time,
                                        y = !!.x),
                                    colour = "red",
                                    size = 1) +
                         scale_y_continuous(limits = c(0, NA),
                                            name = .y,
                                            expand = expansion(add = c(0, 
                                                                       0.01))) +
                         scale_x_discrete(name = "Test time",
                                          labels = NULL) + 
                         theme_strength() %+replace%
                         theme(axis.ticks.x = element_blank())
                       
                       return(current_plot)
                     }) %>% 
    set_names(c("force",
                "root_clino", "trunk_clino",
                "trunk_elasto", "root_elasto"))
  
  return(pull_plots)
}

plot_all_pulls <- map(.x = seq(1:30),
                       .f = plot_pull_test)

agg_png("fig-output/pull-trajectories.png",
        res = 900, scaling = 1, units = "cm",
        width = 25, height = 140) # max acceptable height under ragg defaults
wrap_plots(list_c(plot_all_pulls), ncol = 5)
dev.off()

#
##### 9. plot Rhizophora cutting data #####

# tree_id 15 was one of three Rhizophora trees with dual pulls (before
# and after root cutting)
before <- plot_pull_test(15) 
after <- plot_pull_test(15, "gladstone_22")

before_dat <- before$force$data 
after_dat <- after$force$data 

time_at_025 <- before_dat %>% 
  filter(abs(clino2 - 0.25) == 
           min(abs(clino2 - 0.25))) %>% 
  distinct(clino2, .keep_all = TRUE) %>% 
  pull(time)
# time at 0.25 varies by ~0.0025 sec between before and after datasets

label_at_025 <- tibble(x = time_at_025 - 0.02,
                       y = 0.19,
                       label = "at 0.25°\ninclination",
                       size = 0.5,
                       vjust = "top",
                       hjust = "right")

plot_rs_cut_force <- ggplot(mapping = aes(x = time,
                                          y = force_h)) +
  geom_line(data = before_dat,
            colour = triR_col) +
  geom_line(data = after_dat,
            colour = triM_col) +
  geom_vline(xintercept = time_at_025,
             linetype = "dashed") +
  annotate("text", label = "at 0.25°\ninclination",
           x = time_at_025 - 0.01, 
           y = 0.195, 
           size = 2.5,
           lineheight = 0.75,
           fontface = "italic",
           vjust = "top",
           hjust = "right") +
  scale_y_continuous(limits = c(0, 0.20),
                     name = "Horizontal force (kN)",
                     expand = expansion(add = c(0, 
                                                0.01))) +
  scale_x_discrete(name = "Test time",
                   labels = NULL,
                   expand = expansion()) + 
  theme_strength() %+replace%
  theme(axis.title.x = element_text(margin = margin(5, 0, 0, 0)),
        axis.ticks.x = element_blank())

plot_rs_cut_clino <- ggplot(mapping = aes(x = time,
                                          y = clino2)) +
  geom_line(data = before_dat,
            colour = triR_col) +
  geom_line(data = after_dat,
            colour = triM_col) +
  scale_y_continuous(limits = c(0, 0.6),
                     labels = function(x) format(x, nsmall = 2),
                     name = "Inclination at root plate (°)",
                     expand = expansion(add = c(0, 
                                                0.01))) +
  scale_x_discrete(name = "Test time",
                   labels = NULL,
                   expand = expansion()) + 
  theme_strength() %+replace%
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

agg_png("fig-output/figure-4.png",
        res = 1080, scaling = 0.8, units = "cm",
        width = 8, height = 12)
plot_rs_cut_clino / plot_rs_cut_force
dev.off()

#