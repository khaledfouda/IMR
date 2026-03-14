# copy first part from before:>
require(tidyverse)
devtools::load_all()
require(magrittr)
source("article/movielens/preprocess.R")
out <- prepare_results_for_analysis()
dat = out$dat; fits=out$fits; res.df=dat$res


#f$lambda_gamma
#f <- fit.imr3;f$fit <- NULL; f
# --- movie→genre map for the 5 genres ------------------------------------
genres <- c("Documentary", "Musical", "Drama", "Fantasy", "Children's",
            "War", "Action", "Sci-Fi", "Horror", "Animation")
Z_sel <- dat$Z[,genres , drop = FALSE]
arr   <- which(Z_sel != 0, arr.ind = TRUE)


mapz <- data.frame(
  genre = colnames(Z_sel)[arr[, 2]],
  movie =  as.character(arr[, 1]),
  stringsAsFactors = FALSE
)  %>%
  mutate(movie = as.numeric(movie))
mapx <- as.data.frame(dat$X) %>%
  mutate(user = seq(1:nrow(dat$X))) %>%
  mutate(group = case_when(
    .G == 1 & .A1 == 1 ~ "Male (25-34)",
    .G == 0 & .A1 == 1 ~ "Female (25-34)",
    .G == 1 & .A2 == 1 ~ "Male (35-49)",
    .G == 0 & .A2 == 1 ~ "Female (35-49)",
    .G == 1 & .A3 == 1 ~ "Male (50+)",
    .G == 0 & .A3 == 1 ~ "Female (50+)",
    .G == 1  ~ "Male (0-24)",
    .G == 0  ~ "Female (0-24)",
  )) %>%
  dplyr::select(user, group)



mapz %<>%
  group_by(movie) %>%
  mutate(n = n()) %>%
  ungroup %>%
  filter(n == 1) %>%
  dplyr::select(-n)




movies <- unique(arr[, 1])
subYh <- (out$out[[3]]$xbeta+out$out[[3]]$gammaz)[,movies]
subYh <- out$out[[3]]$estimates[,movies] - subYh
subYh  <- out$out[[3]]$estimates[,movies]
colnames(subYh) = movies
as.data.frame(subYh) %>%
  rownames_to_column("row") %>%
  pivot_longer(-row, names_to="col", values_to="value") %>%
  mutate(across(everything(), as.numeric)) %>%
  rename(movie=col, user=row, estimate = value) %>%
  inner_join(mapz, "movie", relationship="many-to-one") %>%
  inner_join(mapx, "user", relationship ="many-to-one") ->
  ume


ume %>%
  #count(movie, genre, group, name="n_cell") %>%
  group_by(genre, movie, group) %>%
  mutate(median_estim =  median(estimate)) %>% head()

q7 <- function(x, p) stats::quantile(x, probs = p, na.rm = TRUE, type = 7, names = FALSE)

ume %>%
  summarise(
    q1 = q7(estimate, 0.25),
    q3 = q7(estimate, 0.75),
    pooled_IQR = q3-q1,
    .by = c(movie, genre)
  ) %>%
  dplyr::select(movie, genre, pooled_IQR) ->
  pooled_iqr_tbl

ume %>%
  summarise(
    med = median(estimate, na.rm=TRUE),
    .by = c(movie, genre, group)
  ) %>%
  summarise(
    S_diff = diff(range(med)),
    .by = c(movie, genre)
  ) ->
  median_gap_tbl

# get movie titles:
data.table::fread(
  file = "article_results/movielens/data/movies_Z.dat",
  sep = NULL,
  encoding = "Latin-1",
  header = FALSE
) %>%
  tidyr::separate(
    V1,
    into = c("movie", "title", "genres"),
    sep = "::"
  ) %>%
  as.data.frame() %>%
  mutate(movie = as.numeric(movie)) %>%
  dplyr::select(-genres) %>%
  filter(movie %in% mapz$movie) ->
  titles




pooled_iqr_tbl %>%
  inner_join(median_gap_tbl, by = c("movie", "genre")) %>%
  mutate(
    .by = genre,
    z_iqr = as.numeric(scale(pooled_IQR)),
    z_diff = as.numeric(scale(S_diff)),
    S = 0.5 * z_diff + 0.5 * z_iqr
  ) %>%
  arrange(genre, desc(S), desc(S_diff), desc(pooled_IQR), movie) %>%
  group_by(genre) %>%
  slice_max(order_by = S, n=3, with_ties = FALSE) %>%
  ungroup() %>%
  filter(genre %in% c( "Children's")) %>%
  arrange(genre, desc(S)) ->
  selected_movies

movies_picked <- c(586, 1367, 1592)
  #c("Home Alone (1990)","101 Dalmatians (1996)", "Air Bud (1997)" )
pooled_iqr_tbl %>%
  inner_join(median_gap_tbl, by = c("movie", "genre")) %>%
  mutate(
    .by = genre,
    z_iqr = as.numeric(scale(pooled_IQR)),
    z_diff = as.numeric(scale(S_diff)),
    S = 0.5 * z_diff + 0.5 * z_iqr
  ) %>%
  arrange(genre, desc(S), desc(S_diff), desc(pooled_IQR), movie) %>%
  filter(movie %in% movies_picked) %>%
  filter(genre %in% c( "Children's")) %>%
  arrange(genre, desc(S)) ->
  selected_movies

titles %>%
  filter(movie %in% selected_movies$movie) %>%
  inner_join(selected_movies, "movie", relationship = "one-to-one") %>%
  dplyr::select(-z_iqr, -z_diff) ->
  selected_movies


ume |>
  semi_join(selected_movies, by = c("movie", "genre")) %>%
  inner_join(titles, "movie", relationship = "many-to-one") ->
  ume_sel

gender_lv <- c("Female","Male")
age_lv    <- c("0-24","25-34","35-49","50+")
# group_lv  <- as.vector(outer(gender_lv, age_lv, \(g,a) paste0(g, " (", a, ")")))
group_lv <- c(
  "Female (0-24)",
  "Female (25-34)",
  "Female (35-49)",
  "Female (50+)",
  "Male (0-24)",
  "Male (25-34)",
  "Male (35-49)",
  "Male (50+)"
)

ume_sel <- ume_sel %>%
  mutate(group = factor(group, levels = group_lv),
         age   = str_extract(group, "(?<=\\().+(?=\\))") |> factor(levels = age_lv),
         gender = ifelse(stringr::str_detect(group, "^Female.*"), "Female", "Male")|>
           factor())

# Order movies within genre by selection score (if available)
if (exists("selected_movies")) {
  movie_order <- selected_movies %>% arrange(genre, desc(S)) %>%
    distinct(genre, movie) %>% group_by(genre) %>% mutate(order = row_number()) %>% ungroup()
  ume_sel <- ume_sel %>% left_join(movie_order, by = c("movie","genre"))
}
ume_sel$movie <- factor(ume_sel$movie, levels = unique(ume_sel$movie[order(ume_sel$genre, ume_sel$order)]))


library(ggh4x)
library(grid)

ggplot(ume_sel, aes(x = age, y = estimate, fill = age)) +
  geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 0.7, alpha = 0.9) +
  ggh4x::facet_nested(
    cols = vars(title, gender),
    scales = "fixed",
    # switch = "x",
    strip = ggh4x::strip_nested(
      text_x = ggh4x::elem_list_text(face = c("bold"))
    )
  ) +

  scale_fill_brewer(palette = "Set2", name = "Age") +
  scale_y_continuous("Estimated rating", limits = c(0.5, 5), breaks = seq(0.5, 5, 0.5),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_discrete(name = "Age group") +
  theme_bw(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.background.x = element_blank(),
    panel.spacing.x = unit(0, "pt"),
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    #axis.text.x = element_text(angle = 0, hjust = 0),
    legend.position = "none",
    ggh4x.facet.nestline = element_line(colour = "grey70")
  ) +
  ggtitle("Fitted Movie Ratings (Full Model)",
          subtitle = "Selected children's movies by gender and age group") -> g;g
#
# ggsave("./article_results/movielens/data/plot_intercept_model.png",
#        g, width = 320/25.4, height = 150/25.4, dpi = 600)


ggsave("./article/movielens/data/plot_full_model.png",
       g, width = 9, height = 4, scale=1.2, dpi = 300)


#===============================================================================
# generate a table of xbeta contribution for those movies across age groups
