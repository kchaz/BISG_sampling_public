#'
#' Calculate statistics reported on in the appendix on validating
#' the obituary data names
#'
library(tidyverse)
library(readr)

#TODO: deal with unusual characters

###############################
# Load Data
###############################
df = read_csv("../../data/ObitData/obit_list_092424-fulldehyph-check.csv")
wikidf_old = read_csv("../../data/Wikipedia/WikipediaSurnames20241125.csv")
wikidf = read_csv("../../data/Wikipedia/WikipediaSurnames20251105.csv")

# cleaned first name dataset
dffirst = read_csv("../../data/ObitData/clean_first_names-v2.csv")
dffirst_merge = c(dffirst$first, df$first)

wikifirst = read_csv("../../data/Wikipedia/WikipediaFirstNames20251105.csv")

dfratio = read_csv("../../data/cluster_output/firstnames_jratio-pretruncation.csv")
dfratio2 = read_csv("../../data/cluster_output/firstnames_jratio.csv")


# Names from Jewish data that were not in voter file
j_unused = read_csv("../../data/cluster_output/unused-jewish-first-names.csv")

###############################
# Basic Stats
###############################
total_unique = length(unique(dffirst_merge))
total_unique

# Stats on first names not used
length(j_unused$jnames_leftover)
sample(j_unused$jnames_leftover, size = 100)
sum(grepl("-", j_unused$jnames_leftover))


###############################
# Wikipedia Check - Last Names
###############################

# Check for things that shouldn't be there
sum(grepl("SURNAME", wikidf$surname))
wikidf[grepl("ARAB", wikidf$surname),]
wikidf[grepl("ORIGIN", wikidf$surname),]
wikidf[grepl("ASHKENAZ", wikidf$surname),]
wikidf$surname[nchar(wikidf$surname)>15]
wikidf$surname[nchar(wikidf$surname)<3]
stopifnot(sum(table(wikidf$surname)>1) == 0)

# How much overlap with previous scrape - mostly
mean(wikidf_old$surname %in% wikidf$surname)

# non-overlap - Rubin weird to remove?
wikidf_old$surname[!wikidf_old$surname %in% wikidf$surname]


length(unique(wikidf$surname))
unique_lasts = unique(df$last)

overlap = wikidf$surname[wikidf$surname %in% unique_lasts]
length(overlap)
length(overlap)/length(wikidf$surname)

only_obit = unique_lasts[!unique_lasts %in% wikidf$surname]
length(only_obit)

only_wiki = wikidf$surname[!wikidf$surname %in% unique_lasts]
length(only_wiki)

#only_wiki = data.frame("surname" = only_wiki)
#dim(only_wiki)
#write.csv(only_wiki, file.path("../../data/WikipediaSurnames/","WikipediaSurnamesNotInObit.csv"))


# stats on counts comparison
counts = table(df$last)

mean(counts)
median(counts)
max(counts)

median(counts[overlap])
mean(counts[overlap])
max(counts[overlap])





######################################################
# First Names - Wikipedia check
#######################################################
length(wikifirst$first)

unique_firsts = unique(dffirst_merge)
length(unique_firsts)

overlap = wikifirst$first[wikifirst$first %in% dffirst_merge]
length(overlap)
length(overlap)/length(wikifirst$first)

only_wiki = wikifirst$first[!wikifirst$first %in% dffirst_merge]
only_wiki

counts = table(dffirst_merge)
counts[overlap]

mean(counts)
median(counts)
max(counts)

median(counts[overlap])
mean(counts[overlap])
max(counts[overlap])

counts[counts>10000]

# Comparing wikipedia names in terms of ratio
mean(dfratio$ratio)
median(dfratio$ratio)

subset = dfratio[dfratio$first %in% wikifirst$first,]
mean(subset$ratio)
median(subset$ratio)
min(subset$ratio)
summary(subset$ratio)


#########################################
# How many additional first names from adding the first name data
#########################################

dead_people_firsts = unique(df$first)
mentioned_firsts = unique(dffirst$first)

length(setdiff(dead_people_firsts, mentioned_firsts))
length(setdiff(mentioned_firsts, dead_people_firsts))





##############
# Ratios in general
##############
# Requires name to appear in both voter file and the obituary scraped data

# How many are > 1 vs < 1
mean(dfratio2$ratio >= 1, na.rm = T)
mean(dfratio2$ratio < 1, na.rm = T)

k = 15 #

top_rows <- dfratio %>%
  filter(count_vf >= 10) %>%
  arrange(desc(ratio)) %>%   # sort descending by ratio
  slice_head(n = k)          # take the first k rows
top_rows

bottom_rows <- dfratio %>%
  filter(count_vf > 10) %>%
  arrange(ratio) %>%   # sort descending by ratio
  slice_head(n = k)          # take the first k rows
bottom_rows


#names with >10 ratios
j_trunc = read_csv("../../data/cluster_output/first-names-that-were-truncated.csv")
dim(j_trunc)



library(ggplot2)

p <- ggplot(dfratio, aes(x = log10(ratio))) +
  geom_histogram(
    bins = 30,
    fill = "black",
    color = "white",
    alpha = 0.8
  ) +
  scale_x_continuous(
    name = "Ratio",
    breaks = seq(floor(min(log10(dfratio$ratio), na.rm = TRUE)),
                 ceiling(max(log10(dfratio$ratio), na.rm = TRUE))),
    labels = function(x) parse(text = paste0("10^", x))
  ) +
  labs(
    title = "Distribution of First Name Ratios",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.3),
    panel.grid.major.y = element_line(linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.ticks.x = element_line(linewidth = 0.8),
    axis.ticks.y = element_line(linewidth = 0.8),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    panel.spacing.x = unit(0.8, "lines"),
    panel.background = element_rect(fill = "white", color = "gray85")
  )

print(p)

ggsave("../../paper/Figures/first_ratio_histogram.pdf", p,
       width = 8, height = 6, dpi = 320)









