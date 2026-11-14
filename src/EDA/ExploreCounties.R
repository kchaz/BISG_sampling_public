#
# Explore the AJPP data from Brandeis on geographic distribution of Jewish
# people at the county level
#
# Author: Kyla Chasalow

library(tidyverse)
filename <- "../../data/JewishCountyCounts/brandeis_fips_pop_counts-proc.csv"
df <- read_csv(filename)

# Counts
#--------------
summary(df$all_adults)
summary(df$jewish_adults)
summary(df$jewish_all) #include children? but there's no "all" column?
df[which.max(df$jewish_adults),]

adults_v1 = df$all_adults
adults_v2 = (df$all_adults - mean(df$all_adults))/sd(df$all_adults)
adultjews_v1 = df$jewish_adults
adultjews_v2 = (df$jewish_adults - mean(df$jewish_adults))/sd(df$jewish_adults)

par(mfrow = c(1,2))
hist(adults_v1, col = scales::alpha("red",.4), breaks = 30,
     main = "Geographic Concentration of All Adults",
     xlab = "County Counts",
     cex.lab = 1.5,
     freq = F)
hist(adultjews_v1,
     col = scales::alpha("blue",.5),
     main = "Geographic Concentration of Jewish Adults",
     breaks = 30,
     freq = F,
     cex.lab = 1.5,
     xlab = "County Counts")
legend("topright", legend = c("All Adults", "Jewish Adults"),
       col = c( scales::alpha("red",.4),scales::alpha("blue",.5)),
       pch = 20,
       cex = 2)



# entries to left of line are those for which represents higher faction of Jewish
# population than of US population. The biggest one is Los Angeles.
# NY is smaller because broken up into Manhattan, Brooklyn etc.
par(mfrow = c(1,1))
plot(df$county_given_jewish_adult, df$adult_prop)
abline(a=0,b=1, col = "gray")
#plot(df$jewish_adult_given_county, df$adult_prop)
#abline(a=0,b=1, col ="gray")


#P(R=1)
sum(df$jewish_adults)/sum(df$all_adults)

# relationship between # of FIPs and # of Jewish people
# mostly only those with almost no Jews that have been grouped
# so probably pretty safe to distribute those probs to all fips within
plot(df$num_fips, df$jewish_adults,
     xlab = "Number of FIPS codes Assigned",
     ylab = "Number of Jewish Adults")


# get counts by state
state_counts = group_by(df, state) %>%
  summarize(jewish_adults = sum(jewish_adults),
            all_adults = sum(all_adults))
state_counts$jewish_adult_given_state = state_counts$jewish_adults/state_counts$all_adults
state_counts$state_given_jewish_adult = state_counts$jewish_adults/sum(state_counts$jewish_adults)

par(mfrow=c(1,2))
hist(state_counts$jewish_adult_given_state, xlim = c(0,1))
hist(state_counts$state_given_jewish_adult, xlim = c(0,1))


# how much variation is there in probability by state?
propsum = df %>% group_by(state) %>%
  summarize(meanprop = mean(jewish_adult_given_county),
            sdprop = sd(jewish_adult_given_county),
            count = n(),
  )
hist(propsum$meanprop)
hist(propsum$sdprop)

propsum %>% filter(state == "FL")
df %>% filter(state == "FL")  %>% select(jewish_adult_given_county, county_given_jewish_adult)
