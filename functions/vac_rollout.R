if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}

rollout_children <- data.frame(
                      coronavac = 0,
                      astrazeneca = 0,
                      pfizer = 0,
                      date = seq(from = as.Date("2022-01-13"), to = as.Date("2022-03-31"), by= "days"))

rollout_feb_mar <- 20000000 - 1248000*3
rollout_jan <- 1248000

remaining_days_count <- length(rollout_children[rollout_children$date >= "2022-02-03" & rollout_children$date < "2022-04-01",2])
rollout_children[rollout_children$date >= "2022-01-13" & rollout_children$date < "2022-01-20","pfizer"] <- round(rollout_jan/7)
rollout_children[rollout_children$date >= "2022-01-20" & rollout_children$date < "2022-01-27","pfizer"] <- round(rollout_jan/7)
rollout_children[rollout_children$date >= "2022-01-27" & rollout_children$date < "2022-02-03","pfizer"] <- round(rollout_jan/7)
rollout_children[rollout_children$date >= "2022-02-03" & rollout_children$date < "2022-04-01","pfizer"] <- round(rollout_feb_mar/remaining_days_count)

write.csv(rollout_children, file = "DATA/rollout_children.csv", row.names = FALSE)

## Plot
rollout_children %>%
  ggplot(aes(x = date, y = pfizer)) +
  #geom_line(size = 1, color = "red") +
  geom_point(size = 2, color = "red") +
  ylim(0, max(rollout_children$pfizer)) +
  xlab("") + ylab("Daily vaccine rollout/n") +
  labs(title = "Children vaccination rollout (Pfizer)/n") +
  scale_x_date(date_breaks = "week", date_labels = "%d/%m", limits = as.Date(c("2022-01-01","2022-03-31"))) +
  theme_classic()
