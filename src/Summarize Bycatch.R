require(dplyr)
require(reshape2)


setwd("./data")

data <- read.csv("WCMarineMammals.csv", header = T)
effort.data <- read.csv("ObsEffort.csv")
setwd("..")


interaction.types <- as.character(unique(data$Interaction))

interactions.used <-
  interaction.types[c(3, 4, 6, 7, 8)] # only looking at entanglements, killed, lethal removals

gear.2.use <- "Bottom Trawl"
by.species <- data %>%
  filter(Interaction %in% interactions.used, Gear == gear.2.use) %>%
  group_by(Species, Year, Gear) %>%
  summarise(sigma_n = sum(No.individuals),
            mean.coverage = mean(Coverage.rate))

trawl.effort <- effort.data %>%
  filter(Gear == gear.2.use) %>%
  select(Year, N.Obs.Halus)

species.2.use <- "California Sea Lion"
# Cycle through years, get

tmp.byspecies <- by.species %>%
  filter(Species == species.2.use)

yearlist <- as.vector(trawl.effort$Year)

extract.interaction <- function(year, x) {
  index <- which(x$Year == year)
  if (!length(index)) {
    return(0)
  } else {
    return(as.numeric(x[index, 4]))
  }
}

interaction.list <-
  as.vector(sapply(yearlist, FUN = extract.interaction, x = tmp.byspecies))

# Function for likelihood profile
lik.profile <- function(k, N) {
  phat <- k / N
  NLML <-
    -dbinom(k, N, phat, log = TRUE) # negative log maximum likelihood
  # moving upwards
  inc.test <- 0.00001
  NLL.prior <- NLML
  p.prior <- phat
  p.test <- phat + inc.test
  
  NLL.test <- -dbinom(k, N, p.test, log = TRUE)
  if (NLL.test > (NLML + 1.92)) {
    #this means that the increment already overshot, try a smaller increment
    p.prior <- phat
    p.test <- phat
    NLL.prior <- NLML
    NLL.test <- NLML
    inc.test <- 0.000001
  }
  while (NLL.test <= (NLML + 1.92)) {
    p.prior <- p.test
    NLL.prior <- NLL.test
    p.test <- p.prior + inc.test
    NLL.test <- -dbinom(k, N, p.test, log = TRUE)
  }
  ub <-
    approx(c(NLL.prior, NLL.test), c(p.prior, p.test), xout = (NLML + 1.92))$y
  # Now search for the lower bound, tricker because of 0
  if (phat == 0)
    lb = 0
  if (phat > 0) {
    # here I'll multiply to reduce p.test so that it cannot go negative
    p.test <- phat
    NLL.prior <- NLML
    NLL.test <- NLML
    inc.test <- 0.9 # change by 10% each iteration
    
    while (NLL.test <= (NLML + 1.92)) {
      p.prior <- p.test
      NLL.prior <- NLL.test
      p.test <- p.prior * inc.test
      NLL.test <- -dbinom(k, N, p.test, log = TRUE)
    }
    lb <-
      approx(c(NLL.prior, NLL.test),
             c(p.prior, p.test),
             xout = (NLML + 1.92))$y
    
  }
  return(c(lb, phat, ub))
}


# Loop through years
ci.output <- matrix(NA, nrow = length(yearlist), ncol = 3)
for (i in 1:length(yearlist)) {
  ci.output[i, 1:3] <-
    lik.profile(interaction.list[i], trawl.effort$N.Obs.Halus[i])
}
# simple plot of results, set up y axis (adjust as needed)
ylim <- c(0, 0.008)
par(xpd = TRUE, las = 1)
plot(
  yearlist,
  ci.output[, 2],
  type = "p",
  pch = 21,
  bg = "black",
  yaxs = "i",
  xlab = "Year",
  las = 1,
  ylab = "Interaction Probability",
  ylim = ylim
)

for (i in 1:length(yearlist)) {
  lines(rep(yearlist[i], 2), ci.output[i, c(1, 3)], col = "black")
}
