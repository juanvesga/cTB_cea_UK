# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root                 <- here::here()
	infile0              <- file.path(root, "data", "raw", "TBdeath_linelist_v2.csv")
	infile               <- file.path(root, "data", "raw", "p_cfr_ets.csv")
	infile_parameters    <- file.path(root, "output", "parameters.qs2")
	outfile              <- file.path(root, "output", "p_cfr_ets.qs2")

	# Packages
	source(file.path(root, "R", "modify_attach.R"))
	modify_attach(qs2,     include.only = c("qs_read", "qs_save"))

	library(dplyr)
	library(ggplot2)
	library(sjPlot)

} else {

	# Input
	args    <- commandArgs(trailingOnly = TRUE)
	infile  <- args[1]
	outfile <- args[2]

	# Packages
	library(qs2, include.only = "qs_save")
}


# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------
dat              <- read.csv(infile)
df               <- read.csv(infile0)
parameters       <- qs_read(infile_parameters)

unique(df$startedtreat)
unique(df$tb_type2)
unique(df$tbrelateddeath)

# Select a clean dataset and get binary death outcome
df2<-df %>%
	filter(startedtreat=="started") %>%
	filter(!is.na(tb_type2)) %>%
	filter(!is.na(age)) %>%
	filter(age>15) %>%
	filter(age<=100) %>%
	mutate(tbdeath = case_when(tbrelateddeath=="TB related death" ~ 1,
							   tbrelateddeath=="Other outcome" ~ 0) ) %>%
	mutate(origin = case_when(ukborn=="non_uk_born" ~ "non_uk_born",
							  ukborn=="uk-born" ~ "uk_born",
							  ukborn=="missing" ~ "non_uk_born",) ) %>%
	mutate(age.cat = cut(age, breaks = seq(15,100,5))) %>%
	select(c(age.cat,origin,tb_type2,tbdeath))

head(df2,10)

df2$tbdeath<-as.factor(df2$tbdeath)
#df2$ukborn<-as.factor(df2$ukborn)
df2$origin<-as.factor(df2$origin)
df2$tb_type2<-as.factor(df2$tb_type2)
df2$age.cat<-as.factor(df2$age.cat)



#plot
ggplot(df2, aes(age.cat, fill = tbdeath)) +
	geom_bar() +
	coord_flip()

ggplot(df2[df2$tbdeath==1,], aes(age.cat, fill = tbdeath)) +
	geom_bar() +
	coord_flip()

ggplot(df2[df2$tbdeath==1,], aes(age.cat, fill = tbdeath)) +
	geom_bar() +
	coord_flip() +
	facet_wrap(~tb_type2)

ggplot(df2[df2$tbdeath==1,], aes(age.cat, fill = tbdeath)) +
	geom_bar() +
	coord_flip() +
	facet_wrap(~origin)



# logistic regression ----------------------------------------------------

model <- glm(tbdeath ~ .,
		  data = df2,
		  family = "binomial"
)

# print results
# summary(model)


# Predict CFR by age -----------------------------------------------------------------

data_ptb <- data.frame( age.cat = levels(df2$age.cat) )
data_ptb$origin<-"non_uk_born"
data_ptb$tb_type2<-"pulmonary"

data_eptb <- data.frame( age.cat = levels(df2$age.cat) )
data_eptb$origin<-"non_uk_born"
data_eptb$tb_type2<-"extrapulmonary"

# predict probability of TB death
pred_ptb <- predict(model,
				newdata = data_ptb,
				type = "response")

pred_eptb <- predict(model,
				newdata = data_eptb,
				type = "response")

# Interpolate to get per year CFRs
cfr_ptb <- approx(
	x = dat$mid_age[1:17],
	y = pred_ptb,
	xout = seq(16,100),
	method = "linear"
)

cfr_eptb <- approx(
	x = dat$mid_age[1:17],
	y = pred_eptb,
	xout = seq(16,100),
	method = "linear"
)


# extract standard error of the intercept
st_err<-summary(model)$coefficients[1, 2]

# # Sample from trubcated normal (to avoid negatives)
# rtruncnorm <- function(n, mu, sigma, low, high) {
# 	# find quantiles that correspond the the given low and high levels.
# 	p_low <- pnorm(low, mu, sigma)
# 	p_high <- pnorm(high, mu, sigma)
#
# 	# draw quantiles uniformly between the limits and pass these
# 	# to the relevant quantile function.
# 	if(n>1){
# 		outres<-qnorm(runif(n, p_low, p_high), mu, sigma)
# 	}else{
# 		outres<-1
# 	}
# 	return(outres)
#
# }
# rands <- rtruncnorm(100, 1, st_err, low = 0, high = Inf)
#
#
#
# ptb<-cfr_ptb$y%*%t(rands)
# eptb<-cfr_eptb$y%*%t(rands)
#
# ages<-dat$mid_age[1:17]
#
# matplot(cfr_ptb$x,
# 		ptb,
# 		type="l",
# 		col="grey45",
# 		main="Pulmonary TB",
# 		ylab = "EPTB CFR",
# 		xlab = "Age",
# 		ylim = c(0,0.5))
# points(ages,
# 	   dat$prop[1:17],
# 	   col="blue",
# 	   pch =19)
#
#
#
# matplot(cfr_eptb$x,
# 		eptb,
# 		type="l",
# 		col="grey45",
# 		main="Extra-Pulmonary TB",
# 		ylab = "EPTB CFR",
# 		xlab = "Age",
# 		ylim = c(0,0.5))
# points(ages,
# 	   dat$prop[20:36],
# 	   col="blue",
# 	   pch =19)

# 1. age, sex and chest pain on prob of disease
# plot_model(model,
# 		   type = "pred",
# 		   terms = c("age.cat", "origin", "tb_type2")#,
# 		   #ci.lvl = NA # remove confidence bands
# ) +
# 	labs(y = "Prob(tb death)")


dat<-list(
	ptb   = cfr_ptb,
	eptb   = cfr_eptb,
	st_err= st_err)



# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(dat, outfile)
