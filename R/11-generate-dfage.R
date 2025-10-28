# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root <- here::here()
	infile  <- file.path(root, "output", "parameters.qs2")
	outfile <- file.path(root, "output", "dfage.qs2")

	# Packages
	source(file.path(root, "R", "modify_attach.R"))
	modify_attach(qs2,          include.only = c("qs_read", "qs_save"))
	modify_attach(fitdistrplus, include.only = "fitdist")
	modify_attach(truncdist,    include.only = "rtrunc")

} else {

	# Input
	args <- commandArgs(trailingOnly = TRUE)
	infile  <- args[1]
	outfile  <- args[2]

	# Packages
	library(qs2,          include.only = c("qs_read", "qs_save"))
	library(fitdistrplus, include.only = "fitdist")
	library(truncdist,    include.only = "rtrunc")
}

# -------------------------------------------------------------------------                                   │    │
# Load parameters ---------------------------------------------------------                                   │    │
# -------------------------------------------------------------------------                                   │    │
parameters <- qs_read(infile)


# -------------------------------------------------------------------------
# Generate dfage ----------------------------------------------------------
# -------------------------------------------------------------------------
age_dist    <- parameters$age_dist
cohort_size <- parameters$cohort_size

dfage <- with(age_dist, {
	ages <- rep(age_point, freq * 100)
	fit.gamma <- fitdist(ages, "gamma")
	v <- round(rtrunc(
		n = cohort_size,
		spec = "gamma",
		a=15,
		b=100,
		shape=fit.gamma$estimate[["shape"]],
		rate=fit.gamma$estimate[["rate"]]
	))
	agecounts <- hist(v, breaks=85, plot = FALSE)
	out <- data.frame(age = agecounts$breaks[-1], n = agecounts$counts)
	out$value <- out$n / cohort_size

	#if (mode=="test"){

	# Set all population at 25 years old (mid point betwen 16-35)
	# TODO - Check this with JV. For the age point above we were actually using
	#        26 (not 25) due to the rounding.
	# TODO - is this staying. If so let's delete the rest.
	out$value <- 0
	out$value[round((35-16)/2)] <- 1
	#}
	out
})

# -------------------------------------------------------------------------
# Save dfage --------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(dfage, outfile)
