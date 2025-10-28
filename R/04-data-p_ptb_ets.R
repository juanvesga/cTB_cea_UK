# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root    <- here::here()
	infile  <- file.path(root, "data", "raw", "ets_ptb_raw.csv")
	outfile <- file.path(root, "output", "p_ptb_ets.qs2")

	# Packages
	source(file.path(root, "R", "modify_attach.R"))
	modify_attach(qs2, include.only = "qs_save")

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
dat <- read.csv(infile)

# interpolation by year of age
p_ptb <- approx(
	x = dat$mid_age,
	y = dat$prop_ptb,
	xout = seq(16,100),
	method = "linear"
)


# Calculate the ratios
ratio <- c(1,p_ptb$y[2:length(p_ptb$y)] / p_ptb$y[1:(length(p_ptb$y) - 1)])

# Print the result



dat<-list(
	p_ptb   = p_ptb,
	st_err  = dat$se[1],
	ratio   = ratio
	)



# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(dat, outfile)
