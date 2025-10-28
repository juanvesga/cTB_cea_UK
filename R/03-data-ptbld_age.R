# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root    <- here::here()
	infile  <- file.path(root, "data", "raw", "ptbld_age.csv")
	outfile <- file.path(root, "output", "ptbld_age.qs2")

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

# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(dat, outfile)
