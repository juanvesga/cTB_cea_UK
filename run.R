# Check packages ----------------------------------------------------------
required <- c(
	"here",
	"qs2",
	"readxl",
	"freedom",
	"fitdistrplus",
	"truncdist",
	"data.table",
	"gridExtra",
	"odin",
	"tools",
	"matrixStats"
)

installed <- rownames(installed.packages())
not_installed <- required[!required %in% installed]

if (length(not_installed)) {
	msg <- "Some required packages are missing. Please install the following:\n"
	pkgs <- paste("\t", not_installed, collapse = "\n")
	msg <- sprintf("%s%s", msg, pkgs)
	stop(msg)
} else {
	# Set the root directory --------------------------------------------------
	root <- here::here()

	# List the fils we need to run in order
	files <- file.path(
		root, "R",
		c(
			"01-parameters.R",
			"02-data-conmat.R",
			"03-data-ptbld_age.R",
			"04-data-p_ptb_ets.R",
			"05-data-AE_rates_NICE.R",
			"06-data-age_rates.R",
			"07-data-p_cfr_ets.R",
			"08-data-qaly_input.R",
			"09-data-pars_cohort_discrete2_new.R",
			"10-generate-samples.R",
			"11-generate-dfage.R",
			"12-generate-qale.R",
			"13-generate-model_parameters.R",
			"14-run-incremental-case.R" # ,
			# "15-run-incremental-case-vary.R"
		)
	)

	# Loop over the files
	for (f in files) {
		cat(sprintf("Running %s\n", f))
		local(source(f, local = TRUE))
	}

	message("Analysis complete. Results are in the 'results/' directory.")
}

rm(list = ls())