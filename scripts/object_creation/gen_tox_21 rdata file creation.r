gen_tox_data <- data.frame(
	casn=character(),
	compound=character(),
	assay=character(),
	endpoint=character(),
	s9_positive=logical(),
	concentration=double(),
	concentration_unit=character(),
	response=double(),
	response_type=character(),
	stringsAsFactors=FALSE
);

bmd_data <- data.frame(
	casn=character(),
	compound=character(),
	assay=character(),
	endpoint=character(),
	s9_positive=logical(),
	model_type=character(),
	aic=double(),
	log_likelihood=double(),
	bmr=double(),
	bmdl=double(),
	bmd=double(),
	bmdu=double(),
	stringsAsFactors=FALSE
);

trend_data <- data.frame(
  casn=character(),
  compound=character(),
  assay=character(),
  endpoint=character(),
  s9_positive=logical(),
  trend=logical(),
  stringsAsFactors=FALSE
);

assay_description <- data.frame(
  Assay = c("Ames II", "TGR", "MicroFlow", "CometChip", "MultiFlow"),
  Description = c(
    "Bacterial reverse mutation assay",
    "Mammalian cell mutagenicity assay",
    "Mammalian cell chromosome abnormality",
    "Mammalian cell DNA damage",
    "Mammalian cell genotoxicity biomarker"
  ),
  stringsAsFactors = FALSE
);

save(gen_tox_data, bmd_data, trend_data, assay_description, file="data/GenTox21.RData");