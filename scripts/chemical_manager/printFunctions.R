#Functions to print chemical information outside of a GUI

#v0.7

#print chemical overview
printChemicalOverview = function( chemical = NULL, listpos = 1 ) {
  print( 
    paste0("[", listpos, "] ", chemical$chem_info$casn, " / ", chemical$chem_info$name )
  );
  print(
    paste0("    # of assay endpoints: ", nrow(chemical$assay_info))
  );
  print(
    paste0("    available model data: ", paste0( names(chemical$model_data), collapse=", ") )
  );
  print(
    paste0("    available model results: ", paste0( names(chemical$model_results), collapse = ", ") )
  );
  print(
    paste0("    available analysis results: ", paste0( names(chemical$analysis_results), collapse=", ") )
  );
  print("----------------------------------------------");
}

printChemicalOverviewList = function( chemlist = NULL ) {
  listpos <- 1;
  for (chem in chemlist) {
    printChemicalOverview ( chem, listpos );
    listpos <- listpos + 1;
  }
}


#obtain basic chemical information
printChemInfo = function( chem, css_model_name, oed_model_name, linesep = ", " ){
  #obtain css value and units
  css <- chem$model_results[[css_model_name]]$value;
  css_units <- chem$model_results[[css_model_name]]$units;
  css_assumptions <- paste0( paste0(chem$model_results[[css_model_name]]$assumptions, collapse = linesep));
  if(is.na(css)){
    css_units <- "";
  }
  
  #obtain minimum oed value and units
  if ( !is.null(chem$analysis_results$min_oed) ) {
    min_oed <- chem$analysis_results$min_oed$value;
  } else {
    min_oed <- NA;
  }
  min_oed_units <- chem$analysis_results$min_oed$units;
  
  #obtain minimum ac50 value and units
  if ( !is.null(chem$analysis_results$min_ac50) ) {
    min_ac50 <- chem$analysis_results$min_ac50$value;
  } else {
    min_ac50 <- NA;
  }
  min_ac50_units <- chem$analysis_results$min_ac50$units;
  
  aei <- chem$analysis_results$min_ac50$aeid;
  if ( !is.null(aei) && !length(aei) == 0){
    chemassayrow <- filter(chem$assay_info, aeid == aei);
  } else {
    chemassayrow <- list(assay_component_endpoint_name = NA,
                         organism = NA,
                         tissue = NA,
                         cell_short_name = NA,
                         biological_process_target = NA,
                         intended_target_family = NA,
                         intended_target_family_sub = NA,
                         gene_name = NA)
  }
  
  adj_casn <- chem$chem_info$casn;
  adj_name <- chem$chem_info$name;
  if (startsWith(adj_casn, "*") ){
    adj_casn <- gsub(pattern = "[*]", replacement = "", adj_casn);
    adj_name <- gsub(pattern = "[*]", replacement = "", adj_name);
  }
  
  datarow <- list(casn = adj_casn,
                  name = adj_name,
                  cytotoxic = chem$chem_info$cytotoxic,
                  cytotoxicity_um = chem$chem_info$cytotoxicity_um,
                  css = css,
                  css_units = css_units,
                  css_model = css_model_name,
                  css_assumptions = css_assumptions,
                  min_ac50 = min_ac50,
                  min_ac50_units = min_ac50_units,
                  assay_component_endpoint_name = chemassayrow$assay_component_endpoint_name,
                  organism = chemassayrow$organism,
                  tissue = chemassayrow$tissue,
                  cell_short_name = chemassayrow$cell_short_name,
                  biological_process_target = chemassayrow$biological_process_target,
                  intended_target_family = chemassayrow$intended_target_family,
                  intended_target_family_sub = chemassayrow$intended_target_family_sub,
                  gene_name = chemassayrow$gene_name,
                  min_oed = min_oed,
                  min_oed_units = min_oed_units,
                  oed_model = oed_model_name);
  
  print( 
    paste0(chem$chem_info$casn, " / ", chem$chem_info$name )
  );
  for (name in names(datarow)){
    print( paste0("  ", name, ": ", datarow[[name]])  );
  }
  
  invisible ( datarow );
}

#print chemical assays
printChemicalAssays = function ( chem ){
  assays <- distinct(chem$assay_info, assay_name) %>% arrange();
  print( 
    paste0(chem$chem_info$casn, " / ", chem$chem_info$name )
  );
  assaytree <- list();
  alistpos <- 1;
  for (assay in assays$assay_name){
    assaytree[[assay]] <- list();
    print( paste0("  [", alistpos, "] ", assay) );
    alistpos <- alistpos + 1;
    assay_components <- filter(chem$assay_info, assay_name == assay ) %>%
      select(assay_component_endpoint_name) %>% arrange();
    aclistpos <- 1;
    for (assay_component in assay_components$assay_component_endpoint_name){
      assaytree[[assay]][[assay_component]] <- assay_component;
      print( paste0("    [", aclistpos, "] ", assay_component) );
      aclistpos <- aclistpos + 1;
    }
  }
  
  invisible ( assaytree );
}


printAssayInfo = function(chem, assay_component_endpoint, oed_model_name ){

  #get assay info for selected endpoint
  datarow <- filter(chem$assay_info, assay_component_endpoint_name == assay_component_endpoint);
  
  #in case of duplicates, choose the first element only. this also de-lists the items
  organism <- datarow[["organism"]][[1]];
  tissue <- datarow[["tissue"]][[1]];
  cell <- datarow[["cell_short_name"]][[1]];
  bpt <- datarow[["biological_process_target"]][[1]];
  itf <- datarow[["intended_target_family"]][[1]];
  itfs <- datarow[["intended_target_family_sub"]][[1]];
  gn <- datarow[["gene_name"]][[1]];
  ac50 <- datarow[["ac50"]][[1]];
  acc <- datarow[["ac_cutoff"]][[1]];
  act <- datarow[["ac_top"]][[1]];
  aeid <- datarow[["aeid"]][[1]];
  
  oed <- chem$model_results[[oed_model_name]]$oed[ chem$model_results[[oed_model_name]]$aeid == aeid ];
  
  if(!is.na(oed)){
    oed_units <- chem$model_results[[oed_model_name]]$oed_units[ chem$model_results[[oed_model_name]]$aeid == aeid ];
  } else {
    oed_units <- "";
  }
  
  printlist <- list(assay_endpoint_name = assay_component_endpoint,
                    organism = organism,
                    tissue = tissue,
                    cell = cell,
                    bio_process_target = bpt,
                    int_target_family = itf,
                    int_target_subfamily = itfs,
                    gene = gn,
                    ac50 = signif(ac50, digits = 4),
                    ac_cutoff = signif(acc, digits = 4),
                    ac_top = signif(act, digits = 4),
                    oed = signif(oed, digits = 4),
                    oed_units = oed_units,
                    oed_model = oed_model_name);
  
  print( 
    paste0(chem$chem_info$casn, " / ", chem$chem_info$name, " / ", printlist$assay_endpoint_name )
  );
  for (name in names(printlist)){
    print( paste0("  ", name, ": ", printlist[[name]])  );
  }
  
  invisible (printlist);
  
}
