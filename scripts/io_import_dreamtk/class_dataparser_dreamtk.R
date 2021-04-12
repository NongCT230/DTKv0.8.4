

# Data Parser class for dreamtk ----------------------------------------------


# v0.7

# creates a Chemical object, or a list of them



Class.DataParser.dreamtk <- R6Class("Class.DataParser.dreamtk",
  
  cloneable = FALSE,
  lock_class = TRUE,
  
  #self variables and functions
  private = list(
    
    chemtablevars = NULL,
    assaytablevars = NULL
    
  ),
  
  #public variables and functions
  public = list(
    
    #constructor
    initialize = function() {
      
      private$chemtablevars <- Class.DataImporter.dreamtk.tableinfo$chemtablevars;
      
      private$assaytablevars <- c(
        unlist( Class.DataImporter.dreamtk.tableinfo$ac50tablevars ),
        unlist( Class.DataImporter.dreamtk.tableinfo$assaytablevars )
      ) %>% unique ();
       
    #     list ("casn", "aeid", "ac50", "ac_top", "ac_cutoff", "assay_name", "assay_component_name", 
    #                           "assay_component_endpoint_name", "organism", "tissue", "cell_short_name", 
    #                           "biological_process_target", "intended_target_family", "intended_target_family_sub", "gene_name");
    },
    
    #finalizer
    finalize = function() {},
    
    createChemical = function( cheminfo = NULL, assayinfo = NULL ){
      
      if (is.null(cheminfo) || is.null(assayinfo)){
        e <- "DataParser.dreamtk$createChemical(): Chemical or assay info missing (null)!";
        logerror(e);
        stop(e);
      }

      if (is.atomic(cheminfo) || is.atomic(assayinfo) || 
          !all(private$private$chemtablevars %in% names(cheminfo)) || 
          !all(private$assaytablevars %in% names(assayinfo))){
        e <- "DataParser.dreamtk$createChemical(): Invalid chemical or assay info!";
        logerror(e);
        stop(e);
      }
      
      if (length(cheminfo$casn) > 1){
        e <- "DataParser.dreamtk$createChemical(): Table has too many chemicals. Use $createListOfChemicals() instead.";
        logerror(e);
        stop(e);
      } else if ( length(cheminfo$casn) == 0 ){
        e <- "DataParser.dreamtk$createChemical(): Table has no chemicals. Nothing to parse.";
        logerror(e);
        stop(e);
      }
      
      chemical <- Class.Chemical$new();
          
      chemical$chem_info <- list(casn = cheminfo$casn[1], 
                                 name = cheminfo$name[1], 
                                 cytotoxicity_um = cheminfo$cytotoxicity_um[1],
                                 cytotoxic = cheminfo$cytotoxic[1]);

      if (is.na(chemical$chem_info$name)){ #chemical without a name or that has not been found in the database marked with *
        chemical$chem_info$name <- paste0("*", chemical$chem_info$casn);
        chemical$chem_info$casn <- paste0("*", chemical$chem_info$casn);
        warning(paste0(chemical$chem_info$casn, " not found in database."));
      }
      
      chemical$assay_info <- filter(assayinfo, casn == chemical$chem_info$casn);   
      return( chemical );
          
    },
      
    createListOfChemicals = function( cheminfo = NULL, assayinfo = NULL ){
      
      logdebug("Parsing chemical list...");
      
      if (is.null(cheminfo) || is.null(assayinfo)){
        e <- "DataParser.dreamtk$createChemical(): Chemical or assay info missing (null)!";
        logerror(e);
        stop(e);
      }

      if (is.atomic(cheminfo) || is.atomic(assayinfo) || 
          !all(private$chemtablevars %in% names(cheminfo)) || 
          !all(private$assaytablevars %in% names(assayinfo))){
        e <- "DataParser.dreamtk$createChemical(): Invalid chemical or assay info!";
        logerror(e);
        stop(e);
      }
      
      if ( length(cheminfo$casn) == 0 ){
        e <- "DataParser.dreamtk$createChemical(): Table has no chemicals. Nothing to parse.";
        logerror(e);
        stop(e);
      }
      
      #create chemlist using for. Same efficiency, better readability
      chemlist <- vector("list", length(cheminfo$casn));
      for (i in 1:length(cheminfo$casn)){
        chemlist[[i]] <- Class.Chemical$new();

        chemlist[[i]]$chem_info <- list(casn = cheminfo$casn[i],
                                        name = cheminfo$name[i],
                                        cytotoxicity_um = cheminfo$cytotoxicity_um[i],
                                        cytotoxic = cheminfo$cytotoxic[i]);

        if (is.na(chemlist[[i]]$chem_info$name)){ #chemical without a name or that has not been found in the database marked with *
          chemlist[[i]]$chem_info$name <- paste0("*", chemlist[[i]]$chem_info$casn);
          chemlist[[i]]$chem_info$casn <- paste0("*", chemlist[[i]]$chem_info$casn);
          warning(paste0(chemlist[[i]]$chem_info$casn, " not found in database."));
        }

        chemlist[[i]]$assay_info <- filter(assayinfo, casn == chemlist[[i]]$chem_info$casn);
      }
      
      return( chemlist );
      
    }

  )
  
)
