
# Data structure: Chemical -------------------------------------------

# v0.7

#General parameters we want:

#   chem_info: casn, name, target, cytotoxic concentration, toxicity flag

#
#   assay_info: aeid, name, component name, component endpoints, organism, tissue,
# cell, biological process target, intended target family, ac50, ac cutoff, ac top asymptote
#
#   model_data: list with some data objects for each model keyed by model name
#
#   model_results: list with some data objects for each model results keyed by model name
#
#   analysis_results: list with some data objects for each analysis operation keyed by analysis name
 

Class.Chemical <- R6Class(
  "Class.Chemical",
  
  lock_class = TRUE,
  
  #self variables and functions
  private = list(),
  
  #public variables and functions
  public = list(
    classname = NULL, #we want to use this to id class instead of typeof(object) which returns "environment" for all R6 classes
    
    chem_info = NULL,
    assay_info = NULL,
    model_data = NULL,
    model_results = NULL,
    analysis_results = NULL,


    #constructor
    initialize = function(){
      self$classname = "Class.Chemical";
      
      #not initializing tables to default values saves a whole lot of computation time for long chem lists!
      #self$assay_info = tibble(aeid = 0, assay_name = "", assay_component_name = "", assay_component_endpoint_name ="",
      #                         organism = "", tissue = "", cell_short_name = "", biological_process_target = "", 
      #                         intended_target_family = "", ac50 = 0, ac_top = 0, ac_cutoff = 0);
      self$chem_info = list(casn = "", name = "", cytotoxicity_um = NA, cytotoxic = NA);
      self$assay_info = tibble();
      self$model_data = list();
      self$model_results = list();
      self$analysis_results = list();

    },
    
    #finalizer
    finalize = function() {},

	#calculates the amount of hits for a Chemical object, saves analysis results to Chemical object
	calculateHitCount = function(){
		 if (!is.null(self$assay_info) && length(self$assay_info$aeid) > 0){
			#shortcut being used here. we don't need to sum we only need the length as all hitc = 1
			hits = select(self$assay_info, hitc) %>% filter(hitc == 1);
			
			self$analysis_results$hitcount = length(hits$hitc);
			
		 }else{
			self$analysis_results$hitcount = 0;
		 }
		 
	},
	#gabriel maybe should remake the whole chemical info box so that it filters better.
    #calculates minimum Ac50 for a Chemical object, saves analysis results to Chemical object
    calculateMinAc50 = function(back,hit){
	
      if (!is.null(self$assay_info) && length(self$assay_info$aeid) > 0){
        withCallingHandlers( #catch warnings that clutter up the console
          {
            valid_ac50 <- self$assay_info;
            valid_ac50$ac50 <- 10^(valid_ac50$ac50); #ac50 are log10 values, convert to regular value in uM units
			if(!back){
				valid_ac50 = subset(valid_ac50, intended_target_family != "background measurement");
			}
			if(hit){
				valid_ac50 = subset(valid_ac50, hitc == 1);
			}
            min_ac50 <- min(valid_ac50$ac50, na.rm = TRUE)[[1]]; #grab minimum ac50 values
            
            #find aeid for this value and set it to self
            valid_ac50 <- filter(valid_ac50, ac50 == min_ac50);
            if(length(valid_ac50$aeid)>0){
              min_aeid <- valid_ac50$aeid[[1]];
            } else {
              min_aeid <- NA;
            }
          },
          warning=function(w) { 
            logdebug(str_c("calculateMinAc50() warning:", w));
            invokeRestart("muffleWarning");
          }
        )
      } else {
        min_ac50 <- NA;
        min_aeid <- NA;
      }
      if( !is.finite(min_ac50) ){
        min_ac50 <- NA;
        min_aeid <- NA;
      }
      if(!is.na(min_ac50)){
        min_ac50_units <- "uM";
      } else {
        min_ac50_units <- "";
      }
      
      self$analysis_results$min_ac50 <- list(value = min_ac50, aeid = min_aeid, units = min_ac50_units);
      invisible (min_ac50);
    },
	
	 calculateAVGAc50 = function(back,hit){
      if (!is.null(self$assay_info) && length(self$assay_info$aeid) > 0){
        withCallingHandlers( #catch warnings that clutter up the console
          {
		  
			
            valid_ac50 <- self$assay_info;
            valid_ac50$ac50 <- 10^(valid_ac50$ac50); #ac50 are log10 values, convert to regular value in uM units
			if(!back){
				valid_ac50 = subset(valid_ac50, intended_target_family != "background measurement");
			}
			if(hit){
				valid_ac50 = subset(valid_ac50, hitc == 1);
			}

            avg_ac50 <- mean(valid_ac50$ac50, na.rm = TRUE)[[1]]; #grab average ac50 values
            
          },
          warning=function(w) { 
            logdebug(str_c("calculateAVGAc50() warning:", w));
            invokeRestart("muffleWarning");
          }
        )
      } else {
        avg_ac50 <- NA;
      }
      if( !is.finite(avg_ac50) ){
        avg_ac50 <- NA;
      }
      if(!is.na(avg_ac50)){
        avg_ac50_units <- "uM";
      } else {
        avg_ac50_units <- "";
      }
      
      self$analysis_results$avg_ac50 <- list(value = avg_ac50, units = avg_ac50_units);
      invisible (avg_ac50);
    },
	
    #calculate minimum OED for a Chemical object given a certain name for the OED model,saves analysis results to Chemical object
    calculateMinOED = function (oed_model_name ){
      if (!is.null(self$model_results[[oed_model_name]]) && length(self$model_results[[oed_model_name]]$oed) > 0){
        withCallingHandlers(
          {
		  #thinking Gabriel
            valid_oed <- filter(self$model_results[[oed_model_name]], background != "Y")$oed;
            min_oed <- min(valid_oed, na.rm = TRUE)[[1]]; #grab first minimum oed value 
          },
          warning=function(w) { 
            logdebug(str_c("Chemical list selection warning:", w));
            invokeRestart("muffleWarning");
          }
        )
      } else {
        min_oed <- NA;
      }
      if(!is.finite(min_oed)){
        min_oed <- NA;
      }
      if(!is.na(min_oed)){
        min_oed_units <- self$model_results[[oed_model_name]]$oed_units[[ which(self$model_results[[oed_model_name]]$oed == min_oed) ]];
      } else {
        min_oed_units <- "";
      }
      
      self$analysis_results$min_oed <- list(value = min_oed, units = min_oed_units);
      invisible (min_oed);
    },
	
	#calculate mean OED for a Chemical object given a certain name for the OED model,saves analysis results to Chemical object
	calculateMeanOED = function (oed_model_name ){
	  if (!is.null(self$model_results[[oed_model_name]]) && length(self$model_results[[oed_model_name]]$oed) > 0){
	    withCallingHandlers(
	      {
	        #thinking Gabriel
	        valid_oed <- filter(self$model_results[[oed_model_name]], background != "Y")$oed;
	        mean_oed <- mean(valid_oed, na.rm = TRUE)[[1]]; 
	      },
	      warning=function(w) { 
	        logdebug(str_c("Chemical list selection warning:", w));
	        invokeRestart("muffleWarning");
	      }
	    )
	  } else {
	    mean_oed <- NA;
	  }
	  if(!is.finite(min_oed)){
	    mean_oed <- NA;
	  }
	  if(!is.na(min_oed)){
	    mean_oed_units <- self$model_results[[oed_model_name]]$oed_units[[ which(self$model_results[[oed_model_name]]$oed == mean_oed) ]];
	  } else {
	    mean_oed_units <- "";
	  }
	  
	  self$analysis_results$mean_oed <- list(value = mean_oed, units = mean_oed_units);
	  invisible (mean_oed);
	}
	
	

  )
  
)
