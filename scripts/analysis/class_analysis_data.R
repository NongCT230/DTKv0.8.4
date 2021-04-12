

# Basic statistical data  --------------------------------------------


#v0.7

#class contains basic statistical data and data calculating functions
# - minac50, minoed, target family counts, min and avg ac50 per target family, scalar top

#the design goal with this is that this is part of a composition where the plotting classes need their data as an attribute.
#this will be faster because it erases the need to copy paste. This does not mean that a class can't have their specific data.

Class.Analysis.Data <- R6Class("Class.Analysis.Data",
  
  #private variables and functions
	#those are private because the data should only ever be calculated with the functions. This allows us to never change this. but you have to use getters.
  private = list(
  
	ttl_chem = NULL, #chemicals in the basic stats table
	hitless_chem_name = NULL, #chemicals without any hits (or that don't have enough data to show some stuff.
	hitless_chem_casn = NULL, #chemicals without any hits (or that don't have enough data to show some stuff.
	
    basic_stat_tbl = NULL, #analysis table,
    target_family_counts = NULL, #function for data filtering,
    target_family_ac_tbl = NULL, #target family counts table,
    scalar_top_tbl = NULL, #target family chem activity table – min, avg values per family,
	  
	  #TODO check this part
	  gentox_tbl = NULL, #table for GenToxData
	  bmd_tbl = NULL #table for BMD Data
  ),
  
  #public variables and functions
  public = list(
    
    #constructor
    initialize = function( chemicals = NULL) {	
    },
    
    #finalizer
    finalize = function() {
    },
    
    
    #builds a stat dataframe from a list of chemical objects
    buildBasicStatsTable = function( chemicals, analyse_background = TRUE ){

	  
	  private$ttl_chem = length(chemicals);
	  
	  
      basic_stat_tbl <- tribble(~casn, ~aeid, ~assay_component_endpoint_name, ~intended_target_family,
                                ~intended_target_family_sub, ~ac50, ~ac_cutoff, ~ac_top, ~name, ~above_cutoff, ~cytotoxicity_um, ~cytotoxic, ~hitc, ~gene_name, ~gene_symbol, ~burst_assay);
      for (chem in chemicals){
        if (length(chem$assay_info$aeid)> 0){
          assaydata <- select(chem$assay_info, casn, aeid, assay_component_endpoint_name, 
                              intended_target_family, intended_target_family_sub, ac50, ac_cutoff, ac_top, gene_name,gene_symbol, hitc, burst_assay) %>% filter(hitc == 1); #this is the only spot we need to filter on. for stats. unless....
		  if(!analyse_background){
			assaydata  <-  subset(assaydata, intended_target_family != "background measurement");
		  }
		  if(nrow(assaydata)>0){			
			  assaydata <- add_column(assaydata, name = chem$chem_info$name);
			  assaydata <- add_column(assaydata, cytotoxicity_um = chem$chem_info$cytotoxicity_um);
			  ab_c <- ifelse(assaydata$ac50 > assaydata$ac_cutoff, "Y", "N");
			  ab_cito <- ifelse(10^assaydata$ac50 > chem$chem_info$cytotoxicity_um, "Y", "N");
			  assaydata <- add_column(assaydata, above_cutoff = ab_c, cytotoxic = ab_cito);
			  basic_stat_tbl <- bind_rows(basic_stat_tbl, assaydata);
		  }else{
			logdebug(paste0("No assay information found for ", chem$chem_info$casn));
			self$addHitlessChemInfo(chem$chem_info);
		  }
        } else {
          logdebug(paste0("No assay information found for ", chem$chem_info$casn));
          self$addHitlessChemInfo(chem$chem_info);
        }
		
      }
      private$basic_stat_tbl <- basic_stat_tbl;
      invisible (basic_stat_tbl);
    },
    
	
	
    #add OED results to basic stats table
    addOEDtoBasicStatsTable = function( chemicals, oed_model_name ){
      oed_tbl_sig <- tribble(~casn, ~aeid, ~oed);
      for (chem in chemicals){
        
        if (!is.null(chem$model_results[[oed_model_name]]) && length(chem$model_results[[oed_model_name]]$oed) > 0){
          oeddata <- select(chem$model_results[[oed_model_name]], aeid, oed);
          oeddata <- add_column(oeddata, casn = chem$chem_info$casn);
          oed_tbl_sig <- bind_rows(oed_tbl_sig, oeddata);
        }else{
          logdebug(paste0("No OED results found for ", chem$chem_info$casn, " for model ", oed_model_name));
          next;
		  #gabriel where we might add the httk stuff?
        }
        
      }
	  if(nrow(private$basic_stat_tbl) > 0){
		  private$basic_stat_tbl <- left_join(private$basic_stat_tbl, oed_tbl_sig, by=c("casn", "aeid"));
	  }
    },
    
    #target family counts
    computeTargetFamilyCounts = function(){
      
      target_family_count_tbl <- select(private$basic_stat_tbl, casn, name, intended_target_family) %>% drop_na(intended_target_family);
      
      target_family_counts <- plyr::count(target_family_count_tbl, vars="intended_target_family"); #colnames: intended_target_family, n
      
      private$target_family_counts <- target_family_counts;
      
      invisible ( target_family_counts );
    },
    
    #target family minimum and mean ac50
    computeTargetFamilyMinAndMeanAc50 = function( target_families ){
      
      target_family_ac_tbl <- tribble(~intended_target_family, ~avg_ac50, ~min_ac50);
      for ( target_family in target_families ){
        #Gabriel not how drop_na really works. gotta test if it works here.
        df <- filter(private$basic_stat_tbl, intended_target_family == target_family) %>% 
          drop_na(ac50);
        
        df$ac50 <- 10^(df$ac50); #convert from log10 form to uM
        avg_ac50 <- mean( df$ac50, na.rm=TRUE); #arithmetic mean
        min_ac50 <- min( df$ac50, na.rm=TRUE);
        
		
        target_family_ac_tbl <- add_row(target_family_ac_tbl, 
                                        intended_target_family = target_family, avg_ac50 = avg_ac50, min_ac50 = min_ac50);
        
      }
      if (is.null(private$target_family_counts)){
        self$computeTargetFamilyCounts();
      }

	if(nrow(target_family_ac_tbl) > 0 && nrow(private$target_family_counts > 0)){
      
	  target_family_ac_tbl <- full_join(private$target_family_counts, target_family_ac_tbl, by = "intended_target_family");
      
      private$target_family_ac_tbl <- target_family_ac_tbl;
	  
    }else{
		private$target_family_ac_tbl <- target_family_ac_tbl;
	}
    },
    #calculate scalar top values
    computeScalarTop = function (){
      scalar_top_tbl <- select(private$basic_stat_tbl, casn, name, aeid, 
                               assay_component_endpoint_name, ac50, ac_cutoff, 
                               ac_top, above_cutoff,cytotoxic, gene_name, gene_symbol, contains("oed"));
      scalar_top <- scalar_top_tbl$ac_top / 10^(scalar_top_tbl$ac50); #ac_top needs further verification, as per documentation it is not log10, but looking at data who knows
      #Gabriel
      scalar_top_tbl <- add_column(scalar_top_tbl, scalar_top = scalar_top);
      private$scalar_top_tbl <- scalar_top_tbl;
      
      invisible ( scalar_top_tbl );
    },
	
	#build the GenTox Data Table
	#TODO will need to be checked
	computeGenToxTable = function(chemicals)
	{
	  load("data/GenTox21.RData");
	  
	  gentox_tbl <- tribble(~casn, ~compound, ~assay, ~endpoint, ~s9, 
	                            ~concentration, ~concentration_unit, ~response, ~response_type);
	  for (chem in chemicals){
	    if (length(chem$assay_info$aeid)> 0){
	      chemGenToxData <- gen_tox_data[gen_tox_data$casn == chem$chem_info$casn,] %>%
	        dplyr::mutate(s9 = ifelse(s9_positive, "+", "-"));
	      chemGenToxData <- chemGenToxData[, names(chemGenToxData) != "s9_positive"];
	      if(nrow(chemGenToxData)>0){			
	        gentox_tbl <- bind_rows(gentox_tbl, chemGenToxData);
	      }else{
	        logdebug(paste0("No GenTox information found for ", chem$chem_info$casn));
	        self$addHitlessChemInfo(chem$chem_info);
	      }
	    } else {
	      logdebug(paste0("No GenTox information found for ", chem$chem_info$casn));
	      self$addHitlessChemInfo(chem$chem_info);
	    }
	    
	  }
	  private$gentox_tbl <- gentox_tbl;
	  invisible (gentox_tbl);
	},
	
	#build the BMD Data Table
	#TODO will need to be checked
	computeBMDTable = function(chemicals)
	{
	  load("data/GenTox21.RData");
	  
	  bmd_tbl <- tribble(~casn, ~compound, ~assay, ~endpoint, ~s9, 
	                        ~bmr, ~bmdl, ~bmd, ~bmdu);
	  for (chem in chemicals){
	    if (length(chem$assay_info$aeid)> 0){
	      chemBMDData <- bmd_data[bmd_data$casn == chem$chem_info$casn,] %>%
	        dplyr::mutate(s9 = ifelse(s9_positive, "+", "-"));
	      chemBMDData <- chemBMDData[, names(chemBMDData) != "s9_positive"];
	      if(nrow(chemBMDData)>0){			
	        bmd_tbl <- bind_rows(bmd_tbl, chemBMDData);
	      }else{
	        logdebug(paste0("No BMD information found for ", chem$chem_info$casn));
	        self$addHitlessChemInfo(chem$chem_info);
	      }
	    } else {
	      logdebug(paste0("No BMD information found for ", chem$chem_info$casn));
	      self$addHitlessChemInfo(chem$chem_info);
	    }
	    
	  }
	  private$bmd_tbl <- bmd_tbl;
	  invisible (bmd_tbl);
	},
	
	addHitlessChemInfo = function(chem){
		if(!(chem$casn %in% private$hitless_chem_casn)){
			private$hitless_chem_name = c(private$hitless_chem_name, chem$name);
			private$hitless_chem_casn = c(private$hitless_chem_casn, chem$casn);
		}
	},
    
	getHitlessChemInfoByName = function(){
		return (private$hitless_chem_name);
	},
	getHitlessChemInfoByCasn = function(){
		return (private$hitless_chem_casn);
	},
	
    getBasicStatsTable = function (){
      return (private$basic_stat_tbl);
    },
    
    getTargetFamilyCountsTable = function (){
      return (private$target_family_counts);
    },
    
    getTargetFamilies = function (){
      if (!is.null(private$target_family_counts)){
        return (private$target_family_counts$intended_target_family);
      } else {
        return (list());
      }
      
    },
    
    getTargetFamilyAcTable = function (){
      return (private$target_family_ac_tbl);
    },
    
    getScalarTopTable = function (){
      return (private$scalar_top_tbl);
    },
	
	  getGenToxTable = function () {
	    return (private$gentox_tbl);
	  },
	
	  getBMDTable = function () {
	    return (private$bmd_tbl);
	  },
	
    basicStatsDataExists = function (){
      if (is.null(private$basic_stat_tbl) || length(private$basic_stat_tbl$casn) <= 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
    
    targetFamilyCountsExist = function (){
      if (is.null(private$target_family_counts) || length(private$target_family_counts$freq) <= 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
    
    targetFamilyAcDataExists = function (){
      if (is.null(private$target_family_ac_tbl) || length(private$target_family_ac_tbl$avg_ac50) <= 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
    
    scalarTopDataExists = function (){
      if (is.null(private$scalar_top_tbl) || length(private$scalar_top_tbl$scalar_top) <= 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
    
    oedValuesExist = function(){
      if(!self$basicStatsDataExists() ||  is.null(private$basic_stat_tbl$oed) || all(is.na(private$basic_stat_tbl$oed))){
        return (FALSE);
      } else {
        return (TRUE);
      }
    }
    
  )
)
