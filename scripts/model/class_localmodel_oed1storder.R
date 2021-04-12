

# Calculate OED (oral equivalent dose) ------------------------------------

#v0.7

# model based on:
#   Wetmore et al. (2015). Incorporating high-throughput exposure predictions with dosimetry-adjusted in vitro 
#     bioactivity to inform chemical toxicity testing. Toxicological Sciences, 148(1), 121-136.
# model assumptions:
#   constant chemical exposure rate at 1mg/kgBW/day, with assumption of linear relationship to chemical dose
#   uses in virto concentration at 50% of maximum activity (ac50) assays to estimate
#   eliminates ac50 values that are below assay cutoff threshold (noise threshold)

Class.LocalModel.OED1storder <- R6Class("Class.LocalModel.OED1storder",
                            
  inherit = Class.LocalModel,

  #self variables and functions
  private = list(
    cssmodel = NULL,
    activity_variable = NULL,
    
    #acquires, verifies and  returns Css data from a given model
    getCssDataFromChemical = function ( chemical=NULL, modelname = NULL ){
      
      cssdata <- private$getModelResultsFromChemical( chemical, modelname );

      if(is.null(cssdata)){
        e <- paste0("LocalModel.OED1storder$getCssFromChemical(): Css model result for chemical ", chemical$chem_info$casn, 
                    " does not exist. Could not find result for model: ", modelname, 
                    ". Run the model on the Chemical or assign results manually.");
        logerror(e);
        stop(e);
      };
      
      if ( is.atomic(cssdata) || 
           is.null( cssdata$parameter ) || 
           is.null( cssdata$value ) ||
           is.null( cssdata$unit ) ||
           is.na( cssdata$parameter ) ||
           cssdata$parameter != "css"){
        e <- "LocalModel.OED1storder$getCssFromChemical(): The specified css model result exists, 
                but does not have expected structure. Expected to have at least $parameter='css', $value, and $unit.";
        logerror(e);
        stop(e);
      }
      
      if (is.na( cssdata$unit ) || 
          cssdata$unit != "uM"){
        logwarn("LocalModel.OED1storder$getCssFromChemical(): This model only accepts Css in uM units. Returning NA.");
        return (NA);
      }

      return ( cssdata$value );
    },
    #Gabriel this could be done by the chemical itself.
    #acquires, verifies and returns ac50 data from a given chemical
    #NOTE: ac50 from data are in log10 form, we need to return original uM units
    getAc50DataFromChemical = function ( chemical=NULL ){
      
      if(is.null(chemical$assay_info)){
        e <- str_c("LocalModel.OED1storder$getAc50DataFromChemical(): Assay data for chemical ", chemical$chem_info$casn, 
                   " does not exist. DataImporter...$importAssayInfo...() first or assign assay table to Chemical manually");
        logerror(e);
        stop(e);
      };
      
      if ( is.atomic(chemical$assay_info) || 
           is.null( chemical$assay_info$aeid ) || 
           is.null( chemical$assay_info$ac50 ) ||
           is.null( chemical$assay_info$ac_cutoff) ){
        e <- "LocalModel.OED1storder$getAc50FromChemical(): Chemical assay data exists, 
                but does not have expected structure. Expected to have at least $aeid, $ac50, and $ac_cutoff.";
        logerror(e);
        stop(e);
      }
      
      trimmed_ac50 <- filter ( chemical$assay_info, ac50 > ac_cutoff) %>%
        select( aeid, ac50, intended_target_family ); #why intended_target_family too? because id background measurement!
      
      #ac50 is in log10 form, so it needs to be transformed back into uM units
      trimmed_ac50$ac50 <- 10^( trimmed_ac50$ac50 );

      if( length(trimmed_ac50$aeid) < length(chemical$assay_info$aeid) ) {
        logwarn("LocalModel.OED1storder$getAc50FromChemical(): Ac50 values below cutoff were ignored. 
                OED for these assays will be NA.");
      }
      
      return ( trimmed_ac50 );
      
      
    }
  ),
  
  #public variables and functions
  public = list(

    #constructor
    initialize = function() {
      super$initialize();
      private$model_name <- "oed.1storder";
      private$activity_variable <- "ac50";
    },
    
    #finalizer
    finalize = function() {},
    
    #sets model parameters
    setParameters = function ( cssmodel = NULL, overwrite_result = FALSE, activity_variable = "ac50" ){
      if (is.null(cssmodel) || !is.character(cssmodel)){
        e <- "LocalModel.OED1storder$setParameters(): Css model identifier is a character string. 
                Use cssmodel = MyCssModel$getModelName() for the LocalModel.Css... you have used";
        logerror(e);
        stop(e);
      }
      
      if (is.null(activity_variable) || !is.character(activity_variable)){
        e <- "LocalModel.OED1storder$setParameters(): activity_variable parameter must be a character string. Default is 'ac50'";
        logerror(e);
        stop(e);
      }
      
      private$cssmodel <- cssmodel;
      private$activity_variable <- activity_variable;
      
      if(is.logical(overwrite_result)){
        private$overwrite_result <- overwrite_result;
      }
      
    },
    
    #runs the  model
    compute = function (){
      
      if(is.null(private$chemicals)){
        e <- "LocalModel.OED1storder$compute(): Can not compute model without a Chemical list. $setChemicals() first";
        logerror(e);
        stop(e);
      };
      
      if(is.null(private$cssmodel)){
        e <- "LocalModel.OED1storder$compute(): Css model identifier not given. 
                Use LocalModel.OED...$setParameters() to set the model identifier";
        logerror(e);
        stop(e);
      };
      
      #compute model
      for ( chem in private$chemicals){
        
        if (private$modelResultsExist(chem, private$model_name) && !private$overwrite_result){
          loginfo( paste0("OED results for ", chem$chem_info$casn, " already exist."));
          next;
        }
        
        loginfo( str_c("Calculating Oral Equivalent Doses for ", chem$chem_info$casn) );
        

        #obtain activity data -> should be a table with variables: aeid, <activity_variable>, intended_target_family
        if (private$activity_variable == "ac50"){
          activitydata <- private$getAc50DataFromChemical( chem );
        } else {
          activitydata <- NULL;
        }
        
        #OUTPUT RESULT FORMAT, set it to default values that are returned if model is not successful
        modelresult <- tibble(aeid = chem$assay_info$aeid, oed = NA, oed_units = NA, background = NA, activity = private$activity_variable);
        
        #CHECKPOINT: ensure we have activity data
        if ( is.null(activitydata) ){
          logwarn( str_c("No function to import activity variable '", private$activity_variable, "' exists. Results are set to NA.") );
          private$setModelResultsToChemical ( chem, modelresult ); #set NA result
          next;
        }
        
        #CHECKPOINT: ensure the activity data table is not empty, otherwise set oed results to NA and exit
        if (length(activitydata$aeid) == 0){
          logwarn( str_c("No above-threshold activity data is available for ", chem$chem_info$casn, ". Results are NA.") );
          private$setModelResultsToChemical ( chem, modelresult ); #set NA result
          next;
        }
        
        #CHECKPOINT: obtain css that is not NA, otherwise set result to NA and exit 
        Css <- private$getCssDataFromChemical( chem, private$cssmodel );
        if (is.na(Css)){
          logwarn( str_c("Css value for ", chem$chem_info$casn, " is NA. OED results are NA.") );
          private$setModelResultsToChemical ( chem, modelresult ); #set NA result
          next;
        }
        
        #compute OED: oed = activity(uM) * 1(mg/kgBW/day) / Css(uM);
        oed <- activitydata[[private$activity_variable]] / Css;
        
        #aeid vector to assist assigning oed values into correct modelresult table spots
        aei <- activitydata$aeid;
        
        #background measurement flagging
        bg <- activitydata$intended_target_family;
        
        #save oed to computed data table
        for (rowindex in 1:length(oed)){
          modelresult$oed[ (modelresult$aeid == aei[rowindex]) ] <- oed[rowindex];
          modelresult$oed_units[ (modelresult$aeid == aei[rowindex]) ] <- "mg/kg/day";
          modelresult$background[ (modelresult$aeid == aei[rowindex]) ] <- ifelse( bg[rowindex] == "background measurement", "Y", "N" );
        }
        private$setModelResultsToChemical ( chem, modelresult );
        loginfo(str_c( chem$chem_info$casn, ": OED values calculated and saved to ", private$model_name, ". ",length(oed), " values saved."));

      }
      
      private$chemicals <- NULL;
      invisible(self);
    }
    
  )
);
