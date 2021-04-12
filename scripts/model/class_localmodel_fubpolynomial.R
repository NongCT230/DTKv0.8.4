

# Model Fubs based on polynomial fit for experimental data -------------------------


#v0.7

# model based on:
# MATLAB fit() with poly33 model, computing 3rd order surface fit of pKa and logKow values to human Fub values from 660 chemicals in httk v1.7

# model assumptions used:
#  Fub is only dependent on pKa and logKow
#  organism = human
#  where pKa data is not available, assume neutral pKa of tissue pH
#  steady state -> d/dt = 0 for all values
#  regression based on 660 chemicals with Fub data available from httk v1.7. Assuming the sample is representative of all chemicals to be looked at

#model:
# polynomial and factors, 95% CI in comments:
# Fub <- p00 + p10*pka + p01*logkow + p20*pka^2 + p11*pka*logkow + p02*logkow^2 + p30*pka^3 + p21*pka^2*logkow + p12*pka*logkow^2 + p03*logkow^3;
  # p00 = 724.647500685581e-3; #(0.5483, 0.901)
  # p10 = -136.731999695322e-3; #(-0.2251, -0.04832)
  # p01 = -156.795688064424e-3; #(-0.2251, -0.08852)
  # p20 = 25.6262119918870e-3; #(0.01119, 0.04007)
  # p11 = 11.6344257993534e-3; #(-0.008332, 0.0316)
  # p02 = -486.802594042371e-6; #(-0.01242, 0.01144)
  # p30 = -1.14381437411737e-3; #(-0.00189, -0.0003977)
  # p21 = -923.357055364979e-6; #(-0.002438, 0.0005911)
  # p12 = -513.530793023847e-6; #(-0.001966, 0.0009392)
  # p03 = 1.10259381028834e-3; #(0.0005038, 0.001701)
# Fit quality:
#   sse: 13.9782
#   rsquare: 0.5587
#   dfe: 257
#   adjrsquare: 0.5433
#   rmse: 0.2332


Class.LocalModel.FubPolynomial <- R6Class("Class.LocalModel.FubPolynomial",
  
  inherit = Class.LocalModel,
  
  #private variables and functions
  private = list(
    polynomial_coefficients = NULL,
    pH_plasma = NULL,
    import.dreamtk = NULL,
    import.dreamtk.fub = NULL
  ),
  
  #public variables and functions
  public = list(
    
    #constructor
    initialize = function() {
      super$initialize();
      private$model_name <- "fub.polynomial";
      
      private$polynomial_coefficients <- list();
      private$polynomial_coefficients$p00 = 724.647500685581e-3; #(0.5483, 0.901)
      private$polynomial_coefficients$p10 = -136.731999695322e-3; #(-0.2251, -0.04832)
      private$polynomial_coefficients$p01 = -156.795688064424e-3; #(-0.2251, -0.08852)
      private$polynomial_coefficients$p20 = 25.6262119918870e-3; #(0.01119, 0.04007)
      private$polynomial_coefficients$p11 = 11.6344257993534e-3; #(-0.008332, 0.0316)
      private$polynomial_coefficients$p02 = -486.802594042371e-6; #(-0.01242, 0.01144)
      private$polynomial_coefficients$p30 = -1.14381437411737e-3; #(-0.00189, -0.0003977)
      private$polynomial_coefficients$p21 = -923.357055364979e-6; #(-0.002438, 0.0005911)
      private$polynomial_coefficients$p12 = -513.530793023847e-6; #(-0.001966, 0.0009392)
      private$polynomial_coefficients$p03 = 1.10259381028834e-3; #(0.0005038, 0.001701)
      
    },
    
    #finalizer
    finalize = function() {
    },
    
    #
    setParameters = function ( overwrite_data = FALSE, overwrite_result = FALSE ){
      if(is.logical(overwrite_data)){
        private$overwrite_data <- overwrite_data;
      }
      if(is.logical(overwrite_result)){
        private$overwrite_result <- overwrite_result;
      }
    },
    
    #
    setImporters = function( generalimport = NULL, modelimport = NULL){
      if (!is.null(generalimport) && !is.null(generalimport$classname) && generalimport$classname == "Class.Dataimporter.dreamtk"){
        private$import.dreamtk <- generalimport;
      }
      if (!is.null(modelimport) && !is.null(modelimport$classname) && modelimport$classname == "Class.DataImporter.dreamtk.FubPolynomial"){
        private$import.dreamtk.fub <- modelimport;
      }
    },
    
    #
    importPhysiologyData = function (){
      
      if(is.null(private$import.dreamtk)){
        e <- "LocalModel.FubPolynomial$importPhysiologyData(): Missing dreamtk data importer object";
        logerror(e);
        stop(e);
      };
      
      #what we need to import
      parameters <- c("pH_plasma");
      
      #import physiology data
      private$import.dreamtk$loadPhysiologyData( params = parameters);
      physiology_data <- private$import.dreamtk$getPhysiologyData();
      
      #obtain constant physiology parameters for the calculation, these are same for all chemicals

      #blood plasma pH, pH_plasma, untiless fraction
      private$pH_plasma <- filter( physiology_data, parameter == "pH_plasma", species == "human")$value;
      
      #verify imported data integrity, if values donot exist in the table they will be numeric(0)
      if( length(private$pH_plasma)==0 ) {
        private$pH_plasma <- NA;
        w <- "LocalModel.FubPolynomial$importPhysiologyData(): Unable to obtain fub physiology values from database.";
        logwarn(w);
        warning(w);
      }
      
      invisible(self);
      
    },
    
    #
    importModelData = function () {

      if(is.null(private$chemicals) || is.null(private$caslist)){
        e <- "LocalModel.FubPolynomial$importModelData(): Can not import model data without a chemical. $setChemicals() first";
        logerror(e);
        stop(e);
      };
      
      if(is.null(private$import.dreamtk.fub)){
        e <- "LocalModel.FubPolynomial$importModelData(): Missing fub model data importer object";
        logerror(e);
        stop(e);
      };

      #do the model data import for the chemical
      private$import.dreamtk.fub$loadModelData( chemicals = private$caslist );
      modeldata <- private$import.dreamtk.fub$getModelData();
      
      #attach model data to chemicals
      for ( chem in private$chemicals ){
        
        if (private$modelDataExists(chem) && !private$overwrite_data){
          loginfo( paste0("Fub polynomial model data for ", chem$chem_info$casn, " already exists."));
          next;
        }
        
        modeldata_chem <- filter(modeldata, casn == chem$chem_info$casn) %>% slice(1); #grab first row of data
        
        self$addModelDataToChemical ( chem, modeldata_chem  );
        
      }

      invisible(self);
      
    },
    
    #runs the  model
    compute = function () {
      
      if(is.null(private$chemicals)){
        e <- "LocalModel.FubPolynomial$compute(): Can not compute model without a chemical. $setChemicals() first";
        logerror(e);
        stop(e);
      };
      
      if(is.null(private$pH_plasma)){
        e <- "LocalModel.FubPolynomial$compute(): Can not compute model without physiology data. $importPhysiologyData() first";
        logerror(e);
        stop(e);
      };

      #compute model
      for ( chem in private$chemicals){
        
        if (private$modelResultsExist(chem, private$model_name) && !private$overwrite_result){
          loginfo( paste0("Fub results for ", chem$chem_info$casn, " already exist."));
          next;
        }
        
        if(!private$modelDataExists(chem)){
          e <- str_c("LocalModel.FubPolynomial$compute(): Model data for chemical ", chem$chem_info$casn, 
                     " does not exist. $importModelData() first");
          logerror(e);
          stop(e);
        };
        
        
        loginfo( str_c("Calculating Fub value for ", chem$chem_info$casn) );
        
        #OUTPUT RESULT FORMAT, set it to default values that are returned if model is not successful
        modelresult <- list(parameter = "fub", value = NA, units = NA, organism = NA);
  
        #model data should not be duplicated normally, but if it happens obtain only the first row of model data
        model.datarow <- slice( private$getModelDataFromChemical( chem ), 1);
  
        #CHECKPOINT: ensure the model data table is not empty, otherwise set result to NA and exit
        if (length(model.datarow$log_kow) == 0){
          
          logwarn( str_c("No fub model data is available for ", chem$chem_info$casn, ". Results are NA.") );
          private$setModelResultsToChemical ( chem, modelresult ); #set NA result
          next;
          
        }
  
        #obtain model parameters
        
        #logkow, octanol:Water partition coefficient at 25*C, unitless fraction log
        logkow <- model.datarow$log_kow;
  
        #first and most neutral ionization constant for the chemical
        pka <- model.datarow$pka;
  
        #CHECKPOINT: parameters are not NA (except pKa = NA assumed neutral compound), otherwise set result to NA and exit
        model.vars <- list(logkow = logkow);
        if( any( is.na(model.vars) ) ){
          
          logwarn( str_c("Some essential fub model data variables for ", chem$chem_info$casn, " are NA. Results are NA.") );
          logwarn( str_c("NA variables are: ", str_c( names( which(is.na(model.vars))), collapse=", " ) ) );
          private$setModelResultsToChemical ( chem, modelresult ); #set NA result
          next;
          
        }
  
        #assume chemical is neutral if pka is NA, set pKa to tissue pH
        if (is.na(pka)) pka <- private$pH_plasma;
  
        #unbound chemical fraction in plasma, Fub, unitless fraction
        Fub <- private$polynomial_coefficients$p00 + 
          private$polynomial_coefficients$p10*pka + 
          private$polynomial_coefficients$p01*logkow + 
          private$polynomial_coefficients$p20*pka^2 + 
          private$polynomial_coefficients$p11*pka*logkow + 
          private$polynomial_coefficients$p02*logkow^2 +
          private$polynomial_coefficients$p30*pka^3 + 
          private$polynomial_coefficients$p21*pka^2*logkow + 
          private$polynomial_coefficients$p12*pka*logkow^2 + 
          private$polynomial_coefficients$p03*logkow^3;
  
        if (Fub > 1) { Fub <- 1 }
        else if (Fub < 0) { Fub <- 0 };
  
        #SUCCESS!
        loginfo("Fub:Success!")
        modelresult <- list(parameter = "fub", value = Fub, units = "fraction", organism = "human");
        private$setModelResultsToChemical ( chem, modelresult );
      
      }
      
      private$chemicals <- NULL;
      invisible(self);
      
    },
    
    #
    getModelPhysiologyData = function(){
      
      if( is.null(private$pH_plasma) ){
        return (NULL);
      } else {
        return(    list(pH_plasma = private$pH_plasma) );
      }
    }
    
  )
)
