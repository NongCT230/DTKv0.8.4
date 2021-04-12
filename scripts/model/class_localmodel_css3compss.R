

# Model Css based on 3 components at steady state -------------------------


#v0.7

# model based on:
#   Pearce et al. (2017) httk: R Package for  High-Throughput Toxicokinetics. Journal of Statistical Software, 79(4), 1-25

# model assumptions:
#  3 compartments participate in establishing Css: gut, liver, rest of body
#  average physiological and tissue values
#  steady state -> d/dt = 0 for all concentration gradients
#  constant chemical exposure rate at 1mg/kgBW/day

# overestimates Css compared to Wetmore 2015 median data.


Class.LocalModel.Css3compss <- R6Class("Class.LocalModel.Css3compss",
  
  inherit = Class.LocalModel,
  
  #private variables and functions
  private = list(
    avgBW = NULL,
    kdose = NULL,
    Qgfr = NULL,
    Qliver = NULL,
    Qgut = NULL,
    avgLW = NULL,
    Rhohep = NULL,
    sig_threshold = NULL,
    import.dreamtk = NULL,
    import.dreamtk.3compss = NULL,
    
    na_fub = NULL,
    na_rb2p = NULL,
    na_clint = NULL,

    #acquires, verifies and  returns Fub data from a given model
    getFubDataFromChemicalModel = function ( chemical = NULL, modelname = NULL ){
      
      fubdata <- private$getModelResultsFromChemical( chemical, modelname );
      
      if(is.null(fubdata)){
        logwarn( paste0("LocalModel.Css3compss$getFubFromChemical(): Fub model result for chemical ", chemical$chem_info$casn, 
                        " does not exist. Could not find result for model: ", modelname, 
                        ". Run the model on the Chemical or enter results manually."));
        return (NA);
      };
      
      if ( is.atomic(fubdata) || 
           is.null( fubdata$parameter ) || 
           is.null( fubdata$value ) ||
           is.na( fubdata$parameter ) ||
           fubdata$parameter != "fub"){
        logwarn("LocalModel.Css3compss$getFubFromChemical(): The specified fub model result exists, 
                 but does not have expected structure. Expected to have at least $parameter='fub', and $value");
        return (NA);
      }
      
      return ( fubdata$value );
    },
    
    computeMissingFub = function(){
      #compute missing values for Fub
      fubPolymodel <- private$getAttachedModel ( "fub.polynomial" );
      if (!is.null(fubPolymodel)){
        #model attached, compute
        fubPolymodel$setChemicals(private$na_fub);
        fubPolymodel$importModelData();
        fubPolymodel$compute();
      }
    }
    
    #compute missing Rblood2plasma -todo
    
    #compute missing Clint -todo
    
  ),
  
  #public variables and functions
  public = list(
    
    #constructor
    initialize = function() {
      super$initialize();
      private$model_name <- "css.3compss";
      private$sig_threshold <- 0.05;
    },
    
    #finalizer
    finalize = function() {
    },

    #
    setImporters = function( generalimport = NULL, modelimport = NULL){

      if (!is.null(generalimport) && !is.null(generalimport$classname) && generalimport$classname == "Class.Dataimporter.dreamtk"){
        private$import.dreamtk <- generalimport;
      }
      if (!is.null(modelimport) && !is.null(modelimport$classname) && modelimport$classname == "Class.Dataimporter.dreamtk.3compss"){
        private$import.dreamtk.3compss <- modelimport;
      }
    },
    
    #
    setParameters = function ( sig.threshold = 0.05, overwrite_data = FALSE, overwrite_result = FALSE ){
     if(!is.null(sig.threshold) && !is.na(sig.threshold) && is.numeric(sig.threshold) ){
        private$sig_threshold <- sig.threshold;
     }
      if(is.logical(overwrite_data)){
        private$overwrite_data <- overwrite_data;
      }
      if(is.logical(overwrite_result)){
        private$overwrite_result <- overwrite_result;
      }
    },
    
    #
    importPhysiologyData = function (){
      
      if(is.null(private$import.dreamtk)){
        logerror("LocalModel.Css3compss$importPhysiologyData(): Missing dreamtk data importer object. $setImporters() first");
        stop();
      };
      
      #what we need to import
      parameters <- c("Average BW", "GFR", "Flow_liver", "Flow_gut");
      
      #import required physiology data
      private$import.dreamtk$loadPhysiologyData( params = parameters);
      physiology_data <- private$import.dreamtk$getPhysiologyData();
      
      #obtain constant physiology parameters for the calculation, these are same for all chemicals
      #average body weight, avgBW, kg
      private$avgBW <- filter( physiology_data, parameter == "Average BW", species == "human")$value;
      
      #dosing constant, kdose, mg/kgBW/day
      private$kdose <- 1;
      
      #glomerular filtration rate, Qgfr, ml/min/kg^3/4, (need L/day)
      private$Qgfr <- filter(physiology_data, parameter == "GFR", species == "human")$value;
      #convert to L/day
      private$Qgfr <- private$Qgfr * private$avgBW^(3/4) * (1/10^3) * 60 * 24;
      
      #blood flow through liver, Qliver, ml/min/kg^3/4, (need L/day)
      private$Qliver <- filter(physiology_data, parameter == "Flow_liver", species == "human")$value;
      #convert to L/day
      private$Qliver <- private$Qliver * private$avgBW^(3/4) * (1/10^3) * 60 * 24;
      
      #blood flow through gut, Qgut, ml/min/kg^3/4, (need L/day)
      private$Qgut <- filter(physiology_data, parameter == "Flow_gut", species == "human")$value;
      #convert to L/day
      private$Qgut <- private$Qgut * private$avgBW^(3/4) * (1/10^3) * 60 * 24;
      
      #average liver weight, avgLW, grams (Johnson et al., 2005)
      private$avgLW <- 1596;
      
      #hepatocyte density, Rhohep, cells/gram liver (Ito & Houston, 2004)
      private$Rhohep <- 1.1*10^8;
      
      #verify imported data integrity, if values donot exist in the table they will be numeric(0)
      if( length(private$avgBW)==0 || length(private$Qgfr)==0 || length(private$Qliver)==0 || length(private$Qgut)==0 ) {
        private$avgBW <- NA;
        private$Qgfr <- NA;
        private$Qliver <- NA;
        private$Qgut <- NA;
        w <- "LocalModel.Css3compss$importPhysiologyData(): Unable to obtain 3css physiology values from database.";
        logwarn(w);
        warning(e);
      }
      
      invisible(self);
      
    },
    
    #
    importModelData = function () {

      if(is.null(private$chemicals) || is.null(private$caslist)){
        e <- "LocalModel.Css3compss$importModelData(): Can not import model data without a list of Chemicals. $setChemicals() first";
        logerror(e);
        stop(e);
      };
      
      if(is.null(private$import.dreamtk.3compss)){
        e <- "LocalModel.Css3compss$importModelData(): Missing 3css model data importer object. $setImporters() first";
        logerror(e);
        stop(e);
      };

      #do the model data import for the chemical
      private$import.dreamtk.3compss$loadModelData( chemicals = private$caslist );
      modeldata <- private$import.dreamtk.3compss$getModelData();

      #attach model data to chemicals and determine which chemicals have missing data, add these to correspondig lists for further processing
      private$na_fub <- list();
      private$na_rb2p <- list();
      private$na_clint <- list();
      
      for ( chem in private$chemicals ){
        
        if (private$modelDataExists(chem) && !private$overwrite_data){
          loginfo( paste0("Css 3css model data for ", chem$chem_info$casn, " already exists."));
          next;
        }
        
        modeldata_chem <- filter(modeldata, casn == chem$chem_info$casn) %>% slice(1); #grab first row of data

        self$addModelDataToChemical ( chem, modeldata_chem  );
        
        if(length(modeldata_chem$casn) > 0){
          if (is.na(modeldata_chem$human_funbound_plasma)){
            private$na_fub <- append(private$na_fub, chem);
          }
          if (is.na(modeldata_chem$human_rblood2plasma)){
            private$na_rb2p <- append(private$na_rb2p, chem);
          }
          if (is.na(modeldata_chem$human_clint)){
            private$na_clint <- append(private$na_clint, chem);
          }
        }

      }
      
      invisible(self);
      
    },
    
    #runs the  model
    compute = function () {
      
      if(is.null(private$chemicals)){
        e <- "LocalModel.Css3compss$compute(): Can not compute model without a Chemical list. $setChemicals() first";
        logerror(e);
        stop(e);
      };
      
      if(is.null(private$avgBW)){
        e <- "LocalModel.Css3compss$compute(): Can not compute model without physiology data. $importPhysiologyData() first";
        logerror(e);
        stop(e);
      };
      
      #compute missing values
      private$computeMissingFub();
      
      #compute model
      for ( chem in private$chemicals){

        if (private$modelResultsExist(chem, private$model_name) && !private$overwrite_result){
          loginfo( paste0("Css 3css results for ", chem$chem_info$casn, " already exist."));
          next;
        }
        
        if(!private$modelDataExists(chem)){
          e <- str_c("LocalModel.Css3compss$compute(): Model data for chemical ", chem$chem_info$casn, 
                     " does not exist. $importModelData() first")
          logerror(e);
          stop(e);
        };
        
        
        loginfo( str_c("Calculating 3 compartment steady state Css values for ", chem$chem_info$casn) );
        
        #OUTPUT RESULT FORMAT, set it to default values that are returned if model is not successful
        modelresult <- list(parameter = "css", value = NA, units = NA, organism = NA, assumptions = "missing model parameters");
        
        #model data should not be duplicated normally, but if it happens obtain only the first row of model data
        model.datarow <- slice( private$getModelDataFromChemical( chem ), 1);
  
        #CHECKPOINT: ensure the model data table is not empty, otherwise set result to NA and exit
        if (length(model.datarow$casn) == 0){
          modelresult$assumptions <- paste0(modelresult$assumptions,
                                            ": MW, Fub, Clint, Rblood2plasma");
          logwarn( str_c("No 3compss model data is available for ", chem$chem_info$casn, ". Results are NA.") );
          private$setModelResultsToChemical ( chem, modelresult ); #set NA result
          next;
          
        }
  
        #obtain model parameters
        
        #create assumption set, starting assumption is experimental values
        assumptions = list(fub = "fub:experimental", rblood2plasma = "rblood2plasma:experimental", clint = "clint:experimental");
        
        #chemical molecular weight, mw, g/mol
        mw <- model.datarow$mw;
  
        #unbound chemical fraction, Fub, unitless
        Fub <- model.datarow$human_funbound_plasma;
        
        if (is.na(Fub)){
          # check if Polynomial model has already been computed
          if ( private$modelResultsExist( chem, "fub.polynomial" ) ){
            assumptions$fub = "fub:3rd order polynomial model estimation.";
          }
          # whether model data exists or not, we retrieve computed value or NA
          Fub <- private$getFubDataFromChemicalModel( chem, "fub.polynomial" );
        }
        assumptions$fub = paste0(assumptions$fub, " Fub: ", Fub);
         
        if ( !is.na(Fub) && Fub == 0){ #Fub values below limit of detection are coded as 0 in the table
          Fub <- 0.005;
          assumptions$fub = paste0(assumptions$fub, ", fub:value below limit of detection, set to 0.005");
        }
  
        #ratio of blood to plasma concentration of chemical, Rblood2plasma, unitless
        Rblood2plasma <- model.datarow$human_rblood2plasma;
        if (is.na(Rblood2plasma)){ #presently not computing NA values, assuming they are unity, which is fair for chemicals we are currently examining
          Rblood2plasma <- 1;
          assumptions$rblood2plasma = "rblood2plasma:assumed unity value";
        }
        assumptions$rblood2plasma = paste0(assumptions$rblood2plasma, " Rb2p: ", Rblood2plasma);
  
        #intrinsic invitro hepatic clearance, Clint, uL/min/10^6 cells
        Clint <- model.datarow$human_clint;
        #verify is statistically significant at threshold, make NA if over threshold
        Clint_p <- model.datarow$human_clint_p;
        if (!is.na(Clint_p) && Clint_p >= private$sig_threshold) Clint <- NA;
        if ( is.na(Clint) ){
          # run another model to compute?
        }
        assumptions$clint = paste0(assumptions$clint, " Clint: ", Clint);
        
        modelresult$assumptions <- assumptions;
  
        #whole liver hepatic clearance, Clmetabolism, L/day
        #units L/day = uL/min/cells * 1/10^6 * cells/gram_liver * gram_liver * 60min/1hr * 24hr/1day * 1L/10^6uL
        Clmetabolism <- Clint * (1/10^6) * private$Rhohep * private$avgLW * 60 * 24 * (1/10^6);
  
        #CHECKPOINT: parameters are not NA, otherwise set result to NA and exit
        model.vars <- list(avgBW = private$avgBW, kdose = private$kdose, Qgfr = private$Qgfr, Qliver = private$Qliver, 
                           Qgut = private$Qgut, avgLW = private$avgLW, Rhohep = private$Rhohep, mw = mw, Fub = Fub, 
                           Rblood2plasma = Rblood2plasma, Clmetabolism = Clmetabolism);
        if( any( is.na(model.vars) ) ){
          missing <- str_c( names( which(is.na(model.vars))), collapse=", " );
          modelresult$assumptions$missing <- paste0("missing: ", missing);
          logwarn( str_c("Some essential 3compss model data variables for ", chem$chem_info$casn, " are NA. Results are NA.") );
          logwarn( str_c("NA variables are: ", missing ) );
          private$setModelResultsToChemical ( chem, modelresult ); #set NA result
          next;
          
        }
  
        #calculate components of the Css formula (Pearce et al., 2017)
        Qratio_nom <- (private$Qliver + private$Qgut)*Fub*Clmetabolism; # ( L/day )^2
        Qratio_denom <- (private$Qliver + private$Qgut) + Fub*Clmetabolism/Rblood2plasma; # L/day
        Qratio <- Qratio_nom/Qratio_denom; # L/day
  
        #calculate Css, mg/kgBW/L
        #units mg/kgBW/L = mg/kgBW/day / (L/day + L/day)
        Css <- private$kdose / (Fub*private$Qgfr + Qratio);
  
        #convert Css to uM concentration units
        #uM = mg/kgBW/L * 1g/1000mg * 1/g/mol * kgBW * 10^6umol/1mol
        Css <- Css * (10^-3) * (1/mw) * private$avgBW * (10^6);
  
        #CHECKPOINT: Css may be infinite if denominator ends up a zero, set result to NA and exit
        if(!is.finite(Css)){
          modelresult$assumptions$cssval <- paste0("css: infinite value obtained");
          logwarn( str_c("An infinite value of Css was obtained for ", chem$chem_info$casn, ". Results are NA.") );
          private$setModelResultsToChemical ( chem, modelresult );
          next;
        }
  
        #SUCCESS!
        loginfo("Css: Success!")
        modelresult <- list(parameter = "css", value = Css, units = "uM", organism = "human", assumptions = assumptions);
        private$setModelResultsToChemical ( chem, modelresult );
      
      }
      
      private$chemicals <- NULL;
      invisible(self);
      
    },
    
    #
    getModelPhysiologyData = function(){
      
      if( is.null(private$avgBW) ){
        return (NULL);
      } else {
        return(    list(avgBW = private$avgBW,
                        kdose = private$kdose,
                        Qgfr = private$Qgfr,
                        Qliver = private$Qliver,
                        Qgut = private$Qgut,
                        avgLW = private$avgLW,
                        Rhohep = private$Rhohep,
                        sig_cutoff = private$sig_threshold) );
      }
    }
    
  )
)
