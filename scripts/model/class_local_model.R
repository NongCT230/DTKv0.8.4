
# Local Model class -------------------------------------------------------

#v0.7

#local model -> models and saves parameters for a single chemical
#this class is a prototype, which contains the basic methods expected to be shared by all local models
#actual local model classes are expected to extend this class

Class.LocalModel <- R6Class("Class.LocalModel",
                            
  cloneable = TRUE,
  lock_class = FALSE,

  #self variables and functions
  private = list(
    
    chemicals = NULL,
    caslist = NULL,
    model_name = NULL,
    parameters = NULL,
    attached.models = NULL,
    overwrite_data = NULL,
    overwrite_result = NULL,
    
    modelDataExists = function ( chemical = NULL ){
      if(is.null(chemical)){
        e <- "LocalModel$modelDataExists(): No chemical given. Set some chemicals and add model data first.";
        logerror(e);
        stop(e);
      }
      
      if ( is.null(chemical$model_data[[private$model_name]]) ){
        return (FALSE);
      }
      
      return( TRUE );
      
    },
    
    getModelDataFromChemical = function ( chemical = NULL ){
      if( private$modelDataExists(chemical) ){
        return( chemical$model_data[[private$model_name]] );
      } else {
        return (NULL);
      }
      
    },
    
    modelResultsExist = function ( chemical = NULL, modelname = NULL ){
      if(is.null(chemical)){
        e <- "LocalModel$modelResultsExist(): Can not get model results. $setChemicals() first, then request results.";
        logerror(e);
        stop(e);
      }
      
      if (!is.character(modelname) || length(modelname) != 1 ){
        e <- "LocalModel$modelResultsExist(): Can not get model results. Model name supplied is not a character string.";
        logerror(e);
        stop(e);
      }
      
      if(is.null(chemical$model_results)){
        e <- "LocalModel$modelResultsExist(): Can not get model results. Chemical has not been initialized with model_results list";
        logerror(e);
        stop(e);
      }
      
      if (is.null(chemical$model_results[[modelname]])){
        return (FALSE);
      }
      
      return( TRUE );
      
    },
    
    getModelResultsFromChemical = function ( chemical = NULL, modelname = NULL ){
      if ( private$modelResultsExist( chemical, modelname ) ){
        return( chemical$model_results[[modelname]] );
      } else {
        return (NULL);
      }
    },
    
    setModelResultsToChemical = function ( chemical = NULL, modelresults = NULL ) {
      
      if(is.null(chemical)){
        e <- "LocalModel$setModelResultsToChemical(): Can set model results. $setChemicals() first.";
        logerror(e);
        stop(e);
      }
      
      chemical$model_results[[private$model_name]] <- modelresults;
      
    },
    
    #
    isValidChemical = function( chemical ){
      if (is.null(chemical) || is.atomic( chemical ) || 
          is.null(chemical$classname) || chemical$classname != "Class.Chemical" ||
          is.null(chemical$chem_info) || is.null(chemical$chem_info$casn) || is.na(chemical$chem_info$casn)){
        return (FALSE);
      }else{
        return (TRUE);
      }
    },
    
    #
    getAttachedModel = function( modelname ){
      if (!is.character(modelname) && length(modelname) != 1){
        e <- "LocalModel$getAttachedModel(): Model name should be single a character string.";
        logerror(e);
        stop(e);
      };
      
      if(is.null(private$attached.models) || !is.list(private$attached.models) || length(private$attached.models) == 0){
        logwarn("LocalModel$getAttachedModel(): No models attached. pass model list(...) to $attachModels(). returning NULL.");
        return (NULL);
      };
      
      for (model in private$attached.models){
        
        if (!is.null(model) && !is.atomic(model) && is.function(model$getModelName) &&
            model$getModelName() == modelname ){
          
          loginfo("LocalModel$getAttachedModel(): Found attached model to return.");
          return (model);
        }
      };
      
      logwarn("LocalModel$getAttachedModel(): No requested models attached. returning NULL.");
      return (NULL);
      
    }
    
  ),
  
  #public variables and functions
  public = list(

    #constructor
    initialize = function() {
      private$model_name <- "model_name";
      private$overwrite_data <- FALSE;
      private$overwrite_result <- FALSE;
    },
    
    #finalizer
    finalize = function() {},
    
    #sets chemicals to model and save data to
    setChemicals = function( chemicals = NULL ){
      
      if( !is.list(chemicals) ){
        e <- "LocalModel$setChemicals(): List of chemicals expected. If setting a single chemical, pass list(chemical)";
        logerror(e);
        stop(e);
      }
      
      validchems <- list();
      caslist <- list();
      
      for ( chem in chemicals ){
        if ( private$isValidChemical(chem) ){
          validchems <- append(validchems, chem);
          caslist <- append(caslist, chem$chem_info$casn);
        } else {
          w <- "LocalModel$setChemical(): Removed invalid Chemical object from model chemlist.";
          logwarn(w);
        }
      }

      private$chemicals <- validchems;
      private$caslist <- caslist;
      

    },
    
    #attaches other models that could be used in computation
    attachModels = function( modelslist = NULL ){
      if (!is.null(modelslist)){
        private$attached.models <- modelslist;
      }
    },
    
    #sets model parameters
    setParameters = function ( params = NULL ){
      if (is.null(params)){
        logwarn("LocalModel$setParameters(): Model parameters do not exist.");
      }
      
      private$parameters <- params;
      
    },
     
    addModelDataToChemical = function ( chemical = NULL, modeldata = NULL ) {
      if(is.null(chemical) || is.null(chemical$model_data) || is.atomic(chemical$model_data)){
        e <- "LocalModel$addModelDataToChemical(): Can not add model data. No valid Chemical given.";
        logerror(e);
        stop(e);
      }
      
      chemical$model_data[[private$model_name]] <- modeldata;
      
    },
    
    #runs the  model
    #
    getModelName = function(){
      return (private$model_name);
    }
    
  )
);
