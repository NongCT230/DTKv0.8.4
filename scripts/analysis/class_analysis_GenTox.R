#class contains GT data and GT analysis calculations.
# - buildGTtable and GTplot

#Please note that this part was named GenTox before being renamed to GeneTox. Some old name could still be in the code

Class.Analysis.GT <- R6Class("Class.Analysis.GT",
  
	private = list(
		chem_data = NULL,
		chem_data_gt = NULL,
		assay_desc = NULL,
		logBase = 10,
		hoverSignificantDigits = 3
	),
	
	public = list(
		
		initialize = function(){
			self$basicData <- Class.Analysis.Data$new();
		},
    
		#finalizer
		finalize = function() {},
		
		getChemData = function(){return (private$chem_data);},
		getChemDataGT = function(){return (private$chem_data_gt);},
		
		basicData = NULL,
		
		#functions
		
		#build the table that will be used for all the grahps
		buildGTTable = function(chemList, label_by = "name", pValue = 1, model = 1, logValues = FALSE, preCalcBMR = TRUE, customBMRVal = -1, showAssayDesc = FALSE, progress = NULL, graph = 1) {
		  loginfo("Building GeneTox table");
		  
		  private$chem_data <- tribble(~chemical, ~endpoint_grp, ~bmr, ~bmdl, ~bmd, ~bmdu, ~trend);
		  
		  load("data/GenTox21.RData");
		  
		  private$chem_data_gt <- gen_tox_data;
		  
		  if (preCalcBMR)
		  {
  		  DataTable <- bmd_data;
  		  
  		  if (model == 1) #best model by loglikelihood
  		  {
  		    bestModelOrder <- c("exponential", "linear", "hill curve", "geometric mean");
  		    
  		    DataTable <- DataTable %>% 
  		      group_by(casn, compound, assay, endpoint, s9_positive) %>%
  		      filter(log_likelihood == max(log_likelihood, na.rm = TRUE)) %>%
  		      filter(aic == min(aic, na.rm = TRUE)) %>%
  		      filter(model_type == bestModelOrder[min(match(model_type, bestModelOrder))]) %>%
  		      ungroup() %>% as.data.frame();
  		  }
  		  else if (model == 2) #best model by AIC
  		  {
  		    bestModelOrder <- c("exponential", "linear", "hill curve", "geometric mean");
  		    
  		    DataTable <- DataTable %>% 
  		      group_by(casn, compound, assay, endpoint, s9_positive) %>%
  		      filter(aic == min(aic, na.rm = TRUE)) %>%
  		      filter(log_likelihood == max(log_likelihood, na.rm = TRUE)) %>%
  		      filter(model_type == bestModelOrder[min(match(model_type, bestModelOrder))]) %>%
  		      ungroup() %>% as.data.frame();
  		  }
  		  else if (model == 3) #linear
  		  {
  		    DataTable <- DataTable %>% filter(., .$model_type == "linear");
  		  }
  		  else if (model == 4) #exponential
  		  {
  		    DataTable <- DataTable %>% filter(., .$model_type == "exponential");
  		  }
  		  else if (model == 5) #hill curve
  		  {
  		    DataTable <- DataTable %>% filter(., .$model_type == "hill curve");
  		  }
  		  else #if (model == 6) #geometric mean
  		  {
  		    DataTable <- DataTable %>% filter(., .$model_type == "geometric mean");
  		  }
  		  
  		  DataTable <- DataTable %>% filter(., !is.na(.$bmd));
  		  
  		  if (nrow(DataTable) < 1)
  		  {
  		    return();
  		  }
  		  
  		  DataTable <- join(as.data.frame(apply(DataTable, 2, tolower), stringsAsFactors = FALSE), 
  		                    as.data.frame(apply(trend_data, 2, tolower), stringsAsFactors = FALSE));
  		  
  		  #Types seems to not be what they should at this point. This fixes it.
  		  DataTable$bmr <- as.numeric(DataTable$bmr);
  		  DataTable$bmdl <- as.numeric(DataTable$bmdl);
  		  DataTable$bmd <- as.numeric(DataTable$bmd);
  		  DataTable$bmdu <- as.numeric(DataTable$bmdu);
  		  DataTable$trend <- as.logical(DataTable$trend);
  		  
  		  if (label_by == "name")
  		  {
  		    DataTable <- mutate(DataTable, chemical = tools::toTitleCase(tolower(as.character(compound)))) %>%
  		      filter(tolower(chemical) %in% tolower(chemList));
  		    private$chem_data_gt <- mutate(private$chem_data_gt, chemical = tools::toTitleCase(tolower(as.character(compound)))) %>%
  		      filter(tolower(chemical) %in% tolower(chemList));
  		  }
  		  else
  		  {
  		    DataTable <- mutate(DataTable, chemical = as.character(casn)) %>%
  		      filter(chemical %in% chemList);
  		    private$chem_data_gt <- mutate(private$chem_data_gt, chemical = as.character(casn)) %>%
  		      filter(chemical %in% chemList);
  		  }
		  }
		  else #Need To run throught Proast
		  {
		    
		    lqt <- attr(Class.Analysis.GT, "lastQuery")
		    
		    if (is.null(lqt) || lqt[1, "bmr"] != customBMRVal)
		    {
		      ToAnalyseTable <- gen_tox_data %>%
		        mutate(., endpoint_grp = tools::toTitleCase(paste(.$assay, .$endpoint, ifelse((.$s9_positive), "+S9", "-S9"), sep = " ")));
		      
		      if (label_by == "name")
		      {
		        ToAnalyseTable <- mutate(ToAnalyseTable, chemical = tools::toTitleCase(tolower(as.character(compound)))) %>%
		          filter(tolower(chemical) %in% tolower(chemList));
		      }
		      else
		      {
		        ToAnalyseTable <- mutate(ToAnalyseTable, chemical = as.character(casn)) %>%
		          filter(chemical %in% chemList);
		      }
		      
		      if(nrow(ToAnalyseTable) < 1)
		      {
		        return();
		      }
		      
		      ToAnalyseTable <- mutate(ToAnalyseTable, analyseGroup = paste(chemical, "|", endpoint_grp, sep = "")) %>%
		        select(analyseGroup, concentration, response);
		      
		      ToAnalyseTable <- ToAnalyseTable[order(ToAnalyseTable[,1]),]
		      
		      source("./scripts/analysis/Modified_Proast65.5.r")
		      
		      subsetStartRow <- 1
		      prevChemical <- ToAnalyseTable[1, "analyseGroup"]
		      currChemical <- ""
		      currDose <- 0
		      currResp <- 0
		      prevRow <- 0
		      currSubset <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("analysisGroup", "Dose", "Resp"))
		      
		      results <- setNames(data.frame(matrix(ncol = 18, nrow = 0)), c("analysisGroup", "Hill_loglik", "Hill_AIC", "Hill_BMDL", "Hill_BMD", "Hill_BMDU", "Exp_loglik", "Exp_AIC", "Exp_BMDL", "Exp_BMD", "Exp_BMDU", "Geo_BMDL", "Geo_BMD", "Geo_BMDU"))
		      
		      for (currRow in 1:nrow(ToAnalyseTable)) {
		        currChemical <- ToAnalyseTable[currRow, "analyseGroup"]
		        currDose <- ToAnalyseTable[currRow, "concentration"]
		        currMN <- ToAnalyseTable[currRow, "response"]
		        
		        currGeoBMD <- "NA"
		        currGeoBMDL <- "NA"
		        currGeoBMDU <- "NA"
		        
		        
		        # encounter a new chemical
		        # use the subset to run PROAST
		        # export PROAST results
		        if (currChemical != prevChemical) {
		          
		          withCallingHandlers(f.proast(currSubset, BMR = customBMRVal),
		                              warning=function(w) {
		                                if (grepl("produc", w$message))
		                                  #remove message "Warning in log(x): NaNs produced"
		                                  invokeRestart("muffleWarning")
		                              } )
		          
		          # get data from proast results
		          currHillBMD <- proastTemp$HILL$CED
		          currHillBMDL<- proastTemp$HILL$conf.int[,1]
		          currHillBMDU<- proastTemp$HILL$conf.int[,2]
		          currHill_loglik <- proastTemp$HILL$loglik
		          currHill_aic <- proastTemp$HILL$aic
		          
		          currExpBMD <- proastTemp$EXP$CED
		          currExpBMDL<- proastTemp$EXP$conf.int[,1]
		          currExpBMDU<- proastTemp$EXP$conf.int[,2]
		          currExp_loglik <- proastTemp$EXP$loglik
		          currExp_aic <- proastTemp$EXP$aic
		          
		          currGeoBMD <- sqrt(currExpBMD * currHillBMD)
		          currGeoBMDL <- sqrt(currExpBMDL * currHillBMDL)
		          currGeoBMDU <- sqrt(currExpBMDU * currHillBMDU)
		          
		          # append data to results data frame
		          results <- rbind(results, data.frame(analysisGroup = prevChemical, Hill_loglik = currHill_loglik, Hill_AIC = currHill_aic, Hill_BMDL = currHillBMDL, Hill_BMD = currHillBMD, Hill_BMDU = currHillBMDU, Exp_loglik = currExp_loglik, Exp_AIC = currExp_aic, Exp_BMDL = currExpBMDL, Exp_BMD = currExpBMD, Exp_BMDU = currExpBMDU, Geo_BMDL = currGeoBMDL, Geo_BMD = currGeoBMD, Geo_BMDU = currGeoBMDU))
		          
		          # reset currSubset
		          rm(currSubset)
		          currSubset <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("analysisGroup", "Dose", "MN"))
		        }
		        
		        # if the end of the file
		        else if (currRow == nrow(ToAnalyseTable)) {
		          
		          # add last line of the file to currSubset
		          currSubset <- rbind(currSubset, data.frame(analysisGroup = currChemical, Dose = currDose, MN = currMN))
		          
		          withCallingHandlers(f.proast(currSubset, BMR = customBMRVal),
		                              warning=function(w) {
		                                if (grepl("produc", w$message))
		                                  #remove message "Warning in log(x): NaNs produced"
		                                  invokeRestart("muffleWarning")
		                              } )
		          
		          # get data from proast results
		          currHillBMD <- proastTemp$HILL$CED
		          currHillBMDL<- proastTemp$HILL$conf.int[,1]
		          currHillBMDU<- proastTemp$HILL$conf.int[,2]
		          currHill_loglik <- proastTemp$HILL$loglik
		          currHill_aic <- proastTemp$HILL$aic
		          
		          currExpBMD <- proastTemp$EXP$CED
		          currExpBMDL<- proastTemp$EXP$conf.int[,1]
		          currExpBMDU<- proastTemp$EXP$conf.int[,2]
		          currExp_loglik <- proastTemp$EXP$loglik
		          currExp_aic <- proastTemp$EXP$aic
		          
		          currGeoBMD <- sqrt(currExpBMD * currHillBMD)
		          currGeoBMDL <- sqrt(currExpBMDL * currHillBMDL)
		          currGeoBMDU <- sqrt(currExpBMDU * currHillBMDU)
		          
		          # append data to results data frame
		          results <- rbind(results, data.frame(analysisGroup = prevChemical, Hill_loglik = currHill_loglik, Hill_AIC = currHill_aic, Hill_BMDL = currHillBMDL, Hill_BMD = currHillBMD, Hill_BMDU = currHillBMDU, Exp_loglik = currExp_loglik, Exp_AIC = currExp_aic, Exp_BMDL = currExpBMDL, Exp_BMD = currExpBMD, Exp_BMDU = currExpBMDU, Geo_BMDL = currGeoBMDL, Geo_BMD = currGeoBMD, Geo_BMDU = currGeoBMDU))
		          
		          # reset currSubset
		          rm(currSubset)
		          currSubset <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("analysisGroup", "Dose", "MN"))
		          
		        }
		        #currChemical equals prevChemical, add data to the currSubset
		        #currSubsetRow <- c(currChemical, currDose, currMN)
		        currSubset <- rbind(currSubset, data.frame(analysisGroup = currChemical, Dose = currDose, MN = currMN))
		        
		        # set previous to current
		        prevChemical <- currChemical
		        prevRow <- currRow
		        
		        progress$inc(3/(10*nrow(ToAnalyseTable)), detail = "Generating GeneTox Table")
		      }
		      
		      DataTable <- tribble(~chemical, ~endpoint_grp, ~model_type, ~aic, ~log_likelihood, ~bmr, ~bmdl, ~bmd, ~bmdu, ~analysisGroup);
		      
		      for (i in 1:nrow(results))
		      {
		        DataTable <- rbind(DataTable, data.frame(chemical = str_extract(results[i, "analysisGroup"], regex(".+?(?=\\|)")),
		                                                 endpoint_grp = str_extract(results[i, "analysisGroup"], regex("(?<=\\|).+")),
		                                                 bmr = customBMRVal,
		                                                 model_type = "hill curve",
		                                                 aic = results[i, "Hill_AIC"],
		                                                 log_likelihood = results[i, "Hill_loglik"],
		                                                 bmdl = results[i, "Hill_BMDL"],
		                                                 bmd = results[i, "Hill_BMD"],
		                                                 bmdu = results[i, "Hill_BMDU"],
		                                                 analysisGroup = results[i, "analysisGroup"],
		                                                 stringsAsFactors = FALSE
		        )) %>%
		          rbind(., data.frame(chemical = str_extract(results[i, "analysisGroup"], regex(".+?(?=\\|)")),
		                              endpoint_grp = str_extract(results[i, "analysisGroup"], regex("(?<=\\|).+")),
		                              bmr = customBMRVal,
		                              model_type = "exponential",
		                              aic = results[i, "Exp_AIC"],
		                              log_likelihood = results[i, "Exp_loglik"],
		                              bmdl = results[i, "Exp_BMDL"],
		                              bmd = results[i, "Exp_BMD"],
		                              bmdu = results[i, "Exp_BMDU"],
		                              analysisGroup = results[i, "analysisGroup"],
		                              stringsAsFactors = FALSE
		          )) %>%
		          rbind(., data.frame(chemical = str_extract(results[i, "analysisGroup"], regex(".+?(?=\\|)")),
		                              endpoint_grp = str_extract(results[i, "analysisGroup"], regex("(?<=\\|).+")),
		                              bmr = customBMRVal,
		                              model_type = "geometric mean",
		                              aic = NA_real_,
		                              log_likelihood = NA_real_,
		                              bmdl = results[i, "Geo_BMDL"],
		                              bmd = results[i, "Geo_BMD"],
		                              bmdu = results[i, "Geo_BMDU"],
		                              analysisGroup = results[i, "analysisGroup"],
		                              stringsAsFactors = FALSE
		          ));
		        
		        progress$inc(1/(10*nrow(results)), detail = "Generating GeneTox Table")
		      }
		    }
		    else
		    {
		      DataTable <- lqt;
		    }
		    attr(Class.Analysis.GT, "lastQuery") <<- DataTable;
		    
		    if (model == 1) #best model by loglikelihood
		    {
		      bestModelOrder <- c("exponential", "linear", "hill curve", "geometric mean");
		      
		      DataTable <- DataTable %>% 
		        group_by(analysisGroup) %>%
		        filter(log_likelihood == max(log_likelihood, na.rm = TRUE)) %>%
		        filter(aic == min(aic, na.rm = TRUE)) %>%
		        filter(model_type == bestModelOrder[min(match(model_type, bestModelOrder))]) %>%
		        ungroup() %>% as.data.frame();
		    }
		    else if (model == 2) #best model by AIC
		    {
		      bestModelOrder <- c("exponential", "linear", "hill curve", "geometric mean");
		      
		      DataTable <- DataTable %>% 
		        group_by(analysisGroup) %>%
		        filter(aic == min(aic, na.rm = TRUE)) %>%
		        filter(log_likelihood == max(log_likelihood, na.rm = TRUE)) %>%
		        filter(model_type == bestModelOrder[min(match(model_type, bestModelOrder))]) %>%
		        ungroup() %>% as.data.frame();
		    }
		    else if (model == 3) #linear
		    {
		      DataTable <- DataTable %>% filter(., .$model_type == "linear");
		    }
		    else if (model == 4) #exponential
		    {
		      DataTable <- DataTable %>% filter(., .$model_type == "exponential");
		    }
		    else if (model == 5) #hill curve
		    {
		      DataTable <- DataTable %>% filter(., .$model_type == "hill curve");
		    }
		    else #if (model == 6) #geometric mean
		    {
		      DataTable <- DataTable %>% filter(., .$model_type == "geometric mean");
		    }
		    
		    DataTable <- DataTable %>% filter(., !is.na(.$bmd));
		    
		    if (nrow(DataTable) < 1)
		    {
		      return();
		    }
		    
		    if (label_by == "name")
		    {
		      private$chem_data_gt <- mutate(private$chem_data_gt, chemical = tools::toTitleCase(tolower(as.character(compound)))) %>%
		        filter(tolower(chemical) %in% tolower(chemList));
		    }
		    else
		    {
		      private$chem_data_gt <- mutate(private$chem_data_gt, chemical = as.character(casn)) %>%
		        filter(chemical %in% chemList);
		    }
		    
		    t_table <- trend_data %>%
		      mutate(., endpoint_grp = tools::toTitleCase(paste(.$assay, .$endpoint, ifelse((.$s9_positive), "+S9", "-S9"), sep = " ")));
		    
		    if (label_by == "name")
		    {
		      t_table <- mutate(t_table, chemical = tools::toTitleCase(tolower(as.character(compound)))) %>%
		        filter(tolower(chemical) %in% tolower(chemList));
		    }
		    else
		    {
		      t_table <- mutate(t_table, chemical = as.character(casn)) %>%
		        filter(chemical %in% chemList);
		    }
		    
		    DataTable <- join(as.data.frame(apply(DataTable, 2, tolower), stringsAsFactors = FALSE), 
		                      as.data.frame(apply(t_table, 2, tolower), stringsAsFactors = FALSE));
		    
		    if (label_by == "name")
		    {
		      DataTable <- mutate(DataTable, chemical = tools::toTitleCase(tolower(as.character(compound)))) %>%
		        filter(tolower(chemical) %in% tolower(chemList));
		    }
		    else
		    {
		      DataTable <- mutate(DataTable, chemical = as.character(casn)) %>%
		        filter(chemical %in% chemList);
		    }
		    
		    #Types seems to not be what they should at this point. This fixes it.
		    DataTable$bmr <- as.numeric(DataTable$bmr);
		    DataTable$bmdl <- as.numeric(DataTable$bmdl);
		    DataTable$bmd <- as.numeric(DataTable$bmd);
		    DataTable$bmdu <- as.numeric(DataTable$bmdu);
		    DataTable$trend <- as.logical(DataTable$trend);
		  }
		  
		  DataTable <- DataTable %>%
		    mutate(., endpoint_grp = tools::toTitleCase(paste(.$assay, .$endpoint, ifelse((.$s9_positive), "+S9", "-S9"), sep = " "))) %>%
		    mutate(., model_type = tools::toTitleCase(tolower(as.character(model_type))));
		  
		  private$chem_data_gt <- private$chem_data_gt %>%
		    mutate(., endpoint_grp = tools::toTitleCase(paste(.$assay, .$endpoint, ifelse((.$s9_positive), "+S9", "-S9"), sep = " ")));
		  
		  if (graph != 7)
		  {
		      private$chem_data_gt <- private$chem_data_gt %>%
		          select(chemical, endpoint_grp, concentration, concentration_unit, response, response_type);
		  }
		  else
		  {
		      private$chem_data_gt <- private$chem_data_gt %>%
		          select(chemical, endpoint_grp, concentration, concentration_unit, response, response_type, assay);
		  }
		  
		  if (logValues)
		  {
		    if (graph == 7)
		    {
		        DataTable[, "bmdl"] <- DataTable[, "bmdl"] + 1
		        DataTable[, "bmd"] <- DataTable[, "bmd"] + 1
		        DataTable[, "bmdu"] <- DataTable[, "bmdu"] + 1
		    }
		    DataTable[,"bmdl"] <- log(DataTable[,"bmdl"], private$logBase);
		    DataTable[,"bmd"] <- log(DataTable[,"bmd"], private$logBase);
		    DataTable[,"bmdu"] <- log(DataTable[,"bmdu"], private$logBase);
		  }
		  
		  correctEndpointGrp <- function(x) {
		    x[,"endpoint_grp"] <- x[, "endpoint_grp"] %>%
		      gsub("Tgr ", "TGR ", ., fixed = TRUE) %>%
		      gsub(" Fe1 ", " FE1 ", ., fixed = TRUE) %>%
		      gsub(" +s9", " +S9", ., fixed = TRUE) %>%
		      gsub("flow ", "Flow ", ., fixed = TRUE) %>%
		      gsub(" P53 ", " p53 ", ., fixed = TRUE) %>%
		      gsub(" %cytotoxicity ", " %Cytotoxicity ", ., fixed = TRUE) %>%
		      gsub(" % ", " %", ., fixed = TRUE) %>%
		      gsub(" Ii ", " II ", ., fixed = TRUE) %>%
		      gsub(" Tamix ", " TAMIX ", ., fixed = TRUE) %>%
		      gsub(" Ta98 ", " TA98 ", ., fixed = TRUE) %>%
		      gsub("Dna ", "DNA ", ., fixed = TRUE) %>%
		      gsub(" g-H2ax ", " g-H2aX ", ., fixed = TRUE) %>%
		      gsub("Cometchip ", "CometChip ", ., fixed = TRUE);
		    return(x);
		  }
		  
		  DataTable <- correctEndpointGrp(DataTable);
		  private$chem_data_gt <- correctEndpointGrp(private$chem_data_gt);
		  
		  if (graph != 7)
		  {
		      private$chem_data <- select(DataTable, chemical, endpoint_grp, bmr, bmdl, bmd, bmdu, trend, model_type);
		  }
		  else
		  {
		      private$chem_data <- select(DataTable, chemical, endpoint_grp, bmr, bmdl, bmd, bmdu, trend, model_type, assay);
		  }
		  
		  if (showAssayDesc){private$assay_desc <- assay_description}
		},
		
		#create and show the plots for single BMR preset values
		plotGeneToxSingleBMR = function(input, output, bmrValue = 5, progress, logVal = FALSE, showPlotData = FALSE, forceScale = FALSE, forceOrder = FALSE) {
		  loginfo("Creating GeneTox plots");
		  private$chem_data <- private$chem_data %>% filter(.$bmr == bmrValue);
		  chemicals <- private$chem_data['chemical'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  endpoints <- private$chem_data['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  
		  if (length(input$select_chemical_gt) > length(chemicals))
		  {
		    logwarn(paste("Missing chemical info for: ", paste(input$select_chemical_gt[which(!input$select_chemical_gt %in% chemicals)], collapse = ", "), sep = ""));
		    showNotification("Some chemicals might be missing because no data was associated with them.", type = "warning", duration = 5);
		  }
		  
		  chemDataTable <- list();
		  
		  xmin <- min(min(private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != -Inf,'bmdl'], na.rm = TRUE), min(private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != -Inf,'bmd'], na.rm = TRUE), min(private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != -Inf,'bmdu'], na.rm = TRUE), na.rm = TRUE)
		  xmax <- max(max(private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != Inf,'bmdl'], na.rm = TRUE), max(private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != Inf,'bmd'], na.rm = TRUE), max(private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != Inf,'bmdu'], na.rm = TRUE), na.rm = TRUE)
		  
		  if (showPlotData)
		  {
  		  insertUI(
  		    selector = "#ui_gt_plotContainerCenter",
  		    where = "afterBegin",
  		    ui = box(status = "primary", title = "Plotting Data Table", collapsible = TRUE, width = 12,
            column(width = 8, offset = 2,
                   wellPanel(
                     h4("Plotting Data Table"),
                     DT::dataTableOutput(outputId = "table_GTPlottingData1")
                   )
            )
  		    )
  		  )
  		  output$table_GTPlottingData1 <- DT::renderDataTable(server = FALSE, {
  		    datatable(eval(parse(text=paste("private$chem_data %>% select(chemical, endpoint_grp, bmr, bmdl, bmd, bmdu, trend, model_type) %>% ",
  		                                    "dplyr::rename(\"Chemical\" = chemical, \"Endpoint Group\" = endpoint_grp, \"BMR\" = bmr, ", ifelse(logVal, "\"log(BMDL)\"", "\"BMDL\""),
  		                                    " = bmdl, ", ifelse(logVal, "\"log(BMD)\"", "\"BMD\""), " = bmd, ", ifelse(logVal, "\"log(BMDU)\"", "\"BMDU\""), " = bmdu, \"Trend\" = trend, ",
  		                                    "\"Model Type\" = model_type)", sep=""))),
  		              selection="none", filter="bottom", extensions = "Buttons",
  		              options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
  		                           pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
  		                           scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
  		      formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
  		    
  		  })
		  }
		  
		  if (!is.null(private$assay_desc))
		  {
		    insertUI(
		      selector = "#ui_gt_plotContainerCenter",
		      where = "beforeEnd",
		      ui = box(status = "primary", title = "Assay Description", collapsible = TRUE, width = 12,
		               column(width = 8, offset = 2,
		                      wellPanel(
		                        h4("Assay Description"),
		                        DT::dataTableOutput(outputId = "table_GTAssayDescription1")
		                      )
		               )
		      )
		    )
		    output$table_GTAssayDescription1 <- DT::renderDataTable(server = FALSE, {
		      datatable(private$assay_desc,
		                selection="none", filter="bottom", extensions = "Buttons",
		                options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                             pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                             scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		        formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		      
		    })
		  }
		  
		  insertUI(
		    selector = "#ui_gt_plotContainerCenter",
		    where = "beforeEnd",
		    ui = fluidRow(id = "ui_gt_plots",
                column(12,
                  box(status = "primary", title = NULL, collapsible = FALSE, solidHeader = FALSE, width = 12,
                     dropdownButton(
                       tags$h3("About the plot"),
                       tags$h5("This is a plot showing the confidence interval for each chemicals and endpoints."),
                       tags$h5("Red lines represent a dataset without trend, according to TCPL's MC4 test."),
                       tags$h5("Blue stars represent a value that is tending towards positive or negative infinity."),
                       circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
                       tooltip = NULL)
                   )
                )
		    )
		  );
		  
		  i <- 1;
		  if (input$gt_comparaison == 1)
		  {
  		  for (c in chemicals)
  		  {
  		    chemDataTable[[c]] <- private$chem_data %>% filter(chemical == c);
  		    for (endpt in endpoints)
  		    {
  		      if (nrow(chemDataTable[[c]] %>% filter(endpoint_grp == endpt)) ==0)
  		      {
  		        chemDataTable[[c]] = rbind(chemDataTable[[c]], as.data.frame(list(chemical=c, endpoint_grp = endpt, bmr = NA_real_, bmdl = NA_real_, bmd = NA_real_, bmdu = NA_real_, trend = FALSE, model_type=NA)));
  		      }
  		    }
  		    if ((i-1) %% 2 == 0)
  		    {
  		      insertUI(
  		        selector = "#ui_gt_plotContainerLeft",
  		        where = "beforeEnd",
  		        ui = fluidRow(id = "ui_gt_plots",
  		                      column(12,
  		                             box(status = "primary", title = "", collapsible = TRUE, width = 12,
  		                                 plotlyOutput(outputId = c, height=length(endpoints)*51+95)
  		                             )
  		                      )
  		        )
  		      );
  		    }
  		    else
  		    {
  		      insertUI(
  		        selector = "#ui_gt_plotContainerRight",
  		        where = "beforeEnd",
  		        ui = fluidRow(id = "ui_gt_plots",
  		                      column(12,
  		                             box(status = "primary", title = "", collapsible = TRUE, width = 12,
  		                                 plotlyOutput(outputId = c, height=length(endpoints)*51+95)
  		                             )
  		                      )
  		        )
  		      );
  		    }
  		    local({
  		      myChem <- c;
    		    output[[myChem]] <- renderPlotly(
    		      {
    		        if (logVal)
    		        {
    		          NInfVal <- min(ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] != -Inf, ])['bmdl']) == 0, Inf, min(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] != -Inf, ])['bmdl'], na.rm = TRUE)), 
    		                         ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmd"]) & (chemDataTable[[myChem]])["bmd"] != -Inf, ])['bmd']) == 0, Inf, min(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmd"]) & (chemDataTable[[myChem]])["bmd"] != -Inf, ])['bmd'], na.rm = TRUE)), 
    		                         ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] != -Inf, ])['bmdu']) == 0, Inf, min(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] != -Inf, ])['bmdu'], na.rm = TRUE)),
    		                         ifelse(forceScale, xmin, Inf), na.rm = TRUE)-1;
    		          InfVal <- max(ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
    		                        ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmd"]) & (chemDataTable[[myChem]])["bmd"] != Inf, ])['bmd']) == 0, -Inf, max(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmd"]) & (chemDataTable[[myChem]])["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
    		                        ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)),
    		                        ifelse(forceScale, xmax, -Inf), na.rm = TRUE)+1;
    		        }
    		        else
    		        {
    		          NInfVal <- -1;
    		          InfVal <- max(ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
    		                        ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmd"]) & (chemDataTable[[myChem]])["bmd"] != Inf, ])['bmd']) == 0, -Inf, max(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmd"]) & (chemDataTable[[myChem]])["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
    		                        ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)), 
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)*1.1;
    		        }
    		        
    		        if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] == Inf,]) > 0)
    		          {chemDataTable[[myChem]][!is.na(chemDataTable[[myChem]]["bmdl"]) & chemDataTable[[myChem]]["bmdl"] == Inf,"bmdl"] <- InfVal;}
    		        if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmd"]) & (chemDataTable[[myChem]])["bmd"] == Inf,]) > 0)
    		          {chemDataTable[[myChem]][!is.na(chemDataTable[[myChem]]["bmd"]) & chemDataTable[[myChem]]["bmd"] == Inf,"bmd"] <- InfVal;}
    		        if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] == Inf,]) > 0)
    		          {chemDataTable[[myChem]][!is.na(chemDataTable[[myChem]]["bmdu"]) & chemDataTable[[myChem]]["bmdu"] == Inf,"bmdu"] <- InfVal;}
    		        if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] == -Inf,]) > 0)
    		          {chemDataTable[[myChem]][!is.na(chemDataTable[[myChem]]["bmdl"]) & chemDataTable[[myChem]]["bmdl"] == -Inf,"bmdl"] <- NInfVal;}
    		        if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmd"]) & (chemDataTable[[myChem]])["bmd"] == -Inf,]) > 0)
    		          {chemDataTable[[myChem]][!is.na(chemDataTable[[myChem]]["bmd"]) & chemDataTable[[myChem]]["bmd"] == -Inf,"bmd"] <- NInfVal;}
    		        if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] == -Inf,]) > 0)
    		          {chemDataTable[[myChem]][!is.na(chemDataTable[[myChem]]["bmdu"]) & chemDataTable[[myChem]]["bmdu"] == -Inf,"bmdu"] <- NInfVal;}
    		        
    		        if (forceOrder)
    		        {
    		          p <- ggplot(chemDataTable[[myChem]], aes(x=as.numeric(bmd), y=endpoint_grp, color = trend));
    		        }
    		        else
    		        {
    		          p <- ggplot(chemDataTable[[myChem]], aes(x=as.numeric(bmd), y=reorder(reorder(endpoint_grp, -as.numeric(bmd)), -as.numeric(is.na(bmd))), color = trend));
    		        }
    		          
  		          p <- p + ggtitle(myChem) +
  		            labs(y="", x=ifelse(logVal, "log(BMD) Values [uM]", "BMD Values [uM]")) +
    		          geom_point(size = 2, mapping = aes(text = paste(
    		            paste("Endpoint: ", endpoint_grp, sep=""),
    		            paste(ifelse(logVal, "log(BMDL): ", "BMDL: "), ifelse(bmdl == NInfVal, "-Inf", ifelse(bmdl == InfVal, "Inf", formatC(signif(bmdl), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
    		            paste(ifelse(logVal, "log(BMD): ", "BMD: "), ifelse(bmd == NInfVal, "-Inf", ifelse(bmd == InfVal, "Inf", formatC(signif(bmd), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
    		            paste(ifelse(logVal, "log(BMDU): ", "BMDU: "), ifelse(bmdu == NInfVal, "-Inf", ifelse(bmdu == InfVal, "Inf", formatC(signif(bmdu), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
    		            paste("Trend: ", ifelse(trend, "Yes", "No"), sep=""),
    		            paste("Model: ", model_type, sep=""),
    		            sep="\n"))) +
    		          scale_color_manual(values = setNames(c("black", "red"), c(TRUE, FALSE))) +
    		          theme(plot.title = element_text(size=22, hjust=0.5), legend.position = "none");
  		          
  		          if (forceScale)
  		          {
  		            p <- p +
  		              xlim(xmin-1.5, ifelse(logVal, xmax+1.5, xmax*1.21))
  		            #Allow space for the Inf values and blue dots/stars
  		          }
    		        if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] == NInfVal, ]) > 0)
    		        {
    		          p <- p +
    		            geom_point(data = (chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdl"]) & (chemDataTable[[myChem]])["bmdl"] == NInfVal, ], 
    		                       mapping = aes(x = (NInfVal-0.5), y = endpoint_grp, shape = 8, text = "This BMDL value is tending towards negative infinity."), color = "blue") +
    		            scale_shape_identity();
    		        }
  		          if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] == InfVal, ]) > 0)
  		          {
  		            p <- p +
  		              geom_point(data = (chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])["bmdu"]) & (chemDataTable[[myChem]])["bmdu"] == InfVal, ], 
  		                         mapping = aes(x = (ifelse(logVal, InfVal+0.5, InfVal*1.1)), y = endpoint_grp, shape = 8, text = "This BMDU value is tending towards infinity."), color = "blue") +
  		              scale_shape_identity();
  		          }
    		        if (nrow((chemDataTable[[myChem]])[(!is.na((chemDataTable[[myChem]])["bmdl"]) & !is.na((chemDataTable[[myChem]])["bmdu"])), ]) > 0)
    		        {
    		          p <- p + geom_errorbarh(data = (chemDataTable[[myChem]])[(!is.na((chemDataTable[[myChem]])["bmdl"]) & !is.na((chemDataTable[[myChem]])["bmdu"])), ], 
    		                         aes(y=endpoint_grp, xmax=as.numeric(bmdu), xmin=as.numeric(bmdl), x=as.numeric(bmd), text = "", 
    		                             height=length(endpoints)/9), size=0.2, na.rm = TRUE);
    		        } 
    		        
    		        ggplotly(p, tooltip = "text");
    		      }
    		    )
  		    })
  		    i <- i + 1;
  		    progress$inc(1/(2*length(chemicals)), detail = "Generating GeneTox Plots");
  		  }
		  }
		  else
		  {
		    for (e in endpoints)
		    {
		      chemDataTable[[e]] <- private$chem_data %>% filter(endpoint_grp == e);
		      for (ch in chemicals)
		      {
		        if (nrow(chemDataTable[[e]] %>% filter(chemical == ch)) ==0)
		        {
		          chemDataTable[[e]] = rbind(chemDataTable[[e]], as.data.frame(list(chemical=ch, endpoint_grp = e, bmr = NA_real_, bmdl = NA_real_, bmd = NA_real_, bmdu = NA_real_, trend = FALSE, model_type = NA)));
		        }
		      }
		      if ((i-1) %% 2 == 0)
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerLeft",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotlyOutput(outputId = e, height=length(chemicals)*59+71)
		                               )
		                        )
		          )
		        );
		      }
		      else
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerRight",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotlyOutput(outputId = e, height=length(chemicals)*59+71)
		                               )
		                        )
		          )
		        );
		      }
		      local({
		        myEndpt <- e;
		        output[[myEndpt]] <- renderPlotly(
		          {
		            if (logVal)
		            {
		              NInfVal <- min(ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] != -Inf, ])['bmdl']) == 0, Inf, min(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] != -Inf, ])['bmdl'], na.rm = TRUE)), 
		                             ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmd"]) & (chemDataTable[[myEndpt]])["bmd"] != -Inf, ])['bmd']) == 0, Inf, min(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmd"]) & (chemDataTable[[myEndpt]])["bmd"] != -Inf, ])['bmd'], na.rm = TRUE)), 
		                             ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] != -Inf, ])['bmdu']) == 0, Inf, min(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] != -Inf, ])['bmdu'], na.rm = TRUE)),
		                             ifelse(forceScale, xmin, Inf), na.rm = TRUE)-1;
		              InfVal <- max(ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
		                            ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmd"]) & (chemDataTable[[myEndpt]])["bmd"] != Inf, ])['bmd']) == 0, -Inf, max(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmd"]) & (chemDataTable[[myEndpt]])["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
                                ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)),
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)+1;
		            }
		            else
		            {
		              NInfVal <- -1;
		              InfVal <- max(ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
		                            ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmd"]) & (chemDataTable[[myEndpt]])["bmd"] != Inf, ])['bmd']) == 0, -Inf, max(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmd"]) & (chemDataTable[[myEndpt]])["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
		                            ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)), 
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)*1.1;
		            }
		            
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] == Inf,]) > 0)
		            {chemDataTable[[myEndpt]][!is.na(chemDataTable[[myEndpt]]["bmdl"]) & chemDataTable[[myEndpt]]["bmdl"] == Inf,"bmdl"] <- InfVal;}
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmd"]) & (chemDataTable[[myEndpt]])["bmd"] == Inf,]) > 0)
		            {chemDataTable[[myEndpt]][!is.na(chemDataTable[[myEndpt]]["bmd"]) & chemDataTable[[myEndpt]]["bmd"] == Inf,"bmd"] <- InfVal;}
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] == Inf,]) > 0)
		            {chemDataTable[[myEndpt]][!is.na(chemDataTable[[myEndpt]]["bmdu"]) & chemDataTable[[myEndpt]]["bmdu"] == Inf,"bmdu"] <- InfVal;}
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] == -Inf,]) > 0)
		            {chemDataTable[[myEndpt]][!is.na(chemDataTable[[myEndpt]]["bmdl"]) & chemDataTable[[myEndpt]]["bmdl"] == -Inf,"bmdl"] <- NInfVal;}
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmd"]) & (chemDataTable[[myEndpt]])["bmd"] == -Inf,]) > 0)
		            {chemDataTable[[myEndpt]][!is.na(chemDataTable[[myEndpt]]["bmd"]) & chemDataTable[[myEndpt]]["bmd"] == -Inf,"bmd"] <- NInfVal;}
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] == -Inf,]) > 0)
		            {chemDataTable[[myEndpt]][!is.na(chemDataTable[[myEndpt]]["bmdu"]) & chemDataTable[[myEndpt]]["bmdu"] == -Inf,"bmdu"] <- NInfVal;}
		            
		            if (forceOrder)
		            {
		              p <- ggplot(chemDataTable[[myEndpt]], aes(x=as.numeric(bmd), y=chemical, color = trend));
		            }
		            else
		            {
		              p <- ggplot(chemDataTable[[myEndpt]], aes(x=as.numeric(bmd), y=reorder(reorder(chemical, -as.numeric(bmd)), -as.numeric(is.na(bmd))), color = trend));
		            }
		            
		            p <- p + ggtitle(myEndpt) +
		              labs(y="", x=ifelse(logVal, "log(BMD) Values [uM]", "BMD Values [uM]")) +
		              geom_point(size = 2, mapping = aes(text = paste(
		                paste("Chemical: ", chemical, sep=""),
		                paste(ifelse(logVal, "log(BMDL): ", "BMDL: "), ifelse(bmdl == NInfVal, "-Inf", ifelse(bmdl == InfVal, "Inf", formatC(signif(bmdl), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
		                paste(ifelse(logVal, "log(BMD): ", "BMD: "), ifelse(bmd == NInfVal, "-Inf", ifelse(bmd == InfVal, "Inf", formatC(signif(bmd), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
		                paste(ifelse(logVal, "log(BMDU): ", "BMDU: "), ifelse(bmdu == NInfVal, "-Inf", ifelse(bmdu == InfVal, "Inf", formatC(signif(bmdu), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
		                paste("Trend: ", ifelse(trend, "Yes", "No"), sep=""),
		                paste("Model: ", model_type, sep=""),
		                sep="\n"))) +
		              scale_color_manual(values = setNames(c("red", "black"), c(TRUE, FALSE))) +
		              theme(plot.title = element_text(size=22, hjust=0.5), legend.position = "none")
		            
		            if (forceScale)
		            {
		              p <- p +
		                xlim(xmin-1.5, ifelse(logVal, xmax+1.5, xmax*1.21))
		              #Allow space for the Inf values and blue dots/stars
		            }
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] == NInfVal, ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdl"]) & (chemDataTable[[myEndpt]])["bmdl"] == NInfVal, ], 
		                           mapping = aes(x = NInfVal-0.5, y = chemical, shape = 8, text = "This BMDL value is tending towards negative infinity."), color = "blue") +
		                scale_shape_identity();
		            }
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] == InfVal, ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])["bmdu"]) & (chemDataTable[[myEndpt]])["bmdu"] == InfVal, ], 
		                           mapping = aes(x = ifelse(logVal, InfVal+0.5, InfVal*1.1), y = chemical, shape = 8, text = "This BMDU value is tending towards infinity."), color = "blue") +
		                scale_shape_identity();
		            }
		            if (nrow((chemDataTable[[myEndpt]])[(!is.na((chemDataTable[[myEndpt]])["bmdl"]) & !is.na((chemDataTable[[myEndpt]])["bmdu"])), ]) > 0)
		            {
		              p <- p +geom_errorbarh(data = (chemDataTable[[myEndpt]])[(!is.na((chemDataTable[[myEndpt]])["bmdl"]) & !is.na((chemDataTable[[myEndpt]])["bmdu"])), ],
		                             aes(y=chemical, xmax=as.numeric(bmdu), xmin=as.numeric(bmdl), x=as.numeric(bmd), text = "", 
		                                 height=length(chemicals)/9), size=0.2, na.rm = TRUE);
		            }
		            
		              ggplotly(p, tooltip = "text");
		          }
		        )
		      })
		      i <- i + 1;
		      progress$inc(1/(2*length(endpoints)), detail = "Generating GeneTox Plots");
		      
		    }
		  }
		},
		
		#create and show the plots for all precalculated BMR values
		#Will need updating before re-implementation
		plotGeneToxMultiBMR = function(input, output, progress, logVal = FALSE, showPlotData = FALSE, forceScale = FALSE) {
		  loginfo("Creating GeneTox plots");
		  chemicals <- private$chem_data['chemical'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  endpoints <- private$chem_data['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  BMDres <- as.character(input$gt_bmd_result_chooser);
		  
		  if (length(input$select_chemical_gt) > length(chemicals))
		  {
		    logwarn(paste("Missing chemical info for: ", paste(input$select_chemical_gt[which(!input$select_chemical_gt %in% chemicals)], collapse = ", "), sep = ""));
		    showNotification("Some chemicals might be missing because no data was associated with them.", type = "warning", duration = 5);
		  }
		  
		  chemDataTable <- list();
		  
		  ymin <- min(private$chem_data[BMDres], na.rm = TRUE);
		  ymax <- max(private$chem_data[BMDres], na.rm = TRUE);
		  
		  if (showPlotData)
		  {
  		  insertUI(
  		    selector = "#ui_gt_plotContainerCenter",
  		    where = "afterBegin",
  		    ui = box(status = "primary", title = "Plotting Data Table", collapsible = TRUE, width = 12,
  		             column(width = 8, offset = 2,
  		                    wellPanel(
  		                      h4("Plotting Data Table"),
  		                      DT::dataTableOutput(outputId = "table_GTPlottingData2")
  		                    )
  		             )
  		    )
  		  )
  		  output$table_GTPlottingData2 <- DT::renderDataTable(server = FALSE, {
  		    
  		    datatable(eval(parse(text=paste("private$chem_data %>% select(chemical, endpoint_grp, bmr, ", BMDres, ", trend, model_type) %>%",
  		                                    "dplyr::rename(\"Chemical\" = chemical, \"Endpoint Group\" = endpoint_grp, \"BMR\" = bmr, \"", 
  		                                    ifelse(logVal, paste("log(", toupper(BMDres), ")", sep=""), toupper(BMDres)), 
  		                                    "\" = ", BMDres, ", \"Trend\" = trend, \"Model Type\" = model_type)", sep=""))),
  		              selection="none", filter="bottom", extensions = "Buttons",
  		              options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
  		                           pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
  		                           scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
  		      formatStyle(seq(5), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
  		    
  		  })
		  }
		  
		  if (!is.null(private$assay_desc))
		  {
		    insertUI(
		      selector = "#ui_gt_plotContainerCenter",
		      where = "beforeEnd",
		      ui = box(status = "primary", title = "Assay Description", collapsible = TRUE, width = 12,
		               column(width = 8, offset = 2,
		                      wellPanel(
		                        h4("Assay Description"),
		                        DT::dataTableOutput(outputId = "table_GTAssayDescription2")
		                      )
		               )
		      )
		    )
		    output$table_GTAssayDescription2 <- DT::renderDataTable(server = FALSE, {
		      datatable(private$assay_desc,
		                selection="none", filter="bottom", extensions = "Buttons",
		                options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                             pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                             scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		        formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		      
		    })
		  }
		  
		  i <- 1;
		  if (input$gt_comparaison == 1)
		  {
		    for (c in chemicals)
		    {
		      chemDataTable[[c]] <- private$chem_data %>% filter(chemical == c);
		      for (endpt in endpoints)
		      {
		        if (nrow(chemDataTable[[c]] %>% filter(endpoint_grp == endpt)) ==0)
		        {
		          chemDataTable[[c]] = rbind(chemDataTable[[c]], as.data.frame(list(chemical=c, endpoint_grp = endpt, bmr = NA_real_, bmdl = NA_real_, bmd = NA_real_, bmdu = NA_real_, trend = FALSE, model_type = NA)));
		        }
		      }
		      if ((i-1) %% 2 == 0)
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerLeft",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotOutput(outputId = c, height=405)
		                               )
		                        )
		          )
		        );
		      }
		      else
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerRight",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotOutput(outputId = c, height=405)
		                               )
		                        )
		          )
		        );
		      }
		      local({
		        myChem <- c;
		        output[[myChem]] <- renderPlot(
		          {
		            gg_color_list <- function(n) {
		              hues = seq(15, 375, length = n+1)
		              hcl(h = hues, l = 65, c = 100)[1:n]
		            }
		            
		            graphColors <- gg_color_list(length(endpoints));
		            
		            p <- ggplot(chemDataTable[[myChem]], aes(x=as.numeric(bmr), y=as.numeric(eval(parse(text=BMDres))), color = endpoint_grp, group = endpoint_grp, shape = as.character(trend))) +
		              ggtitle(myChem) +
		              labs(y=paste(ifelse(logVal, paste("log(", toupper(BMDres), ")", sep=""), toupper(BMDres)),"Values",sep=" "), x="BMR Values (%)", color="Endpoints") +
		              geom_line() +
		              geom_point(size = 3) +
		              scale_shape_manual(values = setNames(c(16, 4, 1), c("TRUE", "FALSE", "Circle"))) +
		              scale_size_identity() +
		              scale_color_manual(values = setNames(c(graphColors, "#000000"), c(sort(endpoints), "black")), breaks = sort(endpoints)) +
		              guides(shape = FALSE) +
		              xlim(0, 24) +
		              theme(plot.title = element_text(size=22, hjust=0.5));
		            
		            if (forceScale)
		            {
		              p <- p +
		                ylim(ymin, ymax)
		            }
		            if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == log(private$logZeroVal, private$logBase), ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == log(private$logZeroVal, private$logBase), ], 
		                           mapping = aes(x=as.numeric(bmr), y=as.numeric(eval(parse(text=BMDres))), size = 9, shape = "Circle", stroke = 1.5, color = "black"))+
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            if (logVal & nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == log(private$logInfVal, private$logBase), ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == log(private$logInfVal, private$logBase), ], 
		                           mapping = aes(x=as.numeric(bmr), y=as.numeric(eval(parse(text=BMDres))), size = 9, shape = "Circle", stroke = 1.5, color = "black"))+
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            if (!logVal & nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == private$logInfVal, ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == private$logInfVal, ], 
		                           mapping = aes(x=as.numeric(bmr), y=as.numeric(eval(parse(text=BMDres))), size = 9, shape = "Circle", stroke = 1.5, color = "black"))+
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            p
		          }
		        )
		      })
		      i <- i + 1;
		      progress$inc(1/(2*length(chemicals)), detail = "Generating GeneTox Plots");
		    }
		  }
		  else
		  {
		    for (e in endpoints)
		    {
		      chemDataTable[[e]] <- private$chem_data %>% filter(endpoint_grp == e);
		      for (ch in chemicals)
		      {
		        if (nrow(chemDataTable[[e]] %>% filter(chemical == ch)) ==0)
		        {
		          chemDataTable[[e]] = rbind(chemDataTable[[e]], as.data.frame(list(chemical=ch, endpoint_grp = e, bmr = NA_real_, bmdl = NA_real_, bmd = NA_real_, bmdu = NA_real_, trend = FALSE, model_type = NA)));
		        }
		      }
		      if ((i-1) %% 2 == 0)
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerLeft",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotOutput(outputId = e, height=405)
		                               )
		                        )
		          )
		        );
		      }
		      else
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerRight",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotOutput(outputId = e, height=405)
		                               )
		                        )
		          )
		        );
		      }
		      local({
		        myEndpt <- e;
		        output[[myEndpt]] <- renderPlot(
		          {
		            gg_color_list <- function(n) {
		              hues = seq(15, 375, length = n+1)
		              hcl(h = hues, l = 65, c = 100)[1:n]
		            }
		            
		            graphColors <- gg_color_list(length(chemicals))
		            
		            p <- ggplot(chemDataTable[[myEndpt]], aes(x=as.numeric(bmr), y=as.numeric(eval(parse(text=BMDres))), color = chemical, group = chemical, shape = as.character(trend))) +
		              ggtitle(myEndpt) +
		              labs(y=paste(ifelse(logVal, paste("log(", toupper(BMDres), ")", sep=""), toupper(BMDres)),"Values",sep=" "), x="BMR Values (%)", color = "Chemicals") +
		              geom_line() +
		              geom_point(size = 3) +
		              scale_shape_manual(values = setNames(c(16, 4, 1), c("TRUE", "FALSE", "Circle"))) +
		              scale_size_identity() +
		              scale_color_manual(values = setNames(c(graphColors, "#000000"), c(sort(chemicals), "black")), breaks = sort(chemicals)) +
		              guides(shape = FALSE) +
		              xlim(0, 24) +
		              theme(plot.title = element_text(size=22, hjust=0.5));
		            
		            if (forceScale)
		            {
		              p <- p +
		                ylim(ymin, ymax)
		            }
		            if(nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == log(private$logZeroVal, private$logBase), ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == log(private$logZeroVal, private$logBase), ], 
		                           mapping = aes(x=as.numeric(bmr), y=as.numeric(eval(parse(text=BMDres))), size = 9, shape = "Circle", stroke = 1.5, color = "black")) +
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            if (logVal & nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == log(private$logInfVal, private$logBase), ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == log(private$logInfVal, private$logBase), ], 
		                           mapping = aes(x=as.numeric(bmr), y=as.numeric(eval(parse(text=BMDres))), size = 9, shape = "Circle", stroke = 1.5, color = "black"))+
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            if (!logVal & nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == private$logInfVal, ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == private$logInfVal, ], 
		                           mapping = aes(x=as.numeric(bmr), y=as.numeric(eval(parse(text=BMDres))), size = 9, shape = "Circle", stroke = 1.5, color = "black"))+
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            p
		          }
		        )
		      })
		      i <- i + 1;
		      progress$inc(1/(2*length(endpoints)), detail = "Generating GeneTox Plots");
		    }
		  }
		},
		
		#create and show the plots for all precalculated BMR values (Distribution)
		plotGeneToxBMRDistribution = function(input, output, progress, logVal = FALSE, showPlotData = FALSE, forceScale = FALSE) {
		  loginfo("Creating GeneTox plots");
		  chemicals <- private$chem_data['chemical'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  endpoints <- private$chem_data['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  BMDres <- as.character(input$gt_bmd_result_chooser);
		  
		  if (length(input$select_chemical_gt) > length(chemicals))
		  {
		    logwarn(paste("Missing chemical info for: ", paste(input$select_chemical_gt[which(!input$select_chemical_gt %in% chemicals)], collapse = ", "), sep = ""));
		    showNotification("Some chemicals might be missing because no data was associated with them.", type = "warning", duration = 5);
		  }
		  
		  chemDataTable <- list();
		  
		  xmin <- min(private$chem_data[!is.na(private$chem_data[input$gt_bmd_result_chooser]) & private$chem_data[input$gt_bmd_result_chooser] != -Inf, input$gt_bmd_result_chooser], na.rm = TRUE);
		  xmax <- max(private$chem_data[!is.na(private$chem_data[input$gt_bmd_result_chooser]) & private$chem_data[input$gt_bmd_result_chooser] != Inf, input$gt_bmd_result_chooser], na.rm = TRUE);
		  
		  if (showPlotData)
		  {
  		  insertUI(
  		    selector = "#ui_gt_plotContainerCenter",
  		    where = "afterBegin",
  		    ui = box(status = "primary", title = "Plotting Data Table", collapsible = TRUE, width = 12,
  		             column(width = 8, offset = 2,
  		                    wellPanel(
  		                      h4("Plotting Data Table"),
  		                      DT::dataTableOutput(outputId = "table_GTPlottingData3")
  		                    )
  		             )
  		    )
  		  )
  		  output$table_GTPlottingData3 <- DT::renderDataTable(server = FALSE, {
  		    
  		    datatable(eval(parse(text=paste("private$chem_data %>% select(chemical, endpoint_grp, bmr, ", BMDres, ", trend, model_type) %>%",
  		                                    "dplyr::rename(\"Chemical\" = chemical, \"Endpoint Group\" = endpoint_grp, \"BMR\" = bmr, \"", 
  		                                    ifelse(logVal, paste("log(", toupper(BMDres), ")", sep=""), toupper(BMDres)), 
  		                                    "\" = ", BMDres, ", \"Trend\" = trend, \"Model Type\" = model_type)", sep=""))),
  		              selection="none", filter="bottom", extensions = "Buttons",
  		              options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
  		                           pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
  		                           scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
  		      formatStyle(seq(5), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
  		    
  		  })
		  }
		  
		  if (!is.null(private$assay_desc))
		  {
		    insertUI(
		      selector = "#ui_gt_plotContainerCenter",
		      where = "beforeEnd",
		      ui = box(status = "primary", title = "Assay Description", collapsible = TRUE, width = 12,
		               column(width = 8, offset = 2,
		                      wellPanel(
		                        h4("Assay Description"),
		                        DT::dataTableOutput(outputId = "table_GTAssayDescription3")
		                      )
		               )
		      )
		    )
		    output$table_GTAssayDescription3 <- DT::renderDataTable(server = FALSE, {
		      datatable(private$assay_desc,
		                selection="none", filter="bottom", extensions = "Buttons",
		                options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                             pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                             scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		        formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		      
		    })
		  }
		  insertUI(
		    selector = "#ui_gt_plotContainerCenter",
		    where = "beforeEnd",
		    ui = fluidRow(id = "ui_gt_plots",
		                  column(12,
		                         box(status = "primary", title = NULL, collapsible = FALSE, solidHeader = FALSE, width = 12,
		                             dropdownButton(
		                               tags$h3("About the plot"),
		                               tags$h5("This is a plot showing the distribution of the BMD over different BMR values for each chemicals and endpoints."),
		                               tags$h5("Crosses represent a dataset without trend, according to TCPL's MC4 test."),
		                               tags$h5("Circled points represent a value that is tending towards positive or negative infinity."),
		                               circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
		                               tooltip = NULL)
		                         )
		                  )
		    )
		  );
		  
		  i <- 1;
		  if (input$gt_comparaison == 1)
		  {
		    for (c in chemicals)
		    {
		      chemDataTable[[c]] <- private$chem_data %>% filter(chemical == c);
		      for (endpt in endpoints)
		      {
		        if (nrow(chemDataTable[[c]] %>% filter(endpoint_grp == endpt)) ==0)
		        {
		          chemDataTable[[c]] = rbind(chemDataTable[[c]], as.data.frame(list(chemical=c, endpoint_grp = endpt, bmr = NA_real_, bmdl = NA_real_, bmd = NA_real_, bmdu = NA_real_, trend = FALSE, model_type = NA)));
		        }
		      }
		      if ((i-1) %% 2 == 0)
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerLeft",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotOutput(outputId = c, height=405)
		                               )
		                        )
		          )
		        );
		      }
		      else
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerRight",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotOutput(outputId = c, height=405)
		                               )
		                        )
		          )
		        );
		      }
		      local({
		        myChem <- c;
		        output[[myChem]] <- renderPlot(
		          {
		            gg_color_list <- function(n) {
		              hues = seq(15, 375, length = n+1)
		              hcl(h = hues, l = 65, c = 100)[1:n]
		            }
		            
		            graphColors <- gg_color_list(4);
		            
		            if (logVal)
		            {
		              NInfVal <- min(ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myChem]])[input$gt_bmd_result_chooser] != -Inf, ])[input$gt_bmd_result_chooser]) == 0, Inf, min(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myChem]])[input$gt_bmd_result_chooser] != -Inf, ])[input$gt_bmd_result_chooser], na.rm = TRUE)), 
		                             ifelse(forceScale, xmin, Inf), na.rm = TRUE)-1;
		              InfVal <- max(ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myChem]])[input$gt_bmd_result_chooser] != Inf, ])[input$gt_bmd_result_chooser]) == 0, -Inf, max(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myChem]])[input$gt_bmd_result_chooser] != Inf, ])[input$gt_bmd_result_chooser], na.rm = TRUE)), 
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)+1;
		            }
		            else
		            {
		              NInfVal <- -1;
		              InfVal <- max(ifelse(nrow(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myChem]])[input$gt_bmd_result_chooser] != Inf, ])[input$gt_bmd_result_chooser]) == 0, -Inf, max(((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myChem]])[input$gt_bmd_result_chooser] != Inf, ])[input$gt_bmd_result_chooser], na.rm = TRUE)), 
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)*1.1;
		            }
		            
		            if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myChem]])[input$gt_bmd_result_chooser] == Inf,]) > 0)
		            {chemDataTable[[myChem]][!is.na(chemDataTable[[myChem]][input$gt_bmd_result_chooser]) & chemDataTable[[myChem]][input$gt_bmd_result_chooser] == Inf,input$gt_bmd_result_chooser] <- InfVal;}
		            if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myChem]])[input$gt_bmd_result_chooser] == -Inf,]) > 0)
		            {chemDataTable[[myChem]][!is.na(chemDataTable[[myChem]][input$gt_bmd_result_chooser]) & chemDataTable[[myChem]][input$gt_bmd_result_chooser] == -Inf,input$gt_bmd_result_chooser] <- NInfVal;}
		            
		            p <- ggplot(chemDataTable[[myChem]], aes(x=as.numeric(eval(parse(text=BMDres))), y=endpoint_grp, color = str_pad(as.character(bmr), 2, pad = "0"), shape = as.character(trend))) +
		              ggtitle(myChem) +
		              labs(y="Endpoints", x=paste(ifelse(logVal, paste("log(", toupper(BMDres), ")", sep=""), toupper(BMDres)),"Values [uM]",sep=" "), color="BMR Values [%]") +
		              geom_point(size = 3) +
		              scale_shape_manual(values = setNames(c(16, 4, 1), c("TRUE", "FALSE", "Circle"))) +
		              scale_size_identity() + 
		              scale_color_manual(values = setNames(c(graphColors, "#000000"), c("05", "10", "15", "20", "black")), breaks = c("05", "10", "15", "20")) +
		              guides(shape = FALSE) +
		              theme(plot.title = element_text(size=22, hjust=0.5));
		            
		            if (forceScale)
		            {
		              p <- p +
		                xlim(xmin-1.5, ifelse(logVal, xmax+1.5, xmax*1.21))
		            }
		            if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == NInfVal, ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == NInfVal, ], 
		                           mapping = aes(x=as.numeric(eval(parse(text=BMDres))), y=endpoint_grp, size = 9, shape = "Circle", stroke = 1.5, color = "black")) +
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            if (nrow((chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == InfVal, ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myChem]])[!is.na((chemDataTable[[myChem]])[BMDres]) & (chemDataTable[[myChem]])[BMDres] == InfVal, ], 
		                           mapping = aes(x=as.numeric(eval(parse(text=BMDres))), y=endpoint_grp, size = 9, shape = "Circle", stroke = 1.5, color = "black"))+
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            p
		          }
		        )
		      })
		      i <- i + 1;
		      progress$inc(1/(2*length(chemicals)), detail = "Generating GeneTox Plots");
		    }
		  }
		  else
		  {
		    for (e in endpoints)
		    {
		      chemDataTable[[e]] <- private$chem_data %>% filter(endpoint_grp == e);
		      for (ch in chemicals)
		      {
		        if (nrow(chemDataTable[[e]] %>% filter(chemical == ch)) ==0)
		        {
		          chemDataTable[[e]] = rbind(chemDataTable[[e]], as.data.frame(list(chemical=ch, endpoint_grp = e, bmr = NA_real_, bmdl = NA_real_, bmd = NA_real_, bmdu = NA_real_, trend = FALSE, model_type = NA)));
		        }
		      }
		      if ((i-1) %% 2 == 0)
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerLeft",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotOutput(outputId = e, height=405)
		                               )
		                        )
		          )
		        );
		      }
		      else
		      {
		        insertUI(
		          selector = "#ui_gt_plotContainerRight",
		          where = "beforeEnd",
		          ui = fluidRow(id = "ui_gt_plots",
		                        column(12,
		                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
		                                   plotOutput(outputId = e, height=405)
		                               )
		                        )
		          )
		        );
		      }
		      local({
		        myEndpt <- e;
		        output[[myEndpt]] <- renderPlot(
		          {
		            gg_color_list <- function(n) {
		              hues = seq(15, 375, length = n+1)
		              hcl(h = hues, l = 65, c = 100)[1:n]
		            }
		            
		            graphColors <- gg_color_list(4)
		            
		            if (logVal)
		            {
		              NInfVal <- min(ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser] != -Inf, ])[input$gt_bmd_result_chooser]) == 0, Inf, min(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser] != -Inf, ])[input$gt_bmd_result_chooser], na.rm = TRUE)), 
		                             ifelse(forceScale, xmin, Inf), na.rm = TRUE)-1;
		              InfVal <- max(ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser] != Inf, ])[input$gt_bmd_result_chooser]) == 0, -Inf, max(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser] != Inf, ])[input$gt_bmd_result_chooser], na.rm = TRUE)), 
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)+1;
		            }
		            else
		            {
		              NInfVal <- -1;
		              InfVal <- max(ifelse(nrow(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser] != Inf, ])[input$gt_bmd_result_chooser]) == 0, -Inf, max(((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser] != Inf, ])[input$gt_bmd_result_chooser], na.rm = TRUE)), 
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)*1.1;
		            }
		            
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser] == Inf,]) > 0)
		            {chemDataTable[[myEndpt]][!is.na(chemDataTable[[myEndpt]][input$gt_bmd_result_chooser]) & chemDataTable[[myEndpt]][input$gt_bmd_result_chooser] == Inf,input$gt_bmd_result_chooser] <- InfVal;}
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser]) & (chemDataTable[[myEndpt]])[input$gt_bmd_result_chooser] == -Inf,]) > 0)
		            {chemDataTable[[myEndpt]][!is.na(chemDataTable[[myEndpt]][input$gt_bmd_result_chooser]) & chemDataTable[[myEndpt]][input$gt_bmd_result_chooser] == -Inf,input$gt_bmd_result_chooser] <- NInfVal;}
		            
		            p <- ggplot(chemDataTable[[myEndpt]], aes(x=as.numeric(eval(parse(text=BMDres))), y=chemical, color = str_pad(as.character(bmr), 2, pad = "0"), shape = as.character(trend))) +
		              ggtitle(myEndpt) +
		              labs(y="Chemicals", x=paste(ifelse(logVal, paste("log(", toupper(BMDres), ")", sep=""), toupper(BMDres)),"Values [uM]",sep=" "), color = "BMR Values [%]") +
		              geom_point(size = 3) +
		              scale_shape_manual(values = setNames(c(16, 4, 1), c("TRUE", "FALSE", "Circle"))) +
		              scale_size_identity() +
		              scale_color_manual(values = setNames(c(graphColors, "#000000"), c("05","10","15","20","black")), breaks = c("05", "10", "15", "20")) +
		              guides(shape = FALSE) +
		              theme(plot.title = element_text(size=22, hjust=0.5));
		            
		            if (forceScale)
		            {
		              p <- p +
		                xlim(xmin-1.5, ifelse(logVal, xmax+1.5, xmax*1.21))
		            }
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == NInfVal, ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == NInfVal, ], 
		                           mapping = aes(x=as.numeric(eval(parse(text=BMDres))), y=chemical, size = 9, shape = "Circle", stroke = 1.5, color = "black")) +
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            if (nrow((chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == InfVal, ]) > 0)
		            {
		              p <- p +
		                geom_point(data = (chemDataTable[[myEndpt]])[!is.na((chemDataTable[[myEndpt]])[BMDres]) & (chemDataTable[[myEndpt]])[BMDres] == InfVal, ], 
		                           mapping = aes(x=as.numeric(eval(parse(text=BMDres))), y=chemical, size = 9, shape = "Circle", stroke = 1.5, color = "black"))+
		                labs(caption = "Circled values are tending towards positive or negative infinity.");
		            }
		            p
		          }
		        )
		      })
		      i <- i + 1;
		      progress$inc(1/(2*length(endpoints)), detail = "Generating GeneTox Plots");
		      
		    }
		  }
		},
		
		#Heatmap with dendrogram generation
		plotGeneToxHeatMap = function (input, output, bmrValue = 5, progress, logVal = FALSE, showPlotData = FALSE)
		{
		  loginfo("Creating GeneTox plots");
		  private$chem_data <- private$chem_data %>% filter(.$bmr == bmrValue);
		  BMDres <- input$gt_bmd_result_chooser
		  
		  private$chem_data_gt <- aggregate(private$chem_data_gt$concentration, by = list(private$chem_data_gt$chemical, private$chem_data_gt$endpoint_grp), max);
		  colnames(private$chem_data_gt) <- c("chemical", "endpoint_grp", "concentration");
		  if (logVal)
		  {
		    private$chem_data_gt[, "concentration"] <- log(private$chem_data_gt[, "concentration"], private$logBase);
		  }

		  chemicals <- private$chem_data['chemical'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  endpoints <- private$chem_data['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  
		  if (length(input$select_chemical_gt) > length(chemicals))
		  {
		    logwarn(paste("Missing chemical info for: ", paste(input$select_chemical_gt[which(!input$select_chemical_gt %in% chemicals)], collapse = ", "), sep = ""));
		    showNotification("Some chemicals might be missing because no data was associated with them.", type = "warning", duration = 5);
		  }
		  
		  heatMapData <- matrix(NA_real_, nrow = length(endpoints), ncol = length(chemicals));
		  trendData <- matrix(FALSE, nrow = length(endpoints), ncol = length(chemicals));
		  spareInfo <- matrix(paste(ifelse(logVal, paste("log(", toupper(BMDres),")",sep=""), toupper(BMDres)), ": NA</br>P-Value: NA</br>Model: NA<extra></extra>",sep=""), nrow = length(endpoints), ncol = length(chemicals));
		  
		  
		  #paste(ifelse(logVal, paste("log(", toupper(BMDres),")",sep=""), toupper(BMDres)), ": ", ifelse(heatMapData == log(private$logZeroVal, private$logBase), "-Inf", "%{z}"), sep=""),
		  colnames(heatMapData) <- colnames(trendData) <- colnames(spareInfo)<- chemicals;
		  rownames(heatMapData) <- rownames(trendData) <- rownames(spareInfo) <- endpoints;
		  		  
		  if (nrow(private$chem_data) > 0)
		  {
		    if (logVal)
		    {
		      NInfVal <- ifelse(nrow((private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] != -Inf, ])[BMDres]) == 0, Inf, min((private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] != -Inf, ])[BMDres], na.rm = TRUE))-1;
		      InfVal <- ifelse(nrow((private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] != Inf, ])[BMDres]) == 0, -Inf, max((private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] != Inf, ])[BMDres], na.rm = TRUE))+1;
		    }
		    else
		    {
		      NInfVal <- -1;
		      InfVal <- ifelse(nrow((private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] != Inf, ])[BMDres]) == 0, -Inf, max((private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] != Inf, ])[BMDres], na.rm = TRUE))*1.1;
		    }
		    
		    if (nrow(private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] == Inf,]) > 0)
		    {private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] == Inf, BMDres] <- InfVal;}
		    if (nrow(private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] == -Inf,]) > 0)
		    {private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] == -Inf, BMDres] <- NInfVal;}
		    
  		  for (i in seq(nrow(private$chem_data)))
  		  {
  		    heatMapData[private$chem_data[i, "endpoint_grp"], private$chem_data[i, "chemical"]] <- private$chem_data[i, BMDres];
  		    trendData[private$chem_data[i, "endpoint_grp"], private$chem_data[i, "chemical"]] <- private$chem_data[i, "trend"];
  		    spareInfo[private$chem_data[i, "endpoint_grp"], private$chem_data[i, "chemical"]] <- paste(ifelse(logVal, paste("log(", toupper(BMDres),")",sep=""), toupper(BMDres)), ": ", 
  		                                                                                                  ifelse(private$chem_data[i, BMDres] == NInfVal, "-Inf", 
  		                                                                                                         ifelse(private$chem_data[i, BMDres] == InfVal, "Inf", 
  		                                                                                                                formatC(signif(private$chem_data[i, BMDres]), digits=private$hoverSignificantDigits, flag="#"))),
  		                                                                                               "</br>Trend: ", ifelse(private$chem_data[i, "trend"], "Yes", "No"), 
  		                                                                                               "</br>Model: ", private$chem_data[i, "model_type"],
  		                                                                                               ifelse(private$chem_data[i, BMDres] == NInfVal, 
  		                                                                                                      paste("</br>This", toupper(BMDres), "value is tending towards negative infinity."), 
  		                                                                                                      ifelse(private$chem_data[i, BMDres] == InfVal, 
  		                                                                                                             paste("</br>This", toupper(BMDres), "value is tending towards infinity."), "")),
  		                                                                                               "<extra></extra>", sep="");
  		  }
		  }
		  else
		  {
		    return(NULL);
		  }
		  
		  maxConc <- heatMapData;
		  maxConc.NaRow <- (which(is.na(maxConc))-1) %% nrow(maxConc) + 1;
		  maxConc.NaCol <- (which(is.na(maxConc))-1) %/% nrow(maxConc) + 1;
		  
		  if (length(maxConc.NaRow) > 0)
		  {
		    for (i in seq(maxConc.NaRow))
		    {
		      r = maxConc.NaRow[i];
		      c = maxConc.NaCol[i];
		      maxConc[rownames(maxConc)[r], colnames(maxConc)[c]] <- as.numeric((private$chem_data_gt %>% filter(chemical == colnames(maxConc)[c], endpoint_grp == rownames(maxConc)[r]))["concentration"]);
		    }
		  }
		  
		  progress$inc(1/4, detail = "Generating GeneTox Plots");
		  
		  if (length(chemicals) > 1)
		  {
  		  dd.col <- as.dendrogram(hclust(dist(t(maxConc))));
  		  col.order <- order.dendrogram(dd.col);
		  }
		  else
		  {
		    col.order = c(1);
		  }
		  
		  if (length(endpoints) > 1)
		  {
  		  dd.row <- as.dendrogram(hclust(dist(maxConc)));
  		  row.order <- order.dendrogram(dd.row);
		  }
		  else
		  {
		    row.order = c(1);
	    }
		  
		  if (!(length(chemicals) > 1 && length(endpoints) > 1))
		  {
		    heatMapData <- matrix(heatMapData[row.order, col.order], nrow = length(endpoints), ncol = length(chemicals));
		    trendData <- matrix(trendData[row.order, col.order], nrow = length(endpoints), ncol = length(chemicals));
		    spareInfo <- matrix(spareInfo[row.order, col.order], nrow = length(endpoints), ncol = length(chemicals));
		    
		    colnames(heatMapData) <- colnames(trendData) <- colnames(spareInfo) <- chemicals[col.order];
		    rownames(heatMapData) <- rownames(trendData) <- rownames(spareInfo) <- endpoints[row.order];
		  }
		  else 
		  {
		    heatMapData <- heatMapData[row.order, col.order];
		    trendData <- trendData[row.order, col.order];
		    spareInfo <- spareInfo[row.order, col.order];
		  }
		  
		  valToInf <- trendAnnotation <- matrix("", nrow = nrow(trendData), ncol= ncol(trendData))
		  
		  colnames(trendAnnotation) <- colnames(valToInf) <- colnames(trendData);
		  rownames(trendAnnotation) <- rownames(valToInf) <- rownames(trendData);
		  
		  for (r in seq(nrow(trendData)))
		  {
		    for (c in seq(ncol(trendData)))
		    {
		      trendAnnotation[r, c] <- ifelse(!is.na(trendData[r, c]) & trendData[r, c], "", "\ud7");
		      valToInf[r, c] <- ifelse(!is.na(heatMapData[r, c]) & (heatMapData[r, c] == NInfVal | heatMapData[r, c] == InfVal), "O", "");
		    }
		  }
		  
		  if (showPlotData)
		  {
  		  insertUI(
  		    selector = "#ui_gt_plotContainerCenter",
  		    where = "afterBegin",
  		    ui = box(status = "primary", title = "Plotting Data Table", collapsible = TRUE, width = 12,
  		             column(width = 8, offset = 2,
  		                    wellPanel(
  		                      h4("Plotting Data Table"),
  		                      DT::dataTableOutput(outputId = "table_GTPlottingData4")
  		                    )
  		             )
  		    )
  		  )
  		  output$table_GTPlottingData4 <- DT::renderDataTable(server = FALSE, {
  		    
  		    datatable(eval(parse(text=paste("private$chem_data %>% select(chemical, endpoint_grp, bmr, ", BMDres, ", trend, model_type) %>%",
  		                                    "dplyr::rename(\"Chemical\" = chemical, \"Endpoint Group\" = endpoint_grp, \"BMR\" = bmr, \"", 
  		                                    ifelse(logVal, paste("log(", toupper(BMDres), ")", sep=""), toupper(BMDres)), 
  		                                    "\" = ", BMDres, ", \"Trend\" = trend, \"Model Type\" = model_type)", sep=""))),
  		              selection="none", filter="bottom", extensions = "Buttons",
  		              options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
  		                           pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
  		                           scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
  		      formatStyle(seq(5), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
  		    
  		  })
		  }
		  
		  if (!is.null(private$assay_desc))
		  {
		    insertUI(
		      selector = "#ui_gt_plotContainerCenter",
		      where = "beforeEnd",
		      ui = box(status = "primary", title = "Assay Description", collapsible = TRUE, width = 12,
		               column(width = 8, offset = 2,
		                      wellPanel(
		                        h4("Assay Description"),
		                        DT::dataTableOutput(outputId = "table_GTAssayDescription4")
		                      )
		               )
		      )
		    )
		    output$table_GTAssayDescription4 <- DT::renderDataTable(server = FALSE, {
		      datatable(private$assay_desc,
		                selection="none", filter="bottom", extensions = "Buttons",
		                options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                             pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                             scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		        formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		      
		    })
		  }
		  
		  insertUI(
		    selector = "#ui_gt_plotContainerCenter",
		    where = "beforeEnd",
		    ui = fluidRow(id = "ui_gt_plots",
		                  column(12,
		                         box(status = "primary", title = NULL, collapsible = FALSE, solidHeader = FALSE, width = 12,
		                             dropdownButton(
		                               tags$h3("About the plot"),
		                               tags$h5(paste("This is a heatmap of", toupper(BMDres), "with dendrograms. The color of each",
		                                             "square represents the", toupper(BMDres), "Value, the reddest being the lowest value.", sep = " ")),
		                               tags$h5("The \ud7 symbol represents a dataset without trend, according to TCPL's MC4 test."),
		                               tags$h5("The O symbol represents a value that is tending towards positive or negative infinity."),
		                               circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
		                               tooltip = NULL),
		                             plotlyOutput(outputId = "plot", height=810)
		                         )
		                  )
		    )
		  );
		  progress$inc(1/4, detail = "Generating GeneTox Plots");
		  
		  output$plot <- renderPlotly({
		    empty <- ggplot() + theme_void();
		    
		    if (length(endpoints) > 1)
		    {
  		    dgx <- (ggplot(segment(dendro_data(reorder(dd.row, row.order), type = "rectangle"))) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  		      coord_flip() + labs(x = "", y = "") + theme_minimal() + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())) %>%
  		      ggplotly(tooltip = "");
		    }
		    else
		    {
		      dgx <- ggplot() + theme_void();
		    }
		    
		    if (length(chemicals) > 1)
		    {
  		    dgy <- (ggplot(segment(dendro_data(reorder(dd.col, col.order), type = "rectangle"))) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  		      labs(x = "", y = "") + theme_minimal() + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())) %>%
  		      ggplotly(tooltip = "");
		    }
		    else
		    {
		      dgy <- ggplot() + theme_void();
		    }

		    hmdz <- heatMapData;
		    hmdz[is.na(hmdz)] <- 0;
		    
		    hm <- plot_ly() %>%
		      add_trace(x = colnames(heatMapData), y = rownames(heatMapData), z = hmdz, type = "heatmap", colors = c("gray40", "gray40"), showscale = FALSE, colorscale = list(c(0, "rgb(0, 0, 0)"), list(1, "rgb(0 , 0, 0)")), hoverinfo = "text", hovertext="") %>%
		      add_trace(x = colnames(heatMapData), y = rownames(heatMapData), z = heatMapData, text=spareInfo, type = "heatmap", colors=c("red", "white"), colorscale = list(c(0, "red"), c(1, "white")), 
		                hovertemplate=paste(
		                  "Chemical: %{x}</br>",
		                  "Endpoint: %{y}",
		                  "%{text}",
		                  sep="</br>"))%>%
		      layout(xaxis=list(title="Chemicals"), yaxis=list(title="Endpoints")) %>%
		      add_annotations(x = rep(colnames(heatMapData), each = nrow(heatMapData)), y = rep(rownames(heatMapData), times = ncol(heatMapData)), text = trendAnnotation, showarrow = FALSE, ax = 20, ay = -20) %>%
		      add_annotations(x = rep(colnames(heatMapData), each = nrow(heatMapData)), y = rep(rownames(heatMapData), times = ncol(heatMapData)), text = valToInf, showarrow = FALSE, ax = 20, ay = -20) %>%
		      layout(annotations = list(x=1, y=-0.4, text = "The \ud7 symbol represents data without a trend.\nThe O symbol represents a value that is tending towards positive or negative infinity.", 
		                                showarrow = FALSE, xref='paper', yref='paper', xanchor='right', yanchor='auto'));
		    
		    if(length(endpoints) > 1)
		    {
  		     plot <- subplot(dgy, empty, hm, dgx, nrows = 2, margin = 0, widths = c(0.8, 0.2), heights = c(0.2, 0.8),
  		             shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE);
		    }
		    else
		    {
		      plot <- subplot(dgy, hm, nrows = 2, margin = 0, widths = c(1), heights = c(0.2, 0.8),
		              shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE);
		    }
		      
		    plot <- plot %>% layout( margin=list(b = 200, l = 200, t = 200, r = 200), hovermode = "compare", title="Hierachically Clustered BMD Heatmap" );
		  });
		},
		
		#BMD Range Graph (Multiple Chemicals)
		plotGeneToxSingleBMRAll = function(input, output, bmrValue = 5, progress, logVal = FALSE, showPlotData = FALSE) {
		  loginfo("Creating GeneTox plots");
		  private$chem_data <- private$chem_data %>% filter(., .$bmr == bmrValue) %>%
		    mutate(., chemEndpt = paste(.$chemical, " - ", .$endpoint_grp, sep=""));
		  
		  #Check to make sure every chemical has at least 1 part of data to it.
		  chemicals <- private$chem_data['chemical'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  if (length(input$select_chemical_gt) > length(chemicals))
		  {
		    logwarn(paste("Missing chemical info for: ", paste(input$select_chemical_gt[which(!input$select_chemical_gt %in% chemicals)], collapse = ", "), sep = ""));
		    showNotification("Some chemicals might be missing because no data was associated with them.", type = "warning", duration = 5);
		  }
		  
		  if (showPlotData)
		  {
  		  insertUI(
  		    selector = "#ui_gt_plotContainerCenter",
  		    where = "afterBegin",
  		    ui = box(status = "primary", title = "Plotting Data Table", collapsible = TRUE, width = 12,
  		             column(width = 8, offset = 2,
  		                    wellPanel(
  		                      h4("Plotting Data Table"),
  		                      DT::dataTableOutput(outputId = "table_GTPlottingData5")
  		                    )
  		             )
  		    )
  		  )
  		  output$table_GTPlottingData5 <- DT::renderDataTable(server = FALSE, {
  		    datatable(eval(parse(text=paste("private$chem_data %>% select(chemical, endpoint_grp, bmr, bmdl, bmd, bmdu, trend, model_type) %>% ",
  		                                    "dplyr::rename(\"Chemical\" = chemical, \"Endpoint Group\" = endpoint_grp, \"BMR\" = bmr, ", ifelse(logVal, "\"log(BMDL)\"", "\"BMDL\""),
  		                                    " = bmdl, ", ifelse(logVal, "\"log(BMD)\"", "\"BMD\""), " = bmd, ", ifelse(logVal, "\"log(BMDU)\"", "\"BMDU\""), " = bmdu, \"Trend\" = trend, ",
  		                                    "\"Model Type\" = model_type)", sep=""))),
  		              selection="none", filter="bottom", extensions = "Buttons",
  		              options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
  		                           pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
  		                           scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
  		      formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
  		    
  		  })
		  }
		  
		  if (!is.null(private$assay_desc))
		  {
		    insertUI(
		      selector = "#ui_gt_plotContainerCenter",
		      where = "beforeEnd",
		      ui = box(status = "primary", title = "Assay Description", collapsible = TRUE, width = 12,
		               column(width = 8, offset = 2,
		                      wellPanel(
		                        h4("Assay Description"),
		                        DT::dataTableOutput(outputId = "table_GTAssayDescription5")
		                      )
		               )
		      )
		    )
		    output$table_GTAssayDescription5 <- DT::renderDataTable(server = FALSE, {
		      datatable(private$assay_desc,
		                selection="none", filter="bottom", extensions = "Buttons",
		                options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                             pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                             scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		        formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		      
		    })
		  }
		  
		  insertUI(
		    selector = "#ui_gt_plotContainerCenter",
		    where = "beforeEnd",
		    ui = fluidRow(id = "ui_gt_plots",
		                  column(12,
		                         box(status = "primary", title = NULL, collapsible = FALSE, solidHeader = FALSE, width = 12,
		                             dropdownButton(
		                               tags$h3("About the plot"),
		                               tags$h5("This is a plot showing the grouped confidence interval for all chemicals and endpoints."),
		                               tags$h5("Red lines represent a dataset without trend, according to TCPL's MC4 test."),
		                               tags$h5("Blue stars represent a value that is tending towards positive or negative infinity."),
		                               circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
		                               tooltip = NULL),
		                             plotlyOutput(outputId = "BMRRange", height=nrow(private$chem_data)*59+71)
		                         )
		                  )
		    )
		  );
		  
      output$BMRRange <- renderPlotly(
        {
          if (logVal)
          {
            NInfVal <- min(ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != -Inf, ])['bmdl']) == 0, Inf, min((private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != -Inf, ])['bmdl'], na.rm = TRUE)), 
                           ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != -Inf, ])['bmd']) == 0, Inf, min((private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != -Inf, ])['bmd'], na.rm = TRUE)), 
                           ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != -Inf, ])['bmdu']) == 0, Inf, min((private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != -Inf, ])['bmdu'], na.rm = TRUE)), 
                           na.rm = TRUE)-1;
            InfVal <- max(ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max((private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
                          ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != Inf, ])['bmd']) == 0, -Inf, max((private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
                          ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max((private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)),
                          na.rm = TRUE)+1;
          }
          else
          {
            NInfVal <- -1;
            InfVal <- max(ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max((private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
                          ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != Inf, ])['bmd']) == 0, -Inf, max((private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
                          ifelse(nrow((private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max((private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)), 
                          na.rm = TRUE)*1.1;
          }
          
          if (nrow(private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] == Inf,]) > 0)
          {private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] == Inf,"bmdl"] <- InfVal;}
          if (nrow(private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] == Inf,]) > 0)
          {private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] == Inf,"bmd"] <- InfVal;}
          if (nrow(private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] == Inf,]) > 0)
          {private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] == Inf,"bmdu"] <- InfVal;}
          if (nrow(private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] == -Inf,]) > 0)
          {private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] == -Inf,"bmdl"] <- NInfVal;}
          if (nrow(private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] == -Inf,]) > 0)
          {private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] == -Inf,"bmd"] <- NInfVal;}
          if (nrow(private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] == -Inf,]) > 0)
          {private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] == -Inf,"bmdu"] <- NInfVal;}
          
          p <- ggplot(private$chem_data, aes(x=as.numeric(bmd), y=reorder(chemEndpt, -as.numeric(bmd)), color = trend)) +
            ggtitle("BMD Grouped Confidence Interval Graph") +
            labs(y="", x=ifelse(logVal, "log(BMD) Values [uM]", "BMD Values [uM]")) +
            geom_point(size = 2, mapping = aes(text = paste(
              paste("Chemical: ", chemical, sep=""),
              paste("Endpoint: ", endpoint_grp, sep=""),
              paste(ifelse(logVal, "log(BMDL): ", "BMDL: "), ifelse(bmdl == NInfVal, "-Inf", ifelse(bmdl == InfVal, "Inf", formatC(signif(bmdl), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
              paste(ifelse(logVal, "log(BMD): ", "BMD: "), ifelse(bmd == NInfVal, "-Inf", ifelse(bmd == InfVal, "Inf", formatC(signif(bmd), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
              paste(ifelse(logVal, "log(BMDU): ", "BMDU: "), ifelse(bmdu == NInfVal, "-Inf", ifelse(bmdu == InfVal, "Inf", formatC(signif(bmdu), digits=private$hoverSignificantDigits, flag="#"))), sep=""),
              paste("Trend: ", ifelse(trend, "Yes", "No"), sep=""),
              paste("Model: ", model_type, sep=""),
              sep="\n"))) +
            scale_color_manual(values = setNames(c("black", "red"), c(TRUE, FALSE))) +
            theme(plot.title = element_text(size=22, hjust=0.5), legend.position = "none");
          if (nrow(private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] == NInfVal, ]))
          {
            p <- p +
              geom_point(data = (private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] == NInfVal, ]),
                         mapping = aes(x = NInfVal-0.5, y = chemEndpt, shape = 8, text = "This BMDL value is tending towards negative infinity."), color = "blue") +
              scale_shape_identity();
          }
          if (nrow(private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] == InfVal, ]) > 0)
          {
            p <- p +
              geom_point(data = (private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] == InfVal, ]), 
                         mapping = aes(x = ifelse(logVal, InfVal+0.5, InfVal*1.1), y = chemEndpt, shape = 8, text = "This BMDU value is tending towards infinity."), color = "blue") +
              scale_shape_identity();
          }
          if (nrow(private$chem_data[(!is.na(private$chem_data["bmdl"]) & !is.na(private$chem_data["bmdu"])), ]) > 0)
          {
            p <- p + geom_errorbarh(data = private$chem_data[(!is.na(private$chem_data["bmdl"]) & !is.na(private$chem_data["bmdu"])), ], 
                                    aes(y=chemEndpt, xmax=as.numeric(bmdu), xmin=as.numeric(bmdl), x=as.numeric(bmd), text = "", 
                                        height=nrow(private$chem_data)/9), size=0.2, na.rm = TRUE);
          } 
          
          ggplotly(p, tooltip = "text");
        }
      )
      progress$inc(1/2, detail = "Generating GeneTox Plots");
		},
		
		#Show images of the Curve plot
		plotGeneToxCurvePlot = function(input, output, progress, bmrValue=5) {
		  loginfo("Creating GeneTox plots");
		  chemicals <- private$chem_data['chemical'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  endpoints <- private$chem_data['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  
		  if (length(input$select_chemical_gt) > length(chemicals))
		  {
		    logwarn(paste("Missing chemical info for: ", paste(input$select_chemical_gt[which(!input$select_chemical_gt %in% chemicals)], collapse = ", "), sep = ""));
		    showNotification("Some chemicals might be missing because no data was associated with them.", type = "warning", duration = 5);
		  }
		  
		  chemDataTable <- list();
		  
		  if (!is.null(private$assay_desc))
		  {
		    insertUI(
		      selector = "#ui_gt_plotContainerCenter",
		      where = "beforeEnd",
		      ui = box(status = "primary", title = "Assay Description", collapsible = TRUE, width = 12,
		               column(width = 8, offset = 2,
		                      wellPanel(
		                        h4("Assay Description"),
		                        DT::dataTableOutput(outputId = "table_GTAssayDescription6")
		                      )
		               )
		      )
		    )
		    output$table_GTAssayDescription6 <- DT::renderDataTable(server = FALSE, {
		      datatable(private$assay_desc,
		                selection="none", filter="bottom", extensions = "Buttons",
		                options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                             pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                             scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		        formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		      
		    })
		  }
		  
		  insertUI(
		    selector = "#ui_gt_plotContainerCenter",
		    where = "beforeEnd",
		    ui = fluidRow(id = "ui_gt_plots",
		                  column(12,
		                         box(status = "primary", title = NULL, collapsible = FALSE, solidHeader = FALSE, width = 12,
		                             dropdownButton(
		                               tags$h3("About the plot"),
		                               tags$h5("This is a plot showing something."),
		                               circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
		                               tooltip = NULL)
		                         )
		                  )
		    )
		  );
		  
		  i <- 1;
		  if (input$gt_comparaison == 1)
		  {
		    for (c in chemicals)
		    {
		      for (endpt in endpoints)
		      {
		        if ((i-1) %% 2 == 0)
		        {
		          insertUI(
		            selector = "#ui_gt_plotContainerLeft",
		            where = "beforeEnd",
		            ui = fluidRow(id = "ui_gt_plots",
		                          column(12,
		                                 box(status = "primary", title = paste(c, " - ", endpt, " - ", bmrValue, sep=""), collapsible = TRUE, width = 12,
		                                     imageOutput(outputId = paste(c, " - ", endpt, " - ", bmrValue, sep=""), height=310)
		                                 )
		                          )
		            )
		          );
		        }
		        else
		        {
		          insertUI(
		            selector = "#ui_gt_plotContainerRight",
		            where = "beforeEnd",
		            ui = fluidRow(id = "ui_gt_plots",
		                          column(12,
		                                 box(status = "primary", title = paste(c, " - ", endpt, " - ", bmrValue, sep=""), collapsible = TRUE, width = 12,
		                                     imageOutput(outputId = paste(c, " - ", endpt, " - ", bmrValue, sep=""), height=310)
		                                 )
		                          )
		            )
		          );
		        }
		        i <- i + 1;
		        local({
		          myVal <- paste(c, " - ", endpt, " - ", bmrValue, sep="");
		          output[[myVal]] <- renderImage(
		            {
		              return(
		                list(
		                  src = "data\\GeneTox_curve_fit_plots\\Sample_Curve_Fit_Plots.png",
		                  filetype = "image/png",
		                  alt = "curve plot",
		                  height = 300
		                )
		              )
		            }
		          )
		        })
		        progress$inc(1/(2*length(chemicals)*length(endpoints)), detail = "Generating GeneTox Plots");
		      }
		    }
		  
		  }
		  else
		  {
		    for (endpt in endpoints)
		    {
		      for (c in chemicals)
		      {
		        if ((i-1) %% 2 == 0)
		        {
		          insertUI(
		            selector = "#ui_gt_plotContainerLeft",
		            where = "beforeEnd",
		            ui = fluidRow(id = "ui_gt_plots",
		                          column(12,
		                                 box(status = "primary", title = paste(c, " - ", endpt, " - ", bmrValue, sep=""), collapsible = TRUE, width = 12,
		                                     imageOutput(outputId = paste(c, " - ", endpt, " - ", bmrValue, sep=""), height=310)
		                                 )
		                          )
		            )
		          );
		        }
		        else
		        {
		          insertUI(
		            selector = "#ui_gt_plotContainerRight",
		            where = "beforeEnd",
		            ui = fluidRow(id = "ui_gt_plots",
		                          column(12,
		                                 box(status = "primary", title = paste(c, " - ", endpt, " - ", bmrValue, sep=""), collapsible = TRUE, width = 12,
		                                     imageOutput(outputId = paste(c, " - ", endpt, " - ", bmrValue, sep=""), height=310)
		                                 )
		                          )
		            )
		          );
		        }
		        i <- i + 1;
		        local({
		          myVal <- paste(c, " - ", endpt, " - ", bmrValue, sep="");
		          output[[myVal]] <- renderImage(
		            {
		              return(
		                list(
		                  src = "data\\GeneTox_curve_fit_plots\\Sample_Curve_Fit_Plots.png",
		                  filetype = "image/png",
		                  alt = "curve plot",
		                  height = 300
		                )
		              )
		            }
		          )
		        })
		        progress$inc(1/(2*length(chemicals)*length(endpoints)), detail = "Generating GeneTox Plots");
		      }
		    }
		    
		  }
		},
		
		#create and show the plots the GeneToxPi Graphs
		plotGeneToxPi = function(input, output, bmrValue = 5, progress, logVal = FALSE, showPlotData = FALSE, forceScale = FALSE, forceOrder = FALSE) {
		  loginfo("Creating GeneTox plots");
		  private$chem_data <- private$chem_data %>% filter(.$bmr == bmrValue);
		  
		  print(private$chem_data)
		  
		  chemicals <- private$chem_data['chemical'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  endpoints <- private$chem_data['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  endpointsMain <- (private$chem_data[private$chem_data[,"assay"] != "multiflow",])['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  endpointsMulti <- (private$chem_data[private$chem_data[,"assay"] == "multiflow",])['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		  BMDres <- as.character(input$gt_bmd_result_chooser);
		  
		  print(private$chem_data)
		  
		  if (length(input$select_chemical_gt) > length(chemicals))
		  {
		    logwarn(paste("Missing chemical info for: ", paste(input$select_chemical_gt[which(!input$select_chemical_gt %in% chemicals)], collapse = ", "), sep = ""));
		    showNotification("Some chemicals might be missing because no data was associated with them.", type = "warning", duration = 5);
		  }
		  
		  chemDataTableGeneral <- list();
		  chemDataTableMultiFlow <- list();
		  
		  xmin <- min(min(private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != -Inf,'bmdl'], na.rm = TRUE), min(private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != -Inf,'bmd'], na.rm = TRUE), min(private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != -Inf,'bmdu'], na.rm = TRUE), na.rm = TRUE)
		  xmax <- max(max(private$chem_data[!is.na(private$chem_data["bmdl"]) & private$chem_data["bmdl"] != Inf,'bmdl'], na.rm = TRUE), max(private$chem_data[!is.na(private$chem_data["bmd"]) & private$chem_data["bmd"] != Inf,'bmd'], na.rm = TRUE), max(private$chem_data[!is.na(private$chem_data["bmdu"]) & private$chem_data["bmdu"] != Inf,'bmdu'], na.rm = TRUE), na.rm = TRUE)
		  ymin <- min(private$chem_data[!is.na(private$chem_data[BMDres]) & private$chem_data[BMDres] > 0 & private$chem_data[BMDres] != -Inf, BMDres], na.rm = TRUE)
		  
		  if (showPlotData)
		  {
		    insertUI(
		      selector = "#ui_gt_plotContainerCenter",
		      where = "afterBegin",
		      ui = box(status = "primary", title = "Plotting Data Table", collapsible = TRUE, width = 12,
		               column(width = 8, offset = 2,
		                      wellPanel(
		                        h4("Plotting Data Table"),
		                        DT::dataTableOutput(outputId = "table_GTPlottingData7")
		                      )
		               )
		      )
		    )
		    output$table_GTPlottingData7 <- DT::renderDataTable(server = FALSE, {
		      datatable(eval(parse(text=paste("private$chem_data %>% select(chemical, endpoint_grp, bmr, ", BMDres, ", trend, model_type) %>% ",
		                                      "dplyr::rename(\"Chemical\" = chemical, \"Endpoint Group\" = endpoint_grp, \"BMR\" = bmr, ", ifelse(logVal, paste("\"log(", toupper(BMDres)," +1)\"", sep =""), paste("\"", toupper(BMDres), "\"", sep="")),
		                                      " = ", BMDres, ", \"Trend\" = trend, ",
		                                      "\"Model Type\" = model_type)", sep=""))),
		                selection="none", filter="bottom", extensions = "Buttons",
		                options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                             pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                             scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		        formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		      
		    })
		  }
		  
		  if (!is.null(private$assay_desc))
		  {
		    insertUI(
		      selector = "#ui_gt_plotContainerCenter",
		      where = "beforeEnd",
		      ui = box(status = "primary", title = "Assay Description", collapsible = TRUE, width = 12,
		               column(width = 8, offset = 2,
		                      wellPanel(
		                        h4("Assay Description"),
		                        DT::dataTableOutput(outputId = "table_GTAssayDescription7")
		                      )
		               )
		      )
		    )
		    output$table_GTAssayDescription7 <- DT::renderDataTable(server = FALSE, {
		      datatable(private$assay_desc,
		                selection="none", filter="bottom", extensions = "Buttons",
		                options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                             pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                             scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		        formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		      
		    })
		  }
		  
		  insertUI(
		    selector = "#ui_gt_plotContainerCenter",
		    where = "beforeEnd",
		    ui = fluidRow(id = "ui_gt_plots",
		                  column(12,
		                         box(status = "primary", title = NULL, collapsible = FALSE, solidHeader = FALSE, width = 12,
		                             dropdownButton(
		                               tags$h3("About the plot"),
		                               tags$h5("This is NOT a plot showing the confidence interval for each chemicals and endpoints."),
		                               tags$h5("Red lines represent a dataset without trend, according to TCPL's MC4 test."),
		                               tags$h5("Blue stars represent a value that is tending towards positive or negative infinity."),
		                               circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
		                               tooltip = NULL)
		                         )
		                  )
		    )
		  );
		  
		  i <- 1;
		    for (c in chemicals)
		    {
		      chemDataTableGeneral[[c]] <- private$chem_data %>% filter(chemical == c) %>% filter (assay != "multiflow");
		      chemDataTableMultiFlow[[c]] <- private$chem_data %>% filter(chemical == c) %>% filter (assay == "multiflow");
		      #for (endpt in endpoints)
		      #{
		      #  if (nrow(chemDataTableGeneral[[c]] %>% filter(endpoint_grp == endpt)) ==0)
		      #  {
		      #    chemDataTableGeneral[[c]] = rbind(chemDataTableGeneral[[c]], as.data.frame(list(chemical=c, endpoint_grp = endpt, bmr = NA_real_, bmdl = NA_real_, bmd = NA_real_, bmdu = NA_real_, trend = FALSE, model_type=NA, assay="")));
		      #  }
	           # if (nrow(chemDataTableMultiFlow[[c]] %>% filter(endpoint_grp == endpt)) ==0)
	           # {
	           #   chemDataTableMultiFlow[[c]] = rbind(chemDataTableMultiFlow[[c]], as.data.frame(list(chemical=c, endpoint_grp = endpt, bmr = NA_real_, bmdl = NA_real_, bmd = NA_real_, bmdu = NA_real_, trend = FALSE, model_type=NA, assay="")));
	           # }
		      #}
		      
    	        insertUI(
    	          selector = "#ui_gt_plotContainerLeft",
    	          where = "beforeEnd",
    	          ui = fluidRow(id = "ui_gt_plots",
    	                        column(12,
    	                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
    	                                   plotOutput(outputId = paste(c," - Main",sep=""), height=350)
    	                               )
    	                        )
    	          )
    	        );
		       
    	        insertUI(
    	          selector = "#ui_gt_plotContainerRight",
    	          where = "beforeEnd",
    	          ui = fluidRow(id = "ui_gt_plots",
    	                        column(12,
    	                               box(status = "primary", title = "", collapsible = TRUE, width = 12,
    	                                   plotOutput(outputId = paste(c," - Multiflow", sep=""), height=350)
    	                               )
    	                        )
    	          )
    	        );
		        
		      local({
		        myChem <- c;
		        output[[paste(myChem," - Main",sep="")]] <- renderPlot(
		          {
		            if (logVal)
		            {
		              NInfVal <- min(ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdl"]) & (chemDataTableGeneral[[myChem]])["bmdl"] != -Inf, ])['bmdl']) == 0, Inf, min(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdl"]) & (chemDataTableGeneral[[myChem]])["bmdl"] != -Inf, ])['bmdl'], na.rm = TRUE)), 
		                             ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmd"]) & (chemDataTableGeneral[[myChem]])["bmd"] != -Inf, ])['bmd']) == 0, Inf, min(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmd"]) & (chemDataTableGeneral[[myChem]])["bmd"] != -Inf, ])['bmd'], na.rm = TRUE)), 
		                             ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdu"]) & (chemDataTableGeneral[[myChem]])["bmdu"] != -Inf, ])['bmdu']) == 0, Inf, min(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdu"]) & (chemDataTableGeneral[[myChem]])["bmdu"] != -Inf, ])['bmdu'], na.rm = TRUE)),
		                             ifelse(forceScale, xmin, Inf), na.rm = TRUE)-1;
		              InfVal <- max(ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdl"]) & (chemDataTableGeneral[[myChem]])["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdl"]) & (chemDataTableGeneral[[myChem]])["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
		                            ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmd"]) & (chemDataTableGeneral[[myChem]])["bmd"] != Inf, ])['bmd']) == 0, -Inf, max(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmd"]) & (chemDataTableGeneral[[myChem]])["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
		                            ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdu"]) & (chemDataTableGeneral[[myChem]])["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdu"]) & (chemDataTableGeneral[[myChem]])["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)),
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)+1;
		            }
		            else
		            {
		              NInfVal <- -1;
		              InfVal <- max(ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdl"]) & (chemDataTableGeneral[[myChem]])["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdl"]) & (chemDataTableGeneral[[myChem]])["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
		                            ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmd"]) & (chemDataTableGeneral[[myChem]])["bmd"] != Inf, ])['bmd']) == 0, -Inf, max(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmd"]) & (chemDataTableGeneral[[myChem]])["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
		                            ifelse(nrow(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdu"]) & (chemDataTableGeneral[[myChem]])["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max(((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdu"]) & (chemDataTableGeneral[[myChem]])["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)), 
		                            ifelse(forceScale, xmax, -Inf), na.rm = TRUE)*1.1;
		            }
		            
		            if (nrow((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdl"]) & (chemDataTableGeneral[[myChem]])["bmdl"] == Inf,]) > 0)
		            {chemDataTableGeneral[[myChem]][!is.na(chemDataTableGeneral[[myChem]]["bmdl"]) & chemDataTableGeneral[[myChem]]["bmdl"] == Inf,"bmdl"] <- InfVal;}
		            if (nrow((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmd"]) & (chemDataTableGeneral[[myChem]])["bmd"] == Inf,]) > 0)
		            {chemDataTableGeneral[[myChem]][!is.na(chemDataTableGeneral[[myChem]]["bmd"]) & chemDataTableGeneral[[myChem]]["bmd"] == Inf,"bmd"] <- InfVal;}
		            if (nrow((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdu"]) & (chemDataTableGeneral[[myChem]])["bmdu"] == Inf,]) > 0)
		            {chemDataTableGeneral[[myChem]][!is.na(chemDataTableGeneral[[myChem]]["bmdu"]) & chemDataTableGeneral[[myChem]]["bmdu"] == Inf,"bmdu"] <- InfVal;}
		            if (nrow((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdl"]) & (chemDataTableGeneral[[myChem]])["bmdl"] == -Inf,]) > 0)
		            {chemDataTableGeneral[[myChem]][!is.na(chemDataTableGeneral[[myChem]]["bmdl"]) & chemDataTableGeneral[[myChem]]["bmdl"] == -Inf,"bmdl"] <- NInfVal;}
		            if (nrow((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmd"]) & (chemDataTableGeneral[[myChem]])["bmd"] == -Inf,]) > 0)
		            {chemDataTableGeneral[[myChem]][!is.na(chemDataTableGeneral[[myChem]]["bmd"]) & chemDataTableGeneral[[myChem]]["bmd"] == -Inf,"bmd"] <- NInfVal;}
		            if (nrow((chemDataTableGeneral[[myChem]])[!is.na((chemDataTableGeneral[[myChem]])["bmdu"]) & (chemDataTableGeneral[[myChem]])["bmdu"] == -Inf,]) > 0)
		            {chemDataTableGeneral[[myChem]][!is.na(chemDataTableGeneral[[myChem]]["bmdu"]) & chemDataTableGeneral[[myChem]]["bmdu"] == -Inf,"bmdu"] <- NInfVal;}
		          
		            p <- ggplot(chemDataTableGeneral[[myChem]], aes(x=match(chemDataTableGeneral[[myChem]]$endpoint_grp, chemDataTableGeneral[[myChem]]$endpoint_grp), y=1/eval(parse(text=paste("chemDataTableGeneral[[myChem]]$",BMDres,sep=""))), fill=endpoint_grp)) + ggtitle(paste(myChem, " - Main", sep="")) +
		              geom_bar(width = 1, stat='identity') +
		              coord_polar(theta="x", start=0)+
		              theme(plot.title = element_text(size=22, hjust=0.5))+
		              theme(axis.title.x = element_blank())+
		              theme(axis.text.x = element_blank())+
		              ylab(ifelse(logVal, paste("1/log(", BMDres, " +1)", sep=""), paste("1/", BMDres, sep="")))+
		              labs(fill = "Endpoint");
		            
		            if (forceScale)
		            {
		              p <- p +
		                ylim(0, ifelse(logVal, (1/ymin)+1.5, (1/ymin)*1.21))
		              #Allow space for the Inf values and blue dots/stars
		            }
		            
		            p
		          }
		        )
		        
		        output[[paste(myChem," - Multiflow",sep="")]] <- renderPlot(
		            {
		                if (logVal)
		                {
		                    NInfVal <- min(ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdl"]) & (chemDataTableMultiFlow[[myChem]])["bmdl"] != -Inf, ])['bmdl']) == 0, Inf, min(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdl"]) & (chemDataTableMultiFlow[[myChem]])["bmdl"] != -Inf, ])['bmdl'], na.rm = TRUE)), 
		                                   ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmd"]) & (chemDataTableMultiFlow[[myChem]])["bmd"] != -Inf, ])['bmd']) == 0, Inf, min(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmd"]) & (chemDataTableMultiFlow[[myChem]])["bmd"] != -Inf, ])['bmd'], na.rm = TRUE)), 
		                                   ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdu"]) & (chemDataTableMultiFlow[[myChem]])["bmdu"] != -Inf, ])['bmdu']) == 0, Inf, min(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdu"]) & (chemDataTableMultiFlow[[myChem]])["bmdu"] != -Inf, ])['bmdu'], na.rm = TRUE)),
		                                   ifelse(forceScale, xmin, Inf), na.rm = TRUE)-1;
		                    InfVal <- max(ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdl"]) & (chemDataTableMultiFlow[[myChem]])["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdl"]) & (chemDataTableMultiFlow[[myChem]])["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
		                                  ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmd"]) & (chemDataTableMultiFlow[[myChem]])["bmd"] != Inf, ])['bmd']) == 0, -Inf, max(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmd"]) & (chemDataTableMultiFlow[[myChem]])["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
		                                  ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdu"]) & (chemDataTableMultiFlow[[myChem]])["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdu"]) & (chemDataTableMultiFlow[[myChem]])["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)),
		                                  ifelse(forceScale, xmax, -Inf), na.rm = TRUE)+1;
		                }
		                else
		                {
		                    NInfVal <- -1;
		                    InfVal <- max(ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdl"]) & (chemDataTableMultiFlow[[myChem]])["bmdl"] != Inf, ])['bmdl']) == 0, -Inf, max(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdl"]) & (chemDataTableMultiFlow[[myChem]])["bmdl"] != Inf, ])['bmdl'], na.rm = TRUE)), 
		                                  ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmd"]) & (chemDataTableMultiFlow[[myChem]])["bmd"] != Inf, ])['bmd']) == 0, -Inf, max(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmd"]) & (chemDataTableMultiFlow[[myChem]])["bmd"] != Inf, ])['bmd'], na.rm = TRUE)), 
		                                  ifelse(nrow(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdu"]) & (chemDataTableMultiFlow[[myChem]])["bmdu"] != Inf, ])['bmdu']) == 0, -Inf, max(((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdu"]) & (chemDataTableMultiFlow[[myChem]])["bmdu"] != Inf, ])['bmdu'], na.rm = TRUE)), 
		                                  ifelse(forceScale, xmax, -Inf), na.rm = TRUE)*1.1;
		                }
		                
		                if (nrow((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdl"]) & (chemDataTableMultiFlow[[myChem]])["bmdl"] == Inf,]) > 0)
		                {chemDataTableMultiFlow[[myChem]][!is.na(chemDataTableMultiFlow[[myChem]]["bmdl"]) & chemDataTableMultiFlow[[myChem]]["bmdl"] == Inf,"bmdl"] <- InfVal;}
		                if (nrow((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmd"]) & (chemDataTableMultiFlow[[myChem]])["bmd"] == Inf,]) > 0)
		                {chemDataTableMultiFlow[[myChem]][!is.na(chemDataTableMultiFlow[[myChem]]["bmd"]) & chemDataTableMultiFlow[[myChem]]["bmd"] == Inf,"bmd"] <- InfVal;}
		                if (nrow((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdu"]) & (chemDataTableMultiFlow[[myChem]])["bmdu"] == Inf,]) > 0)
		                {chemDataTableMultiFlow[[myChem]][!is.na(chemDataTableMultiFlow[[myChem]]["bmdu"]) & chemDataTableMultiFlow[[myChem]]["bmdu"] == Inf,"bmdu"] <- InfVal;}
		                if (nrow((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdl"]) & (chemDataTableMultiFlow[[myChem]])["bmdl"] == -Inf,]) > 0)
		                {chemDataTableMultiFlow[[myChem]][!is.na(chemDataTableMultiFlow[[myChem]]["bmdl"]) & chemDataTableMultiFlow[[myChem]]["bmdl"] == -Inf,"bmdl"] <- NInfVal;}
		                if (nrow((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmd"]) & (chemDataTableMultiFlow[[myChem]])["bmd"] == -Inf,]) > 0)
		                {chemDataTableMultiFlow[[myChem]][!is.na(chemDataTableMultiFlow[[myChem]]["bmd"]) & chemDataTableMultiFlow[[myChem]]["bmd"] == -Inf,"bmd"] <- NInfVal;}
		                if (nrow((chemDataTableMultiFlow[[myChem]])[!is.na((chemDataTableMultiFlow[[myChem]])["bmdu"]) & (chemDataTableMultiFlow[[myChem]])["bmdu"] == -Inf,]) > 0)
		                {chemDataTableMultiFlow[[myChem]][!is.na(chemDataTableMultiFlow[[myChem]]["bmdu"]) & chemDataTableMultiFlow[[myChem]]["bmdu"] == -Inf,"bmdu"] <- NInfVal;}
		                
		                p <- ggplot(chemDataTableMultiFlow[[myChem]], aes(x=match(chemDataTableMultiFlow[[myChem]]$endpoint_grp, chemDataTableMultiFlow[[myChem]]$endpoint_grp), y=1/eval(parse(text=paste("chemDataTableMultiFlow[[myChem]]$",BMDres,sep=""))), fill=endpoint_grp)) + ggtitle(paste(myChem, " - MultiFlow", sep="")) +
		                    geom_bar(width = 1, stat='identity') +
		                    coord_polar(theta="x", start=0)+
		                    theme(plot.title = element_text(size=22, hjust=0.5))+
		                    theme(axis.title.x = element_blank())+
		                    theme(axis.text.x = element_blank())+
		                    ylab(ifelse(logVal, paste("1/log(", BMDres, " +1)", sep=""), paste("1/", BMDres, sep="")))+
		                    labs(fill = "Endpoint");
		                
		                if (forceScale)
		                {
		                    p <- p +
		                        ylim(0, ifelse(logVal, (1/ymin)+1.5, (1/ymin)*1.21))
		                    #Allow space for the Inf values and blue dots/stars
		                }
		                
		                p
		            }
		        )
		      })
		      i <- i + 1;
		      progress$inc(1/(2*length(chemicals)), detail = "Generating GeneTox Plots");
		    }
		},
		#Table of Endpoints with Threshold
		plotEndpointsTable = function(input, output, progress, threshold, showPlotData = FALSE) {
		    loginfo("Creating GeneTox plots");
		    chemicals <- private$chem_data['chemical'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		    endpoints <- private$chem_data['endpoint_grp'] %>% unique() %>% unlist() %>% unname() %>% as.character();
		    BMDres <- as.character(input$gt_bmd_result_chooser);
		    
		    if (length(input$select_chemical_gt) > length(chemicals))
		    {
		        logwarn(paste("Missing chemical info for: ", paste(input$select_chemical_gt[which(!input$select_chemical_gt %in% chemicals)], collapse = ", "), sep = ""));
		        showNotification("Some chemicals might be missing because no data was associated with them.", type = "warning", duration = 5);
		    }
		    
		    chemDataTable <- matrix(nrow=length(chemicals), ncol=length(endpoints), dimnames = list(chemicals, endpoints));
		    
		    xmin <- min(private$chem_data[!is.na(private$chem_data[input$gt_bmd_result_chooser]) & private$chem_data[input$gt_bmd_result_chooser] != -Inf, input$gt_bmd_result_chooser], na.rm = TRUE);
		    xmax <- max(private$chem_data[!is.na(private$chem_data[input$gt_bmd_result_chooser]) & private$chem_data[input$gt_bmd_result_chooser] != Inf, input$gt_bmd_result_chooser], na.rm = TRUE);
		    
		    if (showPlotData)
		    {
		        insertUI(
		            selector = "#ui_gt_plotContainerCenter",
		            where = "afterBegin",
		            ui = box(status = "primary", title = "Plotting Data Table", collapsible = TRUE, width = 12,
		                     column(width = 8, offset = 2,
		                            wellPanel(
		                                h4("Plotting Data Table"),
		                                DT::dataTableOutput(outputId = "table_GTPlottingData8")
		                            )
		                     )
		            )
		        )
		        output$table_GTPlottingData8 <- DT::renderDataTable(server = FALSE, {
		            
		            datatable(eval(parse(text=paste("private$chem_data %>% select(chemical, endpoint_grp, bmr, ", BMDres, ", trend, model_type) %>%",
		                                            "dplyr::rename(\"Chemical\" = chemical, \"Endpoint Group\" = endpoint_grp, \"BMR\" = bmr, \"", 
		                                            toupper(BMDres), 
		                                            "\" = ", BMDres, ", \"Trend\" = trend, \"Model Type\" = model_type)", sep=""))),
		                      selection="none", filter="bottom", extensions = "Buttons",
		                      options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                                   pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                                   scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		                formatStyle(seq(5), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		            
		        })
		    }
		    
		    if (!is.null(private$assay_desc))
		    {
		        insertUI(
		            selector = "#ui_gt_plotContainerCenter",
		            where = "beforeEnd",
		            ui = box(status = "primary", title = "Assay Description", collapsible = TRUE, width = 12,
		                     column(width = 8, offset = 2,
		                            wellPanel(
		                                h4("Assay Description"),
		                                DT::dataTableOutput(outputId = "table_GTAssayDescription8")
		                            )
		                     )
		            )
		        )
		        output$table_GTAssayDescription8 <- DT::renderDataTable(server = FALSE, {
		            datatable(private$assay_desc,
		                      selection="none", filter="bottom", extensions = "Buttons",
		                      options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
		                                   pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
		                                   scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
		                formatStyle(seq(7), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
		            
		        })
		    }
		    insertUI(
		        selector = "#ui_gt_plotContainerCenter",
		        where = "beforeEnd",
		        ui = fluidRow(id = "ui_gt_plots",
		                      column(12,
		                             box(status = "primary", title = NULL, collapsible = FALSE, solidHeader = FALSE, width = 12,
		                                 dropdownButton(
		                                     tags$h3("About the plot"),
		                                     tags$h5("Something, something."), #Change this description
		                                     circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
		                                     tooltip = NULL)
		                             )
		                      )
		        )
		    );
		    
	        for (c in chemicals)
	        {
	            for (endpt in endpoints)
	            {
	                if (nrow(private$chem_data %>% filter(chemical == c) %>% filter(endpoint_grp == endpt)) == 0)
	                {
	                    chemDataTable[c, endpt] = "NA";
	                }
	                else if (((private$chem_data %>% filter(chemical == c) %>% filter(endpoint_grp == endpt))[1, BMDres]) < threshold)
	                {
	                    chemDataTable[c, endpt] = "-"
	                }
	                else if (((private$chem_data %>% filter(chemical == c) %>% filter(endpoint_grp == endpt))[1, BMDres]) > threshold)
	                {
	                    chemDataTable[c, endpt] = "+"
	                }
	                else
	                {
	                    chemDataTable[c, endpt] = "0"
	                }
	            }
	            
                insertUI(
                    selector = "#ui_gt_plotContainerCenter",
                    where = "beforeEnd",
                    ui = fluidRow(id = "ui_gt_plots",
                                  column(12,
                                         box(status = "primary", title = "", collapsible = TRUE, width = 12,
                                             tableOutput(outputId = "table8")
                                         )
                                  )
                    )
                );
	            output$table8 <- renderTable(chemDataTable, rownames = TRUE, spacing="xs")
	            progress$inc(1/(2*length(chemicals)), detail = "Generating GeneTox Plots");
	        }
		}
	)
)