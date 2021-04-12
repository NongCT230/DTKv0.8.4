# Basic statistical data  --------------------------------------------


#v0.7

#class contains basic statistical data and data calculating functions
# - minac50, minoed, target family counts, min and avg ac50 per target family, scalar top

#the design goal with this is that this is part of a composition where the basicstats data is all managed by another class and all the plotting is done by this class.
#The reason I am doing this is that it allows for easy expension. Since all you need is to include the data. 
Class.Analysis.BERAnalysis <- R6Class("Class.Analysis.BERAnalysis",
  #private variables and functions
  private = list(
  
	#BER plot variables
	pltBER = NULL,
	tblBER = NULL,
	tblMeanBER = NULL,
	pltBERvsAc50 = NULL,
	pltOEDvsAc50 = NULL
),

#public variables and functions
  public = list(
	BERData =  NULL,
	basicData = NULL,
    #constructor
    initialize = function(  ) {
		self$BERData <- Class.Analysis.BERData$new(); #Done here because it will be shared across all instances of the object otherwise.
		self$basicData <- Class.Analysis.Data$new(  ); #Done here because it will be shared across all instances of the object otherwise.
	},
    #finalizer
    finalize = function() {
    },
	
	
	
	#BER Analysis Region
	
	#computes the Agregated table for BER analysis
    computerMeanTableBER = function(){
      if (self$BERData$calcBERStatsDataExists() && self$basicData$basicStatsDataExists()) {
		private$tblMeanBER <- tribble(~casn,~name,~average_of_oral_consumer_products_exposure, ~minAC50, ~minOED,
		                              ~min_oed_over_oral_consumer_products_exposure, ~mean_oed_over_oral_consumer_products_exposure);
		#getting our private data
        calc_BER_stat_tbl <- self$BERData$getCalcBERStatsTable();
		chemical_casn_list <- unique(calc_BER_stat_tbl[["casn"]]);
		Ac50_stat_tbl <- self$basicData$getBasicStatsTable();
		Ac50_stat_tbl<- filter(Ac50_stat_tbl, above_cutoff == "Y", ac50 >= -2, ac50 <= 10000, cytotoxicity_um > 0.01) %>% drop_na(ac50);
		
		for(chemical_casn in chemical_casn_list){
			casn_tbl <- filter(calc_BER_stat_tbl, casn == chemical_casn);
			ac50_tbl <- filter(Ac50_stat_tbl, casn == chemical_casn);
			#ac50_tbl <- filter(ac50_tbl, ac50 > ac_cutoff);
			avg_mean <- mean(casn_tbl$direct_ingestion + casn_tbl$direct_vapor + casn_tbl$direct_aerosol);
			oral_ber <- (casn_tbl$direct_ingestion + casn_tbl$direct_vapor + casn_tbl$direct_aerosol);
			avg_ac50 <- mean(10^ac50_tbl$ac50);
			min_ac50 <- min(10^ac50_tbl$ac50);
			
			avg_oed <- mean(ac50_tbl$oed);
			min_oed <- min(ac50_tbl$oed);
			
			
			theorical_95th_percentile <- quantile(oral_ber,0.95);
			closest_to_95th <- absolute.min(oral_ber - theorical_95th_percentile) + theorical_95th_percentile ; #gets the closest value to the 95th percentile.
			
			
			private$tblMeanBER <- add_row(private$tblMeanBER, casn = chemical_casn, 
						name = as.character(filter(calc_BER_stat_tbl, casn == chemical_casn) %>% select(name) %>% distinct(name)),
						average_of_oral_consumer_products_exposure = signif(avg_mean, digits = 5), minAC50 = min_ac50, minOED = min_oed,
						min_oed_over_oral_consumer_products_exposure = signif(ifelse(is.nan(min_oed) || is.infinite(min_oed),0,min_ac50 / closest_to_95th), digits = 5),
						mean_oed_over_oral_consumer_products_exposure = signif(ifelse(is.nan(avg_oed) || is.infinite(avg_oed),0,avg_ac50 / closest_to_95th), digits = 5));
				
			
		}
		private$tblMeanBER <- mutate(private$tblMeanBER, ac50Rank = rank(private$tblMeanBER$minAC50), oedRank = rank(private$tblMeanBER$minOED));
		
		
		invisible(private$tblMeanBER);
     }   
    },
	
	# plots the box chart comparing the values of the oral exposure vs the AC50 of chemicals
	# October, 2019: Convert horizontal Box Plot to Vertical to accomodate different monitor styles. Plot now sorte by OED.
	plotBERvsAc50 = function( label_by = "casn"){
	  if (self$BERData$calcBERStatsDataExists() && self$basicData$basicStatsDataExists()) {
	    BER_data <- self$BERData$getCalcBERStatsTable();
	    BER_data$oral_ber <- BER_data$direct_ingestion + BER_data$direct_vapor + BER_data$direct_aerosol;
	    BER_data <- filter(BER_data, oral_ber > 0.00000001) %>%  drop_na(oral_ber);
	    
	    basic_data <- self$basicData$getBasicStatsTable();
	    basic_data <- filter(basic_data, above_cutoff == "Y", ac50 >= -2, ac50 <= 10000, cytotoxicity_um > 0.01) %>% drop_na(ac50);
	    
	    # October, 2019: sort the basic_data dataframe by OED and create a new sorted_basic_data dataframe
	    # the newly sorted casn and name columns will be used plot the box plots
	    sorted_basic_data <- basic_data[order(as.double(basic_data$oed), decreasing = TRUE),];
	    sorted_casn <- sorted_basic_data$casn;
	    sorted_name <- sorted_basic_data$name;
	    
	    
	    if(label_by == "casn"){
	      private$pltBERvsAc50 <- plot_ly() %>%
	        #add_trace(data = BER_data, x = ~casn, y = ~signif(oral_ber, digits = 5), type = 'box', name = 'Daily Intake') %>% # horizontal plot
	        #add_trace(data = basic_data, x = ~casn, y = ~signif(oed, digits = 5), type = 'box', name = 'Equivalent Dose') %>% # horizontal plot
	        add_trace(data = sorted_basic_data, x = ~signif(oed, digits = 5), y = ~casn, type = 'box', name = 'Equivalent Dose') %>% # vertical plot
	        add_trace(data = BER_data, x = ~signif(oral_ber, digits = 5), y = ~casn, type = 'box', name = 'Daily Intake') %>% # vertical plot
	        layout(title = 'Consumer Product Exposure vs  In Vitro Tox Equivalent Dose',
	               xaxis = list(title = "Dose - mg/kg/day", type = "log"),
	               yaxis = list(title = "", categoryarray = ~sorted_casn, categoryorder = "array"), # use sorted_casn to plot ordred boxes
	               autosize = T);
	      
	    }else{
	      private$pltBERvsAc50 <- plot_ly() %>%
	        #add_trace(data = BER_data, x = ~name, y = ~signif(oral_ber, digits = 5), type = 'box', name = 'Daily intake') %>% # horizontal plot
	        #add_trace(data = basic_data, x = ~name, y =~signif(oed, digits = 5), type = 'box', name = 'Equivalent Dose') %>% # horizontal plot
	        add_trace(data = sorted_basic_data, x = ~signif(oed, digits = 5), y = ~name, type = 'box', name = 'Equivalent Dose') %>% # vertical plot
	        add_trace(data = BER_data, x = ~signif(oral_ber, digits = 5), y = ~name, type = 'box', name = 'Daily intake') %>% # vertical plot
	        layout(title = 'Consumer Product Exposure vs In Vitro Tox Equivalent Dose',
	               xaxis = list(title = "Dose - mg/kg/day", type = "log"),
	               yaxis = list(title = "", categoryarray = ~sorted_name, categoryorder = "array"), # use sorted_name to plot ordred boxes
	               autosize = T);
	      
	    }
	    
	    invisible(private$pltBERvsAc50);
	  }   
	},

	#plots the box chart comparing the values of the oral exposure vs the AC50 of chemicals
	# October, 2019: Convert horizontal Box Plot to Vertical to accomodate different monitor styles. Plot now sorte by AC50.
	plotOEDvsAc50 = function( label_by = "casn"){
	  if (self$BERData$calcBERStatsDataExists() && self$basicData$basicStatsDataExists()) {
	    basic_data <- self$basicData$getBasicStatsTable();
	    basic_data <- filter(basic_data, above_cutoff == "Y", ac50 >= -2, ac50 <= 10000, cytotoxicity_um > 0.01) %>% drop_na(ac50);
	    
	    # October, 2019: sort the basic_data dataframe by AC50 and create a new sorted_basic_data dataframe
	    # the newly sorted casn and name columns will be used plot the box plots
	    sorted_basic_data <- basic_data[order(as.double(basic_data$ac50), decreasing = TRUE),];
	    sorted_name <- sorted_basic_data$name;
	    sorted_casn <- sorted_basic_data$casn;
	    sorted_ac50 <- sorted_basic_data$ac50;
	    
	    if(label_by == "casn"){
	      private$pltOEDvsAc50 <- plot_ly() %>%
	        #add_trace(data = basic_data, x = ~casn, y = ~signif(10^ac50, digits = 5), type = 'box', name = 'AC50') %>% # horizontal plot
	        #add_trace(data = basic_data, x = ~casn, y = ~signif(oed, digits = 5), type = 'box', name = 'Equivalent Dose') %>% # horizontal plot
	        add_trace(data = sorted_basic_data, x = ~signif(oed, digits = 5), y = ~casn, type = 'box', name = 'Equivalent Dose') %>% # vertical plot
	        add_trace(data = sorted_basic_data, x = ~signif(10^ac50, digits = 5), y = ~casn, type = 'box', name = 'AC50') %>% # vertical plot
	        layout(title = 'In Vitro Tox Concentrations vs Equivalent Doses',
	               xaxis = list(title = "Dose - mg/kg/day", type = "log"),
	               yaxis = list(title = "", categoryarray = ~sorted_casn, categoryorder = "array"), # use sorted_casn to plot ordred boxes
	               autosize = T);
	      
	    }else{
	      private$pltOEDvsAc50 <- plot_ly() %>%
	        #add_trace(data = basic_data, x = ~name, y = ~signif(10^ac50, digits = 5), type = 'box', name = 'AC50') %>% # horizontal plot
	        #add_trace(data = basic_data, x = ~name, y =~signif(oed, digits = 5), type = 'box', name = 'Equivalent Dose') %>% # horizontal plot
	        add_trace(data = sorted_basic_data, x = ~signif(oed, digits = 5), y = ~name, type = 'box', name = 'Equivalent Dose') %>% # vertical plot
	        add_trace(data = sorted_basic_data, x = ~signif(10^ac50, digits = 5), y = ~name, type = 'box', name = 'AC50') %>% # vertical plot
	        layout(title = 'In Vitro Tox Concentrations vs Equivalent Doses',
	               xaxis = list(title = "Dose - mg/kg/day", type = "log"),
	               yaxis = list(title = "", categoryarray = ~sorted_name, categoryorder = "array"), # use sorted_name to plot ordred boxes
	               autosize = T);
	      
	    }
	    
	    invisible(private$pltOEDvsAc50);
	  }   
	},
	
	
	
	#plots the BER box chart
    plotBER = function( label_by = "casn" ){
      if (self$BERData$calcBERStatsDataExists() && self$basicData$basicStatsDataExists()) {
        private$tblBER = tribble(~casn,~name,~product_use_category, ~average_BER);
		#getting our private data
    calc_BER_stat_tbl <- self$BERData$getCalcBERStatsTable();
		chemical_casn_list <- unique(calc_BER_stat_tbl[["casn"]]);
		Ac50_stat_tbl <- self$basicData$getBasicStatsTable();
		Ac50_stat_tbl<- filter(Ac50_stat_tbl, above_cutoff == "Y", ac50 >= -2, ac50 <= 10000, cytotoxicity_um > 0.01) %>% drop_na(ac50);
		
		for(chemical_casn in chemical_casn_list){
			
			  casn_tbl <- filter(calc_BER_stat_tbl, casn == chemical_casn);
			  ac50_tbl <- filter(Ac50_stat_tbl, casn == chemical_casn);
			  avg_mean <- mean(casn_tbl$direct_ingestion + casn_tbl$direct_vapor + casn_tbl$direct_aerosol);
			  oral_ber <- (casn_tbl$direct_ingestion + casn_tbl$direct_vapor + casn_tbl$direct_aerosol);
			  avg_ac50 <- mean(10^ac50_tbl$ac50);
			  min_ac50 <- min(10^ac50_tbl$ac50);
			  
			  avg_oed <- mean(ac50_tbl$oed);
			  min_oed <- min(ac50_tbl$oed);
			  
			  theorical_95th_percentile <- quantile(oral_ber,0.95);
			  closest_to_95th <- absolute.min(oral_ber - theorical_95th_percentile) + theorical_95th_percentile ; #gets the closest value to the 95th percentile.
			  
			  
			  private$tblBER <- add_row(private$tblBER, casn = chemical_casn, 
			                             name = as.character(filter(calc_BER_stat_tbl, casn == chemical_casn) %>% select(name) %>% distinct(name)),
			                             average_BER = signif(ifelse(is.nan(min_oed) || is.infinite(min_oed),0,min_ac50 / closest_to_95th), digits = 5));
			  
			  
			

		}
		
		if(label_by == "casn"){
			private$pltBER <- plot_ly() %>%
			  add_trace(data = private$tblBER, x = ~casn, y = ~average_BER, type = 'box') %>%
			  layout(title = 'BER',
					 xaxis = list(title = ""),
					 yaxis = list(title = "")); #, type = "log"));
	
		}else{
			private$pltBER <- plot_ly() %>%
			  add_trace(data = private$tblBER, x = ~name, y = ~average_BER, type = 'box') %>%
			  layout(title = 'BER',
					 xaxis = list(title = ""),
					 yaxis = list(title = "")); #, type = "log"));
		
		}
		
		invisible(private$pltBER);
     }   
    },

	getMeanBERTable = function (){
		return (private$tblMeanBER);
    },
	getBERvsAc50plot = function(){
		return (private$pltBERvsAc50);
	},
	getOEDvsAc50plot = function(){
	  return (private$pltOEDvsAc50);
	},
	getBERplot = function(){
		return (private$pltBER);
	}
 
	)
)