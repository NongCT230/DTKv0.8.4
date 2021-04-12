# Basic statistical data  --------------------------------------------


#v0.7

#class contains basic statistical data and data calculating functions
# - minac50, minoed, target family counts, min and avg ac50 per target family, scalar top

#the design goal with this is that this is part of a composition where the basicstats data is all managed by another class and all the plotting is done by this class.
#The reason I am doing this is that it allows for easy expension. Since all you need is to include the data. 
Class.Analysis.BasicAnalysis <- R6Class("Class.Analysis.BasicAnalysis",
  #private variables and functions
  private = list(
  
	#basic Analysis Variables
	plttfc = NULL, #target family counts plot,
    pltminac50 = NULL, #target family minimum ac50 plot,
    pltavgac50 = NULL, #target family average ac50 plot,
    pltst_ac50 = NULL, #ac50 vs scalar top scatterplot,
    pltst_oed = NULL, #oed vs scalar top scatterplot,
	pltAc50Box = NULL, #ac50 box chart,
	
	#toxpi Analysis Variables
	target_family_ac_list = NULL,
    target_family_ac_list_tox = NULL,
    target_family_ac_list_group = NULL,
    pltpi = NULL,
    pltpi_group = NULL,
	
	#clustermap Analysis Variables
	
	target_family_matrix = NULL,
    assay_endpoint_matrix = NULL,
    plttf_dendroheatmap = NULL,
    pltassay_dendroheatmap = NULL,
	
	
    
    # helper function for creating dendograms - makes a ggplot
    ggdend = function(df) {
      ggplot() +
        geom_segment(data = df, aes_string(x="x", y="y", xend="xend", yend="yend")) +
        labs(x = "", y = "") + theme_minimal() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank())
    }
	
	
),

#public variables and functions
  public = list(
  
	
	#functions
	
    #constructor
    initialize = function(  ) {
		self$basicData <- Class.Analysis.Data$new(  ); #Done here because it will be shared across all instances of the object otherwise.
    },
    
    #finalizer
    finalize = function() {
    },
	
	basicData = NULL,
	
	
	#Basic Analysis Region
	
	 #plots pie chart of target family counts
	 #preconditions: target_family_counts is calculated
    plotTargetFamilyCounts = function(){
      
      if ( self$basicData$targetFamilyCountsExist()) {
	  
		#getting our private data
        target_family_counts_local <- self$basicData$getTargetFamilyCountsTable() %>% arrange(intended_target_family);
		
        private$plttfc <- plot_ly(target_family_counts_local, labels = ~intended_target_family, values = ~freq,
                                  type = "pie",
                                  hole = 0.5,
                                  pull = 0,
                                  sort = FALSE,
                                  direction = "clockwise",
                                  textposition = 'inside',
                                  textinfo = 'label+percent',
                                  insidetextfont = list(size = 14, color = "#FFFFFF"),
                                  outsidetextfont = list(size = 14, color = "#000000"),
                                  hoverinfo = 'text',
                                  text = ~paste("</br> Target family: ", intended_target_family,
                                                "</br> Count: ", freq),
                                  marker = list(colors = rainbow(length(target_family_counts_local$n), s = 1/2, v = 3/4),
                                                line = list(color = '#FFFFFF', width = 1)),
                                  showlegend = FALSE) %>%
          
          layout(title = "Target Family Proportions",
                 xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                 yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                 margin = list(
                   l = 20,
                   r = 10,
                   t = 30,
                   b = 10
                 ));
        
      } else {
        
        private$plttfc <- NULL;

      }
      
      invisible ( private$plttfc );
      
    },
	#preconditions : basicData is calculated.
	#plots bar chart family targets vs minac50
    plotTargetFamilyMinAc50 = function(){
      
      if ( self$basicData$targetFamilyAcDataExists() && length(self$basicData$getTargetFamilyAcTable()$min_ac50)>0 ) {
        
		target_family_ac_tbl_local <- self$basicData$getTargetFamilyAcTable() %>% arrange(intended_target_family);
		
        private$pltminac50 <- plot_ly(target_family_ac_tbl_local, y =target_family_ac_tbl_local$intended_target_family  , x = target_family_ac_tbl_local$min_ac50, 
                                      type = 'bar', orientation = 'h',
									  color = target_family_ac_tbl_local$intended_target_family,
                                      colors = rainbow(nrow(target_family_ac_tbl_local), s = 1/2, v = 3/4), 
                                      marker = list(line = list(color = 'rgba(0, 0, 0, 1.0)', width = 1)),
									  hoverinfo = "text",
                                      text = ~paste( "</br> Target family: ", intended_target_family,
                                                     "</br> Min ac50: ", min_ac50)) %>% 
          
          layout(title = "Min Ac50 vs Target Family",
                 xaxis = list(title = "Minimum ac50 (uM)",
                              titlefont = list(size = 14)),
                 yaxis = list(title = "Target Family",
                              titlefont = list(size = 14)),
                 showlegend = FALSE,
                 margin = list(
                   l = 200
                 )) %>% 
          
          add_annotations(xref = 'x', yref = 'y',
                          x = target_family_ac_tbl_local$min_ac50 + 0.07,  y = target_family_ac_tbl_local$intended_target_family,
                          text = signif(target_family_ac_tbl_local$min_ac50, 4),
                          font = list(family = 'Arial', size = 12),
                          showarrow = FALSE)
        
      } else {
        
        private$pltminac50 <- NULL;
        
      }
      
      invisible ( private$pltminac50 );
      
    },
    
    #plots bar chart family targets vs avgac50
    plotTargetFamilyAvgAc50 = function(){
      
      if ( self$basicData$targetFamilyAcDataExists()) {
      
		target_family_ac_tbl_local <- self$basicData$getTargetFamilyAcTable() %>% arrange(intended_target_family);
	  

        
        private$pltavgac50 <- plot_ly(target_family_ac_tbl_local, y = target_family_ac_tbl_local$intended_target_family, x = target_family_ac_tbl_local$avg_ac50, 
                                      type = 'bar', orientation = 'h',
                                      color = target_family_ac_tbl_local$intended_target_family,
                                      colors = rainbow(nrow(target_family_ac_tbl_local), s = 1/2, v = 3/4), 
                                      marker = list(line = list(color = 'rgba(0, 0, 0, 1.0)', width = 1)),
                                      hoverinfo = "text",
                                      text = ~paste( "</br> Target family: ", intended_target_family,
                                                     "</br> Avg ac50: ", avg_ac50)) %>% 
          
          layout(title = "Average Ac50 vs Target Family",
                 xaxis = list(title = "Average ac50 (uM)",
                              titlefont = list(size = 14)),
                 yaxis = list(title = "Target Family",
                              titlefont = list(size = 14)),
                 showlegend = FALSE,
                 margin = list(
                   l = 200
                 )) %>% 
          
          add_annotations(xref = 'x', yref = 'y',
                          x = target_family_ac_tbl_local$avg_ac50 + 0.07,  y = target_family_ac_tbl_local$intended_target_family,
                          text = signif(target_family_ac_tbl_local$avg_ac50, 4),
                          font = list(family = 'Arial', size = 12),
                          showarrow = FALSE)
        
      } else {
        
        private$pltavgac50 <- NULL;
        
      }
      
      return ( private$pltavgac50 );
      
    },
	 
    #scatterplot of ac50 vs scalar top
    plotAc50VsScalarTop = function( label_by ){
      if ( self$basicData$scalarTopDataExists() ) {
        
		scalar_top_tbl_local <- self$basicData$getScalarTopTable();
		
        above_cutoff_labels <- sapply(scalar_top_tbl_local$above_cutoff, function(x){
          if(x == "Y"){
            return ("AboveCutoff");
          } else {
            return ("BelowCutoff");
          }
        })
        
        #depending on the number of chemicals we may or may not want to show legend
        if( length( unique(scalar_top_tbl_local$casn) )>20 ){
          sl <- FALSE;
        } else {
          sl <- TRUE;
        }
		

		
        
        #note: ac50 is log10 value, scalar_top is not
        private$pltst_ac50 <- plot_ly(scalar_top_tbl_local, x = ~scalar_top, y = ~ac50, color = scalar_top_tbl_local[[label_by]], alpha = 0.5,
                           colors = rainbow(10, s = 3/4, v=3/4),
                           type = "scatter", mode = "markers",
                           symbol = above_cutoff_labels, symbols = c("circle", "x"),
                           hoverinfo = "text",
						   
                           text = ~paste( "</br> Name: ", name,
                                         "</br> Casn: ", casn,
                                         "</br> Assay endpoint: ", assay_component_endpoint_name,
                                         "</br> Above cutoff: ", above_cutoff,
                                         "</br> Ac50 (uM): ", signif(10^ac50, digits=5),
                                         "</br> Scalar top: ", signif(scalar_top, digits=5))  ) %>% 
                          
          layout(title = "Log Ac50 vs Log scalar top",
                 xaxis = list(title = "Log scalar top",
                              titlefont = list(size = 14),
                              type = "log"),
                 yaxis = list(title = "Log ac50",
                              titlefont = list(size = 14)),
                 showlegend = sl,
                 legend = list(orientation = "h")
				 );
        
      } else {
        
        private$pltst_ac50 <- NULL;
        
      }
      
      return ( private$pltst_ac50 );
    },
    
	  #scatterplot of oed vs scalar top
    plotOEDVsScalarTop = function( label_by ){
      
      if ( self$basicData$scalarTopDataExists()) {
        
		scalar_top_tbl_local <- self$basicData$getScalarTopTable();
		
        above_cutoff_labels <- sapply(scalar_top_tbl_local$above_cutoff, function(x){
          if(x == "Y"){
            return ("AboveCutoff");
          } else {
            return ("BelowCutoff");
          }
        })
        
        #depending on the number of chemicals we may or may not want to show legend
        if( length( unique(scalar_top_tbl_local$casn) )>20 ){
          sl <- FALSE;
        } else {
          sl <- TRUE;
        }
        
        private$pltst_oed <- plot_ly(scalar_top_tbl_local, x = ~scalar_top, y = ~oed, color = scalar_top_tbl_local[[label_by]], alpha = 0.5,
                                 colors = rainbow(10, s = 3/4, v=3/4),
                                 type = "scatter", mode = "markers",
                                 symbol = above_cutoff_labels, symbols = c("circle", "x"),
                                 hoverinfo = "text",
                                 text = ~paste( "</br> Name: ", name,
                                                "</br> Casn: ", casn,
                                                "</br> Assay endpoint: ", assay_component_endpoint_name,
                                                "</br> Above cutoff: ", above_cutoff,
                                                "</br> OED (mg/kg/day): ", signif(oed, digits=5),
                                                "</br> Scalar top: ", signif(scalar_top, digits=5))  ) %>% 
          
          layout(title = "OED vs Log scalar top",
                 xaxis = list(title = "Log scalar top",
                              titlefont = list(size = 14),
                              type = "log"),
                 yaxis = list(title = "OED (mg/kg/day)",
                              titlefont = list(size = 14)),
                 showlegend = sl,
                 legend = list(orientation = "h")   );
        
      } else {
        
        private$pltst_oed <- NULL;
        
      }
      
      return ( private$pltst_oed );
    },
	#box plot of the Ac50 values by chemical.
	plotAc50Box = function( label_by = "casn"){
      if (self$basicData$basicStatsDataExists()) {
	
		basic_data_burst = filter(self$basicData$getBasicStatsTable(), burst_assay == 1);
		basic_data_burst_no = filter(self$basicData$getBasicStatsTable(), burst_assay == 0);
		if(label_by == "casn"){
			private$pltAc50Box <- plot_ly() %>%
			  add_trace(data = basic_data_burst, x = ~casn, y = ~signif(ac50, digits = 5), type = 'box', name = 'Burst Assay') %>%
			  add_trace(data = basic_data_burst_no, x = ~casn, y = ~signif(ac50, digits = 5), type = 'box', name = 'Not Burst Assay') %>%
			  layout(title = 'Burst Assay Analysis Ac50 Box Chart',
					 xaxis = list(title = ""),
					 yaxis = list(title = ""));
	
		}else{
			private$pltAc50Box <- plot_ly() %>%
				add_trace(data = basic_data_burst, x = ~name, y = ~signif(ac50, digits = 5), type = 'box', name = 'Burst Assay') %>%
				add_trace(data = basic_data_burst_no, x = ~name, y = ~signif(ac50, digits = 5), type = 'box', name = 'Not Burst Assay') %>%
			  layout(title = 'Burst Assay Analysis Ac50 Box Chart',
					 xaxis = list(title = ""),
					 yaxis = list(title = ""));
		
		}
		
		invisible(private$pltAc50Box);
     }   
    },
	
	getPlotAc50Box = function(){
		return (private$pltAc50Box);
	},

	
	getTargetFamilyCountsPlot = function (){
      return (private$plttfc);
    },
    
    getTargetFamilyMinAcPlot = function (){
      return (private$pltminac50);
    },
    
    getTargetFamilyAvgAcPlot = function (){
      return (private$pltavgac50);
    },
    
    getAc50ScalarTopPlot = function (){
      return (private$pltst_ac50);
    },
    
    getOEDScalarTopPlot = function (){
      return (private$pltst_oed);
    },
    
	#Basic Analysis Region END
	
	#Tox Pi Region Start
	
	 #toxPi tables for plotting
	 computeToxPI = function(){
      
      private$target_family_ac_list <- list();
      
      if (!self$basicData$basicStatsDataExists()){
        logwarn("Analysis.BasicAnalysis$computeToxPI():BasicStatsTable stats table does not exist, returning empty list. Make sure to self$basicData$buildStatsTable() first.");
        return(private$target_family_ac_list);
      }
      
	  basic_table <- self$basicData$getBasicStatsTable();
      chemical_casn_list <- unique(basic_table[["casn"]]);
      for (chemical_casn in chemical_casn_list){
        
        target_family_count_tbl <- filter(basic_table, casn == chemical_casn, above_cutoff == "Y") %>% 
          select(casn, name, intended_target_family) %>%
          drop_na(intended_target_family);
        target_family_counts <- plyr::count(target_family_count_tbl, vars = "intended_target_family"); #colnames: intended_target_family, n
        target_family_ac_tbl <- tribble(~intended_target_family, ~avg_ac50);

		tmp = drop_na(target_family_counts);
		target_family_counts = tmp;
		
        for ( target_family in target_family_counts$intended_target_family ){
          
          df <- filter(basic_table, casn == chemical_casn, intended_target_family == target_family, above_cutoff == "Y") %>% 
            drop_na(ac50);         
		 
          df$ac50 <- 10^(df$ac50); #ac50 is originally in log10
          avg_ac50 <- mean( df$ac50, na.rm=TRUE);
          
          target_family_ac_tbl <- add_row(target_family_ac_tbl, 
                                          intended_target_family = target_family, avg_ac50 = avg_ac50);
        }
		tmp = drop_na(target_family_ac_tbl);
		target_family_ac_tbl = tmp;
        #transform data to emphasize smallest values and normalize to log scale from 0 to 1
        #
        #reciprocal of the ac50 is used to emphasize smallest ac50 in the plot
        re_avg_ac50 <- 1 / target_family_ac_tbl$avg_ac50;
        #normalize to range from 1 to 10
        min_val <- min(re_avg_ac50);
        max_val <- max(re_avg_ac50);
        range_val <- max_val - min_val;
        slope <- 9/range_val;
        re_avg_ac50 <- ( re_avg_ac50 - min_val ) * slope + 1;
        #log scale to ensure extreme values do not drown out the regular ones
        re_avg_ac50 <- log10( re_avg_ac50 );
        #add rearranged ac50 into the table
        target_family_ac_tbl <- add_column(target_family_ac_tbl, re_avg_ac50 = re_avg_ac50);
        
        #add the taget family counts and ac50 table to the list, keyed by target family; then sort by re_avg_ac50		
				
		if(nrow(target_family_ac_tbl)> 0 && nrow(target_family_counts)> 0){
			private$target_family_ac_list[[chemical_casn]] <- full_join(target_family_counts, target_family_ac_tbl, 
                                                                    by = "intended_target_family") %>% 
			arrange(intended_target_family);
			
		
		}else{
			chemical = select(basic_table, casn, name) %>% filter(chemical_casn == casn) %>% distinct(casn,name);
			self$basicData$addHitlessChemInfo(chemical);
		}
      }
      return ( private$target_family_ac_list );
    },
    
    
    #can we have a better name for this Gabriel
    #toxPi tables for plotting
    computeToxPI2 = function(){
      
      private$target_family_ac_list_tox <- list();
      
      if (!self$basicData$basicStatsDataExists()){
        logwarn("Analysis.BasicAnalysis$computeToxPI2():ToxPI stats table does not exist, returning empty list. Make sure to $buildStatsTable() first.");
        return(private$target_family_ac_list_tox);
      }
     
	  basic_table <- self$basicData$getBasicStatsTable();
      chemical_casn_list <- unique(basic_table[["casn"]]);
      for (chemical_casn in chemical_casn_list){
           
        dfa <- filter(basic_table, casn == chemical_casn,above_cutoff == "Y", ac50 >= -2, ac50 <= 1000, cytotoxicity_um > 0.01) %>% drop_na(ac50);
          
        dfa$ac50 <- 10^(dfa$ac50); #ac50 is originally in log10
        avg_ac50 <- mean(dfa$ac50, na.rm = TRUE)
        avg_cyto <- min(dfa$cytotoxicity_um, na.rm=TRUE);
        avg_top <- mean(dfa$ac_top / dfa$ac50, na.rm = TRUE);
        avg_top[avg_top > 1000] <- 1000;
		
	
        #range_tp <- diff(range(dfa$ac_top / dfa$ac50));
		#this is fine because it's linear scaling not brute forcing values to be positive.
		val_ac = (log10(avg_ac50) + 2) * 0.8;
		val_ct = (log10(avg_cyto) + 2) * 0.8;
		val_tp = (log10(avg_top) + 2) * 0.8;
        #val_tp = min(max(-log10(avg_top/range_tp),0),4);
		#val_tp = -log10(avg_top/range_tp);
        #add the taget family counts and ac50 table to the list, keyed by target family; then sort by re_avg_ac50
		if(nrow(dfa)> 0){
		
			private$target_family_ac_list_tox[[chemical_casn]] <- tibble(val_ac, val_ct, val_tp);
		
		}else{
				chemical = select(basic_table, casn, name) %>% filter(chemical_casn == casn) %>% distinct(casn,name);
				self$basicData$addHitlessChemInfo(chemical);
			}
	  }
      return ( private$target_family_ac_list_tox);
    },
    
    #toxPi tables for plotting
    computeToxPIGroup = function(){
      
      private$target_family_ac_list_group <- list();
      
      if (!self$basicData$basicStatsDataExists()){
        logwarn("Analysis.BasicAnalysis$computeToxPIGroup():ToxPI stats table does not exist, returning empty list. Make sure to $buildStatsTable() first.");
        return(private$target_family_ac_list_group);
      }
	  
	  basic_table <- self$basicData$getBasicStatsTable();
	  
      #create a common count table
      all_chem_tfac_table <- select(basic_table, casn, name, intended_target_family, ac50, above_cutoff) %>% filter(above_cutoff == "Y");
      all_chem_tfac_table$ac50 <- 10^(all_chem_tfac_table$ac50);
      all_chem_tfac_table <- group_by(all_chem_tfac_table, casn, intended_target_family) %>% 
        dplyr::summarize(avg_ac50 = mean(ac50, na.rm=TRUE)) %>% ungroup() %>% 
        distinct(casn, name, intended_target_family, avg_ac50) # %>% drop_na(all_chem_tfac_table,intended_target_family);
	  tmp = drop_na(all_chem_tfac_table);
	  all_chem_tfac_table = tmp;

	  
	 
	  
	  
	  
      #transform data to emphasize smallest values and normalize to log scale from 0 to 1
      #
      #reciprocal of the ac50 is used to emphasize smallest ac50 in the plot
      re_avg_ac50 <- 1 / all_chem_tfac_table$avg_ac50;
      #normalize to range from 1 to 10
      min_val <- min(re_avg_ac50);
      max_val <- max(re_avg_ac50);
      range_val <- max_val - min_val;
      slope <- 9/range_val;
      re_avg_ac50 <- ( re_avg_ac50 - min_val ) * slope + 1;
      #log scale to ensure extreme values do not drown out the regular ones
      re_avg_ac50 <- log10( re_avg_ac50 );
      #add rearranged ac50 into the table
      all_chem_tfac_table <- add_column(all_chem_tfac_table, re_avg_ac50 = re_avg_ac50);
      
      
      #create a counts table, which should have all equal counts to space pie slices out evenly
      all_chem_family_counts <- select(all_chem_tfac_table, intended_target_family) %>% distinct(intended_target_family);
      all_chem_family_counts <- add_column(all_chem_family_counts, n = 1);

	  
	  tmp = drop_na(all_chem_family_counts);
	  all_chem_family_counts = tmp;
      
      chemical_casn_list <- unique(basic_table[["casn"]]);
      for (chemical_casn in chemical_casn_list){
        target_family_ac_tbl <- filter(all_chem_tfac_table, casn == chemical_casn) %>% 
          select(intended_target_family, avg_ac50, re_avg_ac50);
		  
        if(nrow(target_family_ac_tbl)> 0 && nrow(all_chem_family_counts)>0){
        #add the taget family counts and ac50 table to the list, keyed by target family; then fill NA values with 0s
			private$target_family_ac_list_group[[chemical_casn]] <- full_join(all_chem_family_counts, target_family_ac_tbl, 
																			  by = "intended_target_family") %>% 
			  replace_na( list(avg_ac50 = 0, re_avg_ac50 = 0) ) %>% arrange(intended_target_family);
		}else{
			chemical = select(basic_table, casn, name) %>% filter(chemical_casn == casn) %>% distinct(casn,name);
			self$basicData$addHitlessChemInfo(chemical);
		}
      }

      return ( private$target_family_ac_list_group );
    },
    
    #plots pie chart of target family counts
    plotToxPI = function( label_by = "name", grouped = FALSE){
      if ( (!grouped & self$toxPIDataExists()) | ( grouped & self$toxPIGroupDataExists()) ) {
        #plotly does not have a suitable polar plotting function, so slower and non-interactive ggplot2 is used
        #toxpi representation of data uses a box plot that is converted to polar coordinates
        
        #plot theme
        toxpi_theme <- theme_minimal()+
          theme(
            legend.title=element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            title = element_text(size=16),
            axis.ticks.x=element_line(size=5),
            axis.text.x = element_blank(),
            legend.position = "bottom",
            legend.text=element_text(size=14),
            axis.text.y = element_blank(),
            plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
          )
        
        #create subplots
        subplots <- list();
        
        if ( grouped ){
          all_casn <- names(private$target_family_ac_list_group);
        } else {
          all_casn <- names(private$target_family_ac_list);
        }
        
        for ( chem_casn in all_casn ){
          
          if (label_by == "name"){
            chem_label <- filter(self$basicData$getBasicStatsTable(), casn == chem_casn)$name[[1]];
          } else {
            chem_label <- chem_casn;
          }
          
          if ( grouped ){
            plotdata <- private$target_family_ac_list_group[[chem_casn]];
          } else {
            plotdata <- private$target_family_ac_list[[chem_casn]];
          }

          #first, need to create the boxes, with left and right coordinates for each
          left <- 0;
          right <- ifelse(grouped, plotdata$n[[1]], plotdata$freq[[1]]);
          center <- (left+right)/2;
		  if ((grouped && length(plotdata$n) >= 2) || (!grouped && length(plotdata$freq) >= 2)){ #checking if there's more than 1 intended_target_family
			  for ( i in 2:ifelse(grouped, length(plotdata$n), length(plotdata$freq))) {
				left <- c(left, (left[[(i-1)]]+ifelse(grouped, plotdata$n[[(i-1)]], plotdata$freq[[(i-1)]])));
				right <- c(right, (right[[(i-1)]]+ifelse(grouped, plotdata$n[[i]], plotdata$freq[[i]])));
				center <- c(center, ((left[[i]]+right[[i]])/2));
			  }
          }
          #also key some colors for boxes and corresponding labels
          colorpalette <- rainbow(length(plotdata$intended_target_family), end = 0.7, s = 3/4, v = 3/4);
          
          #add box positions and colors to the plot
          #since data is sorted, red rainbow color will fall onto the highest re_avg_ac50 value
          plotdata <- add_column(plotdata, left=left, right=right, center = center, colour = colorpalette)
          
          #save data for printing plots later
          if ( grouped ){
            private$target_family_ac_list_group[[chem_casn]] <- plotdata;
          } else {
            private$target_family_ac_list[[chem_casn]] <- plotdata;
          }
          
          #create plot -> secondary geom_rect with show.legend=TRUE is a fake one to force legend mapping with aes() and scale_fill_identity
          pltpi <- ggplot(data=plotdata)+
            geom_rect(aes_string(xmin = "left", xmax = "right", ymin = 0, ymax = "re_avg_ac50"), 
                      size = 0.5, color = "black", fill = plotdata$colour)+
            geom_rect(aes_string(xmin = "left", xmax = "right", ymin = 0, ymax = 0, fill = "colour"), 
                      size = 0.5, color = "black", show.legend = TRUE)+
            scale_x_continuous(breaks=plotdata$center,
                               labels=plotdata$intended_target_family)+
            geom_text_repel(aes_string(x = "center", y = "re_avg_ac50",
                                       label = "intended_target_family", size = 7),
                            nudge_y = max(plotdata$re_avg_ac50)-(plotdata$re_avg_ac50)+1,
                            min.segment.length = 0.5,
                            show.legend = FALSE,
                            color = plotdata$colour) +
            geom_text_repel(aes_string(x = "center", y = "re_avg_ac50",
                                       label = signif(plotdata$re_avg_ac50,3), size = 7),
                            nudge_y = max(plotdata$re_avg_ac50)-(plotdata$re_avg_ac50)+0.1, #was 0.1 for just above the slice
                            min.segment.length = 30,
                            show.legend = FALSE,
                            color = plotdata$colour)+
            toxpi_theme+
            ylim(0, max(plotdata$re_avg_ac50)+0.7)+
            scale_fill_identity(guide=guide_legend(override.aes = list(fill = plotdata$colour)),
                                labels = plotdata$intended_target_family)+
            guides(size = FALSE)+
            coord_polar(theta="x", start=0)+
            ggtitle(chem_label);
          
          #add to subplot
          subplots[[chem_casn]] <- pltpi;
        }
        
        #save all subplots
        if ( grouped ){
          private$pltpi_group <- subplots;
        } else {
          private$pltpi <- subplots;
        }
        
      } else {
        logwarn("Analysis.ToxPI$plotToxPI():ToxPI plot data does not exist. No plots produced. Make sure to $computeToxPI() or $computeToxPIGroup() first.")
        if ( grouped ){
          private$pltpi_group <- NULL;
        } else {
          private$pltpi <- NULL;
        }
        
      }
      if ( grouped ){
        return (private$pltpi_group);
      } else {
        return (private$pltpi);
      }
      
    },
    #Gabriel Better name?
    plotToxPI2 = function( label_by = "name"){
      
	#Gabriel Add Checks for the data we are pulling. 
	  
     #plot theme
        toxpi_theme <- theme_minimal()+
          theme(
            legend.title=element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            title = element_text(size=16),
            axis.ticks.x=element_line(size=5),
            axis.text.x = element_blank(),
            legend.position = "bottom",
            legend.text=element_text(size=14),
            axis.text.y = element_blank(),
            plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
          )

        
        #create subplots
        subplots <- list();
        
        all_casn <- names(private$target_family_ac_list_tox);
        
        for ( chem_casn in all_casn ){
          
          if (label_by == "name"){
            chem_label <- filter(self$basicData$getBasicStatsTable(), casn == chem_casn)$name[[1]];
          } else {
            chem_label <- chem_casn;
          }
		
		   plotdata <- private$target_family_ac_list_tox[[chem_casn]];
		   
           #  create mean of each row

          
          

		   
		   #val_ac <- min(max( -log10(1/plotdata$avg_ac50/plotdata$range_ac*10),0),4);
           #val_ct <- min(max(-log10(1/plotdata$avg_cyto/plotdata$range_ac*10),0),4);
           #val_tp <- min(max(-log10(1/plotdata$avg_top/plotdata$range_tp*10),0),4);
		   
		 #NAN somewhere here.
		  #DUMMY CHECK TO CHANGE WHEN WE FIGURE OUT THE CAUSE OF THE NAN VALUES AND THE SOLUTION. THIS CAUSES THE TOXPIES TO NOT CRASH. GABRIEL 
		   plotdata$val_ac <- if(is.nan(plotdata$val_ac)) 0 else plotdata$val_ac;
		   plotdata$val_ct <- if(is.nan(plotdata$val_ct)) 0 else plotdata$val_ct;
		   plotdata$val_tp <- if(is.nan(plotdata$val_tp)) 0 else plotdata$val_tp;
           plot_label <- c('Ac50', 'Cytotoxity', 'Top Response');
           plot_values <- c(plotdata$val_ac, plotdata$val_ct, plotdata$val_tp);
           colorpalette <- rainbow(3);
           #gabe
			#the limit is set at 4 apparently...

           dfb <- data.frame(plot_label, plot_values);
           pltpi <-  ggplot(dfb, aes(x=plot_label, y=plot_values, fill=plot_label)) +
             geom_col(width = 1, colour = colorpalette) +
             coord_polar(theta="x", start=0) +
             toxpi_theme +
             ylim(0,4) +
             ggtitle(chem_label, subtitle = chem_casn);
        
        #add to subplot
        subplots[[chem_casn]] <- pltpi;
      }
        
        private$pltpi <- subplots;
        return (private$pltpi);
        
    },
	
	    getToxPIPlotData = function (){
      return (private$target_family_ac_list);
    },
    
    getToxPI2PlotData = function (){
      return (private$target_family_ac_list_tox);
    },
    
    
    getToxPIGroupPlotData = function (){
      return (private$target_family_ac_list_group);
    },
    
    getToxPIPlots = function (){
      return (private$pltpi);
    },
    
    getToxPI2Plots = function (){
      return (private$pltpi);
    },
    
    getToxPIGroupPlots = function (){
      return (private$pltpi_group);
    },
    
    toxPIDataExists = function (){
      if (is.null(private$target_family_ac_list) || length(private$target_family_ac_list) == 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
    
    toxPI2DataExists = function (){
      if (is.null(private$target_family_ac_list_tox) || length(private$target_family_ac_list_tox) == 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
    
    toxPIGroupDataExists = function (){
      if (is.null(private$target_family_ac_list_group) || length(private$target_family_ac_list_group) == 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
	
	#TOXPI region END
	
	#CLUSTER Region START
	 #target family matrix for plotting
    computeTargetFamilyMatrix = function(){
      
      loginfo ("Clusterheatmap computing TargetFamilyMatrix");
      
      #CHECKPOINT: empty stat table, presently (Nov 2017) dyplyr::spread crashes on empty df
      if (!self$basicData$basicStatsDataExists()) {
        
        private$target_family_matrix <- matrix(nrow = 0, ncol = 0); #will return empty matrix later
        
      } else {
      
        data3 <- unite(self$basicData$getBasicStatsTable(), intended_target_family, intended_target_family, intended_target_family_sub, sep = " > ") %>% 
          select(casn, ac50, intended_target_family);
        data3$ac50 <- log10( 1/10^(data3$ac50) );
        data3 <- group_by(data3, casn, intended_target_family) %>% dplyr::summarise(avgac50 = mean(ac50, na.rm=TRUE)) %>% ungroup()%>%
          distinct(casn, intended_target_family, avgac50) %>% 
          spread(key = intended_target_family, value = avgac50, fill = -100); #assign -100 to missing values, this is on a LOG10 scale, where no data reaches that far
        rownames <- data3$casn;
        data3 <- select(data3, -casn) %>% as.matrix();
        row.names(data3) <- rownames;
        private$target_family_matrix <- data3;

      }
      
      invisible ( private$target_family_matrix );
    },
    
	  #assay endpoint matrix for plotting
    computeAssayEndpointMatrix = function(){
      
      loginfo ("Clusterheatmap computing AssayEndpointMatrix");
      
      #CHECKPOINT: empty stat table, presently (Nov 2017) dyplyr::spread crashes on empty df
      if (!self$basicData$basicStatsDataExists()) {
        
        private$assay_endpoint_matrix <- matrix(nrow = 0, ncol = 0); #will return empty matrix later
        
      } else {
      
        data3 <- select(self$basicData$getBasicStatsTable(), casn, ac50, assay_component_endpoint_name);
        data3$ac50 <- log10( 1/10^(data3$ac50) );
        data3 <-  spread(data3, key = assay_component_endpoint_name, value = ac50, fill = -100); #assign -100 to missing values, this is on a LOG10 scale, where no data reaches that far
        rownames <- data3$casn;
        data3 <- select(data3, -casn) %>% as.matrix();
        row.names(data3) <- rownames;
        private$assay_endpoint_matrix <- data3;
      
      }
      
      invisible ( private$assay_endpoint_matrix );
    },
    
    #plots dendro-heatmap of target families/subfamilies
    plotTargetFamilyDendroHeatmap = function( label_by = "name"){
      
      if ( self$targetFamilyDataExists() ) {
      
        x <- private$target_family_matrix;
        #labeling
        rownames <- rownames(x);
        labeldf <- filter(self$basicData$getBasicStatsTable(), casn %in% rownames) %>% select(casn, name) %>% distinct(casn, name);

        orderedlabeldf <- data.frame(casn = as.character(rownames)) %>% inner_join(labeldf, by = "casn");
        if (label_by == "name"){
          rownames <- orderedlabeldf$name;
          rownames(x) <- rownames;
        }
        
        #column clustering
        if(nrow(x) > 1){
          dendro.col <- hclust(dist(x));
          #column dendrograms
          dd.col <- reorder(as.dendrogram(dendro.col), rowMeans(x, na.rm = TRUE)) #convert clustering to dendrogram
          dy <- dendro_data(dd.col) #ggdendro function to extract plotting info from dendrogram
          #leaf ordering for plot labels
          col.ord <- order.dendrogram(dd.col);
          # y dendrogram - ggplot
          py <- private$ggdend(dy$segments) + coord_flip();
        } else { #trivial case which crashes hclust(), essentially a single observation that can be considered a single cluster
          col.ord <- 1L;
          py <- ggplot(data.frame(x1 = 0, x2 = 0)) + geom_blank();
        }
          

        #row clustering
        if(ncol(x) > 1){
          dendro.row <- hclust(dist(t(x)));
          #row dendrograms
          dd.row <- reorder(as.dendrogram(dendro.row), colMeans(x, na.rm = TRUE))
          dx <- dendro_data(dd.row) #ggdendro function to extract plotting info from dendrogram
          #leaf ordering for plot labels
          row.ord <- order.dendrogram(dd.row)
          # x dendograms - ggplot
          px <- private$ggdend(dx$segments);
        } else {
          row.ord <- 1L;
          px <- ggplot(data.frame(x1 = 0, x2 = 0)) + geom_blank();
        }
        if(ncol(x) > 2){
        # heatmap
        #rearrange original data by label order and make dataframe
        df <- as.data.frame( x )[col.ord, row.ord];
        df$casn <- row.names(df);
        df$casn <- with(df, factor(casn, levels=casn, ordered=TRUE))
        mdf <- reshape2::melt(df, id.vars="casn") #reshape by casn; creates {casn, variable, value} table
        colnames(mdf) <- c("chem", "target", "activity");
        lower_limit <- mdf %>% filter(activity != -100);
        lower_limit <- min(lower_limit$activity);
        p <- ggplot(mdf, aes_string(x = "target", y = "chem")) + geom_tile(aes_string(fill = "activity"), color = "white") +
          scale_fill_gradient(limits = c(lower_limit, NA),low = "#330000", high = "red",
                              na.value = "green");
        
        p_empty <- plot_ly(type = "scatter", mode = "markers") %>%
          # note that margin applies to entire plot, so we can
          # add it here to make tick labels more readable
          layout(margin = list(l = 200),
                 xaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
                 yaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
                 plot_bgcolor="white",
                 paper_bgcolor="white");
        
        pxly <- ggplotly(px, tooltip = c("y")) %>%
          layout(plot_bgcolor="white",
                 paper_bgcolor="white");
        pyly <- ggplotly(py, tooltip = c("y")) %>%
          layout(plot_bgcolor="white",
                 paper_bgcolor="white");
        ply <- ggplotly(p, tooltip = c("chem", "target", "activity")) %>% 
          layout(xaxis = list(title = "Target Family > Subfamily", tickangle = 90,
                              tickfont = list(size=12)),
                 yaxis = list(title = "Chemical", tickfont = list(size=14)),
                 plot_bgcolor="white",
                 paper_bgcolor="white");

        private$plttf_dendroheatmap <- subplot(pxly,
                                               p_empty,
                                               ply,
                                               pyly,
                                               nrows = 2, margin = 0,
                widths = c(0.8, 0.2), heights = c(0.3, 0.7), shareX = TRUE, shareY = TRUE) %>%
          layout( margin=list(b = 270, l = 200), hovermode = "compare" );

		}
      } else {
        logwarn("Analysis.ClusterHeatmap$plotTargetFamilyDendroHeatmap():Plot data does not exist. No plots produced. Make sure to $computeTargetFamilyMatrix() first.")
        private$plttf_dendroheatmap <- NULL;

      }
      
      invisible ( private$plttf_dendroheatmap );
      
    },
    
    #plots dendro-heatmap of target families/subfamilies
    plotAssayEndpointDendroHeatmap = function( label_by = "name"){
      
      if ( self$assayEndpointDataExists() ) {
        
        x <- private$assay_endpoint_matrix;
        
        #labeling
        rownames <- rownames(x);
        labeldf <- filter(self$basicData$getBasicStatsTable(), casn %in% rownames) %>% select(casn, name) %>% distinct(casn, name);
        orderedlabeldf <- data.frame(casn = as.character(rownames)) %>% inner_join(labeldf, by = "casn");
        if (label_by == "name"){
          rownames <- orderedlabeldf$name;
          rownames(x) <- rownames;
        }
        
        #column clustering
        if(nrow(x) > 1){
          dendro.col <- hclust(dist(x));
          #column dendrograms
          dd.col <- reorder(as.dendrogram(dendro.col), rowMeans(x, na.rm = TRUE)) #convert clustering to dendrogram
          dy <- dendro_data(dd.col) #ggdendro function to extract plotting info from dendrogram
          #leaf ordering for plot labels
          col.ord <- order.dendrogram(dd.col);
          # y dendrogram - ggplot
          py <- private$ggdend(dy$segments) + coord_flip();
        } else {
          col.ord <- 1L;
          py <- ggplot(data.frame(x1 = 0, x2 = 0)) + geom_blank();
        }
        
        
        #row clustering
        if(ncol(x) > 1){
          dendro.row <- hclust(dist(t(x)));
          #row dendrograms
          dd.row <- reorder(as.dendrogram(dendro.row), colMeans(x, na.rm = TRUE))
          dx <- dendro_data(dd.row) #ggdendro function to extract plotting info from dendrogram
          #leaf ordering for plot labels
          row.ord <- order.dendrogram(dd.row)
          # x dendograms - ggplot
          px <- private$ggdend(dx$segments);
        } else {
          row.ord <- 1L;
          px <- ggplot(data.frame(x1 = 0, x2 = 0)) + geom_blank();
        }
        if(ncol(x) > 2){
        # heatmap
        #rearrange original data by label order and make dataframe
        df <- as.data.frame( x )[col.ord, row.ord];
        df$casn <- row.names(df);
        df$casn <- with(df, factor(casn, levels=casn, ordered=TRUE))
        mdf <- reshape2::melt(df, id.vars="casn") #reshape by casn; creates {casn, variable, value} table
        colnames(mdf) <- c("chem", "assay_endpoint", "activity");
        lower_limit <- mdf %>% filter(activity != -100);
        lower_limit <- min(lower_limit$activity);
        p <- ggplot(mdf, aes_string(x = "assay_endpoint", y = "chem")) + geom_tile(aes_string(fill = "activity"), color = "white") +
          scale_fill_gradient(limits = c(lower_limit, NA),low = "#330000", high = "red",
                              na.value = "green");
        
        p_empty <- plot_ly(type = "scatter", mode = "markers") %>%
          # note that margin applies to entire plot, so we can
          # add it here to make tick labels more readable
          layout(margin = list(l = 200),
                 xaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
                 yaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
                 plot_bgcolor="white",
                 paper_bgcolor="white");
        
        pxly <- ggplotly(px, tooltip = c("y")) %>%
          layout(plot_bgcolor="white",
                 paper_bgcolor="white");
        pyly <- ggplotly(py, tooltip = c("y")) %>%
          layout(plot_bgcolor="white",
                 paper_bgcolor="white");
        ply <- ggplotly(p, tooltip = c("chem", "assay_endpoint", "activity")) %>% 
          layout(xaxis = list(title = "Assay Endpoints", tickangle = 90,
                              tickfont = list(size=12)),
                 yaxis = list(title = "Chemical", tickfont = list(size=14)),
                 plot_bgcolor="white",
                 paper_bgcolor="white");
        
        private$pltassay_dendroheatmap <- subplot(pxly,
                                               p_empty,
                                               ply,
                                               pyly,
                                               nrows = 2, margin = 0,
                                               widths = c(0.8, 0.2), heights = c(0.3, 0.7), shareX = TRUE, shareY = TRUE) %>%
          layout( margin=list(b = 200, l = 200), hovermode = "compare" );
       } 
      } else {
        logwarn("Analysis.ClusterHeatmap$plotAssayEndpointDendroHeatmap():Plot data does not exist. No plots produced. Make sure to $computeAssayEndpointMatrix() first.")
        private$pltassay_dendroheatmap <- NULL;
        
      }
      
      invisible ( private$pltassay_dendroheatmap );
      
    },
	
	  getTargetFamilyMatrix = function (){
      return (private$target_family_matrix);
    },
    
    getAssayEndpointMatrix = function (){
      return (private$assay_endpoint_matrix);
    },
    
    getTargetFamilyPlot = function (){
      return (private$plttf_dendroheatmap);
    },
    
    getAssayEndpointPlot = function (){
      return (private$pltassay_dendroheatmap);
    },
    
    targetFamilyDataExists = function (){
      if (is.null(private$target_family_matrix) || 
          nrow(private$target_family_matrix) == 0 || ncol(private$target_family_matrix) == 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
    
    assayEndpointDataExists = function (){
      if (is.null(private$assay_endpoint_matrix) || 
          nrow(private$assay_endpoint_matrix) == 0 || ncol(private$assay_endpoint_matrix) == 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    }
    
 )
)