# SERVER SUPPORT FUNCTIONS

#v0.8.4

# Version number and revision
app.version.major <- "0";
app.version.minor <- "8";
app.version.revision <- "4";

#num columns 
num_plot_cols <- 1;

# Model-related functions --------------------------------------------------

#sets up the required models with the supplied database and importer objects
setupModels = function ( input, output, session,
                         dtk.db, dtk.import,
                         model.3css, dtk.import.3css,
                         model.fubPoly, dtk.import.fubPoly,
                         model.oed1){
  
  #set up importer objects
  dtk.import$setSource(dtk.db);
  dtk.import.3css$setSource(dtk.db);
  dtk.import.fubPoly$setSource(dtk.db);
  #what is this line? I guess it does nothing
  #dtk.import$loadPhysiologyData
  
  #set up model objects
  model.fubPoly$setImporters(dtk.import, dtk.import.fubPoly);
  model.fubPoly$importPhysiologyData();
  
  model.3css$setImporters(dtk.import, dtk.import.3css);
  model.3css$attachModels( list(model.fubPoly) );
  model.3css$setParameters ( sig.threshold = 0.1); #as per [Pearce et al 2017] httk-1.7 paper
  model.3css$importPhysiologyData();
  
  #model.3css$importModelData();
  
  model.oed1$setParameters( cssmodel = model.3css$getModelName() ); #OED model will use Css from 3css model calculations
  
  
}

#runs the required models
runModels = function( chemlist,
                      input, output, session,
                      model.3css, model.oed1,
                      progress=NULL, n=NULL){
  
  #
  if(!is.null(progress) && !is.null(n)) progress$inc(1/n, detail = "Computing Css models");
	model.3css$setChemicals(chemlist);
	model.3css$importModelData();
	model.3css$compute();
	#
	if(!is.null(progress) && !is.null(n)) progress$inc(1/n, detail = "Computing OED models");
	model.oed1$setChemicals(chemlist);
	model.oed1$compute();
	#
	if(!is.null(progress) && !is.null(n)) progress$inc(1/n, detail = "Computing basic stats");
	for (chem in chemlist){
		chem$calculateMinAc50(input$select_background, input$select_hit);
		chem$calculateMinOED(model.oed1$getModelName() );
	  #chem$calculateMeanOED(model.oed1$getModelName() );
	}


}
# Analysis-related functions ----------------------------------------------------
computeAnalyses = function ( input, output, session,
	                           chemicals, stat_choices, 
                             gt_choices, basicAnalysis,
                             oedmodelname, 
                             progress, n ){
							 
	if(any(c("tfcounts","scalartop_ac50","scalartop_oed","ac50_box","tfhm", "assayhm", "toxpi", "toxpi2",  "toxpigroup") %in% stat_choices) ||
	   any(c("GTDataTable", "BMDDataTable") %in% gt_choices)){
		progress$inc(1/n, detail = "Building main analysis tables.");
		loginfo("Stats: Building main analysis tables.");
		basicAnalysis$basicData$buildBasicStatsTable(chemicals, input$analyse_background);
	}

	if("scalartop_oed" %in% stat_choices){
		progress$inc(1/n, detail = "Adding OED analysis to tables.");
		loginfo("Stats: Adding OED analysis to tables.");
		basicAnalysis$basicData$addOEDtoBasicStatsTable (chemicals, oedmodelname );
	}
	if(any("tfcounts" %in% stat_choices)){
		progress$inc(1/n, detail = "Computing target family stats.");
		loginfo("Stats: Computing target family stats.");
		#count target family occurrences
		basicAnalysis$basicData$computeTargetFamilyCounts();
		#find minimum and mean ac50 values per target
		basicAnalysis$basicData$computeTargetFamilyMinAndMeanAc50( basicAnalysis$basicData$getTargetFamilies() );
	}
	if (any(c("tfhm", "assayhm") %in% stat_choices)){
		progress$inc(1/n, detail = "Building ClusterHeatmap analysis tables.");
		loginfo("Stats: Building ClusterHeatmap analysis tables.");
		basicAnalysis$computeTargetFamilyMatrix();
		basicAnalysis$computeAssayEndpointMatrix();
	}
	
	if (any(c("toxpi", "toxpi2", "toxpigroup") %in% stat_choices)) {
		progress$inc(1/n, detail = "Building ToxPI analysis tables.");
		loginfo("Stats: Building ToxPI analysis tables.");
		#calculate toxpi plot data
		if ("toxpi" %in% stat_choices) {
			basicAnalysis$computeToxPI();
		}
		if ("toxpi2" %in% stat_choices) {
			basicAnalysis$computeToxPI2();
		}
		if ("toxpigroup" %in% stat_choices) {
			basicAnalysis$computeToxPIGroup();
		}
	}
	if(any(c("scalartop_ac50","scalartop_oed") %in% stat_choices)){
		progress$inc(1/n, detail = "Computing scalar top.");
		loginfo("Stats: Computing scalar top.");
		basicAnalysis$basicData$computeScalarTop();
	}
  
  if (any(c("GTDataTable") %in% gt_choices))
  {
    progress$inc(1/n, detail = "Computing GenTox Table");
    loginfo("Stats: Computing GenTox Table.");
    basicAnalysis$basicData$computeGenToxTable(chemicals);
  }
  
  if (any(c("BMDDataTable") %in% gt_choices))
  {
    progress$inc(1/n, detail = "Computing BMD Table");
    loginfo("Stats: Computing BMD Table.");
    basicAnalysis$basicData$computeBMDTable(chemicals);
  }
	
}
# Search-related functions -------------------------------------------------------------

#update autocomplete list
updateSearchAutocomplete = function(dtk.import, chem_autocomplete){
	dtk.import$loadChemicalData("*");
	chem_autocomplete$data <- select(dtk.import$getChemicalData(), casn, name);
}

#parses search string for multiple chemicals
parseCasnSearch = function( search_string ){
	if (is.null(search_string) || length(search_string)==0){
		return ("");
	}
	#split string, and remove any characters that are not allowed in sql queries (namely '\\ crash the scripts) 
	chems = unlist(strsplit(search_string, split="\t|\r| +| *,+ +"));
	chems = gsub(pattern = "['\\]+|[,]{2,}", replacement = "", chems);
	chems = sub(pattern = ",$", replacement = "", chems);
	chems = setdiff(chems, c(","));
	chems = unique(chems[nchar(chems)>0]);
	
	if (length(chems)>0){
		return (chems);
	} else {
		return ("");
	}
};
parseNameSearch = function( import.dreamtk, search_string ){
	if (is.null(search_string) || length(search_string)==0){
		return ("");
	}
	#split string, and remove any characters that are not allowed in sql queries (namely '\\ crash the scripts) 
	#chems = unlist(strsplit(search_string, split="\t|\r| +| *,+ +"));
	chems = gsub(pattern = "['\\]+|[,]{2,}", replacement = "", search_string);
	chems = sub(pattern = ",$", replacement = "", chems);
	chems = setdiff(chems, c(","));
	chems = unique(chems[nchar(chems)>0]);
	
	if (length(chems)>0){
		chems = import.dreamtk$findCasnFromChemname( chems );
		return (chems);
	} else {
		return ("");
	}
};

#first pass of duplicate chemical name removal before loading data and computing models
removeChemsThatExistInChemicalManager = function( session, input, output, chems, chem_manager){
	existing <- chem_manager$getListOfChemicalCasn();
	duplicates_firstpass <- which(chems %in% existing);
	if ( length(duplicates_firstpass) > 0) { 
		chems <- chems[ -duplicates_firstpass ];
		showNotification( "Duplicate entries removed.", type = "warning", duration = 5 );
		if ( length(chems) == 0 ) chems <- "";
	}
	return ( chems );
}

# Info getters -------------------------------------------------------

#obtain Hits, calculate if unavailable
getHitCount = function(chem){
	if ( !is.null(chem$analysis_results$hitcount) ) {
		hitcount <- chem$analysis_results$hitcount;
	} else {
		hitcount <- NULL;
	}
	
	if ( is.null(hitcount) ){
		hitcount <- chem$calculateHitCount();
	}
	return ( hitcount );
}

#gabriel this should probably be deleted there's no use for it.
#obtain minimum Ac50, Always calculate.
getMinAc50 = function(chem, back, hit){
	min_ac50 <- chem$calculateMinAc50(back,hit);
	return ( min_ac50 );
}
#gabriel this should probably be deleted there's no use for it.
#obtain average Ac50, Always Calculate.
getAVGAc50 = function(chem, back, hit){

	avg_ac50 <- chem$calculateAVGAc50(back,hit);
	return ( avg_ac50 );
}


#obtain minimum OED, calculate if unavailable
getMinOED = function(chem, oed_model_name ){
	if ( !is.null(chem$analysis_results$min_oed)) {
		min_oed <- chem$analysis_results$min_oed$value;
	} else {
		min_oed <- NULL;
	}
	
	if ( is.null(min_oed) ){
		min_oed <- chem$calculateMinOED(oed_model_name);
	}
	return ( min_oed );
}


#obtain basic chemical information
getBasicChemInfo = function( chem, css_model_name, oed_model_name, linesep = "<br>", back, hit ){
	#obtain css value and units
	css <- chem$model_results[[css_model_name]]$value;
	css_units <- chem$model_results[[css_model_name]]$units;
	css_assumptions <- paste0( linesep, paste0(chem$model_results[[css_model_name]]$assumptions, collapse = linesep));
	if(is.na(css)){
		css_units <- "";
	}
	
	#obtain minimum oed value and units
	min_oed <- getMinOED(chem, oed_model_name);
	min_oed_units <- chem$analysis_results$min_oed$units;
	
	#obtain minimum ac50 value and units
	min_ac50 <- getMinAc50(chem, back, hit);
	min_ac50_units <- chem$analysis_results$min_ac50$units;
	
	avg_ac50 <- getAVGAc50(chem, back, hit);
	avg_ac50_units <- chem$analysis_results$avg_ac50$units;
	
	hitcount <- getHitCount(chem);
	
	aei <- chem$analysis_results$min_ac50$aeid;
	if ( !is.null(aei) && !length(aei) == 0 ){
		chemassayrow <- filter(chem$assay_info, aeid == aei);
	} else {
		chemassayrow <- list(assay_component_endpoint_name = NA,
							organism = NA,
							tissue = NA,
							cell_short_name = NA,
							biological_process_target = NA,
							intended_target_family = NA,
							intended_target_family_sub = NA,
							gene_name = NA)
	}
	
	adj_casn <- chem$chem_info$casn;
	adj_name <- chem$chem_info$name;
	if (startsWith(adj_casn, "*") ){
		adj_casn <- gsub(pattern = "[*]", replacement = "", adj_casn);
		adj_name <- gsub(pattern = "[*]", replacement = "", adj_name);
	}
	datarow <- list(casn = adj_casn,
					name = adj_name,
					cytotoxic = chem$chem_info$cytotoxic,
					cytotoxicity_um = chem$chem_info$cytotoxicity_um,
					css = css,
					css_units = css_units,
					css_model = css_model_name,
					css_assumptions = css_assumptions,
					min_ac50 = min_ac50,
					min_ac50_units = min_ac50_units,
					avg_ac50 = avg_ac50,
					avg_ac50_units = avg_ac50_units,
					hitcount = hitcount,
					assay_component_endpoint_name = chemassayrow$assay_component_endpoint_name,
					organism = chemassayrow$organism,
					tissue = chemassayrow$tissue,
					cell_short_name = chemassayrow$cell_short_name,
					biological_process_target = chemassayrow$biological_process_target,
					intended_target_family = chemassayrow$intended_target_family,
					intended_target_family_sub = chemassayrow$intended_target_family_sub,
					gene_name = chemassayrow$gene_name,
					min_oed = min_oed,
					min_oed_units = min_oed_units,
					oed_model = oed_model_name);
	
	return ( datarow );
}

# GUI-related functions ----------------------------------------------------

#creates a list of chemical listbox choices
createChemicalCasnChoices = function( chem_manager ){
	chem_manager$sortChemicalsByCasn();
	return ( chem_manager$getListOfChemicalCasn() );
}
createChemicalNameChoices = function( chem_manager ){
	chem_manager$sortChemicalsByName();
	return ( chem_manager$getListOfChemicalNames() );
}
#updates chemical list box based on choices
updateChemListByRadioButtonSelection = function( input, output, session, selected_chem, chem_manager ){
	if (input$radio_listtype == "casn") {
		choices <- createChemicalCasnChoices( chem_manager );
		selected_chem <- selected_chem$casn;
	} else {
		choices <- createChemicalNameChoices( chem_manager );
		selected_chem <- selected_chem$name;
	}
	#main chem list
	numchems <- chem_manager$getNumChemicals();
	updateSelectInput(session, inputId = "select_chemical", selected = selected_chem,
					choices = choices, label = paste0("Searched chemical list: ", numchems, " chemicals")
	);
	#stats chem list
	updateSelectizeInput(session, inputId = "select_chemical_stats",
					choices = choices
	);
	#pca chem list
	updateSelectizeInput(session, inputId = "select_chemical_mfa",
					choices = choices
	);
	
	updateSelectizeInput(session, inputId = "select_chemical_ber",
					choices = choices
	);
	#GenTox chem list
	updateSelectizeInput(session, inputId = "select_chemical_gt",
	        choices = choices
  )
}
#updates selected chemical variable
updateSelectedChemical = function( input, output, session, selected_chem, chem_manager){
	if (input$radio_listtype == "casn") {
		selected_chem$casn <- input$select_chemical; #save selection value
		selected_chem$chemical <- chem_manager$getChemical( selected_chem$casn );
		selected_chem$name <- selected_chem$chemical$chem_info$name;
	} else {
		selected_chem$name <- input$select_chemical; #save selection value
		selected_chem$chemical <- chem_manager$getChemicalByName( selected_chem$name );
		selected_chem$casn <- selected_chem$chemical$chem_info$casn;
	}
	selected_chem$chemical <- chem_manager$getChemical( selected_chem$casn );
}

CreateSelectedChemicalUI = function(input,output,session,selected_chem, chem_manager, selected_assay, selected_assay_endpoint, model.3css,model.oed1){
				updateSelectedChemical(input, output, session, selected_chem, chem_manager);
                 #get the selected chemical
                 chem = selected_chem$chemical;
                 selected_assay$name <- input$select_assay;
                 selected_assay_endpoint$name <- input$select_assay_comp;
                 #make a grouped list for the listbox
                 
				 #here we filter based on what the user want
				 tmp = chem$assay_info;
				 
				 if(!input$select_background){
					tmp = subset(tmp, intended_target_family != "background measurement");
				 }
				 if(input$select_hit){
					tmp = subset(tmp, hitc == 1);
				 }
				 assays <- distinct(tmp, assay_name);
				 
                 #clean up the view to remove extraneous labels
                 assays <- sapply(assays, function(x){ return( unname(x) ); },
                 USE.NAMES = TRUE
                 );
                 names(assays) <- NULL;
                 x<<-assays;
                 if(!is.atomic(assays) || is.null(assays) || length(assays[[1]]) == 0 ){
                   assays <- NA;
                   selected_assay$name <- NULL;
                   selected_assay_endpoint$name <- NULL;
                 } else {
                   assays <- sort(assays);
                 }
                 if (!all(selected_assay$name %in% assays)) {
                   selected_assay$name <- NULL;
                 }
                  
                 #update assay list box
                 updateSelectInput(session, inputId = "select_assay", selected = selected_assay$name,
                                   choices = assays
                 );
                 #update chem info html text field
                 output$html_chemicalinfo <- renderUI(
                   {
                     #obtain model names for which values will be taken
                     css_model_name <- model.3css$getModelName();
                     oed_model_name <- model.oed1$getModelName();
                     
                     datarow <- getBasicChemInfo( chem, css_model_name, oed_model_name, linesep="<br>", input$select_background, input$select_hit  );
                     
                     HTML( 
                       str_c( paste0("<strong>Chemical Name: </strong>", datarow$name), 
                              paste0("<strong>CASN: </strong>", datarow$casn), 
                              paste0("<strong>Cytotoxic (<100uM cytotoxicity): </strong>", datarow$cytotoxic),
                              paste0("<strong>Cytotoxicity (uM): </strong>", datarow$cytotoxicity_um ),
							  paste0("<strong>Hit Count: </strong>", datarow$hitcount),
							  paste0("<strong>Ac50(avg): </strong>", signif(datarow$avg_ac50, digits = 4), " ", datarow$avg_ac50_units),
                              paste0("<strong>Ac50(min): </strong>", signif(datarow$min_ac50, digits = 4), " ", datarow$min_ac50_units),
                              paste0("<strong>&emsp; For: </strong>"),
                              paste0("<strong>&emsp; &emsp; Assay endpoint: </strong>", datarow$assay_component_endpoint_name),
                              paste0("<strong>&emsp; &emsp; Organism: </strong>", datarow$organism),
                              paste0("<strong>&emsp; &emsp; Tissue: </strong>", datarow$tissue),
                              paste0("<strong>&emsp; &emsp; Cell: </strong>", datarow$cell_short_name),
                              paste0("<strong>&emsp; &emsp; Biological process target: </strong>", datarow$biological_process_target),
                              paste0("<strong>&emsp; &emsp; Intended target family: </strong>", datarow$intended_target_family, ": ", datarow$intended_target_family_sub),
                              paste0("<strong>&emsp; &emsp; Gene: </strong>", datarow$gene_name),
                              paste0("<strong>Css: </strong>", signif(datarow$css, digits = 4), " ", datarow$css_units),
                              paste0("<strong>Css model: </strong>", datarow$css_model),
                              paste0("<strong>Css assumptions: </strong>", datarow$css_assumptions),
							                #paste0("<strong>OED: </strong>", signif(datarow$mean_oed, digits = 4), " ", datarow$mean_oed_units),
                              paste0("<strong>OED(min): </strong>", signif(datarow$min_oed, digits = 4), " ", datarow$min_oed_units),
                              paste0("<strong>OED model: </strong>", datarow$oed_model),
                              
                              sep = '<br>') 
                     );
                   }
                   
                 );
}

UpdateAssayComponentList = function(input,output,session,selected_chem,selected_assay,selected_assay_endpoint){

 selected_assay$name <- input$select_assay;
                  chem <- selected_chem$chemical;
                  
                  if ( is.null(chem) ||
                       is.null(selected_assay$name) ){
                    assay_components <- NA;
                    selected_assay_endpoint$name <- NULL;
                  } else {
                    #obtain component names for the chosen assay and the right endpoints that hit and aren't background compared to what the user asked
					
				 assay_components <- filter(chem$assay_info, assay_name == selected_assay$name )
				 if(!input$select_background){
					assay_components = subset(assay_components, intended_target_family != "background measurement");
				 }
				 if(input$select_hit){
					assay_components = subset(assay_components, hitc == 1);
				 }
                     
                    assay_components =  select(assay_components, assay_component_endpoint_name);
                    #clean up the view to remove extraneous labels
                    assay_components <- sapply(assay_components, function(x){ return( unname(x) ); },
                    USE.NAMES = TRUE
                    );
                    names(assay_components) <- NULL;
                    if(!is.atomic(assay_components) || is.null(assay_components) || length(assay_components[[1]]) == 0){
                      assay_components <- NA;
                      selected_assay_endpoint$name <- NULL;
                    } else {
                      assay_components <- sort(assay_components);
                    }
                  }
                  
                  if (!all(selected_assay_endpoint$name %in% assay_components)){
                    selected_assay_endpoint$name <- NULL;
                  }
                  #update assay component list box
                  updateSelectInput(session, inputId = "select_assay_comp", selected = selected_assay_endpoint$name,
                                    choices = assay_components
                  );
} 


#disable/enable searching and custom	chemical loading
disableSearching = function( app_flags ){
	app_flags$enable_searching <- FALSE;
	shinyjs::disable(id = "button_search");
	shinyjs::disable(id = "button_load_customchems");
	shinyjs::disable(id = "field_search");
}
enableSearching = function( app_flags ){
	app_flags$enable_searching <- TRUE;
	shinyjs::enable(id = "button_search");
	shinyjs::enable(id = "button_load_customchems");
	shinyjs::enable(id = "field_search");
}

#tab - clear UI outputs
clearTabUI = function( input, output, session, uilist ){
	 
	
	#clear dynamic UI stats output if exists
	for (ui in uilist){
		removeUI(
			selector = ui, immediate = TRUE
		);
	}
}

#Analysis tab - create target family count output UI
createTargetFamilyCountUI = function( input, output, session,
									basicAnalysis, id ){
  if (basicAnalysis$basicData$targetFamilyAcDataExists()){
		insertUI(
			selector = "#ui_stats",
			where = "beforeEnd",
			ui = fluidRow(id = id,
							
				box(status = "primary", title = "Target Family Analysis Plots", collapsible = TRUE, width = 12,
					column(4,
						plotlyOutput(outputId = "plot_stats_tf_pie", height=600)
					),
					column(4,
						plotlyOutput(outputId = "plot_stats_tf_minac50", height=600)
					),
					column(4,
						plotlyOutput(outputId = "plot_stats_tf_avgac50", height=600)
						)
					),
					br(),
				box(status = "primary", title = "Target Family Analysis Table", collapsible = TRUE, width = 12,
					column(8, offset = 2,
						wellPanel(
							h4("Target Family Statistics Table"),
								DT::dataTableOutput(outputId = "table_stats_tf")
							)
						)
					),
					br()
			)
		)
			
	}else{
		insertUI(
			selector = "#ui_stats",
			where = "beforeEnd",
			ui = fluidRow(id = id,
				box(status = "primary", title = "Target Family Analysis", collapsible = TRUE, width = 12,
					column(10, offset = 1,
						h3("No assay data exists. Can not make target family plots.", 
							id = "plot_stats_tf_pie", style = "color:red")
					)
				)
			)
		)
}
	
	#plot pie of target family counts
	plttfc <- basicAnalysis$plotTargetFamilyCounts ();
	
	output$plot_stats_tf_pie <- renderPlotly({
		plttfc;
	})
	
	#plot bar charts of target family vs avgac50 and minac50
	pltmin <- basicAnalysis$plotTargetFamilyMinAc50 ();
	
	output$plot_stats_tf_minac50 <- renderPlotly({
		pltmin;
	})
	
	pltavg <- basicAnalysis$plotTargetFamilyAvgAc50 ();
	
	output$plot_stats_tf_avgac50 <- renderPlotly({
		pltavg;
	})
	
	#always output data table even if it is empty, process table on client side
	output$table_stats_tf <- DT::renderDataTable(server = FALSE, {
		
		datatable(basicAnalysis$basicData$getTargetFamilyAcTable(), 
							selection="none", filter="bottom", extensions = "Buttons",
							options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
										pageLength = 10, searchHighlight = TRUE,
										scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
										formatRound(columns=c("avg_ac50", "min_ac50"), digits=5);
		
	});
	
	
}

#Analysis tab - create target family heatmap output UI
createTargetFamilyHeatmapUI = function( input, output, session,	clusterheatmap, id ){

	#target family heatmap plot
	if (input$radio_listtype == "casn"){
		plttfhm <- clusterheatmap$plotTargetFamilyDendroHeatmap(label_by="casn");
	} else if(input$radio_listtype == "name"){
		plttfhm <- clusterheatmap$plotTargetFamilyDendroHeatmap(label_by="name");
	}
	
	insertUI(
		selector = "#ui_stats",
		where = "beforeEnd",
		ui = fluidRow(id = id,
					if( clusterheatmap$targetFamilyDataExists() && !is.null(plttfhm)){
						box(status = "primary", title = "Hierarchically Clustered Heatmap of Target Subfamily activity ( log10(1/ac50) )", collapsible = TRUE, width = 12,
							fluidRow(
								column(10, offset = 1,
									plotlyOutput(outputId = "plot_stats_tfhm", height = 800)
								)
							)
						)
					} else {
						box(status = "primary", title = "Hierarchically Clustered Heatmap of Target Subfamily activity ( log10(1/ac50) )", collapsible = TRUE, width = 12,
							column(10, offset = 1,
								h3("No target family data exists. Can not make target subfamily heatmaps.", 
									id = "plot_stats_tfhm", style = "color:red")
							)
						)
					}
		)
	)
	
	
	
	if(!is.null(plttfhm)){
		output$plot_stats_tfhm <- renderPlotly({
			plttfhm;
		})
	}
}

#Analysis tab - create assay heatmap output UI
createAssayHeatmapUI = function( input, output, session,
								clusterheatmap, id ){
	
	if (input$radio_listtype == "casn"){
		pltassayhm <- clusterheatmap$plotAssayEndpointDendroHeatmap(label_by="casn");
	} else if(input$radio_listtype == "name"){
		pltassayhm <- clusterheatmap$plotAssayEndpointDendroHeatmap(label_by="name");
	}
	
	insertUI(
		selector = "#ui_stats",
		where = "beforeEnd",
		ui = fluidRow(id = id,
					if( clusterheatmap$assayEndpointDataExists() && !is.null(pltassayhm) ){
						box(status = "primary", title = "Hierarchically Clustered Heatmap of Assay Endpoint activity ( log10(1/ac50) )", collapsible = TRUE, width = 15,
							fluidRow(
								column(10, offset = 1,
										plotlyOutput(outputId = "plot_stats_assayhm", height = 800)
								)
							)
						)
					} else {
		                box(status = "primary", title = "Hierarchically Clustered Heatmap of Assay Endpoint activity ( log10(1/ac50) )", collapsible = TRUE, width = 12,
							column(10, offset = 1,
								h3("No assay data exists. Can not make assay endpoint heatmaps.", 
									id = "plot_stats_tfhm", style = "color:red")
							)
						)
					}
		)
	)
	
	
	#assay endpoint heatmap plot
	if (input$radio_listtype == "casn"){
		pltassayhm <- clusterheatmap$plotAssayEndpointDendroHeatmap(label_by="casn");
	} else if(input$radio_listtype == "name"){
		pltassayhm <- clusterheatmap$plotAssayEndpointDendroHeatmap(label_by="name");
	}
	if(!is.null(pltassayhm)){
		output$plot_stats_assayhm <- renderPlotly({
			pltassayhm;
		})
	}
}

#Analysis tab - create ac50/oed scalartop scatterplot output UI
createScalarTopUI = function( input, output, session,
							basicAnalysis, stat_choices, id ){
	
	insertUI(
		selector = "#ui_stats",
		where = "beforeEnd",
		ui = fluidRow(id = id,
						if( basicAnalysis$basicData$scalarTopDataExists() && any("scalartop_ac50" %in% stat_choices)){
							box(status = "primary", title = "Ac50 vs ScalarTop Analysis Plot", collapsible = TRUE, width = 12,
								column(10, offset = 1,
									plotlyOutput(outputId = "plot_stats_scalartop_ac50", height=600)
								)
							)
						} else if (!any("scalartop_ac50" %in% stat_choices)) {
							# do nothing
						} else {
							box(status = "primary", title = "Ac50 vs ScalarTop Analysis Plot", collapsible = TRUE, width = 12,
								column(10, offset = 1,
									h3("No assay data exists. No ac50 vs ScalarTop data to plot.", 
										id = "plot_stats_scalartop_ac50", style = "color:red")
								)
							)
						},
						br(),
						if( basicAnalysis$basicData$scalarTopDataExists() && any("scalartop_oed" %in% stat_choices) &&
							basicAnalysis$basicData$oedValuesExist()){
							box(status = "primary", title = "OED vs ScalarTop Analysis Plot", collapsible = TRUE, width = 12,
								column(10, offset = 1,
									plotlyOutput(outputId = "plot_stats_scalartop_oed", height=600)
								)
							)
						} else if (!any("scalartop_oed" %in% stat_choices)) {
						# do nothing
						} else {
							box(status = "primary", title = "OED vs ScalarTop Analysis Plot", collapsible = TRUE, width = 12,
								column(10, offset = 1,
									h3("No OED values for assays above cutoff exist. No OED vs ScalarTop data to plot.", 
									id = "plot_stats_scalartop_oed", style = "color:red")
								)
							)
						},
						br(),
						box(status = "primary", title = "ScalarTop Analysis Table", collapsible = TRUE, width = 12,
							column(8, offset = 2,
								wellPanel(
									h4("Chemical Activity Table"),
									DT::dataTableOutput(outputId = "table_stats_scalartop")
								)
							)
						),
						br()	
		)
	)
	

	if (any("scalartop_ac50" %in% stat_choices)){
		#scatterplot of scalartop vs ac50
		if (input$radio_listtype == "casn"){
			pltst_ac50 <- basicAnalysis$plotAc50VsScalarTop(label_by="casn");
		} else if(input$radio_listtype == "name"){
			pltst_ac50 <- basicAnalysis$plotAc50VsScalarTop(label_by="name");
		}
		
		output$plot_stats_scalartop_ac50 <- renderPlotly({
			pltst_ac50;
		})
	}
	
	if (any("scalartop_oed" %in% stat_choices)){
		#scatterplot of scalartop vs oed
		if (input$radio_listtype == "casn"){
			pltst_oed <- basicAnalysis$plotOEDVsScalarTop(label_by="casn");
		} else if(input$radio_listtype == "name"){
			pltst_oed <- basicAnalysis$plotOEDVsScalarTop(label_by="name");
		}
		
		output$plot_stats_scalartop_oed <- renderPlotly({
			pltst_oed;
		})
	}
	
	#always output data table even if it is empty, process table on client side
	output$table_stats_scalartop <- DT::renderDataTable(server = FALSE,  {
		
		round_cols <- c("ac50", "ac_cutoff", "ac_top", "scalar_top");
		if (any("scalartop_oed" %in% stat_choices)){
			round_cols <- c(round_cols, "oed");
		}
		
		sctable <- basicAnalysis$basicData$getScalarTopTable();
		sctable$ac50 <- 10^(sctable$ac50); #for consistency, convert all values to uM units
		sctable$ac_cutoff <- 10^(sctable$ac_cutoff);
		#sctable$ac_top <- 10^(sctable$ac_top);
		datatable( sctable, 
					selection="none", filter = "bottom", extensions = "Buttons",
					options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
					pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
					scrollX=TRUE, scrollCollapse=TRUE), rownames= FALSE) %>%
					formatRound(columns = round_cols, digits=5);
		
	})
	
}
#Analysis tab - create ac50 box output ui
createAc50BoxUI = function(input,output,session,basicAnalysis,id){
	insertUI(
		selector = "#ui_stats",
		where = "beforeEnd",
		ui = fluidRow(id = id,

							if( basicAnalysis$basicData$basicStatsDataExists()){
								box(status = "primary", title = "Burst Assay Analysis", collapsible = TRUE, width = 12,
							
									fluidRow(
										column(10, offset = 1,
										plotlyOutput(outputId = "plot_Ac50_box", height = 800)
										)
									)
								)
							}else{
								box(status = "primary", title = "Burst Assay Analysis", collapsible = TRUE, width = 12,
									column(12, offset = 0,
										h3("No Ac50 Data available for the chemical(s). No graph to plot.", 
										id = "plot_Ac50_box", style = "color:red")
									)
								)
							}	
		)
	)
	
	basicAnalysis$plotAc50Box(label_by = input$radio_listtype);
	output$plot_Ac50_box <- renderPlotly({basicAnalysis$getPlotAc50Box()});
}

#Analysis tab - create toxpi output UI
createToxPIUI = function( input, output, session,
						basicAnalysis, numChems, id){
	insertUI(
		selector = "#ui_stats",
		where = "beforeEnd",
		ui = fluidRow(id = id,
		if( basicAnalysis$toxPIDataExists() ){
			box(status = "primary", title = "ToxPI Plots", collapsible = TRUE, width = 12,
				dropdownButton(
					tags$h3("About plot"),
							HTML( 
							str_c( "<p>Individual ToxPI Plots represent the relative size of each target family and their relative activity levels.</p>",
								   "<p><strong>Slice angle:</strong> relative number of target family endpoints within the chemical.</p>",
								   "<p><strong>Slice height:</strong> target activity expressed as LOG10( 1/ac50 ) normalized within chemical to a range between 0 and 1.</p>"
                            )
                          ),
                          circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
                          tooltip = tooltipOptions(title = "About plot")
                        ),
                        fluidRow(
                          column(10, offset = 1,
                                 plotOutput(outputId = "plot_stats_toxpi", height = ceiling(numChems/num_plot_cols)*450)
                          )
                        )
                    )} else {
                      box(status = "primary", title = "Individual ToxPI Plots", collapsible = TRUE, width = 12,
                          column(10, offset = 1,
                                 h3("No target family data exists. Can not make toxpi plots.", 
                                    id = "plot_stats_toxpi", style = "color:red")
                          )
                      )
                    }
    )
  )

  
  #toxPI plots
  if (input$radio_listtype == "casn"){
    pltpi <- basicAnalysis$plotToxPI(label_by="casn", grouped = FALSE);
  } else if(input$radio_listtype == "name"){
    pltpi <- basicAnalysis$plotToxPI(label_by="name", grouped = FALSE);
  }

  output$plot_stats_toxpi <- renderPlot({
	multiplot2(plotlist=pltpi, cols = num_plot_cols, plot_data = basicAnalysis$getToxPIPlotData());
  });

  
}

#Analysis tab - create toxpi output UI
createToxPI2UI = function( input, output, session,
                          basicAnalysis, numChems, id ){
  
  insertUI(
    selector = "#ui_stats",
    where = "beforeEnd",
    ui = fluidRow(id = id,
                  if( basicAnalysis$toxPI2DataExists() ){
                    box(status = "primary", title = "ToxPI Plots", collapsible = TRUE, width = 12,
                        dropdownButton(
                          tags$h3("About plot"),
                          HTML( 
                            str_c( "<p>Individual ToxPI Plots represent the relative activity levels of response, top scaled response and cytotoxicity concentration for each chemical.</p>",
                                   "<p><strong>Slice height:</strong> target activity expressed as LOG10( ac50 ) normalized within chemical to a range between 0 and 4.</p>"
                            )
                          ),
                          circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
                          tooltip = tooltipOptions(title = "About plot")
                        ),
                        fluidRow(
                          column(10, offset = 1,
                                 plotOutput(outputId = "plot_stats_toxpi2", height = ceiling(numChems/num_plot_cols)*450)
                          )
                        )
                    )} else {
                      box(status = "primary", title = "Individual ToxPI Plots", collapsible = TRUE, width = 12,
                          column(10, offset = 1,
                                 h3("No target family data exists. Can not make toxpi plots.", 
                                    id = "plot_stats_toxpi2", style = "color:red")
                          )
                      )
                    }
    )
  )
  
  
  #toxPI plots
  if (input$radio_listtype == "casn"){
    pltpi <- basicAnalysis$plotToxPI2(label_by="casn");
  } else if(input$radio_listtype == "name"){
    pltpi <- basicAnalysis$plotToxPI2(label_by="name");
  }
  
  output$plot_stats_toxpi2 <- renderPlot({
    multiplot2(plotlist=pltpi, cols = num_plot_cols, plot_data = basicAnalysis$getToxPI2PlotData());
  })
  
}


#Analysis tab - create grouped toxpi output UI
createToxPIGroupUI = function( input, output, session,
                               basicAnalysis, numChems, id ){
  
  num_plot_cols <- 1;
  insertUI(
    selector = "#ui_stats",
    where = "beforeEnd",
    ui = fluidRow(id = id,
                  if( basicAnalysis$toxPIGroupDataExists() ){
                    box(status = "primary", title = "Grouped ToxPI Plots", collapsible = TRUE, width = 12,
                        dropdownButton(
                          tags$h3("About plot"),
                          HTML( 
                            str_c( "<p>Grouped ToxPI Plots represent the target families and their relative activity levels across multiple chemicals.</p>",
                                   "<p><strong>Slice angle:</strong> denotes a target family positioned and coloured the same across all chemicals of interest.</p>",
                                   "<p><strong>Slice height:</strong> target activity expressed as LOG10( 1/ac50 ) normalized for <strong>all</strong> chemicals to a range between 0 and 1.</p>"
                            )
                          ),
                          circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
                          tooltip = NULL
                        ),
                        fluidRow(
                          column(10, offset = 1,
                                 plotOutput(outputId = "plot_stats_toxpigrp", height = ceiling(numChems/num_plot_cols)*450)
                          )
                        )
                    )} else {
                      box(status = "primary", title = "Grouped ToxPI Plots", collapsible = TRUE, width = 12,
                          column(10, offset = 1,
                                 h3("No target family data exists. Can not make toxpi plots.", 
                                    id = "plot_stats_toxpigrp", style = "color:red")
                          )
                      )
                    }
    )
  )
  
  
  #toxPI plots
  if (input$radio_listtype == "casn"){
    pltpigrp <- basicAnalysis$plotToxPI(label_by="casn", grouped = TRUE);
  } else if(input$radio_listtype == "name"){
    pltpigrp <- basicAnalysis$plotToxPI(label_by="name", grouped = TRUE);
  }
  
  output$plot_stats_toxpigrp <- renderPlot({
    multiplot2(plotlist=pltpigrp, cols = num_plot_cols, plot_data = basicAnalysis$getToxPIGroupPlotData());
  })
}

#Analysis Tab - create GenTox Table
createGenToxDataTable = function(input, output, session,
                                 basicAnalysis, id)
{
  insertUI(
    selector = "#ui_stats",
    where = "beforeEnd",
    ui = fluidRow(id = id,
                  
                  if( basicAnalysis$basicData$basicStatsDataExists()){
                    box(status = "primary", title = "Gen Tox Data Table", collapsible = TRUE, width = 12,
                        column(width = 8, offset = 2,
                          wellPanel(
                            h4("GenTox Table"),
                            DT::dataTableOutput(outputId = "table_gentox_genToxDataTable")
                          )
                        )
                    )
                  }else{
                    box(status = "primary", title = "Gen Tox Data Table", collapsible = TRUE, width = 12,
                        column(8, offset = 2,
                               h3("No GenTox Data available for the chemical(s). No graph to plot.", 
                                  id = "plot_genToxData_tbl", style = "color:red")
                        )
                    )
                  }	
    )
  )
  
  #always output data table even if it is empty, process table on client side
  output$table_gentox_genToxDataTable <- DT::renderDataTable(server = FALSE, {
    
    datatable(basicAnalysis$basicData$getGenToxTable(), 
              selection="none", filter="bottom", extensions = "Buttons",
              options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
                           pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
                           scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
      formatStyle(seq(ncol(basicAnalysis$basicData$getGenToxTable())-1), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
    
  })
}

#Analysis Tab - create BMD Table
createBMDDataTable = function(input, output, session,
                                 basicAnalysis, id)
{
  insertUI(
    selector = "#ui_stats",
    where = "beforeEnd",
    ui = fluidRow(id = id,
                  
                  if( basicAnalysis$basicData$basicStatsDataExists()){
                    box(status = "primary", title = "BMD Data Table", collapsible = TRUE, width = 12,
                        column(width = 8, offset = 2,
                          wellPanel(
                            h4("BMD Table"),
                            DT::dataTableOutput(outputId = "table_gentox_bmdDataTable")
                          )
                        )
                    )
                  }else{
                    box(status = "primary", title = "BMD Data Table", collapsible = TRUE, width = 12,
                        column(8, offset = 2,
                               h3("No BMD Data available for the chemical(s). No graph to plot.", 
                                  id = "plot_bmd_tbl", style = "color:red")
                        )
                    )
                  }	
    )
  )
  
  #always output data table even if it is empty, process table on client side
  output$table_gentox_bmdDataTable <- DT::renderDataTable(server = FALSE, {
    
    datatable(basicAnalysis$basicData$getBMDTable(), 
              selection="none", filter="bottom", extensions = "Buttons",
              options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
                           pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
                           scrollX=TRUE, scrollCollapse=TRUE), rownames = FALSE) %>%
      formatStyle(seq(ncol(basicAnalysis$basicData$getBMDTable())-1), "border-right" = "solid 1px", "border-right-color" =  "rgba(221, 221, 221, 0.2)");
    
  })
}


#MFA tab - create MFA output UI
createMFAUI = function( input, output, session,
                        chem_data, chem_data_mfa, id ){
  
  insertUI(
    selector = "#ui_mfa",
    where = "beforeEnd",
    ui = fluidRow(id = id,
                  box(status = "primary", title = "Individual factor map: Target Family group centers.", collapsible = TRUE, width = 12,
                      plotOutput(outputId = "plot_mfa_ind", height=800) #%>% withSpinner(color="#0dc5c1")
                  ),
                  box(status = "primary", title = "Individual factor map: inertia of Target Family group centers.", collapsible = TRUE, width = 12,
                      plotOutput(outputId = "plot_mfa_ind2", height=800) #%>% withSpinner(color="#0dc5c1")
                  ),
                  box(status = "primary", title = "Individual factor map: individual chemical assays and partial inertias.", collapsible = TRUE, width = 12,
                      plotOutput(outputId = "plot_mfa_ind3", height=800) #%>% withSpinner(color="#0dc5c1")
                  ),
                  box(status = "primary", title = "Correlation circle: numeric variables.", collapsible = TRUE, width = 12,
                      plotOutput(outputId = "plot_mfa_var", height=800) #%>% withSpinner(color="#0dc5c1")
                  ),
                  box(status = "primary", title = "Grouping representation plot by variables.", collapsible = TRUE, width = 12,
                      plotOutput(outputId = "plot_mfa_group", height=800) #%>% withSpinner(color="#0dc5c1")
                  )
                  
    )
  );
  
  output$plot_mfa_ind <- renderPlot(
    {
      tryCatch({
        #cpalette <- palette(c("#1468DC", "#FF5C01", "#A0B0AF", "#3B494A", "#5BB2F7", "#C32C01",
        #              "#738786", "#A1DBFF", "#FFA865", "#C3D9E4"));
        #cpalette <- palette(rainbow(30, s = 3/4, v=3/4));
        cpalette <- palette(topo.colors(length(unique(chem_data$intended_target_family))*3));
        plot(chem_data_mfa, choix="ind", habillage="intended_target_family", lab.ind = FALSE,
             palette=cpalette);
      },
      error = function(e){
        logerror("Could not run MFA. FactoMineR probably did not like something about missing data");
        logerror(e);
        showNotification("Could not run MFA. FactoMineR probably did not like something about missing data",
                         type = "error", duration = 10);
        showNotification(paste0("Error detail: ", e),
                         type = "error", duration = 10);
        return (NULL);
      }
      )
    }
  )
  
  output$plot_mfa_ind2 <- renderPlot(
    {
      tryCatch({
        #cpalette <- palette(c("#1468DC", "#FF5C01", "#A0B0AF", "#3B494A", "#5BB2F7", "#C32C01",
        #              "#738786", "#A1DBFF", "#FFA865", "#C3D9E4"));
        #cpalette <- palette(rainbow(30, s = 3/4, v=3/4));
        cpalette <- palette(topo.colors(10));
        plot(chem_data_mfa, choix="ind", habillage="group", invisible = "ind", lab.var = TRUE,
             partial=chem_data_mfa$plotlist$ind1$partial, axes = chem_data_mfa$plotlist$ind1$axes,
             palette=cpalette);
      },
      error = function(e){
        logerror("Could not run MFA. FactoMineR probably did not like something about missing data");
        logerror(e);
        showNotification("Could not run MFA. FactoMineR probably did not like something about missing data",
                         type = "error", duration = 10);
        showNotification(paste0("Error detail: ", e),
                         type = "error", duration = 10);
        return (NULL);
      }
      )
    }
  )
  
  output$plot_mfa_ind3 <- renderPlot(
    {
      tryCatch({
        #cpalette <- palette(c("#1468DC", "#FF5C01", "#A0B0AF", "#3B494A", "#5BB2F7", "#C32C01",
        #              "#738786", "#A1DBFF", "#FFA865", "#C3D9E4"));
        #cpalette <- palette(rainbow(30, s = 3/4, v=3/4));
        cpalette <- palette(topo.colors(10));
        plot(chem_data_mfa, choix="ind", habillage="group", invisible = "quali", lab.ind = TRUE,
             partial=chem_data_mfa$plotlist$ind2$partial, axes = chem_data_mfa$plotlist$ind2$axes,
             palette=cpalette);
      },
      error = function(e){
        logerror("Could not run MFA. FactoMineR probably did not like something about missing data");
        logerror(e);
        showNotification("Could not run MFA. FactoMineR probably did not like something about missing data",
                         type = "error", duration = 10);
        showNotification(paste0("Error detail: ", e),
                         type = "error", duration = 10);
        return (NULL);
      }
      )
    }
  )
  
  output$plot_mfa_var <- renderPlot(
    {
      tryCatch({
        plot(chem_data_mfa, choix="var", palette=palette(topo.colors(10)));
      },
      error = function(e){
        logerror("Could not run MFA. FactoMineR probably did not like something about missing data");
        logerror(e);
        showNotification("Could not run MFA. FactoMineR probably did not like something about missing data",
                         type = "error", duration = 10);
        showNotification(paste0("Error detail: ", e),
                         type = "error", duration = 10);
        return (NULL);
      }
      )
    }
  )
  
  output$plot_mfa_group <- renderPlot(
    {
      tryCatch({
        plot(chem_data_mfa, choix="group", palette=palette(topo.colors(10)));
      },
      error = function(e){
        logerror("Could not run MFA. FactoMineR probably did not like something about missing data");
        logerror(e);
        showNotification("Could not run MFA. FactoMineR probably did not like something about missing data",
                         type = "error", duration = 10);
        showNotification(paste0("Error detail: ", e),
                         type = "error", duration = 10);
        return (NULL);
      }
      )
    }
  )
}
  
#BER tab - create BER output UI
createBERUI = function( input, output, session, BERAnalysis, id ){
 
  insertUI(
    selector = "#ui_ber",
    where = "beforeEnd",
	ui = fluidRow(id = id,
            if( BERAnalysis$BERData$calcBERStatsDataExists()){
              box(status = "primary", title = "In Vitro Tox Concentrations vs Equivalent Doses", collapsible = TRUE, width = 12,
                  dropdownButton(
                    tags$h3("About plot"),
                    HTML( 
                      str_c( 
                        "<p><strong>Consumer Products Exposure  BE:</strong>Consumer Products Exposure  BER of a product use category for a chemical is the summation of the direct incidental ingestion, direct aerosol ingestion and direct vapor ingestion. The value is shown in ug/day.</p>",
                        "<p><strong>Ac50:</strong> This value is shown in uM.</p>"
                      )
                    ),
                    circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
                    tooltip = NULL
                  ),
                  br(),
                  
                  fluidRow(
                    column(10, offset = 1,
                           plotlyOutput(outputId = "plot_OED_vs_Ac50", height = 800)
                    )
                  )
              )
            }else{
              box(status = "primary", title = "In Vitro Tox Concentrations vs Equivalent Doses", collapsible = TRUE, width = 12,
                  column(12, offset = 0,
                         h3("No Consumer Products Exposure Data available for the chemical(s). No graph to plot.", 
                            id = "plot_OED_vs_Ac50", style = "color:red")
                  )
              )
            },
            
            br(),						
	              
	              
	          if( BERAnalysis$BERData$calcBERStatsDataExists()){
							box(status = "primary", title = "Consumer Products Exposure VS In Vitro Tox Equivalent Dose", collapsible = TRUE, width = 12,
								dropdownButton(
									  tags$h3("About plot"),
									  HTML( 
										str_c( 
												"<p><strong>Consumer Products Exposure  BE:</strong>Consumer Products Exposure  BER of a product use category for a chemical is the summation of the direct incidental ingestion, direct aerosol ingestion and direct vapor ingestion. The value is shown in ug/day.</p>",
												"<p><strong>Ac50:</strong> This value is shown in uM.</p>"
										)
									  ),
									  circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
									  tooltip = NULL
								),
								br(),
								
								fluidRow(
									column(10, offset = 1,
									plotlyOutput(outputId = "plot_BER_vs_Ac50", height = 800)
									)
								)
							)
						}else{
							box(status = "primary", title = "Consumer Products Exposure VS In Vitro Tox Equivalent Dose", collapsible = TRUE, width = 12,
								column(12, offset = 0,
									h3("No Consumer Products Exposure Data available for the chemical(s). No graph to plot.", 
									id = "plot_BER_vs_Ac50", style = "color:red")
								)
							)
						},

						br(),
						
						if( BERAnalysis$BERData$calcBERStatsDataExists()){
							box(status = "primary", title = "Consumer Products Exposure Analysis Table", collapsible = TRUE, width = 12,
								dropdownButton(
									  tags$h3("About plot"),
									  HTML( 
										str_c( "<p>This shows the Consumer Products Exposure values in ug/day for each category of exposure.</p>")
										),

									  circle = TRUE, status = "danger", icon = icon("question-circle"), witdth = "300px",
									  tooltip = NULL
								),
								br(),
								column(12, offset = 0,
									wellPanel(
										DT::dataTableOutput(outputId = "table_stats_BER")
										)
									)
								)
						}else{
							box(status = "primary", title = "BE Analysis Table", collapsible = TRUE, width = 12,
								column(12, offset = 0,
									h3("No Consumer Products Exposure Data available for the chemical(s). No Table to plot.", 
									id = "table_stats_BER", style = "color:red")
								)
							)
						},
						br(),
						if( BERAnalysis$BERData$calcBERStatsDataExists()){
							box(status = "primary", title = "BER box chart", collapsible = TRUE, width = 12,
								dropdownButton(
									  tags$h3("About plot"),
									  HTML( 
										str_c( 
												"<p><strong>Oral Consumer Products Exposure:</strong> Oral Consumer Products Exposure of a product use category for a chemical is the summation of the direct incidental ingestion, direct aerosol ingestion and direct vapor ingestion.</p>",
												"<p><strong>Oral BER:</strong> Mean Ac50 of the chemical divided by Oral Consumer Products Exposure</p>"
										)
									  ),
									  circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
									  tooltip = NULL
								),
								br(),
								
								fluidRow(
									column(10, offset = 1,
									plotlyOutput(outputId = "plot_BER", height = 800)
									)
								)
							)
						}else{
							box(status = "primary", title = "BER box chart", collapsible = TRUE, width = 12,
								column(12, offset = 0,
									h3("No BER Data available for the chemical(s). No graph to plot.", 
									id = "plot_BER", style = "color:red")
								)
							)
						},	

						br(),
						if( BERAnalysis$BERData$calcBERStatsDataExists()){
							box(status = "primary", title = "BER Agregated Analysis Table", collapsible = TRUE, width = 12,
								dropdownButton(
									  tags$h3("About plot"),
									  HTML( 
										str_c( "<p><strong>Oral Consumer Products Exposure:</strong> Oral Consumer Products Exposure of a product use category for a chemical is the summation of the direct incidental ingestion, direct aerosol ingestion and direct vapor ingestion.</p>",
												"<p><strong>Mean:</strong> This value is calculating by getting the average of all the product categories of each chemical.</p>",
												"<p><strong>Mean and Min over Oral Consumer Products Exposure:</strong> The Ac50 value is divided by the value closest to the 95th percentile of exposure.</p>"
										)
									  ),
									  circle = TRUE, status = "danger", icon = icon("question-circle"), width = "300px",
									  tooltip = NULL
								),
								br(),
								
								column(12, offset = 0,
									wellPanel(
										DT::dataTableOutput(outputId = "table_stats_mean_BER")
										)
									)
								)
						}else{
							box(status = "primary", title = "Daily Intake Agregated Analysis Table", collapsible = TRUE, width = 12,
								column(12, offset = 0,
									h3("No exposure Data available for the chemical(s). No Table to plot.", 
									id = "table_stats_mean_BER", style = "color:red")
								)
							)
							
						}

						
						
		)
		
	);

	output$table_stats_BER <- DT::renderDataTable(server = FALSE,  {
		datatable( BERAnalysis$BERData$getCalcBERStatsTable(), 
					selection="none", filter = "bottom", extensions = "Buttons",
					options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
					pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
					scrollX=TRUE, scrollCollapse=TRUE), rownames= FALSE) 

		}
	);
	
	output$table_stats_mean_BER <- DT::renderDataTable(server = FALSE,  {
		datatable( BERAnalysis$getMeanBERTable(), 
					selection="none", filter = "bottom", extensions = "Buttons",
					options=list(buttons = c('copy', 'csv', 'excel'), dom = "Blfrtip",
					pageLength = 10, searchHighlight = TRUE, lengthMenu = c(10, 20, 50, 100),
					scrollX=TRUE, scrollCollapse=TRUE), rownames= FALSE) 

		}
	);
	output$plot_BER_vs_Ac50 <- renderPlotly({BERAnalysis$getBERvsAc50plot()});
	output$plot_OED_vs_Ac50 <- renderPlotly({BERAnalysis$getOEDvsAc50plot()});
	output$plot_BER <- renderPlotly({BERAnalysis$getBERplot()});
}

#GenTox21 Importing Function --------------------------------

importTables <- function(tables, tableID = 1) {
  #if there is no table, no need to load everything
  if (length(tables) == 0) 
    return();
  
  load("data/GenTox21.RData");
  
  #They will be used for the end message
  allTable <- TRUE;
  anyTable <- FALSE;
  
  progress <- shiny::Progress$new();
  on.exit({
    progress$close();
  });
  progress$set(message = "", value = 0);
  
  GT <- Class.Analysis.GT$new();
  
  for (t in tables)
  {
    progress$inc(1/length(tables), detail = "Importing Tables.");
    if (ncol(t) < 2)
    {
      allTable <- FALSE;
      next;
    }
    
    if (tableID == 1)
    {
      newRows <- ldply(importTable(t, tableID), data.frame, stringsAsFactors = FALSE)
      gen_tox_data <- rbind(gen_tox_data, newRows);
      trend_checker <- newRows[!duplicated(newRows[,c("compound", "assay", "endpoint", "s9_positive")]),];
      if (nrow(trend_checker) > 0)
      {
        for (i in seq(nrow(trend_checker)))
        {
          if (nrow(trend_data) < 1)
          {
            trd = c();
          }
          else
          {
            trd <- which(trend_data["compound"] == trend_checker[i, "compound"] &
                           trend_data["assay"] == trend_checker[i, "assay"] &
                           trend_data["endpoint"] == trend_checker[i, "endpoint"] &
                           trend_data["s9_positive"] == trend_checker[i, "s9_positive"]);
          }
          
          currChemEndptData <- gen_tox_data[gen_tox_data[,"compound"] == trend_checker[i, "compound"] &
                                              gen_tox_data["assay"] == trend_checker[i, "assay"] &
                                              gen_tox_data["endpoint"] == trend_checker[i, "endpoint"] &
                                              gen_tox_data["s9_positive"] == trend_checker[i, "s9_positive"], ];
          
          print(currChemEndptData)
          
          logDose <- log(currChemEndptData[,"concentration"])
          logMnResp <- log(currChemEndptData[,"response"])
          
          MAD <- median(abs(logMnResp-median(logMnResp)))
          
          paramsMn <- tcplFit(logc=logDose, resp=logMnResp, bmad=MAD, force.fit=TRUE)
          
          aic_hill_mn <- paramsMn$hill_aic
          aic_gnls_mn <- paramsMn$gnls_aic
          
          if (length(trd) == 1)
          {
            trend_data[trd, "trend"] <- (!is.na(aic_hill_mn) || !is.na(aic_gnls_mn))
          }
          else
          {
            trend_data <- rbind(trend_data, data.frame(
              casn = trend_checker[i, "casn"],
              compound = trend_checker[i, "compound"],
              assay = trend_checker[i, "assay"],
              endpoint = trend_checker[i, "endpoint"],
              s9_positive = trend_checker[i, "s9_positive"],
              trend = (!is.na(aic_hill_mn) || !is.na(aic_gnls_mn)), stringsAsFactors = FALSE))
          }
        }
      }
    }
    else if (tableID == 2)
    {
      newBMDs <- ldply(importTable(t, tableID), data.frame, stringsAsFactors = FALSE);
      if (nrow(newBMDs) == 0)
      {
        showNotification("The table was not imported. No valid rows were presented", type = "error", duration = 5);
        return();
      }
      for (i in seq(nrow(newBMDs)))
      {
        if (nrow(bmd_data) < 1)
        {
          bmdToChange = c();
        }
        else
        {
          bmdToChange <- which(bmd_data["compound"] == newBMDs[i, "compound"] &
                                 bmd_data["assay"] == newBMDs[i, "assay"] &
                                 bmd_data["endpoint"] == newBMDs[i, "endpoint"] &
                                 bmd_data["s9_positive"] == newBMDs[i, "s9_positive"] &
                                 bmd_data["model_type"] == newBMDs[i, "model_type"] &
                                 bmd_data["bmr"] == newBMDs[i, "bmr"]);
        }
        
        if (length(bmdToChange) == 1)
        {
          bmd_data[bmdToChange, "aic"] <- newBMDs[i, "aic"];
          bmd_data[bmdToChange, "log_likelihood"] <- newBMDs[i, "log_likelihood"];
          bmd_data[bmdToChange, "bmdl"] <- newBMDs[i, "bmdl"];
          bmd_data[bmdToChange, "bmd"] <- newBMDs[i, "bmd"];
          bmd_data[bmdToChange, "bmdu"] <- newBMDs[i, "bmdu"];
        }
        else
        {
          bmd_data <- rbind(bmd_data, newBMDs[i, ]);
        }
      }
    }
    anyTable <- TRUE;
  } 
  
  save(gen_tox_data, bmd_data, trend_data, assay_description, file="data/GenTox21.RData");
  
  shinyjs::hide("button_importFiles");
  showNotification(ifelse(anyTable, ifelse(allTable, "All tables were imported", "Valid tables were imported."),"No tables were able to be imported"), type = "message", duration = 5);
}

importTable <- function(table, tableID = 1) {
  progressBar <- shiny::Progress$new();
  on.exit({
    progressBar$close();
  });
  progressBar$set(message = "", value = 0);
  if (tableID == 1)
  {
    endpoints <- c(ifelse(!is.na(ept <- str_extract(table[1, grepl("(experiment)s?", tolower(colnames(table)), ignore.case = TRUE)], "(?<=II )(.*)(?= (\\+|-)(s9|S9))")), list(ept),
                          ifelse(!is_empty(ncn <- colnames(table)[!is.na(str_extract(colnames(table), "^[^a-zA-Z\\s%].*"))]), list(ncn),
                                ifelse(!is_empty(per <- colnames(table)[!is.na(str_extract(colnames(table), "^[^a-zA-Z\\s].*"))]), list(per),
                                    ifelse(!is.na(str_extract(colnames(table)[ncol(table)], "^.*(MF|mf).*(.?[^\\w\\n\\s-]).*")), list("Mutation Frequency"),
                                        list(str_extract(table[1, grepl("(experiment)s?", tolower(colnames(table)), ignore.case = TRUE)], "(^(.*)I)|(^.*?(?= (\\+|-)(s9|S9)))")))))));
    rm(ept);
    try(rm(ncn));
    try(rm(per));
  }
  else if (tableID == 2)
  {
    gtAnalysis <- Class.Analysis.GT$new();
    bmrModels <- c("exp", "hill", "geo") 
  }
  
  newRows = list();
  
  casImporter <- Class.DataImporter.dreamtk.basicData$new();
  casImporter$setSource("./data/DreamTKv0.8.RData");
  
  for (j in seq(nrow(table)))
  {
    progressBar$inc(1/nrow(table), detail = "Importing Table.");
    
    if (any(c("cytotoxic") %in% tolower(colnames(table))))
    {
      if (tolower(table[j, grepl("^\\bcytotoxic\\b$", tolower(colnames(table)), ignore.case = TRUE)]) == "yes")
      {
        next;
      }
    }
    
    if (tableID == 1)
    {
      if (any(c("cas", "casn") %in% tolower(colnames(table))))
      {
        cn <- table[j, grepl("(cas)n?", tolower(colnames(table)), ignore.case = TRUE)];
      }
      else
      {
        cn <- casImporter$findCasnFromChemname(table[j, grepl("(compound)s?", tolower(colnames(table)), ignore.case = TRUE)])
      }
      
      for (i in (seq(length(endpoints[[1]]))))
      {
        newRows[[(j-1)*length(endpoints[[1]])+i]] <- list(
          casn=ifelse(!is.na(cn) && !is.null(cn) && cn != table[j, grepl("(compound)s?", tolower(colnames(table)), ignore.case = TRUE)], cn, "123-456-7"),
          compound=table[j, grepl("(compound)s?", tolower(colnames(table)), ignore.case = TRUE)],
          assay=str_extract(table[j, grepl("(experiment)s?", tolower(colnames(table)), ignore.case = TRUE)], "(^(.*)II)|(^.*?(?= (\\+|-)(s9|S9)))"),
          endpoint= (endpoints[[1]])[[i]],
          s9_positive= str_extract(table[j, grepl("(experiment)s?", tolower(colnames(table)), ignore.case = TRUE)], ".?(?=(s9|S9))") == "+",
          concentration=table[j, grepl("(concentration)s?", tolower(colnames(table)), ignore.case = TRUE)],
          concentration_unit=sub("\\).*", "", sub(".*\\(", "", colnames(table)[grepl("(concentration)s?", tolower(colnames(table)), ignore.case = TRUE)])),
          response=ifelse((endpoints[[1]])[[i]] %in% colnames(table),table[j, (endpoints[[1]])[[i]]], table[j, ncol(table)]),
          response_type=ifelse((endpoints[[1]])[[i]] %in% colnames(table), str_extract((endpoints[[1]])[[i]], "^\\W"), str_extract(colnames(table)[ncol(table)], "(^[^\\w\\s]|[(x|X)]?[\\s]?10[\\s]?\\^.*)")));
      }
    }
    else if (tableID == 2)
    {
      cn <- casImporter$findCasnFromChemname(table[j, grepl("(compound)s?", tolower(colnames(table)), ignore.case = TRUE)])
      for (i in seq(length(bmrModels)))
      {
        newRows[[(j-1)*length(bmrModels)+i]] <- list(
          casn=ifelse(!is.na(cn) && !is.null(cn) && cn != table[j, grepl("(compound)s?", tolower(colnames(table)), ignore.case = TRUE)], cn, "123-456-7"),
          compound=table[j, grepl("(compound)s?", tolower(colnames(table)), ignore.case = TRUE)],
          assay=table[j, grepl("(assay)s?", tolower(colnames(table)), ignore.case = TRUE)],
          endpoint= table[j, grepl("(endpoint)s?", tolower(colnames(table)), ignore.case = TRUE)],
          s9_positive= (table[j, grepl("(s9)", tolower(colnames(table)), ignore.case = TRUE)]) == "+",
          model_type= ifelse(i == 1, "exponential", ifelse(i == 2, "hill curve", "geometric mean")),
          aic= ifelse(i == 3, NA, table[j, grepl(paste("(", bmrModels[i], "_aic)s?", sep=""), tolower(colnames(table)), ignore.case = TRUE)]),
          log_likelihood= ifelse(i==3, NA, table[j, grepl(paste("(", bmrModels[i], "_loglik)s?", sep=""), tolower(colnames(table)), ignore.case = TRUE)]),
          bmr = as.numeric(table[j, grepl("(bmr)s?", tolower(colnames(table)), ignore.case = TRUE)])*100,
          bmdl = table[j, grepl(paste("(", bmrModels[i], "_bmdl)s?", sep=""), tolower(colnames(table)), ignore.case = TRUE)],
          bmd = table[j, grepl(paste("(", bmrModels[i], "_bmd)s?$", sep=""), tolower(colnames(table)), ignore.case = TRUE)],
          bmdu = table[j, grepl(paste("(", bmrModels[i], "_bmdu)s?", sep=""), tolower(colnames(table)), ignore.case = TRUE)]);
      }
    }
      
    rm(cn);
  }
  if (tableID == 1) rm(endpoints);
  rm(casImporter);
  return(newRows);
}
