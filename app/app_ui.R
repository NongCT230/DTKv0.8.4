# Web UI for shiny

#v0.7

#dashboard version

# HELPER FUNCTIONS ------------------------------

#adds id parameter to a menuItem. Allows calling the generated list item by selector reference
menuItemAddId = function( menuitem, id){
  #menuitem$attribs$id <- id;
  menuitem$children[[1]]$attribs$id <- id;
  return (menuitem);
} 

# UI COMPONENTS -------------------------------
#header is the top stripe panel
header <- dashboardHeader(title = paste0("DreamTK v",app.version.major,".",app.version.minor, ".", app.version.revision));
RDataFile <- "DreamTKv0.8.RData";
#source()

#side bar is the expandable menu bar on the left side
sidebar <- dashboardSidebar(
  # Custom CSS to hide the default logout panel
  tags$head(tags$style(HTML('.shiny-server-account { display: none; } '))),
  #Some classes might need some special cases.
  tags$head(tags$style(HTML('#gt_custom_choice::placeholder { text-align: left; } 
                            #gt_custom_choice:-ms-input-placeholder { text-align: left !important; }
                            #gt_custom_choice::-ms-input-placeholder { text-align: left; }
                            #gt_custom_choice::-webkit-input-placeholder { text-align: left; }
                            #gt_custom_choice:-moz-placeholder { text-align: left; }
                            #gt_custom_choice::-moz-placeholder { text-align: left; }'))),
  
  # The dynamically-generated user panel
  uiOutput("userpanel"),
  
  div(id = "sidebardiv", 
  sidebarMenu(id = "sidebar",
    menuItemAddId( menuItem("Home", tabName = "hometab", icon = icon("home")), id="li_home" ),
    menuItemAddId( menuItem("Search", tabName = "searchtab", icon = icon("search")), id="li_search" ),
    menuItemAddId( menuItem("Analysis", tabName = "analysistab", icon = icon("pie-chart")), id="li_analysis" ),
    menuItemAddId( menuItem("MFA", tabName = "mfatab", icon = icon("line-chart")), id="li_mfa" ),
	  menuItemAddId( menuItem("BER Analysis", tabName = "bertab", icon = icon("user", lib = "glyphicon")), id="li_ber" ),
    menuItemAddId( menuItem("GeneTox21 Analysis", tabName = "genetoxtab", icon = icon("chart-bar")), id="li_geneTox" ),
    menuItemAddId( menuItem("Save/Load", tabName = "savetab", icon = icon("floppy-o")), id="li_save" ),
	  menuItemAddId( menuItem("Import GeneTox Table", tabName = "importGTTable", icon = icon("file-import")), id="li_importGTTable"), #YANIC
    menuItemAddId( menuItem("Help", icon = icon("question-circle"), startExpanded = FALSE,
                            actionLink(inputId = "link_help_overview", label = "Interface Help"),
                            actionLink(inputId = "link_help_search", label = "Search Help"),
                            actionLink(inputId = "link_help_save", label = "Data Backup Help"),
                            actionLink(inputId = "link_help_analysis", label = "Analysis Help")
                            ),
                   id = "li_help" ),
    menuItemAddId( menuItem("Quick Options", icon = icon("cog"), startExpanded = FALSE,
             radioButtons(inputId = "radio_listtype", label = "Chemical display options:", selected = "name",
                          choices = list("By CASN" = "casn",
                                         "By Name" = "name") ),
             #radioButtons(inputID = "radio_plotorientation", label = "Plot orientation options:", selected = "vertical",
              #            choices = list("Horizontal" ="horizontal",
              #                           "Vertical" = "vertical")),
             
             # This is removed, since SQL DATABASE WILL NOT BE USED, AT LEAST FOR NOW. #YANIC
             hidden(radioButtons(inputId = "radio_dbtype", label = "Database options:", selected = "rdata", width = "100%",
                          choices = list("MySQL: dreamtk_db (v0.7, based on httk-1.7, tcpl-v2)" = "mysql",
                                         "RData: DreamTKv0.8.RData (based on httk-1.7, tcpl-v2)" = "rdata") )),
             hidden(uiOutput(outputId="ui_database_status"))
             ), 
             id="li_quickoptions" )
  ))
);

#body is the main dashboard are with all the required tabs
body <- dashboardBody(
  
  shinyjs::useShinyjs(),
  includeCSS( "./www/progressbar.css"),
  introjsUI(),
  
  tabItems(
    tabItem(tabName = "hometab",
            h2("Home"),
            fluidRow( div(id="welcomeboxdiv",
              box(status = "primary", title = "Welcome", collapsible = TRUE, width = 6,
                  h5("Welcome to DreamTK, an R application which facilitates toxicokinetic analysis of a variety of chemicals."),
                  h5("ToxCast Pipeline for High-Throughput Screening Data (tcpl v2.0) (Filer et al., 2016, US EPA) is the  primary chemical and assay database used by this application."),
                  h5("High-Throughput Toxicokinetics (httk v1.8) (Pearce et al., 2018) database is used to obtain the necessary toxicokinetic constants."),
				  h5("Only Assays that hit their cutoffs are considered for analysis.")
              ))
            ),
            fluidRow(
              box(status = "primary", title = "How To", collapsible = TRUE, width = 6,
                  h5("The following tutorials can help you familiarize with the application. These are also available from the Help tab on the Navigation menu."),
                  actionButton(inputId = "button_tour_overview", label = "Interface Tutorial"),
                  actionButton(inputId = "button_tour_search", label = "Search Tutorial"),
                  actionButton(inputId = "button_tour_save", label = "Data Backup Tutorial"),
                  actionButton(inputId = "button_tour_analysis", label = "Analysis Tutorial")
              )
            ),
            fluidRow(
              box(status = "primary", title = "News and Announcements", collapsible = TRUE, width = 6,
                  includeHTML("./www/news.html")
              )
            )
    ),
    
    tabItem(tabName = "searchtab",
            h3("Search Chemical Database"),
            fluidRow(
              
              tabBox(title = tagList(shiny::icon("compass")), id = "searchtabset",
                     
                     #search panel
                     tabPanel("Search", icon = icon("search"), fluid = TRUE,
                              fluidRow(
                                  column(12,
                                         div(
                                           style = "display: inline-block;vertical-align:baseline; width: 90%;",
                                           selectizeInput(inputId = "field_search", label = "", width = "100%",
                                                          options = list(placeholder = "Enter one or more chemicals separated by spaces",
                                                                         create = TRUE, createOnBlur = TRUE, createFilter = "^[\\w-,()\\[\\]]+$", persist = FALSE,
                                                                         maxOptions = 1000, loadThrottle = 800, openOnFocus = FALSE,
                                                                         delimiter = " ", hideSelected = TRUE, closeAfterSelect = TRUE,
                                                                         plugins = list("restore_on_backspace", "remove_button")),
                                                          multiple = TRUE, choices = NULL)
                                         ),
                                         div(
                                           style = "display: inline-block;vertical-align:70%; width: 9%;",
                                           actionButton(inputId = "button_search", label = "", icon = icon("search"), style="color:#0e76b7"),
                                           bsTooltip(id="button_search", "Search", placement = "bottom", trigger = "hover", options = NULL)
                                         ),
                                         radioButtons(inputId = "radio_searchtype", label = "Search options:", inline = TRUE, selected = "casn",
                                                      choices = list("By CASN" = "casn",
                                                                     "By Name" = "name") )
                                  )
                              )
                     ),
                     
                     #custom file panel
                     tabPanel("Custom File", icon = icon("table"), fluid = TRUE,
                              fluidRow(
                                  column(10, 
                                         fileInput(inputId = "file_customchem", label = "Choose CSV File", accept = c("csv", ".csv"), width = "100%"),
                                         uiOutput(outputId="ui_customchems_status")
                                  ),
                                  column(2,
                                         br(),
                                         actionButton(inputId = "button_load_customchems", label = "Parse file", icon = icon("folder-open-o"), 
                                                      style="padding-left:10px; padding-right:10px; padding-top:10px; padding-bottom:10px; white-space: normal;")
                                  )
                              ),
                              fluidRow(
                                column(12,
                                       br(),
                                       actionLink(inputId = "link_customchem_hint", label = "File requirements hint", icon = icon("hand-o-up")),
                                       hidden(fluidRow(id="panelCustomchemHint",
                                                       column(12,
                                                              h4("Template:"),
                                                              downloadLink(outputId = "button_customchem_template",label = tagList(shiny::icon("table"), "Download CSV Template"), 
                                                                           style="padding-left:10px; padding-right:10px; width:100%; white-space: normal; font-size:20px"),
                                                              HTML( 
                                                                str_c( "<p>Comma-separated value (.csv) file with first row containing column names. ",
                                                                       "<strong>'casn'</strong> column is mandatory. Other columns are optional, and if missing, will be looked up in the database.</p>",
                                                                       "<strong>Table columns:</strong>",
                                                                       "<blockquote style = 'font-size:14px; border-left: 10px solid #fff;'><strong>casn:</strong> CAS Registry Number",
                                                                       "<br><strong>name:</strong> Chemical name",
                                                                       "<br><strong>cytotoxicity_um:</strong> Cytotoxic concentration (<strong>Units:</strong> uM)",
                                                                       "<br><strong>cytotoxic:</strong> is the chemical considered cytotoxic at above concentration? Y/N flag",
                                                                       "<br><strong>mw:</strong> Molecular weight (<strong>Units:</strong> g/mol)",
                                                                       "<br><strong>human_funbound_plasma:</strong> Unbound fraction of chemical in blood (<strong>Units:</strong> none)",
                                                                       "<br><strong>human_clint:</strong> Intrinsic invitro hepatic clearance (<strong>Units:</strong> uL/min/10<sup>6</sup> cells)",
                                                                       "<br><strong>human_rblood2plasma:</strong> Blood to plasma concentration ratio (<strong>Units:</strong> none)",
                                                                       "<br><strong>log_kow:</strong> Octanol:Water partition coefficient at 25*C (<strong>Units:</strong> LOG10 value)",
                                                                       "<br><strong>pka:</strong> Equivalent chemical ionization constant (<strong>Units:</strong> - LOG10 value)</blockquote>"
                                                                )
                                                              )
                                                       )
                                       ))
                                )
                              )
                     )
                     
              )
              
            ),
            fluidRow(
              
              box(status = "primary", title = "Select and examine chemicals and chemical assays", collapsible = TRUE, width = 12, id = "chemlistbox",
                fluidRow(
					column(2,
					#if for some reason we change the colours primary won't work and we will have to figure a way to put colours in those switchtes. Style is not useable. https://rdrr.io/cran/shinyWidgets/src/R/input-pretty.R Gabriel
					prettySwitch(inputId = "select_hit", label = "Only Active Assays", value = TRUE, status = "primary",
									  fill = TRUE, bigger = TRUE, inline = TRUE,
									  width = NULL)
					),
					column(2,
					prettySwitch(inputId = "select_background", label = "Include Background Measurements", value = TRUE, status = "primary",
									  fill = TRUE, bigger = TRUE, inline = TRUE,
									  width = NULL)
					)
					
				),
				fluidRow(
                  column(12,
                         div(
                           style = "display: inline-block;vertical-align:bottom; width: 33%;",
                           selectInput(inputId = "select_chemical", label = "Searched chemical list", selectize = FALSE, size = 10,
                                       choices = NULL )
                         ),	
                         div(
                           style = "display: inline-block;vertical-align:bottom; width: 33%;",
                           selectInput(inputId = "select_assay", label = "Chemical assay list", selectize = FALSE, size = 10,
                                       choices = NULL )
                         ),
                         div(
                           style = "display: inline-block;vertical-align:bottom; width: 33%;",
                           selectInput(inputId = "select_assay_comp", label = "Assay component list", selectize = FALSE, size = 10,
                                       choices = NULL )
                         )
                  )

                ),
                fluidRow(
                  column(1,
                    actionButton(inputId = "button_delete_chem", label = "", icon = icon("trash-o"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;"),
                    bsTooltip(id="button_delete_chem", "Delete selected chemical", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(1,
                    actionButton(inputId = "button_delete_missing", label = "*", icon = icon("trash"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;"),
                    bsTooltip(id="button_delete_missing", "Delete missing chemicals", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(1,
                    actionButton(inputId = "button_list_missing", label = "*", icon = icon("list-ul"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;"),
                    bsTooltip(id="button_list_missing", "List missing chemicals", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(1,
                    actionButton(inputId = "button_clear_list", label = "", icon = icon("times"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;"),
                    bsTooltip(id="button_clear_list", "Clear chemical list", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(2,
                    downloadButton(outputId = "button_savechems", label = "Chem CSV", icon = icon("table"), 
                                   style="padding-left:10px; padding-right:10px; width:95%; white-space: normal;"),
                    bsTooltip(id="button_savechems", "Save information for all listed chemicals as a CSV file", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(2,
                    downloadButton(outputId = "button_saveassays", label = "Assays CSV", icon = icon("table"), style="padding-left:10px; padding-right:10px; width:95%; white-space: normal;"),
                    bsTooltip(id="button_saveassays", "Save assay information for the selected chemical as a CSV file", placement = "bottom", trigger = "hover", options = NULL)
                  )
                )							   
              )
              
          ),
          fluidRow(
            #chemical and assay information panels
            box(status = "primary", title = "Chemical info", collapsible = TRUE, width = 6,
                htmlOutput(outputId = "html_chemicalinfo")
            ),
            box(status = "primary", title = "Assay info", collapsible = TRUE, width = 6,
                htmlOutput(outputId = "html_assayinfo")
            )
          )
    ),
    
    tabItem(tabName = "analysistab",
            h3("Analyse chemicals"),
            #chemical selection panel
            fluidRow(
              box(status = "primary", title = "Select chemicals and desired analysis", collapsible = TRUE, width = 12,
                column(6,
                       wellPanel(id = "stats_selectcontrol",
                         selectizeInput(inputId = "select_chemical_stats", label = "Selected chemicals", 
                                        options = list(placeholder = "Click me to select chemicals",
                                                       maxOptions = 10000, loadThrottle = 800,
                                                       delimiter = " ", hideSelected = TRUE,
                                                       plugins = list("restore_on_backspace", "remove_button")),
                                        multiple = TRUE, choices = NULL),
                         fluidRow(
                           column(3,
                                  actionButton(inputId = "button_stats_selectall", label = "Select All", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;overflow:hidden")
                           ),
                           column(3,
                                  actionButton(inputId = "button_stats_deselectall", label = "Deselect All", icon = icon("ban"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;overflow:hidden")
                           )
                         )
                       ),
					   prettySwitch(inputId = "analyse_background", label = "Include Background Measurements", value = TRUE, status = "primary",
									  fill = TRUE, bigger = TRUE, inline = TRUE,
									  width = NULL),
					   prettySwitch(inputId = "include_GenTox", label = "Include GeneTox21 Data", value = FALSE, status = "primary",
					                fill = TRUE, bigger = TRUE, inline = TRUE,
					                width = NULL)
                ),
                column(6,
                       wellPanel(id = "stats_optionscontrol",
                                 fluidRow(
                                 column(6,
                         checkboxGroupInput(inputId = "checkbox_stats", label = "Select statistics", 
                                            choices = c("Target Family Counts and ac50 values" = "tfcounts",
                                                        "Hierarchical cluster heatmap of Target Subfamily activities" = "tfhm",
                                                        "Hierarchical cluster heatmap of Assay Endpoint activities" = "assayhm",
                                                        "Ac50 vs ScalarTop" = "scalartop_ac50",
                                                        "OED vs ScalarTop" = "scalartop_oed",
														"Burst Assay vs Not Burst Assay" = "ac50_box",
                                                        "Chemical ToxPI Plots (Individual)" = "toxpi",
                                                        "Chemical ToxPI Plots (Cytotoxicity)" = "toxpi2",
                                                        "Chemical ToxPI Plots (Grouped)" = "toxpigroup"),
                                            selected = NULL)
														),
														conditionalPanel("input.include_GenTox == true",
														column(6,
														       checkboxGroupInput(inputId = "checkbox_GenTox", label = "", 
														                          choices = c("GeneTox Data Table" = "GTDataTable",
														                                      "BMD Data Table" = "BMDDataTable"),
														                          selected = NULL))
														)),
						#those are magic numbers. Do not change them.					
						fluidRow(
                           column(4,
								actionButton(inputId = "button_select_all_stats", label = "Select/Deselect all", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; white-space: normal; width:100%; display:inline-block; overflow:hidden")
                                  
                           ),
                           column(3, offset = 0.5,
                                actionButton(inputId = "button_stats_run", label = "Run Stats", icon = icon("bar-chart"), style="padding-left:10px; padding-right:10px; white-space: normal; width:100%; display:inline-block; overflow:hidden")  
                           )
                         ),
					    busyIndicator(text="Working...")
                       )
                )
              )
            ),
            fluidRow(
              box(status = "primary", title = "Analysis results", collapsible = TRUE, width = 12,
                uiOutput(outputId = "ui_stats")
              )
            )

    ),
    
    tabItem(tabName = "mfatab",
            h3("Multiple Factor Analysis"),
            #chemical selection panel
            fluidRow(
              box(status = "primary", title = "Select chemicals and desired analysis", collapsible = TRUE, width = 6,
                  wellPanel(id = "mfa_selectcontrol",
                    selectizeInput(inputId = "select_chemical_mfa", label = "Selected chemicals", 
                                   options = list(placeholder = "Click me to select chemicals",
                                                  maxOptions = 10000, loadThrottle = 800,
                                                  delimiter = " ", hideSelected = TRUE,
                                                  plugins = list("restore_on_backspace", "remove_button")),
                                   multiple = TRUE, choices = NULL),
                    fluidRow(
                      column(3,
                             actionButton(inputId = "button_mfa_selectall", label = "Select All", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(3,
                             actionButton(inputId = "button_mfa_deselectall", label = "Deselect All", icon = icon("ban"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(3,
                             actionButton(inputId = "button_mfa_run", label = "Run MFA", icon = icon("bar-chart"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(1,
                             busyIndicator(text="Working...")
                      )
                    )
                  )
              )
            ),
            fluidRow(
              box(status = "primary", title = "Analysis results", collapsible = TRUE, width = 12,
                  uiOutput(outputId = "ui_mfa")
              )
            )
    ),
	
	tabItem(tabName = "bertab",
            h3("Biological Exposure Ratio Analysis"),
            #chemical selection panel
            fluidRow(
              box(status = "primary", title = "Select chemicals and desired analysis", collapsible = TRUE, width = 6,
                  wellPanel(id = "ber_selectcontrol",
                    selectizeInput(inputId = "select_chemical_ber", label = "Selected chemicals", 
                                   options = list(placeholder = "Click me to select chemicals",
                                                  maxOptions = 10000, loadThrottle = 800,
                                                  delimiter = " ", hideSelected = TRUE,
                                                  plugins = list("restore_on_backspace", "remove_button")),
                                   multiple = TRUE, choices = NULL),
                    fluidRow(
                      column(3,
                             actionButton(inputId = "button_ber_selectall", label = "Select All", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(3,
                             actionButton(inputId = "button_ber_deselectall", label = "Deselect All", icon = icon("ban"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(3,
                             actionButton(inputId = "button_ber_run", label = "Run BER analysis", icon = icon("bar-chart"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(1,
                             busyIndicator(text="Working...")
                      )
                    )
                  )
				  
              ),
			  box(status = "primary", title = "Extra Information and Assumptions", collapsible = TRUE, width = 6, height = 250,
					p("The calculations are based on the", a( href = "https://pubs.acs.org/doi/10.1021/es502513w", "SHEDS-HT", style = "color: blue;", target = "_blank", rel = "noopener noreferrer"), " exposure model."),
					h4("Assumptions"),
					tags$ul(tags$li("Physical activity index  =  1.75"),
						tags$li("Basal Alveolar Ventilation Rate  =  15.7 m",tags$sup("3"),"/day"),
						tags$li("Vapor pressure  =  0.876 Pa")
						),
					
					h4("Warnings"),
					tags$ul(tags$li("Ac50 values are used as a surrogate for OED values as the vast majority of chemicals do not have any values for OED.")))
            ),
            fluidRow(
              box(status = "primary", title = "Analysis results", collapsible = TRUE, width = 12,
                  uiOutput(outputId = "ui_ber")
              )
            )
    ),
    #YANIC
	  tabItem(tabName = "genetoxtab",
	          h3("GeneTox21 BMD Analysis"),
	          fluidRow(
	            box(status = "primary", title = "Select chemicals and desired analysis", collapsible = TRUE, width = 12,
                fluidRow(
                  column(6,
                    wellPanel(id = "gt_selectcontrol",
                      selectizeInput(inputId = "select_chemical_gt", label = "Selected chemicals", 
                                   options = list(placeholder = "Click me to select chemicals",
                                                  maxOptions = 10000, loadThrottle = 800,
                                                  delimiter = " ", hideSelected = TRUE,
                                                  plugins = list("restore_on_backspace", "remove_button")),
                                   multiple = TRUE, choices = NULL),
                      fluidRow(
                        column(3,
                               actionButton(inputId = "button_gt_selectall", label = "Select All", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;overflow:hidden")
                        ),
                        column(3,
                               actionButton(inputId = "button_gt_deselectall", label = "Deselect All", icon = icon("ban"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;overflow:hidden")
                        )
                      ),
                      fillRow(height = 21)
                    )
                  ),
                  column(6,
                    wellPanel(id = "gt_graphType",
                      selectizeInput(inputId = "select_graphType_gt", label = "Selected Graph Type",
                        multiple = FALSE, choices = c(
                          "BMD Confidence Interval Graph" = 1, 
                          #"BMD vs BMR Graph" = 2, #This graph was not deemed useful. YANIC
                          "BMR Distribution Graph" = 3, 
                          "Hierarchically Clustered BMD Heatmap Graph" = 4, 
                          "BMD Grouped Confidence Interval Graph" = 5,
                          "Curve Fit Plots" = 6,
                          "GeneToxPi Plots" = 7,
                          "Endpoint Threshold Table" = 8
                        ), selected = 1),
                      fluidRow(
                        column(6,
                          conditionalPanel("input.select_graphType_gt == 1 || input.select_graphType_gt == 2 || input.select_graphType_gt == 3 || input.select_graphType_gt == 7",
                            strong("Force Scale"),
                            prettySwitch(inputId = "forceScale_switch_gt", label = "", value = FALSE, status = "primary",
                                         fill = TRUE, bigger = TRUE, inline=FALSE,
                                         width = 20),
                            )
                          ),
                        column(6,
                          conditionalPanel("input.select_graphType_gt == 1",
                             strong("Force Ordering"),
                             prettySwitch(inputId = "forceOrder_switch_gt", label = "", value = FALSE, status = "primary",
                                          fill = TRUE, bigger = TRUE, inline = FALSE,
                                          width = 20),
                          ),
                          conditionalPanel("input.select_graphType_gt != 1",
                            fillRow(height = 55)
                          )
                        )
                      )
                    )
                  )
                ),
                fluidRow(
                  column(6,
                    wellPanel(id = "gt_calculation_type",
                      fluidRow(
                        column(6,
                          conditionalPanel("input.select_graphType_gt == 1 || input.select_graphType_gt == 4 || input.select_graphType_gt == 5 || input.select_graphType_gt == 6 || input.select_graphType_gt == 7 || input.select_graphType_gt == 8",
                            conditionalPanel("input.select_graphType_gt != 6",
                               radioButtons(inputId = "gt_bmr_picker", label = "Select BMR Value", choices = c("Precalculated" = 1, "Custom" = 2), selected = 1, inline = TRUE)
                            ),
                            conditionalPanel("input.select_graphType_gt == 6 || input.gt_bmr_picker == 1",
                              fluidRow(
                                column(8,
                                  selectizeInput(inputId = "gt_precalculated_choices", label = "Precalculated Values",
                                    multiple = FALSE, choices = c("5%" = 5, "10%" = 10, "15%" = 15, "20%" = 20), selected = 5
                                  ),
                                  fillRow(height = 32)
                                )
                              ),
                            ),
                            conditionalPanel("input.select_graphType_gt != 6 && input.gt_bmr_picker == 2",
                              h6("Please note that custom values will be calculated on the spot, and might take some time."),
                              fluidRow(
                                column(8,
                                  HTML(
                                    '<div class="form-group shiny-input-container">
                                    <label for="gt_custom_choice">Custom Value</label>
                                    <input id="gt_custom_choice" type="text" style="text-align: right;" class="form-control" placeholder="Enter your BMR value" value="" maxlength="7"/>
                                    </div>'
                                  )
                                ),
                                column(4, style="margin-top: 32px; margin-left: -25px",
                                  p("%")
                                )
                              ),
                            )
                          ),
                          conditionalPanel("input.select_graphType_gt == 2 || input.select_graphType_gt == 3",
                            fillRow(height = 171)
                          ),
                          conditionalPanel("input.select_graphType_gt == 2 || input.select_graphType_gt == 3 || input.select_graphType_gt == 4 || input.select_graphType_gt == 7 || input.select_graphType_gt == 8",
                            fluidRow(
                              column(8,
                                selectizeInput(inputId = "gt_bmd_result_chooser", label = "BMD Result", multiple = FALSE,
                                  choices = c("BMDU" = "bmdu", "BMD" = "bmd", "BMDL" = "bmdl"), selected = "bmd")
                              )
                            )
                          )
                        ),
                        column(6,
                          conditionalPanel("input.select_graphType_gt != 6 && input.select_graphType_gt != 8",
                            strong("Activate Logarithmic Values"),
                            prettySwitch(inputId = "log_switch_gt", label = "", value = FALSE, status = "primary",
                                         fill = TRUE, bigger = TRUE, inline = FALSE,
                                         width = NULL),
                          ),
                          conditionalPanel("input.select_graphType_gt == 6 || input.select_graphType_gt == 8",
                                           fillRow(height = 55)
                          ),
                          conditionalPanel("input.select_graphType_gt != 6",
                            strong("Show Plotting Data"),
                            prettySwitch(inputId = "plottingData_switch_gt", label = "", value = FALSE, status = "primary",
                                         fill = TRUE, bigger = TRUE, inline = FALSE,
                                         width = NULL),
                          ),
                          conditionalPanel("input.select_graphType_gt == 6",
                            fillRow(height = 55)
                          ),
                          strong("Show Assay-Endpoint Description"),
                          prettySwitch(inputId = "description_switch_gt", label = "", value = FALSE, status = "primary",
                                       fill = TRUE, bigger = TRUE, inline = FALSE,
                                       width = NULL),
                          fillRow(height = 88)
                        )
                      )
                    )
                  ),
                  column(6,
                    wellPanel(id = "gt_comparaison_type",
                      conditionalPanel("input.select_graphType_gt == 1 || input.select_graphType_gt == 2 || input.select_graphType_gt == 3 || input.select_graphType_gt == 4 || input.select_graphType_gt == 5 || input.select_graphType_gt == 7 || input.select_graphType_gt == 8",
                        selectizeInput(inputId = "gt_modelComp", label = "Select Comparison Model",
                                     multiple = FALSE, choices = c(
                                       "Best model (by Loglikelihood)" = 1, 
                                       "Best model (by AIC)" = 2, 
                                       #"Linear Model" = 3, 
                                       "Exponential Curve Model" = 4, 
                                       "Hill Curve Model" = 5, 
                                       "Geometric Mean of Hill Curve and Exponential" = 6
                                      ), selected = 1),
                      ),
                      conditionalPanel("input.select_graphType_gt == 6",
                        fillRow(height = 79)
                      ),
                      conditionalPanel("input.select_graphType_gt == 1 || input.select_graphType_gt == 2 || input.select_graphType_gt == 3 || input.select_graphType_gt == 6 || input.select_graphType_gt == 7",
                        selectizeInput(inputId = "gt_comparaison", label = "Select Comparison type",
                          multiple = FALSE, choices = c("By chemicals" = 1, "By endpoints" = 2), selected = 1)
                      ),
                      conditionalPanel("input.select_graphType_gt == 8",
                        numericInput(inputId = "gt_threshold", label = "Threshold", value = 0.05)
                        ),
                      conditionalPanel("input.select_graphType_gt == 4 || input.select_graphType_gt == 5",
                        fillRow(height = 79)
                      )
                    ),
                    wellPanel(id = "gt_run",
                      actionButton(inputId = "button_gt_createPlots", label = "Run Analysis", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;overflow:hidden")
                    )
                  )
                )
              )
	          ),
	          fluidRow(
              box(status = "primary", title = "Analysis results", collapsible = TRUE, width = 12,
                uiOutput(outputId = "ui_gt")
              )
	          )
    ),
	
    tabItem(tabName = "savetab",
            h3("Save/Load workspace"),
            box(status = "primary", title = "Save", width = 3,
                       radioButtons(inputId = "radio_savetype", label = "Save options:", selected = "cm",
                                    choices = list("R Chemical object list" = "cm") ),
                       downloadButton(outputId = "button_savefile", label = "Save", style="padding-left:10px; padding-right:10px; white-space: normal;")
            ),
            box(status = "primary", title = "Load", width = 4,
                       radioButtons(inputId = "radio_loadtype", label = "Load options:", selected = "cm",
                                    choices = list("R Chemical object list" = "cm") ),
                       fileInput(inputId = "file_load", label = "Choose File"),
                       actionButton(inputId = "button_loadfile", label = "Load", style="padding-left:10px; padding-right:10px; white-space: normal;"),
                       busyIndicator(text="Working..."),
                       uiOutput(outputId="ui_load_status")
            )
    ),
  
	tabItem(tabName = "importGTTable",
	  h3("Import GeneTox Table"),
	  fluidRow(
	    box(status = "primary", title = "Type of table", collapsible = TRUE, width = 9,
	        h5("This section allows a change in the table to which this importer will add. All files imported will go to the table shown below."),
	        selectInput(inputId="tableDropDown", label = "", choices=c("GeneTox Table" = 1, "BMD Table" = 2), multiple = FALSE),
	        br(),
	        actionLink(inputId = "link_importGeneTox_hint", label = "File requirements hint", icon = icon("hand-o-up")),
	        hidden(fluidRow(id="panelImportGeneToxHint",
              column(12,
                 #h4("Template:"),
                 # downloadLink(outputId = "button_customchem_template",label = tagList(shiny::icon("table"), "Download CSV Template"), 
                 #              style="padding-left:10px; padding-right:10px; width:100%; white-space: normal; font-size:20px"),
                 HTML(
                   str_c(
                     "<strong>Required format:</strong>",
                     "<blockquote style = 'font-size:14px; border-left: 10px solid #fff;'>.csv (Comma Delimited)",
                     "<br>.txt (Tab Delimited)</blockquote>"
                   )
                 ),
                 conditionalPanel("input.tableDropDown == 1",
                   HTML( 
                     str_c( "<strong>Table columns:</strong>",
                            "<blockquote style = 'font-size:14px; border-left: 10px solid #fff;'><strong>Experiment:</strong> Assay, followed by ?S9 (ex: MultiFlow -S9)",
                            "<br><strong>Compound</strong> Chemical name",
                            "<br><strong>Concentration (uM):</strong> Concentration value in micromolar",
                            "<br><strong>*Any other columns*:</strong> Will be an endpoint if the name doesn't start with a letter",
                            "<br>Last column will be an endpoint if all column names start with a letter</blockquote>"
                      )
                    )
                 ),
                 conditionalPanel("input.tableDropDown == 2",
                   HTML( 
                     str_c( "<strong>Table columns:</strong>",
                            "<blockquote style = 'font-size:14px; border-left: 10px solid #fff;'><strong>Compound:</strong> Chemical name",
                            "<br><strong>Assay:</strong> Assay name",
                            "<br><strong>Endpoint:</strong> Endpoint name",
                            "<br><strong>S9:</strong> + or -, depending on the s9 value",
                            "<br><strong>BMR:</strong> BMR Value, in decimal (ex: 5% = 0.05)",
                            "<br><strong>EXP_loglik:</strong> Loglikelihood of the exponential model",
                            "<br><strong>EXP_AIC:</strong> AIC value of the exponential model",
                            "<br><strong>EXP_BMDL:</strong> BMDL value of the exponential model",
                            "<br><strong>EXP_BMD:</strong> BMD value of the exponential model",
                            "<br><strong>EXP_BMDU:</strong> BMDU value of the exponential model",
                            "<br><strong>HILL_loglik:</strong> Loglikelihood of the hill curve model",
                            "<br><strong>HILL_AIC:</strong> AIC value of the hill curve model",
                            "<br><strong>HILL_BMDL:</strong> BMDL value of the hill curve model",
                            "<br><strong>HILL_BMD:</strong> BMD value of the hill curve model",
                            "<br><strong>HILL_BMDU:</strong> BMDU value of the hill curve model",
                            "<br><strong>GEO_BMDL:</strong> BMDL value of the geometric mean between the exponential and hill curve model",
                            "<br><strong>GEO_BMD:</strong> BMD value of the geometric mean between the exponential and hill curve model",
                            "<br><strong>GEO_BMDU:</strong> BMDL value of the geometric mean between the exponential and hill curve model</blockquote>"
                     )
                   )
                 )
              )
	        ))
	    ),
	  ),
	  fluidRow(
	    box(status = "primary", title = "Single Table Import", collapsible = TRUE, width = 9,
	        h5("This section allows a new file to be imported into DreamTK. Once selected, the table will be seen in the window below."),
	        actionButton(inputId = "button_chooseFile", label="Choose File", style = "padding-left:10px; padding-right:10px; white-space: normal;"),
	        textOutput(outputId = "out_FilePath"),
          fluidRow(
            box(status = "primary", title = "View Table", collapsible = TRUE, width = 12,
              div(style = 'overflow-x: auto; overflow-y: auto; max-height: 80vh',
                tableOutput("table_out")
              )
            )
          ),
	        hidden(actionButton(inputId = "button_importTable", label = "Import", style = "padding-left:10px; padding-right:10px; white-space: normal;"))
      )
	  ),
	  fluidRow(
	    box(status = "primary", title = "Bulk Import", collapsible = TRUE, width = 9,
	        h5("This section allows multiple files to be imported into DreamTK at once. You will only need to select a folder. Please note that all valid files will be imported."),
	        actionButton(inputId = "button_chooseFolder", label = "Choose Folder", style = "padding-left:10px; padding-right:10px; white-space: normal;"),
	        textOutput(outputId = "out_FolderPath"),
	        fluidRow(
	          box(status = "primary", title = "View Files to import", collapsible = TRUE, width = 12,
	              htmlOutput(outputId = "out_filesToLoad")
	          )
	        ),
	        hidden(actionButton(inputId = "button_importFiles", label = "Import", style = "padding-left:10px; padding-right:10px; white-space: normal;"))
      )
	  )
	)
	
	)
);




# main ui function
ui <- dashboardPage( header, sidebar, body, skin = "blue" );