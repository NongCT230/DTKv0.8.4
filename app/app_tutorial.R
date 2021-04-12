steps_basics <- reactive(tribble(~step, ~element, ~intro,
                                 1, NA, "<h4>Welcome to the application Interface overview.</h4> 
                                         <p> This overview will help you familiarize with the basic application interface and may take several minutes to complete. </p>
                                         <p>To navigate the tutorial: <br>- Use Next and  Back buttons, <br>- Left, Right, Enter keyboard keys, or <br>- click on the bulleted steps. 
                                         <br>Press Skip or Esc key to exit.</p>
                                         <i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i> icon indicates an important point.
                                         <br><i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> icon indicates an interactive opportunity.",
                                 
                                 2, NA, "<i class='fa fa-exclamation-triangle' style='color:red; font-size:24px'></i> 
                                         <p style='font-size:16px'><i>Please follow instructions given in each step for proper operation of the tutorial</i></p>",
                                 
                                 3, ".wrapper", "<h4>Application window</h4> 
                                                 <p>The application consists of the header, sidebar, and the body.</p>",
                                 
                                 4, ".main-header", "<h4>Application header</h4> 
                                                     <p>Header shows the application name and version. 
                                                     Please note the version number when you save any work since a different version may have changes in features and output.</p>",
                                 
                                 5, ".sidebar-toggle", "<h4>Sidebar toggle button</h4>
                                                       <p>This button toggles the sidebar (Navigationon Menu).</p>
                                                       <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>                                                    
                                                       Try pressing it now.
                                                       <p> </p>
                                                       <i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i>
                                                       <i>Keep the sidebar open for the next step.</i>",
                                                       
                                 6, ".sidebar-menu", "<h4>Sidebar Navigation Menu</h4>
                                                     <p>This is the sidebar Navigation Menu. It allows you to access all of the application components.</p>
                                                     <p>Home tab is currently active, as indicated by a blue line on the left of the button.</p>
                                                     <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>                                                    
                                                     Try pressing the menu items and observe background changes.",
                                 
                                 7, "#li_home", "<h4>Home tab</h4>
                                                 <p> The first Navigation menu item is the Home tab </p>
                                                 <p> It provides general application information and announcements.</p>",
                                 
                                 8, "#li_search", "<h4>Search tab</h4>
                                                   <p>Search tab allows search of database chemicals or import of custom chemicals.</p> 
                                                   <p>Chemical and assay selection lists allows you to examine chemical and assay information.</p>",
                                 
                                 9, "#li_analysis", "<h4>Analysis tab</h4>
                                                     <p>Analysis tab allows you to compute and plot a variety of statistics for given chemical(s).</p>",
                                 
                                 10, "#li_mfa", "<h4>MFA tab</h4>
                                                 <p>MFA tab is responsible for running multifactor analysis for given chemical(s).</p>",
												 
								 11, "#li_ber", "<h4>BER Analysis tab</h4>
                                                 <p>BER Analysis tab is responsible for running biological exposure ratio analysis for given chemical(s).</p>",
                                 
                                 12, "#li_save", "<h4>Save/Load tab</h4>
                                                   <p>Save the workspace chemicals for loading later using the Save/Load tab.</p> 
                                                   <p>Loading saved chemicals avoids search and model re-computation.</p>",
                                 
                                 13, "#li_help", "<h4>Help links</h4>
                                                  <p>Help facilitates access a variety of help topics and tutorials.</p>",
                                 
                                 14, "#li_quickoptions", "<h4>Quick options</h4>
                                                           <p>Quick options allows for quick access to chemical display and datatabse options from anywhere within the application.</p>",
                                 
                                 15, "#shiny-tab-hometab", "<h4>Application body - Home tab</h4>
                                                           <p>This is the application body where content is displayed.</p>
                                                           <p>Content depends on which sidebar menu item is selected.</p>
                                                           <p>Since Home menu item is selected, the Home tab content is displayed.</p>",
                                 
                                 16, ".col-sm-6", "<h4>Content box</h4>
                                                   <p>This is a content box. Such boxes are used to display a variety of information.",
                                 
                                 17, ".box-tools.pull-right", "<h4>Box minimization</h4>
                                                               <p>Some boxes have a Collapse button to hide the contents. 
                                                               This can help manage the workspace by controlling the display of information.</p>.
                                                               <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>                          
                                                               Try clicking it.",
                                 
                                 18, "#li_quickoptions", "<h4>Quick Options</h4>
                                                           <p>This is the Quick Options rollout.</p>
                                                           <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>                          
                                                           Click on it if it is not already expanded.
                                                           <br>
                                                           <p>The rollout allows adjustment of chemical display and available databases.</p>",
                                 
                                 19, "#li_quickoptions", "<h4>Quick Options</h4>
                                                           <p><i>Chemical display options</i> control how chemicals are displayed in chemical lists and search fields. 
                                                           Listing can be done by chemical name or CAS registry number.</p>
                                                           <p><i>Database options</i> control the source database. 
                                                           Application RData database is always available, but requires a complete initial load before chemicals can be searched.</p>
                                                           <p>MySQL database, if available, loads data on a per-search basis, 
                                                           resulting in faster initial load, but slightly slower individual searches depending on the number of chemicals searched.</p>
                                                           <i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i> RData database is always available
                                                           provided the application has been configured correctly by the administrator. Use it if other databases are unavailable.
                                                           Keep in mind that if you are not using the application from the web, it may not have the most up-to-date database.",
                                                           
                                 20, NA, "<h4>Conclusion</h4> 
                                           <p>This concludes the basics tutorial. Press Done or Esc key on your keyboard to finish.</p>"
));

events_basics <- reactive(list("onchange"=I("switch(this._currentStep){

                                                  case 0:
                                                  //fold all menus on start, goto home tab
                                                          if ($('#li_help').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                              $('#li_help').click();
                                                          };
                                                          if ($('#li_quickoptions').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                              $('#li_quickoptions').click();
                                                          };
                                                          $('#li_home').click();
                                                          break;
                                                  
                                                  case 5:
                                                          $('#li_home').click();
                                                          break;
                                                  
                                                  case 6:
                                                          $('#li_home').click();
                                                          break;
                                                  
                                                  case 7:
                                                          $('#li_search').click();
                                                          break;
                                                  
                                                  case 8:
                                                          $('#li_analysis').click();
                                                          break;
                                                  
                                                  case 9:
                                                          $('#li_mfa').click();
                                                          break;
														  
												  case 10:
                                                          $('#li_ber').click();
                                                          break;
                                                  
                                                  case 11:
                                                          $('#li_save').click();
                                                          break;
                                                  
                                                  case 12:
                                                          //$('#li_help').click();
                                                          break;
                                                  
                                                  case 13:
                                                          //$('#li_quickoptions').click();
                                                          break;
                                                  
                                                  case 14:
                                                          $('#li_home').click();
                                                          break;
                                                  
                                                  case 17:
                                                          if ($('#li_quickoptions').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                              $('#li_quickoptions').click();
                                                          };
                                                          break;
                                                  
                                                  case 18:
                                                          if ($('#li_quickoptions').parent().children('ul').attr('class') == 'treeview-menu'){
                                                              $('#li_quickoptions').click();
                                                          };
                                                          break;
                                                  
                                                  }")
                                 
  ));
  
  steps_search <- reactive(tribble(~step, ~element, ~intro,
                                   1, NA, "<h4>Welcome to the Search overview.</h4> 
                                           <p>This tutorial will help you familiarize with chemical search and information display and may take up to 10 minutes to complete.</p>
                                           <p>Please refer to the Interface tutorial for tutorial and general interface navigation.</p> ",
                                   
                                   2, "#shiny-tab-searchtab", "<h4>Search tab.</h4> 
                                                               <p>Search tab contains all of the chemical and assay search functions.</p>",
                                   
                                   3, ".nav-tabs-custom", "<h4>Search box.</h4> 
                                                           <p>This is the Search tab box which allows searching a chemical within the application database, or inputting custom chemicals for analysis.</p>",
                                   
                                   4, "#searchtabset", "<h4>Search box tabs.</h4> 
                                                        <p>These are the Search tabs which allow you to switch between database search and custom chemical import.</p>
                                                        <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>                                                    
                                                        Try switching tabs.",
                                   
                                   5, ".nav-tabs-custom", "<h4>Database Search tab.</h4> 
                                                            <p>Database Search tab allows you to find chemicals within the database.</p>",
                                   
                                   6, "#radio_searchtype", "<h4>Search options.</h4> 
                                                            <p>These radio buttons control whether search is done by chemical CAS registry number or chemical name.</p>
                                                            <p>Switching search type briefly disables the search in order to update the chemical autocomplete options corresponding to the selection.</p>
                                                            <p>For the purpose of this tutorial, chemicals will be searched by Name.</p> ",
                                   
                                   7, ".nav-tabs-custom > div > div > div > div > div:nth-child(1) > div", "<h4>Search field.</h4> 
                                                                                                            <p>This text field accepts chemical search terms.</p>",
                                   
                                   8, "#button_search", "<h4>Search button.</h4> 
                                                         <p>This button starts retrieval of chemical information and model calculations for chemicals in the search field.</p>",
                                   
                                   9, ".nav-tabs-custom", "<h4>Search tab.</h4> 
                                                           <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>                                                    
                                                           Try searching a chemical by name with 'By Name' selected in the Search options.
                                                           <br>
                                                           <i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i>
                                                           Keep in mind that pressing <i>Enter</i> or <i>Left</i> and <i>Right</i> direction keys during the tutorial advances the tutorial, use mouse to select items.
                                                           <p> </p>
                                                           <p>For example start typing the first 3 letters of <i><strong>ris</strong>peridone</i>. As you type the partial chemical matches will show up in the drop-down list below.</p>
                                                           <p>The list can be scrolled using the scrollbar or by tapping <i>Up</i> and <i>Down</i> keys on your keyboard.</p>
                                                           <p>Clicking on the list item or pressing <i>Enter</i> selects the chemical from the list. Clicking outside of the list or pressing <i>Spacebar</i> enters the chemical as typed without selecting from the list.</p>",
                                   
                                   10, ".nav-tabs-custom", "<h4>Search tab.</h4> 
                                                           <p>Once chemical has been entered or selected it appears as a selectable item.</p>
                                                           <p>This item can be deleted by pressing the <b>x</b> icon next to it, or clicking on the item and pressing <i>Delete</i> or <i>Backspace</i> on your keyboard.</p>
                                                           <p>Clicking on empty space in the search field allows you to continue entering additional chemicals, or to use <i>Backspace</i> key to edit the item.</p>
                                                           <p>The cursor can be positioned between chemical items by using <i>Left</i> and <i>Right</i> keys on your keyboard. 
                                                           (<i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i>This also controls the tutorial!)</p>",
                                   
                                   11, ".nav-tabs-custom", "<h4>Search tab.</h4> 
                                                           <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>                                                    
                                                           You can continue adding chemicals. Press the Search (<i class='fa fa-search' style='color:blue; font-size:14px'></i>)
                                                            button on the right of the search field to begin the chemical information search.",
                                   
                                   12, "#chemlistbox", "<h4>Chemical and assay selection.</h4> 
                                                        <p>Once search is successfully completed, chemicals can be individually examined by selecting items in these list boxes.</p>",
                                   
                                   13, "#select_chemical", "<h4>Searched chemical list.</h4> 
                                                        <p>The chemicals you have searched will appear in this list.</p>
                                                        <p>Any chemical for which absolutely no information was found will be marked by an asterisk <b>*</b></p>
                                                        <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> If your search list contains any chemicals without an asterisk, select one now.",
                                   
                                   14, "#select_assay", "<h4>Chemical assay list</h4> 
                                                        <p>When a chemical is selected from the list, assay selection list will show all of the associated assays.</p>
                                                        <p>This list will show NA if no assays were located in the database.</p>
                                                        <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> If your assay list contains any assays, select one now.",
                                   
                                   15, "#select_assay_comp", "<h4>Chemical assay component selection.</h4> 
                                                        <p>When an assay is selected from the assay list, assay component list will show all of the associated assay components.</p>
                                                        <p>This list will show NA if no assay components for the assay were located in the database.</p>
                                                        <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> If your assay component list contains any assay components, select one now.",
                                   
                                   16, "#html_chemicalinfo", "<h4>Chemical information.</h4> 
                                                        <p>Once a chemical is selected from the Searched chemical list, information about the chemical will appear here.</p>
                                                        <p>Any chemical information that was not found in the database or could not be computed will appear as NA.</p>
                                                        <p>Models used and critical parameters are be shown. Missing model parameters and assumptions are identified as well.</p>",
                                   
                                   17, "#html_assayinfo", "<h4>Assay information.</h4> 
                                                        <p>Once an assay endpoint is selected from the Assay component list, information about the assay component will appear here.</p>
                                                        <p>Any assay information that was not found in the database or could not be computed will appear as NA.</p>",
                                   
                                   18, "#chemlistbox", "<h4>List management and export.</h4> 
                                                        <p>Several list management and export options are available via the buttons under the chemical and assay list boxes.</p>
                                                        <p>Hovering the mouse cursor over the buttons provides tooltips on their functionality</p>",
                                   
                                   19, "#button_delete_chem", "<h4>Delete selected chemical.</h4> 
                                                        <p>Permanently deletes a selected chemical from the list. A new search or import will need to be done to re-add it to the list.</p>",
                                   
                                   20, "#button_delete_missing", "<h4>Delete missing chemicals.</h4> 
                                                        <p>Permanently deletes missing chemicals from the list. These are chemicals for which absolutely no information was found in the database.</p>",
                                   
                                   21, "#button_list_missing", "<h4>List missing.</h4> 
                                                        <p>Creates a comma-separated list of missing chemicals, which can be copied</p>",
                                   
                                   22, "#button_clear_list", "<h4>Clear list.</h4> 
                                                        <p>Permanently deletes ALL chemicals from the list. A new search or import will need to be done to re-add it to the list. Useful for examining independent chemical lists. Don't forget to have the list saved before clearing!</p>",
                                   
                                   23, "#button_savechems", "<h4>Save chemical information.</h4> 
                                                        <p>Save <i>information for all of the chemicals</i> in the list as a comma-separated table file (.csv).</p>",
                                   
                                   24, "#button_saveassays", "<h4>Save assay information.</h4> 
                                                        <p>Save information for <i>all assay endpoints for the selected chemical</i> as a comma-separated table file (.csv).</p>",
                                   
                                   25, ".nav-tabs-custom", "<h4>Search from file.</h4> 
                                                            <p>In addition to searching the database via a text field, a file with the list of chemicals can be supplied.</p>
                                                            <p>Custom File tab provides this functionality.</p>",
                                   
                                   26, ".nav-tabs-custom > div > div:nth-child(2) > div:nth-child(1) > div.col-sm-10 > div.form-group.shiny-input-container > div.input-group", "<h4>Browse file.</h4> 
                                                            <p>Browsing for a .csv file of the correct format will upload the file to the application.</p>",
                                   
                                   27, "#button_load_customchems", "<h4>Parse file.</h4> 
                                                                    <p>Pressing this button will verify the file, fill in any missing information from the database and load chemicals into the Search chemical list.</p>
                                                                    <p>The advantage of using this method of search is that custom model parameters can be specified. 
                                                                    This is particularly useful if only chemical activity data, but no model data is availabe in the database.</p>",
                                   
                                   28, "#link_customchem_hint", "<h4>File requirements hint and template.</h4> 
                                                            <p>File requirements hint provides the expected format for the chemical list file.</p>
                                                            <i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i> <i>casn</i> is the only mandatory column. Other columns do not even have to be present in the file.",
                                   
                                   29, "#button_customchem_template", "<h4>Template.</h4> 
                                                            <p>The template can be obtained by clicking the provided link.</p>",
                                   
                                   30, NA, "<h4>Conclusion.</h4> 
                                           <p>This concludes the search tutorial. Press Done or Esc key on your keyboard to finish.</p>"
  ));
  
  events_search <- reactive(list("onchange"=I("switch(this._currentStep){

                                                    case 0:
                                                    //fold all menus on start, click on search tab
                                                            if ($('#li_help').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                                $('#li_help').click();
                                                            };
                                                            if ($('#li_quickoptions').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                                $('#li_quickoptions').click();
                                                            };
                                                            $('#li_search').click();
                                                            break;
                                                    
                                                    case 4:
                                                            $('#searchtabset').children().eq(0).children('a').click();
                                                            break;

                                                    case 5:
                                                            $('#radio_searchtype').children('div').eq(0).children('label').eq(1).children('input').eq(0).click();
                                                            break;

                                                    case 24:
                                                            $('#searchtabset').children().eq(1).children('a').click();
                                                            break;

                                                    case 26:
                                                            $('#link_customchem_hint').click();
                                                            break;
  }")
  ));
  
  steps_save <- reactive(tribble(~step, ~element, ~intro,
                                 1, NA, "<h4>Welcome to the Data Backup tutorial.</h4> 
                                           <p>This tutorial will help you familiarize with ways to save and load your chemical searches and data. It may take up to 10 minutes to complete.</p>
                                           <p>Please refer to the Interface tutorial for tutorial and general interface navigation.</p>
                                           <p>Please refer to the Search tutorial for searching and managing chemical data</p>",
                                 
                                 2, "#shiny-tab-searchtab", "<h4>Search tab.</h4> 
                                                               <p>Search tab contains all of the chemical and assay search functions.</p>
                                                               <p>It also has the ability to save chemical and assay information for chemicals in the Searched chemical list.</p>
                                                               <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> You can try searching a few chemicals to practice saving and loading.
                                                               <br>
                                                               <i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i>
                                                               Keep in mind that pressing <i>Enter</i> or <i>Left</i> and <i>Right</i> direction keys during the tutorial advances the tutorial, use mouse to select items.",
                                 
                                 3, "#chemlistbox", "<h4>Chemical and assay selection.</h4> 
                                                     <p>Once search is successfully completed, chemicals and associated information can be saved in a comma-separated table file which can be opened by many programs such as Microsoft Excel.</p>",
                                 
                                 4, "#button_savechems", "<h4>Save chemical information.</h4> 
                                                          <p>This button saves <i>information for all of the chemicals</i> in the list as a comma-separated table file (.csv).</p>
                                                          <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> You can try saving the info now. 
                                                          If no chemicals are in the list, an empty .csv file with chemical information headings will be produced.",
                                 
                                 5, ".nav-tabs-custom", "<h4>Load search from .csv file.</h4> 
                                                         <p>As discussed in the Search tutorial, a chemical list can be loaded into the search via the Custom File search tab.</p>
                                                         <p>Because the saved .csv chemical list contains a <i>casn</i> column, it meets the minimum file requirements for repeating the search.</p>
                                                         <p>The search can therefore be easily repeated without retyping the chemicals into the search field</p>",
                                 
                                 6, "#chemlistbox", "<h4>Chemical and assay selection.</h4> 
                                                     <p>To save assay information, a chemical needs to be selected first. Then all of the assay information for that chemical can be saved.</p>
                                                     <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> Select a chemical from the search list if any are available and advance the tutorial.",
                                 
                                 7, "#button_saveassays", "<h4>Save assay information.</h4> 
                                                            <p>This button saves information for <i>all assay endpoints for the selected chemical</i> as a comma-separated table file (.csv).</p>
                                                            <p>The file name is automatically set to chemical CAS registry number, but can be changed by typing a new name.</p>
                                                            <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> You can try saving the info now. 
                                                            If no assays are in the list, an empty .csv file with assay information headings will be produced.",
                                 
                                 8, "#shiny-tab-searchtab", "<h4>Search tab.</h4>
                                                      <p>Saving chemicals and assays to .csv files allows information to be used outside of this application.
                                                      However, re-entering the chemicals or lists of chemicals into the search or custom file may be time consuming.</p>
                                                      <p>This appication provides an ability to directly save the chemical data objects from the application workspace for rapid loading later.
                                                      All chemical, assay and model computation data is saved for each chemical.</p>",
                                 
                                 9, ".sidebar-menu", "<h4>Save/Load menu item</h4>
                                                     <p>Clicking on the Save/Load menu item in the Navigation bar opens the Save/Load tab.</p>",
                                 
                                 10, "#shiny-tab-savetab", "<h4>Save/Load tab</h4>
                                                     <p>The Save/Load tab contains the chemical saving and loading controls, 
                                                     which operate similarly to saving .csv files or loading custom chemical lists in the Search tab.</p>",
                                 
                                 11, "#radio_savetype", "<h4>Save options</h4>
                                                     <p>Save file options can be accessed from this radio button list. 
                                                     Presently there is only one option available (.rds file), but this may expand in future application versions.</p>",
                                 
                                 12, "#button_savefile", "<h4>Save button</h4>
                                                     <p>Selected file option can be saved by pressing this button.</p> 
                                                     <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> Try saving a file now for loading later.",
                                 
                                 13, "#radio_loadtype", "<h4>Load options</h4>
                                                     <p>Load file options can be accessed from this radio button list. 
                                                     Presently there is only one option available (.rds file), but this may expand in future application versions.</p>",
                                 
                                 14, "#shiny-tab-savetab > div.col-sm-4 > div > div.box-body > div:nth-child(2)", "<h4>Load file selection</h4>
                                                     <p>Use the file picker to upload the .rds file.</p> 
                                                     <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> Try uploading a file you saved earlier.",
                                 
                                 15, "#button_loadfile", "<h4>Load button</h4>
                                                     <p>This button loads and verifies the file before adding the containing chemicals to the list.</p> 
                                                     <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> Try loading a file you uploaded earlier.",
                                 
                                 16, "#shiny-tab-savetab", "<h4>File loading</h4>
                                                     <i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i> 
                                                     If a chemical with the same CAS registry number already exists in the chemical list, it will not be overwritten.",
                                 
                                 16, NA, "<h4>Conclusion.</h4>
                                          <p>This concludes the save/load tutorial. Press Done or Esc key on your keyboard to finish.</p>"
                                 
                                 
                                 
  ));
  
  events_save <- reactive(list("onchange"=I("switch(this._currentStep){

                                                    case 0:
                                                    //fold all menus on start
                                                    if ($('#li_help').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                        $('#li_help').click();
                                                    };
                                                    if ($('#li_quickoptions').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                        $('#li_quickoptions').click();
                                                    };
                                                    break;
                                                    
                                                    case 1:
                                                            $('#li_search').click();
                                                            break;

                                                    case 4:
                                                            $('#searchtabset').children().eq(1).children('a').click();
                                                            break;

                                                    case 8:
                                                            $('#li_save').click();
                                                            break;
                                                    
                                                    }")
  ));
  
  
  
  steps_analysis <- reactive(tribble(~step, ~element, ~intro,
                                     1, NA, "<h4>Welcome to the Analysis tutorial.</h4> 
                                           <p>This tutorial will help you familiarize with chemical analysis optinos. It may take up to 10 minutes to complete.</p>
                                           <p>Please refer to the Interface tutorial for tutorial and general interface navigation.</p>
                                           <p>Please refer to the Search tutorial for searching chemicals</p>",
                                     
                                     2, "#shiny-tab-analysistab", "<h4>Analysis tab.</h4> 
                                                                  <p>Analysis tab helps examine and compare several properties of interest for a single or multiple chemicals.</p>
                                                                  <p>Advancing the tutorial will take you to the Search tab to add chemicals to the chemical list for analysis.</p>",
                                     
                                     3, "#shiny-tab-searchtab", "<h4>Search tab.</h4> 
                                                                  <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i> 
                                                                  For the purpose of this tutorial, search by Name for <i>risperione</i>, <i>mevinphos</i> and <i>acrolein</i>. Advance the tutorial when ready.",
                                     
                                     4, "#shiny-tab-analysistab", "<h4>Analysis tab.</h4> 
                                                                  <p>Analysis interface consists of the chemical/analysis selection and results sections.</p>",
                                     
                                     5, "#stats_selectcontrol", "<h4>Chemical selecton.</h4> 
                                                                  <p>Chemicals can be selected by typing the chemicals of interest into the selection field or selecting them from the drop-down list. 
                                                                  The selection field allows selecting multiple chemicals,and operates in the same way as the search field in the Search tab.</p>
                                                                  <p>Only the chemicals in the Search list from the Search tab will be available.</p>",
                                     
                                     6, "#button_stats_selectall", "<h4>Select all.</h4> 
                                                                    <p>This button selects all of the chemicals in the list.</p>
                                                                    <i class='fa fa-exclamation-triangle' style='color:red; font-size:14px'></i>
                                                                    Please be patient when selecting over 500 chemicals. It may take over a minute to populate the selection list.",
                                     
                                     7, "#button_stats_deselectall", "<h4>Deselect all.</h4> 
                                                                    <p>This button deselects all of the chemicals in the list.</p>",
                                     
                                     8, "#stats_selectcontrol", "<h4>Chemical selecton.</h4> 
                                                                <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                Select all of the chemicals now",
                                     
                                     9, "#stats_optionscontrol", "<h4>Analysis options selecton.</h4> 
                                                                  <p>Analysis options are presented in checkboxes. Multiple can be active simulataneously.</p>
                                                                  <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                  Select <i>Target Family counts and ac50 values</i>",
                                     
                                     10, "#button_stats_run", "<h4>Run analysis.</h4> 
                                                                  <p>this button starts analysis and result output.</p>
                                                                  <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                  Press it now and advance the tutorial when you can see the charts.",
                                     
                                     11, "#shiny-tab-analysistab > div:nth-child(3) > div > div", "<h4>Analysis results.</h4> 
                                                                  <p>Analysis results populate the analysis box. They do not appear immediately since plotting functions take time plot, particularly if the data set is large.</p>
                                                                  <p>Use mousewheel, keyboard or right scrollbar to scroll the analysis frame, which may span many pages.</p>",
                                     
                                     12, "#shiny-tab-analysistab > div:nth-child(3) > div > div", "<h4>Analysis plots.</h4> 
                                                                <p>Analysis is in most cases rendered as interactive plots. Hovering over plot items may provide additional information.</p>",
                                     
                                     13, "#shiny-tab-analysistab > div:nth-child(3) > div > div", "<h4>Analysis plots.</h4> 
                                                                <p>Hovering over the plot activates a toolbar at the top left margin of the plot.</p>
                                                                <p>The toolbar buttons provide tooltips as to their functionality when activated. 
                                                                Not all functions may operate, depending on the plot type.</p>
                                                                <p>For example, clicking the <i class='fa fa-camera' style='color:grey; font-size:14px'></i> icon generates a .png image of the plot that can be saved.</p>",
                                     
                                     14, "#shiny-tab-analysistab > div:nth-child(3) > div > div", "<h4>Analysis table.</h4> 
                                                             <p>Analysis table is also generated. It allows direct examination of the computed data.</p>
                                                             <p>The table can be copied to clipboard, or saved as .csv or .xlsx files for processing outside of the app</p>
                                                             <p>A Search field in the upper right corner, and filter fields below each column allow to narrow down the data.</p>",
                                     
                                     15, "#stats_optionscontrol", "<h4>Analysis options selecton.</h4> 
                                                                  <p>Let us look at the plot navigation options for a scatter plot.</p>
                                                                  <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                  Deselect <i>Target Family counts and ac50 values</i>. Select <i>Ac50 vs ScalarTop</i> and re-run the stats.",
                                     

                                     16, "#shiny-tab-analysistab > div:nth-child(3) > div > div", "<h4>Analysis plot navigation.</h4> 
                                                                         <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                         clicking and dragging the mouse pointer on the plot zooms in on that plot region. Try it on  the messy regions of the plot.
                                                                         <br>
                                                                         <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                         doubleclicking on an empty plot space resets the view. Try it out after zooming in",
                                     
                                     17, "#shiny-tab-analysistab > div:nth-child(3) > div > div", "<h4>Analysis plot navigation.</h4> 
                                                                         <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                          Using <i>Compare data on hover</i> tool from the toolbar (next to <i class='fa fa-book' style='color:grey; font-size:14px'></i> 
                                                                          icon) enables examining and comparing the data points. It is particularly helpful in examining overlapping data points.
                                                                         <br>
                                                                         <p>When done, select <i>Show closest data on hover</i> tool to return to default view.</p>",
                                     
                                     18, "#shiny-tab-analysistab > div:nth-child(3) > div > div", "<h4>Analysis plot navigation.</h4> 
                                                                          <p>Clicking on legend items toggles the visibility of corresponding data points.</p>                                                                         
                                                                          <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                          Feel free to experiment with various tools. Advance the tutorial to the next step when ready.",
                                     
                                     19, "#stats_optionscontrol", "<h4>Analysis options selecton.</h4> 
                                                                  <p>Let us look at the ToxPI plots.</p>
                                                                  <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                  Deselect <i>Ac50 vs ScalarTop</i>, select <i>Chemical ToxPI plots</i> and re-run the stats.",
                                     
                                     20, "#shiny-tab-analysistab > div:nth-child(3) > div > div", "<h4>Analysis plots.</h4> 
                                                               <p>Some plots are not interactive. For example these ToxPI plots.</p>
                                                               <p>These plots use custom plotting routines which are not provided with interactive plotting libraries.</p>",
                                     21, "#shiny-tab-bertab", "<h4>Biological Exposure Ratio tab.</h4> 
                                                                  <p>This is the BER tab which performs Biological Exposure Ratio Analysis on the selected chemicals.</p>
                                                                  <p>The control elements of this tab are similar to the Analysis tab.</p>",
                                     
                                     22, "#ber_selectcontrol", "<h4>Chemical selection.</h4> 
                                                                  <p>chemical selection operates in exactly the same way as Analysis tab.</p>
                                                                  <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                  <i>Select all</i> chemicals, then <i>Run MFA</i>.",
                                     
                                     23, "#shiny-tab-bertab > div:nth-child(3) > div > div", "<h4>MFA results.</h4> 
                                                                  <p>BER results appear in the Analysis results box.</p>",
									 
									 
									 
                                     24, "#shiny-tab-mfatab", "<h4>Multiple Factor Analysis tab.</h4> 
                                                                  <p>This is the MFA tab which performs Multiple Factor Analysis on the selected chemicals.</p>
                                                                  <p>The control elements of this tab are similar to the Analysis tab.</p>",
                                     
                                     25, "#mfa_selectcontrol", "<h4>Chemical selection.</h4> 
                                                                  <p>chemical selection operates in exactly the same way as Analysis tab.</p>
                                                                  <i class='fa fa-hand-o-right' style='color:blue; font-size:14px'></i>
                                                                  <i>Select all</i> chemicals, then <i>Run MFA</i>.",
                                     
                                     26, "#shiny-tab-mfatab > div:nth-child(3) > div > div", "<h4>MFA results.</h4> 
                                                                  <p>MFA results appear in the Analysis results box. These plots are non-interactive as well since an R library script generates them.</p>",
                                     
                                     27, NA, "<h4>Conclusion.</h4>
                                              <p>This concludes the Analysis tutorial. Press Done or Esc key on your keyboard to finish.</p>"
  
  ));

  events_analysis <- reactive(list("onchange"=I("switch(this._currentStep){

                                                    case 0:
                                                    //fold all menus on start
                                                    if ($('#li_help').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                        $('#li_help').click();
                                                    };
                                                    if ($('#li_quickoptions').parent().children('ul').attr('class') == 'treeview-menu menu-open'){
                                                        $('#li_quickoptions').click();
                                                    };
                                                    break;
                                                    
                                                    case 1:
                                                            $('#li_analysis').click();
                                                            break;

                                                    case 2:
                                                            $('#li_search').click();
                                                            break;

                                                    case 3:
                                                            $('#li_analysis').click();
                                                            break;

															
													case 20:
                                                            $('#li_ber').click();
                                                            break;
                                                    case 23:
                                                            $('#li_mfa').click();
                                                            break;

                                                    }")
  ));
  