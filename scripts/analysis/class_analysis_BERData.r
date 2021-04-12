

# BER statistical data  --------------------------------------------


#v0.7

#class contains basic statistical data and data calculating functions
#https://pubs.acs.org/doi/suppl/10.1021/es502513w/suppl_file/es502513w_si_001.pdf
#the design goal with this is that this is part of a composition where the plotting classes need their data as an attribute.
#this will be faster because it erases the need to copy paste. This does not mean that a class can't have their specific data.

Class.Analysis.BERData <- R6Class("Class.Analysis.BERData",
  
  #private variables and functions
	#those are private because the data should only ever be calculated with the functions. This allows us to never change this. but you have to use getters.
  private = list(
	
	calc_BER_stat_tbl = NULL,
    BER_stat_tbl = NULL, #analysis table,
	Sheds_stat_tbl = NULL #shedsinformation table for analysis
  ),
  
  #public variables and functions
  public = list(
    
    #constructor
    initialize = function(chemicals = NULL) {	
    },
    
    #finalizer
    finalize = function() {
    },
    
    
    #builds a stat dataframe from a list of chemical objects
    buildBERStatsTable = function(BERData,sheds, chemlist, label_by){
		if(label_by == "casn"){ 
			private$BER_stat_tbl = filter(BERData, casn %in% chemlist);
        } else if (label_by == "name"){
			private$BER_stat_tbl = filter(BERData, name %in% chemlist);
        }
	
		private$Sheds_stat_tbl  = filter(sheds, product_use_category %in% private$BER_stat_tbl$product_use_category);

		invisible (private$BER_stat_tbl);
    },
    
	
 
    
    #computes the exposure quantities.
    computeBER = function(){
      
		  chemical_casn_list <- unique(private$BER_stat_tbl[["casn"]]);
		  BER_tbl = tribble(~casn,~name,~product_use_category,~direct_dermal,~direct_ingestion,~direct_vapor,~direct_aerosol);

		  for (chemical_casn in chemical_casn_list){

			
			
			BER_Calc_tbl <- filter(private$BER_stat_tbl, casn == chemical_casn) %>% 
			  select(casn, name, product_use_category, predicted_weight_fraction_mean );
			  
			BVA = 15.7; # 12.78? or 30.09600??? or 15.7????? or 17.36 if considering males???
			PAI = 1.75;
			vptreatedroomac = 0.876;
			vevap = 5;
			
			for ( puc in BER_Calc_tbl$product_use_category ){
				shedsrow = as.vector(filter(private$Sheds_stat_tbl, product_use_category == puc));
				puc_tbl = as.vector(filter(BER_Calc_tbl, product_use_category == puc));
				#we are doing the calculations in grams for now.
				if(shedsrow$direct_dermal == "1"){
					direct_dermal_val = (shedsrow$frequency/365)*(shedsrow$mass_per_use*puc_tbl$predicted_weight_fraction_mean)*(shedsrow$fret/100) * (shedsrow$fcont/100);
				}else{
					direct_dermal_val = 0;
				}
				#https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/respiratory-minute-volume for 1 of the 2 vars. BVA = 8.64 m3/day McMurray
				#for pai https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3958414/ ? idk what the pai is I'll just set to 1 for now. 1.75 according to Andy.
				if(shedsrow$direct_inhalation_of_vapor == "1"){
					cavap = (shedsrow$mass_per_use*puc_tbl$predicted_weight_fraction_mean) * vptreatedroomac/101325 * shedsrow$duration_of_direct_use * 1/vevap;
					direct_vapor_val =  (shedsrow$frequency/365)*cavap*shedsrow$duration_of_direct_use/1440*BVA*PAI;
				}else{
					direct_vapor_val = 0
				}
				
				if(shedsrow$direct_inhalation_of_aerosol == "1"){
					fairb = 0.95; #to check.
					caaer = (shedsrow$mass_per_use*puc_tbl$predicted_weight_fraction_mean) * 0.0625 * fairb;
					direct_aerosol_val = (shedsrow$frequency/365)*caaer*shedsrow$duration_of_direct_use/1440*BVA*PAI;
				}else{
					direct_aerosol_val = 0
				}
				
				if(shedsrow$direct_incidental_ingestion == "1"){
					direct_ingestion_val = (shedsrow$frequency/365)*(shedsrow$mass_per_use*puc_tbl$predicted_weight_fraction_mean)*(shedsrow$fing/100);
				}else{
					direct_ingestion_val = 0;
				}
				
				#Table of results in mg/kg/day				
				BER_tbl <- add_row(BER_tbl, casn = chemical_casn, 
							name = as.character(filter(private$BER_stat_tbl, casn == chemical_casn) %>% select(name) %>% distinct(name)),
							product_use_category = puc,
      				direct_dermal = (direct_dermal_val) * 10 ^ 6  / (10 ^ 3) / 70 ,  
      				direct_ingestion = (direct_ingestion_val ) * 10 ^ 6  / (10 ^ 3) / 70 ,
      				direct_vapor = (direct_vapor_val ) * 10 ^ 6 / (10 ^ 3) / 70,
      				direct_aerosol = (direct_aerosol_val )* 10 ^ 6 / (10 ^ 3) /70  
							);
  
				
				
				
		  }
		}
		private$calc_BER_stat_tbl = BER_tbl;
		invisible(private$calc_BER_stat_tbl);
	},
    
	getCalcBERStatsTable = function (){
      return (private$calc_BER_stat_tbl);
    },
    
	
    getBERStatsTable = function (){
      return (private$BER_stat_tbl);
    },
    getShedsStatsTable = function (){
      return (private$Sheds_stat_tbl);
    },    
    
	
	calcBERStatsDataExists = function (){
		
		
      if (is.null(private$calc_BER_stat_tbl) || nrow(private$calc_BER_stat_tbl) == 0){
		
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
	
    BERStatsDataExists = function (){
      if (is.null(private$BER_stat_tbl) || nrow(private$BER_stat_tbl) == 0){
        return (FALSE);
      } else {
        return (TRUE);
      }
    }
    
    
  )
)
