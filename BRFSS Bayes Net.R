# Load the bayes net functions
source("BayesianNetworks-template.r")

# Risk factors collected in the survey
cdc_survey_variables = list("income", "exercise", "smoke", "bmi", "bp", "cholesterol", "angina", 
                            "stroke", "attack", "diabetes")
survey_variables_states = c(8, 2, 2, 4, 4, 2, 2, 2, 2, 4)

# Health outcomes and habits 
health_outcomes = c("diabetes", "stroke", "attack", "angina")

habit_variables = c("smoke", "exercise")
bad_habit_values = c("1", "2")
good_habit_values = c("2", "1")

health_variables = c("bp", "cholesterol", "bmi")
bad_health_values = c("1", "1", "3")
good_health_values = c("3", "2", "2")

diabetes_codes = c("1", "2", "3", "4")
diabetes_values = c("Diabetic", "Diabetes only during pregnancy", "Not diabetic", "Pre-diabetic")
stroke_codes = c("1", "2")
stroke_values = c("Will have stroke", "Will not have stroke")
attack_codes = c("1", "2")
attack_values = c("Will have heart attack", "Will not have heart attack")
angina_codes = c("1", "2")
angina_values = c("Will have angina", "Will not have angina")

positive_health_outcomes_codes = c(diabetes_values[1], stroke_values[1], attack_values[1], angina_values[1])

outcome_codes = list(diabetes_codes, stroke_codes, attack_codes, angina_codes)
names(outcome_codes) = health_outcomes
outcome_values = list(diabetes_values, stroke_values, attack_values, angina_values)
names(outcome_values) = health_outcomes

# Load data collected by the CDC in the 2015 Behavioral Risk Factor Surveillance System (BRFSS) survey.
cdc_survey_data = read.csv(file="RiskFactors.csv", header=TRUE, sep=",")

## load libraries
load_libraries = function() {
  library(knitr)
}

## Change table headings for display purposes
## factor_table: factor tables that needs heading changed
## headings_to_change: headings to be changed
## new_headings: new headings
##
## Will return the table with changed headings
change_table_headings = function(factor_table, headings_to_change, new_headings) {
  
  if (!is.null(headings_to_change)) {
    for (heading_index in 1:length(headings_to_change)) {
      names(factor_table)[which(names(factor_table) == headings_to_change[heading_index])] = new_headings[heading_index]
    }
  }
  
  return(factor_table)
  
}  

## Compute number of probabilities that need to be stored for bayes net node 
## variables_in_node: variables that affect the node plus the node itself
##
## Will return number of probabilities that need to be stored for bayes net node 
probabilities_in_node = function(variables_in_node) {
  
  probabilities_to_store = 1.0
  for (survey_variable in variables_in_node) {
    probabilities_to_store = probabilities_to_store * survey_variables_states[which(cdc_survey_variables == survey_variable)]
  }
  return(probabilities_to_store)

}

## Drop columns from table
## factor_table: factor tables that needs columns dropped
## columns_to_drop: columms to be dropped
##
## Will return the table with dropped columns
drop_table_columns = function(factor_table, columns_to_drop) {
  
  if (!is.null(columns_to_drop)) {
    for (column_index in 1:length(columns_to_drop)) {
      factor_table = factor_table[-which(names(factor_table) == columns_to_drop[column_index])]
    }
  }
  
  return(factor_table)
  
}  

## Change column values to be meaningful
## factor_table: factor tables that needs columns dropped
## columns_name: name of the columm
## from_column_codes: codes to change
## to_column_values: values to change codes into
##
## Will return the table with column values changed from codes to meaningful ones
change_to_meaningful_values = function(factor_table, columns_name, from_column_codes, to_column_values) {
  
  if (!is.null(factor_table[[columns_name]])) {
    for (column_value_index in 1:length(from_column_codes)) {
      indexes_to_change = which(factor_table[[columns_name]] == from_column_codes[column_value_index]) 
      factor_table[[columns_name]] = replace(factor_table[[columns_name]], indexes_to_change, rep(c(to_column_values[column_value_index]), each=length(indexes_to_change)))
    }
  }
  
  return(factor_table)
  
}  

## Show formatted table 
## factor_table: factor tables that needs heading changed
## headings_to_change: headings to be changed
## new_headings: new headings
##
## Will print the formatted table
show_formatted_table = function(factor_table, headings_to_change, new_headings) {
  
  factor_table = change_table_headings(factor_table, headings_to_change, new_headings)
  grid.table(head(factor_table), rows=NULL)
  
}  

## Create and return the Bayes Net 
##
## Will return the Bayes Net
create_bayes_net = function() {
  
  risk_factors_bayes_net = list()
  
  # Add income.
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("income"))
  
  # Add smoke status given income
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("smoke","income"))
  
  # Add exercise given income
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("exercise","income"))
  
  # Add bmi given income and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("bmi","income", "exercise"))
  
  # Add blood pressure given exercise, income and smoking
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("bp","exercise", "income", "smoke"))
  
  # Add cholesterol given exercise, income and smoking
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("cholesterol","exercise", "income", "smoke"))
  
  # Add diabetes given bmi
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("diabetes","bmi"))
  
  # Add stroke given bmi, bp and cholesterol
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("stroke","bmi", "bp", "cholesterol"))
  
  # Add attack given bmi, bp and cholesterol
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("attack","bmi", "bp", "cholesterol"))
  
  # Add angina given bmi, bp and cholesterol
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("angina","bmi", "bp", "cholesterol"))
  
  return(risk_factors_bayes_net)

}

## Get health outcomes probability tables
## factor_table: factor tables that needs heading changed
## headings_to_change: headings to be changed
## new_headings: new headings
##
## Will return the health outcomes probability tables
get_health_outcomes = function(bayes_net, health_outcomes, observe_variables, observe_values) {

  health_outcomes_tables = list()

  # Find probability of health outcomes with observed variables like habits
  for (my_outcome in health_outcomes) {
  
    # Do not marginalize the observed variables and the outcome being considered
    outcomes_to_not_marginalize = c(observe_variables, my_outcome)
    outcomes_to_marginalize = cdc_survey_variables[!(cdc_survey_variables %in% outcomes_to_not_marginalize)]

    # Get inference table
    my_outcome_table = infer(bayesnet=bayes_net, 
                             margVars=outcomes_to_marginalize, 
                             obsVars=observe_variables, 
                             obsVals=observe_values)
    
    # Drop columns for observed variables
    my_outcome_table = drop_table_columns(factor_table=my_outcome_table, columns_to_drop=observe_variables)
    
    names(my_outcome_table)[[1]]="Probability"
    
    # Change values from codes to meaningful names
    my_outcome_table = change_to_meaningful_values(factor_table=my_outcome_table, columns_name=my_outcome, from_column_codes=outcome_codes[[my_outcome]], to_column_values=outcome_values[[my_outcome]])
    
    health_outcomes_tables[[length(health_outcomes_tables) + 1]] = my_outcome_table
    
  }
  
  return(health_outcomes_tables)
  
}

## Print the impact of habits on outcomes
## bayes_net: the bayes net
##
## Will print the impact of habits on outcomes
impact_of_habits_on_outcomes = function(bayes_net) {
  
  # What is the probability of the outcome if I have bad habits (smoke and don’t exercise)? 
  outcome_table_1 = get_health_outcomes(bayes_net=bayes_net, 
                                       health_outcomes=health_outcomes, 
                                       observe_variables=habit_variables, 
                                       observe_values=bad_habit_values)

  
  # How about if I have good habits (don’t smoke and do exercise)?
  outcome_table_2 = get_health_outcomes(bayes_net=bayes_net, 
                                       health_outcomes=health_outcomes, 
                                       observe_variables=habit_variables, 
                                       observe_values=good_habit_values)

  print_outcome_probabilities(outcome_table_1=outcome_table_1, 
                              outcome_table_1_observation="Bad Habits", 
                              outcome_table_2=outcome_table_2, 
                              outcome_table_2_observation="Good Habits", 
                              caption_text="Health Outcome Probabilities based on Habits")
  
}

## Print the impact of health on outcomes
## bayes_net: the bayes net
##
## Will print the impact of health on outcomes
impact_of_health_on_outcomes = function(bayes_net) {
  
  # What is the probability of the outcome if I have poor health (high blood pressure, high cholesterol, and overweight)?
  outcome_table_1 = get_health_outcomes(bayes_net=bayes_net, 
                                       health_outcomes=health_outcomes, 
                                       observe_variables=health_variables, 
                                       observe_values=bad_health_values)

  
  # What if I have good health (low blood pressure, low cholesterol, and normal weight)?
  outcome_table_2 = get_health_outcomes(bayes_net=bayes_net, 
                                       health_outcomes=health_outcomes, 
                                       observe_variables=health_variables, 
                                       observe_values=good_health_values)

  print_outcome_probabilities(outcome_table_1=outcome_table_1, 
                              outcome_table_1_observation="Poor Health", 
                              outcome_table_2=outcome_table_2, 
                              outcome_table_2_observation="Good Health", 
                              caption_text="Health Outcome Probabilities based on Health")
  
}

## Print outcome probability table
##
## outcome_table_1: Tables for health outcome 1 
## outcome_table_1_observation: Description for health outcome 1 observed value
## outcome_table_2: Tables for health outcome 2
## outcome_table_2_observation: Description for health outcome 2 observed value
##
## Print formatted table with each outcome on a row with probability under various observed vales
print_outcome_probabilities = function(outcome_table_1, outcome_table_1_observation, outcome_table_2, outcome_table_2_observation, caption_text) {
  
  table_counter = 0
  table_1_probabilities = numeric()
  for (factor_table in outcome_table_1) {
    table_counter = table_counter + 1
    table_1_probabilities[length(table_1_probabilities) + 1] = 
      factor_table[["Probability"]][which(factor_table[[health_outcomes[table_counter]]] == positive_health_outcomes_codes[table_counter])]
  }

  table_counter = 0
  table_2_probabilities = numeric()
  for (factor_table in outcome_table_2) {
    table_counter = table_counter + 1
    table_2_probabilities[length(table_2_probabilities) + 1] = 
      factor_table[["Probability"]][which(factor_table[[health_outcomes[table_counter]]] == positive_health_outcomes_codes[table_counter])]
  }
  
  table_to_print = data.frame(health_outcomes, table_1_probabilities, table_2_probabilities)
  names(table_to_print) = c("Health Outcome", outcome_table_1_observation, outcome_table_2_observation)
  print(kable(table_to_print, caption = caption_text))  

}

## Create a bayes net where habits affect outcomes
##
## Will return a bayes net with edges from smoking to each of the four outcomes and edges from
## exercise to each of the four outcomes.
make_habits_impact_outcomes = function() {
  
  risk_factors_bayes_net = list()
  
  # Add income.
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("income"))
  
  # Add smoke status given income
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("smoke","income"))
  
  # Add exercise given income
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("exercise","income"))
  
  # Add bmi given income and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("bmi","income", "exercise"))
  
  # Add blood pressure given exercise, income and smoking
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("bp","exercise", "income", "smoke"))
  
  # Add cholesterol given exercise, income and smoking
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("cholesterol","exercise", "income", "smoke"))
  
  # Add diabetes given bmi, smoking and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("diabetes","bmi", "smoke", "exercise"))
  
  # Add stroke given bmi, bp, cholesterol, smoking and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("stroke","bmi", "bp", "cholesterol", "smoke", "exercise"))
  
  # Add attack given bmi, bp, cholesterol, smoking and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("attack","bmi", "bp", "cholesterol", "smoke", "exercise"))
  
  # Add angina given bmi, bp, cholesterol, smoking and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("angina","bmi", "bp", "cholesterol", "smoke", "exercise"))
  
  return(risk_factors_bayes_net)
  
}

## Create a bayes net where diabetes affects stroke
##
## Will a bayes net with edges from diabetes to stroke
make_diabetes_impact_stroke_outcome = function() {
  
  risk_factors_bayes_net = list()
  
  # Add income.
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("income"))
  
  # Add smoke status given income
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("smoke","income"))
  
  # Add exercise given income
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("exercise","income"))
  
  # Add bmi given income and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("bmi","income", "exercise"))
  
  # Add blood pressure given exercise, income and smoking
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("bp","exercise", "income", "smoke"))
  
  # Add cholesterol given exercise, income and smoking
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("cholesterol","exercise", "income", "smoke"))
  
  # Add diabetes given bmi, smoking and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("diabetes","bmi", "smoke", "exercise"))
  
  # Add stroke given bmi, bp, cholesterol, smoking, exercise and diabetes
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("stroke","bmi", "bp", "cholesterol", "smoke", "exercise", "diabetes"))
  
  # Add attack given bmi, bp, cholesterol, smoking and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("attack","bmi", "bp", "cholesterol", "smoke", "exercise"))
  
  # Add angina given bmi, bp, cholesterol, smoking and exercise
  risk_factors_bayes_net[[length(risk_factors_bayes_net) + 1]] = createCPTfromData(cdc_survey_data, c("angina","bmi", "bp", "cholesterol", "smoke", "exercise"))
  
  return(risk_factors_bayes_net)
  
}
## Print the impact of diabetes on stroke outcome
## bayes_net: the bayes net
##
## Will print the impact of diabetes on stroke outcome
impact_of_diabetes_on_stroke = function(bayes_net, link_from_diabetes_to_stroke) {
  
  evaluated_outcome_stroke = c("stroke")
  observed_var_diabetes = c("diabetes")
  observed_val_have_diabetes = c("1")
  observed_val_no_diabetes = c("3")
  
  # What is the probability of a stroke if I have diabetes?
  outcome_tables = get_health_outcomes(bayes_net=bayes_net, 
                                       health_outcomes=evaluated_outcome_stroke, 
                                       observe_variables=observed_var_diabetes, 
                                       observe_values=observed_val_have_diabetes)
  table_counter = 0
  
  caption_text_qualifier = ""
  if (isTRUE(link_from_diabetes_to_stroke)) {
    caption_text_qualifier = "- edge from diabetes to stroke"
  } else {
    caption_text_qualifier = "- no edge from diabetes to stroke"
  }
    
  caption_text = "Probability for stroke when diabetes is present"
  
  for (factor_table in outcome_tables) {
    table_counter = table_counter + 1
    print(kable(factor_table, caption = paste(caption_text, caption_text_qualifier)))  
  }
  
  # What is the probability of a stroke if I do not have diabetes?
  outcome_tables = get_health_outcomes(bayes_net=bayes_net, 
                                       health_outcomes=evaluated_outcome_stroke, 
                                       observe_variables=observed_var_diabetes, 
                                       observe_values=observed_val_no_diabetes)
  table_counter = 0
  caption_text = "Probability for stroke when diabetes is not present"
  for (factor_table in outcome_tables) {
    table_counter = table_counter + 1
    print(kable(factor_table, caption = paste(caption_text, caption_text_qualifier)))  
  }
}