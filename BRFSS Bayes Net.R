# Load the bayes net functions
source("BayesianNetworks-template.r")

# Load data collected by the CDC in the 2015 Behavioral Risk Factor Surveillance System (BRFSS) survey.
cdc_survey_data = read.csv(file="RiskFactors.csv", header=TRUE, sep=",")

# Risk factors collected in the survey
cdc_survey_variables = list("income", "exercise", "smoke", "bmi", "bp", "cholesterol", "angina", 
                            "stroke", "attack", "diabetes")

## 1.
# Create the following Bayesian network to analyze the survey results. It will be the product of 
# all the probabilities corresponding to the full joint distribution. 
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

## 2.
# For each of the four health outcomes (diabetes, stroke, heart attack, angina), answer the following by 
# querying the network

## (a)
# What is the probability of the outcome if I have bad habits (smoke and don’t exercise)? How about if I have 
# good habits (don’t smoke and do exercise)?
health_outcomes = c("diabetes", "stroke", "attack", "angina")
habits_variables = c("smoke", "exercise")
bad_habit_values = c("1", "2")

# Find probability of health outcomes with bad habits
for (my_outcome in health_outcomes) {
  
  # Do not marginalize the observed variables and the outcome being considered
  outcomes_to_not_marginalize = c(habits_variables, my_outcome)
  outcomes_to_marginalize = cdc_survey_variables[!(cdc_survey_variables %in% outcomes_to_not_marginalize)]
  print(paste("Outcome: ", my_outcome))
  my_outcome_table = infer(bayesnet=risk_factors_bayes_net, 
                           margVars=outcomes_to_marginalize, 
                           obsVars=habits_variables, 
                           obsVals=bad_habit_values)
  print(my_outcome_table, row.names=FALSE)
}

# Find probability of health outcomes with good habits
good_habit_values = c("2", "1")
for (my_outcome in health_outcomes) {
  
  # Do not marginalize the observed variables and the outcome being considered
  outcomes_to_not_marginalize = c(habits_variables, my_outcome)
  outcomes_to_marginalize = cdc_survey_variables[!(cdc_survey_variables %in% outcomes_to_not_marginalize)]
  print(paste("Outcome: ", my_outcome))
  my_outcome_table = infer(bayesnet=risk_factors_bayes_net, 
                           margVars=outcomes_to_marginalize, 
                           obsVars=habits_variables, 
                           obsVals=good_habit_values)
  print(my_outcome_table, row.names=FALSE)
}

## (b) 
# What is the probability of the outcome if I have poor health (high blood pressure, high cholesterol, and overweight)? 
# What if I have good health (low blood pressure, low cholesterol, and normal weight)?

health_variables = c("bmi", "bp", "cholesterol")
bad_health_values = c("3", "1", "1")

# Find probability of health outcomes with bad health
for (my_outcome in health_outcomes) {
  
  # Do not marginalize the observed variables and the outcome being considered
  outcomes_to_not_marginalize = c(health_variables, my_outcome)
  outcomes_to_marginalize = cdc_survey_variables[!(cdc_survey_variables %in% outcomes_to_not_marginalize)]
  print(paste("Outcome: ", my_outcome))
  my_outcome_table = infer(bayesnet=risk_factors_bayes_net, 
                           margVars=outcomes_to_marginalize, 
                           obsVars=health_variables, 
                           obsVals=bad_health_values)
  print(my_outcome_table, row.names=FALSE)
}

# Find probability of health outcomes with good health
good_health_values = c("2", "3", "2")
for (my_outcome in health_outcomes) {
  
  # Do not marginalize the observed variables and the outcome being considered
  outcomes_to_not_marginalize = c(health_variables, my_outcome)
  outcomes_to_marginalize = cdc_survey_variables[!(cdc_survey_variables %in% outcomes_to_not_marginalize)]
  print(paste("Outcome: ", my_outcome))
  my_outcome_table = infer(bayesnet=risk_factors_bayes_net, 
                           margVars=outcomes_to_marginalize, 
                           obsVars=health_variables, 
                           obsVals=good_health_values)
  print(my_outcome_table, row.names=FALSE)
}

