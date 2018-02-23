# Load data collected by the CDC in the 2015 Behavioral Risk Factor Surveillance System (BRFSS) survey.
cdc_survey_data = read.csv(file="RiskFactors.csv", header=TRUE, sep=",")

# Risk factors collected in the survey
cdc_survey_variables = list("income", "exercise", "smoke", "bmi", "bp", "cholesterol", "angina", 
                            "stroke", "attack", "diabetes")

# Create probability density tables for variables corresponding to the data.
probability_density_tables = list()
table_counter = 0
for (survey_variable in cdc_survey_variables) {
  table_counter = table_counter + 1
  probability_density_tables[[table_counter]] = createCPTfromData(cdc_survey_data, survey_variable)
}

# Create the Bayes Net corresponding to the relationships between variables. It will be the product of 
# all the probabilities corresponding to the full joint distribution. 
table_counter = 0
risk_factors_bayes_net = list()

# Add income.
table_counter = table_counter + 1
risk_factors_bayes_net[[table_counter]] = probability_density_tables[[which(cdc_survey_variables == "income")]]

# Add smoke status given income
table_counter = table_counter + 1
risk_factors_bayes_net[[table_counter]] = productFactor(probability_density_tables[[which(cdc_survey_variables == "smoke")]],
                                                        probability_density_tables[[which(cdc_survey_variables == "income")]])

# Add exercise given income
table_counter = table_counter + 1
risk_factors_bayes_net[[table_counter]] = exercise_given_income =
  productFactor(probability_density_tables[[which(cdc_survey_variables == "exercise")]],
                probability_density_tables[[which(cdc_survey_variables == "income")]])

# Add bmi given income and exercise
table_counter = table_counter + 1
risk_factors_bayes_net[[table_counter]] = bmi_given_exercise_income =
  productFactor(productFactor(probability_density_tables[[which(cdc_survey_variables == "bmi")]],
                probability_density_tables[[which(cdc_survey_variables == "income")]]),
                exercise_given_income)
