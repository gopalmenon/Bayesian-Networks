## Function to create a conditional probability table
## Conditional probability is of the form p(x1 | x2, ..., xk)
## varnames: vector of variable names (strings)
## -- NOTE: first variable listed will be x1, remainder will be parents, x2, ..., xk
## probs: vector of probabilities for the flattened probability table
## levelsList: a list containing a vector of levels (outcomes) for each variable
## See the BayesNetExamples.r file for examples of how this function works
createCPT = function(varnames, probs, levelsList)
{
  ## Check dimensions agree
  if(length(probs) != prod(sapply(levelsList, FUN=length)))
    return(NULL)

  ## Set up table with appropriate dimensions
  m = length(probs)
  n = length(varnames)
  g = matrix(0, m, n)

  ## Convert table to data frame (with column labels)
  g = as.data.frame(g)
  names(g) = varnames

  ## This for loop fills in the entries of the variable values
  k = 1
  for(i in n:1)
  {
    levs = levelsList[[i]]
    g[,i] = rep(levs, each = k, times = m / (k * length(levs)))
    k = k * length(levs)
  }

  return(data.frame(probs = probs, g))
}

## Build a CPT from a data frame
## Constructs a conditional probability table as above, but uses frequencies
## from a data frame of data to generate the probabilities.
createCPTfromData = function(x, varnames)
{
  levelsList = list()

  for(i in 1:length(varnames))
  {
    name = varnames[i]
    levelsList[[i]] = sort(unique(x[,name]))
  }

  m = prod(sapply(levelsList, FUN=length))
  n = length(varnames)
  g = matrix(0, m, n)

  ## Convert table to data frame (with column labels)
  g = as.data.frame(g)
  names(g) = varnames

  ## This for loop fills in the entries of the variable values
  k = 1
  for(i in n:1)
  {
    levs = levelsList[[i]]
    g[,i] = rep(levs, each = k, times = m / (k * length(levs)))
    k = k * length(levs)
  }

  ## This is the conditional probability column
  probs = numeric(m)
  numLevels = length(levelsList[[1]])
  skip = m / numLevels

  ## This chunk of code creates the vector "fact" to index into probs using
  ## matrix multiplication with the data frame x
  fact = numeric(ncol(x))
  lastfact = 1
  for(i in length(varnames):1)
  {
    j = which(names(x) == varnames[i])
    fact[j] = lastfact
    lastfact = lastfact * length(levelsList[[i]])
  }
  ## Compute unnormalized counts of subjects that satisfy all conditions
  a = as.matrix(x - 1) %*% fact + 1
  for(i in 1:m)
    probs[i] = sum(a == i)

  ## Now normalize the conditional probabilities
  for(i in 1:skip)
  {
    denom = 0 ## This is the normalization
    for(j in seq(i, m, skip))
      denom = denom + probs[j]
    for(j in seq(i, m, skip))
    {
      if(denom != 0)
        probs[j] = probs[j] / denom
    }
  }

  return(data.frame(probs = probs, g))
}

## Product of two factors
## A, B: two factor tables
##
## Should return a factor table that is the product of A and B.
## You can assume that the product of A and B is a valid operation.
productFactor = function(table_a, table_b) {

  # Variables in product
  output_columns = union(colnames(table_a)[-1], colnames(table_b)[-1])
  
  # Common variables in input tables
  common_variables = intersect(colnames(table_a)[-1], colnames(table_b)[-1])
  
  # Variables only in Tables A and B
  table_a_only_variables = setdiff(colnames(table_a)[-1], common_variables)
  table_b_only_variables = setdiff(colnames(table_b)[-1], common_variables)
  
  table_a_rows = NROW(table_a)
  table_b_rows = NROW(table_b)
  
  probabilities_column = numeric()
  row_counter = 0
  
  # Define a matrix for the output values with no rows
  output_values_matrix = matrix("", 0, length(output_columns))
  
  # Loop through the tables and find the product
  for (table_a_index in 1:table_a_rows) {
    for (table_b_index in 1:table_b_rows) {
      
      # Join the two tables unconditionally if no common variables  
      if (length(common_variables) == 0) {

        row_counter = row_counter + 1
        column_counter = 0
        
        probabilities_column[row_counter] = table_a$probs[table_a_index] * table_b$probs[table_b_index]
        
        output_row_variables = character(length = length(output_columns))
        
        #Put in values from Table A
        for (column_value in table_a_only_variables) {
          
          column_counter = column_counter + 1
          output_row_variables[column_counter] = table_a[column_value][[1]][table_a_index]
          
        }
        
        #Put in values from Table B
        for (column_value in table_b_only_variables) {
          
          column_counter = column_counter + 1
          output_row_variables[column_counter] = table_b[column_value][[1]][table_b_index]
          
        }

        # Add the output values to the matrix
        output_values_matrix = rbind(output_values_matrix, output_row_variables)
        
      } else {
        
        rows_overlap = TRUE
        
        # Loop through the common variables
        for (common_variable_index in 1:length(common_variables)) {
          
          # Find out if the common variables have the same values in both tables
          if (table_a[common_variables[common_variable_index]][[1]][table_a_index] !=
              table_b[common_variables[common_variable_index]][[1]][table_b_index]) {
            rows_overlap = FALSE
            break
          }

        }
      
        # If the rows overlap, then create a row in the output table
        if (isTRUE(rows_overlap)) {
          
          row_counter = row_counter + 1
          column_counter = 0
          
          probabilities_column[row_counter] = table_a$probs[table_a_index] * table_b$probs[table_b_index]
          
          output_row_variables = character(length = length(output_columns))
          
          #Put in values only in Table A
          for (column_value in table_a_only_variables) {
            
            column_counter = column_counter + 1
            output_row_variables[column_counter] = table_a[column_value][[1]][table_a_index]
            
          }
          
          #Put in values only in Table B
          for (column_value in table_b_only_variables) {
            
            column_counter = column_counter + 1
            output_row_variables[column_counter] = table_b[column_value][[1]][table_b_index]
            
          }
          
          # Put in common values
          for (column_value in common_variables) {
            
            column_counter = column_counter + 1
            output_row_variables[column_counter] = table_a[column_value][[1]][table_a_index]
            
          }
          
          # Add the output values to the matrix
          output_values_matrix = rbind(output_values_matrix, output_row_variables)

        }
      }
      
    }
  }
  
  # Convert the output matrix to a data frame
  output_values_matrix = as.data.frame(output_values_matrix)
  
  # Add column names for the output values 
  names(output_values_matrix) = c(as.character(table_a_only_variables), as.character(table_b_only_variables), as.character(common_variables))
  
  return(data.frame(probs = probabilities_column, output_values_matrix))
  
}

## Marginalize a variable from a factor
## A: a factor table
## margVar: a string of the variable name to marginalize
##
## Should return a factor table that marginalizes margVar out of A.
## You can assume that margVar is on the left side of the conditional.
marginalizeFactor = function(X, margVar)
{
  ## Your code here!
}

## Marginalize a list of variables
## bayesnet: a list of factor tables
## margVars: a vector of variable names (as strings) to be marginalized
##
## Should return a Bayesian network (list of factor tables) that results
## when the list of variables in margVars is marginalized out of bayesnet.
marginalize = function(bayesnet, margVars)
{
  ## Your code here!
}

## Observe values for a set of variables
## bayesnet: a list of factor tables
## obsVars: a vector of variable names (as strings) to be observed
## obsVals: a vector of values for corresponding variables (in the same order)
##
## Set the values of the observed variables. Other values for the variables
## should be removed from the tables. You do not need to normalize the factors
## to be probability mass functions.
observe = function(bayesnet, obsVars, obsVals)
{
  ## Your code here!
}

## Run inference on a Bayesian network
## bayesnet: a list of factor tables
## margVars: a vector of variable names to marginalize
## obsVars: a vector of variable names to observe
## obsVals: a vector of values for corresponding variables (in the same order)
##
## This function should run marginalization and observation of the sets of
## variables. In the end, it should return a single joint probability table. The
## variables that are marginalized should not appear in the table. The variables
## that are observed should appear in the table, but only with the single
## observed value. The variables that are not marginalized or observed should
## appear in the table with all of their possible values. The probabilities
## should be normalized to sum to one.
infer = function(bayesnet, margVars, obsVars, obsVals)
{
  ## Your code here!
}
