# modified from the counting exercise to look at pancreatitis
# to consider spesis outcomes
# started 7-25-20 JLW, rewritten 8-5-20 JLW

library(CohortMethod)
library(SqlRender)

## Check access to Optum??
connectionDetails <- createConnectionDetails(dbms = "postgresql",  server = "server_name.stanford.edu", user="user_name", password="password")
connection = connect(connectionDetails)
source_schema = "source_schema"
results_schema = "results_schema"
results_table = "PancreatitisDrugCombos_JLW"
cdmVersion = 5
outputFolder <- "./pancreatitisDrugCombos"

# analyze relevant tables to expedite study
print("analyzing relevant tables")
sql <- "ANALYZE @source_schema.CONCEPT;"
sql <- render(sql,source_schema=source_schema)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- "ANALYZE @source_schema.CONDITION_OCCURRENCE;"
sql <- render(sql,source_schema=source_schema)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- "ANALYZE @source_schema.observation_period;"
sql <- render(sql,source_schema=source_schema)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

# create table
print("creating relevant tables")
sql <- "IF OBJECT_ID('@results_schema.@results_table', 'U') IS NOT NULL DROP TABLE @results_schema.@results_table;"
sql <- render(sql, results_schema = results_schema,results_table = results_table)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- "CREATE TABLE @results_schema.@results_table (cohort_definition_id INT, cohort_start_date DATE,cohort_end_date DATE,subject_id BIGINT);"
sql <- render(sql, results_schema = results_schema,results_table = results_table)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

# start counting sepsis outcomes and relevant drug exposures
print("Counting single exposures and dme outcomes")
sql <- readSql("count_dmes_outcomes.sql")
cohort_id = 1
dme_id =4192640 
sql <- render(sql, results_schema = results_schema,results_table = results_table,vocabulary_database_schema = source_schema,cohort_id=cohort_id, rel_concept_id= dme_id) # pancreatitis
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

net_drug_str = "19043959,757688,44814542,914335,732309,945286,720810,713823" # pancreatitis-associated drugs with TP53,EDNRA,NFKBIA in their networks
cohort_id = 2
sql <- readSql("count_drug_eras_from_list.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,vocabulary_database_schema = source_schema,cohort_id=cohort_id, drug_list_str =net_drug_str)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql) 

# combo_drug="40220812" # Saililic Acid? 
combo_drug="1112807" # Aspirin
cohort_id = 3
sql <- readSql("count_drug_eras_singleDrug.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,vocabulary_database_schema = source_schema,cohort_id=cohort_id,drug_str=combo_drug)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

nonNet_drug_str = "1317967,713109,1335471,45892531,1333357,1511449,930916,1342001,1363749,797399,1376289,19028106,1503501,789578,901656,1522957,40226742,1318011,43012417,1130585,951279,1331235,1334456,735951,19098548,1103314,1342439,40238052" # pancreatitis drugs without TP53,EDNRA,NFKBIA in their networks
cohort_id = 4
sql <- readSql("count_drug_eras_from_list.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,vocabulary_database_schema = source_schema,cohort_id=cohort_id,drug_list_str=nonNet_drug_str)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

# after individual drugs or groups, look at overlaps
sql <- "ANALYZE @results_schema.@results_table;"
sql <- render(sql,results_schema=results_schema,results_table=results_table)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

print("Looking for overlapping drug eras")
sql <- readSql("look_for_overlaps.sql") # network drug followed by combo drug 
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 5, first_cohort = 2, second_cohort = 3) 
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- readSql("look_for_overlaps.sql") # all other dme drugs followed by combo drug 
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 6, first_cohort =4, second_cohort = 3)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

print("Look for DME outcomes after single or combo drug eras")
# now see if there are any outcomes after these exposures
sql <- readSql("look_for_subsequent_outcomes.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 7, dme_id=1,drug_id=2)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- readSql("look_for_subsequent_outcomes.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 8, dme_id=1,drug_id=3)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- readSql("look_for_subsequent_outcomes.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 9, dme_id=1,drug_id=4)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- readSql("look_for_subsequent_outcomes.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 10, dme_id=1,drug_id=5)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- readSql("look_for_subsequent_outcomes.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 11, dme_id=1,drug_id=6)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

# count cohort numbers
sql <- paste("SELECT cohort_definition_id, COUNT(DISTINCT subject_id) AS count",
"FROM @results_schema.@results_table",
"GROUP BY cohort_definition_id")
sql <- render(sql, results_schema = results_schema,results_table = results_table)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
cohort_counts <- querySql(connection, sql)
print("Finished counting patients in cohort")
print(cohort_counts)

#   COHORT_DEFINITION_ID   COUNT
#1                     1  348519 # pancreatitis 
#2                     2 1648281 # network drugs
#3                     3  418178 # aspirin
#4                     4 8432296 # non-network drugs
#5                     5    6687 # net drug -> aspirin
#6                     6   50380 # non-net drug -> aspirin
#7                     7   28717 # net drugs -> pancreatitis
#8                     8    5408 # aspirin -> pancreatitis
#9                     9  120812 # non-net drugs -> pancreatitis
#10                   10     191 # net drug -> aspirin -> pancreatitis
#11                   11    1252 # non-net drug -> aspirin -> pancreatitis


### FINISH THIS LATER ###
### First compare combo drug + net, compared to combo drug + nonNet
print("Starting Causal Inference, with combo drug")

# now follow new Cohort Method, single studies tutorial, first looking at difference in sepsis outcomes with drug combo
networkDrugs <- c(19043959,757688,44814542,914335,732309,945286,720810,713823,1317967,713109,1335471,45892531,1333357,1511449,930916,1342001,1363749,797399,1376289,19028106,1503501,789578,901656,1522957,40226742,1318011,43012417,1130585,951279,1331235,1334456,735951,19098548,1103314,1342439,40238052)
covSettings <- createDefaultCovariateSettings(excludedCovariateConceptIds = networkDrugs, addDescendantsToExclude = TRUE)

print("Running Cohort Method")
# cohortMethodData <- getDbCohortMethodData(connectionDetails = connectionDetails, cdmDatabaseSchema = source_schema,oracleTempSchema = results_schema, targetId = 5, comparatorId = 6, outcomeIds = 1, studyStartDate = "", studyEndDate = "", exposureDatabaseSchema = results_schema, exposureTable = results_table, outcomeDatabaseSchema = results_schema, outcomeTable = results_table, cdmVersion = cdmVersion, firstExposureOnly = TRUE, removeDuplicateSubjects = TRUE, restrictToCommonPeriod = FALSE, washoutPeriod = 180, covariateSettings = covSettings) # 5.43 hours

#saveCohortMethodData(cohortMethodData, "./pancreatitisDrugCombos/pancreatitisADBR2.zip")
cohortMethodData <- loadCohortMethodData("./pancreatitisDrugCombos/pancreatitisADBR2.zip")

print("Setting up studyPop and PS model")
studyPop <- createStudyPopulation(cohortMethodData = cohortMethodData, outcomeId = 1,
firstExposureOnly = FALSE, restrictToCommonPeriod = FALSE, washoutPeriod = 0, removeDuplicateSubjects = "remove all",
removeSubjectsWithPriorOutcome = FALSE, minDaysAtRisk = 1,riskWindowStart = 0, startAnchor = "cohort start", riskWindowEnd = 30,endAnchor = "cohort end")
ps <- createPs(cohortMethodData = cohortMethodData, population = studyPop)
psauc<-computePsAuc(ps)
print("PS AUC")
print(psauc)
#[1] "PS AUC"
#[1] 0.8988511
outcomeModel <- fitOutcomeModel(population = ps, modelType = "cox",inversePtWeighting = TRUE)
print(outcomeModel)
#
#          Estimate lower .95 upper .95   logRr seLogRr
#treatment  1.56380   1.05921   2.22268 0.44712  0.1891
#Using prior: None
#Using 1 thread(s)
#Fitting outcome model took 0.149 secs
matchedPop <- matchOnPs(ps,caliper = 0.2, caliperScale = "standardized logit", maxRatio= 1)
outcomeModel <- fitOutcomeModel(population = matchedPop, modelType = "cox",) # stratified=TRUE)
print('outcome model with matched pop')
print(outcomeModel)
#Model type: cox
#Stratified: FALSE
#Use covariates: FALSE
#Use inverse probability of treatment weighting: FALSE
#Status: OK 
#
#           Estimate lower .95 upper .95     logRr seLogRr
#treatment 1.0011485 0.5137463 1.9592087 0.0011478  0.3415

### Perform comparison without combo drug
print("Running Cohort Method")
cohortMethodData <- getDbCohortMethodData(connectionDetails = connectionDetails, cdmDatabaseSchema = source_schema,oracleTempSchema = results_schema, targetId = 2, comparatorId = 4, outcomeIds = 1, studyStartDate = "", studyEndDate = "", exposureDatabaseSchema = results_schema, exposureTable = results_table, outcomeDatabaseSchema = results_schema, outcomeTable = results_table, cdmVersion = cdmVersion, firstExposureOnly = TRUE, removeDuplicateSubjects = TRUE, restrictToCommonPeriod = FALSE, washoutPeriod = 180, covariateSettings = covSettings)

saveCohortMethodData(cohortMethodData, "./pancreatitisDrugCombos/pancreatitisNoCombo.zip")
cohortMethodData <- loadCohortMethodData("./pancreatitisDrugCombos/pancreatitisNoCombo.zip")

print("Setting up studyPop and PS model")
studyPop <- createStudyPopulation(cohortMethodData = cohortMethodData, outcomeId = 1,
firstExposureOnly = FALSE, restrictToCommonPeriod = FALSE, washoutPeriod = 0, removeDuplicateSubjects = "remove all",
removeSubjectsWithPriorOutcome = FALSE, minDaysAtRisk = 1,riskWindowStart = 0, startAnchor = "cohort start", riskWindowEnd = 30,endAnchor = "cohort end")
ps <- createPs(cohortMethodData = cohortMethodData, population = studyPop)
psauc<-computePsAuc(ps)
print("PS AUC")
print(psauc)
outcomeModel <- fitOutcomeModel(population = ps, modelType = "cox",inversePtWeighting = TRUE)
print(outcomeModel)
matchedPop <- matchOnPs(ps,caliper = 0.2, caliperScale = "standardized logit", maxRatio= 1)
outcomeModel <- fitOutcomeModel(population = matchedPop, modelType = "cox",) # stratified=TRUE)
print('outcome model with matched pop')
print(outcomeModel)
	#[1] "PS AUC"
	#[1] 0.9598428
	#Using prior: None
	#Using 1 thread(s)
	#Fitting outcome model took 2.18 mins
	#Model type: cox
	#Stratified: FALSE
	#Use covariates: FALSE
	#Use inverse probability of treatment weighting: TRUE
	#Status: OK

	#	  Estimate lower .95 upper .95   logRr seLogRr
	#treatment  1.29712   1.24604   1.34964 0.26015  0.0204
	#Using prior: None
	#Using 1 thread(s)
	#Fitting outcome model took 13.5 secs
	#[1] "outcome model with matched pop"
	#Model type: cox
	#Stratified: FALSE
	#Use covariates: FALSE
	#Use inverse probability of treatment weighting: FALSE
	#Status: OK

	#	  Estimate lower .95 upper .95    logRr seLogRr
	#treatment  0.58017   0.51918   0.64765 -0.54443  0.0564


# investigate the cohort and the PS overlap
#getFollowUpDistribution(population = matchedPop)
#pdf("OptumNeuropath_PS_overlap_72920.pdf") 
#plotPs(ps, scale = "preference", showCountsLabel = TRUE, showAucLabel = TRUE, showEquiposeLabel = TRUE)
#dev.off()
#pdf("OptumNeuropath_FollwUpDist_72920.pdf") 
#plotFollowUpDistribution(population = matchedPop)
#dev.off()
#getAttritionTable(matchedPop)
## show the PS model coefficients that have nonzero values
#getPsModel(ps, cohortMethodData)
