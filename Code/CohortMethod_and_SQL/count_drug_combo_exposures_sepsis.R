# modified from the counting exercise to look at neuropathy
# to consider spesis outcomes
# started 7-25-20 JLW, rewritten 8-5-20 JLW
# modified 8-13-20 JLW to remove drugs with
# any overlapping targets from 2017 DrugBank

library(CohortMethod)
library(SqlRender)

## Check access to Optum??
connectionDetails <- createConnectionDetails(dbms = "postgresql",  server = "server_name.stanford.edu", user="user_name", password="password")
connection = connect(connectionDetails)
source_schema = "source_schema"
results_schema = "results_schema"
results_table = "SepsisDrugCombos_JLW"
cdmVersion = 5
outputFolder <- "./sepsisDrugCombos"

# analyze relevant tables to expedite study
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
sql <- "IF OBJECT_ID('@results_schema.@results_table', 'U') IS NOT NULL DROP TABLE @results_schema.@results_table;"
sql <- render(sql, results_schema = results_schema,results_table = results_table)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- "CREATE TABLE @results_schema.@results_table (cohort_definition_id INT, cohort_start_date DATE,cohort_end_date DATE,subject_id BIGINT);"
sql <- render(sql, results_schema = results_schema,results_table = results_table)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

# start counting sepsis outcomes and relevant drug exposures
sql <- readSql("count_dmes_outcomes.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,vocabulary_database_schema = source_schema,cohort_id=1, rel_concept_id= 132797)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

adrb2_drug_str="914335" # Atropine
sql <- readSql("count_drug_eras_from_list.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,vocabulary_database_schema = source_schema,cohort_id=2, drug_list_str =adrb2_drug_str)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql) 

combo_drug="1154343" # Albuterol, Salbutamol
sql <- readSql("count_drug_eras_singleDrug.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,vocabulary_database_schema = source_schema,cohort_id=3,drug_str=combo_drug)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

nonadrb2_drugs = "1186087,1309770,1105775,40161669,912263,930916,19080458,1151789,1114122,797399,962398,19041065,1366310,937368,901656,720810,964339,1140643"
sql <- readSql("count_drug_eras_from_list.sql")
sql <- render(sql, results_schema = results_schema,results_table = results_table,vocabulary_database_schema = source_schema,cohort_id=4,drug_list_str=nonadrb2_drugs)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

# after individual drugs or groups, look at overlaps
sql <- "ANALYZE @results_schema.@results_table;"
sql <- render(sql,results_schema=results_schema,results_table=results_table)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
executeSql(connection, sql)

# now look for subsequent mirtazapine treatment
sql <- readSql("look_for_overlaps.sql") # paroxetine or atropine followed by albuterol
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 55, first_cohort = 2, second_cohort = 3) 
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

sql <- readSql("look_for_overlaps.sql") # all other sepsis drugs followed by albuterol
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 6, first_cohort =4, second_cohort = 3)
sql <- translate(sql, targetDialect = connectionDetails$dbms)
#executeSql(connection, sql)

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
sql <- render(sql, results_schema = results_schema,results_table = results_table,cohort_id = 10, dme_id=1,drug_id=55)
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
	#1                     1  992596 # sepsis outcomes
	#2                     2  863338 # netwrok drugs
	#3                     3 8090945 # albuterol/ Salbutamol
	#4                     4 4090847 # non-network drugs
	#5                     6  378704 # non-network drugs -> albuterol
	#6                     7   44641 # sepsis after network drugs
	#7                     8  258718 # sepsis after albuterol
	#8                     9  206714 # sepsis after non-network drugs
	#9                    10    1803 # sepsis ater network drugs + albuterol
	#10                   11   46861 # sepsis after non-network drugs + albuterol
	#11                   55   13631 # network drugs -> albuterol


# now follow new Cohort Method, single studies tutorial, first looking at difference in sepsis outcomes with drug combo
networkDrugs <- c(914335,1186087,1309770,1105775,40161669,912263,930916,19080458,1151789,1114122,797399,962398,19041065,1366310,937368,901656,720810,964339,1140643)	
covSettings <- createDefaultCovariateSettings(excludedCovariateConceptIds = networkDrugs, addDescendantsToExclude = TRUE)

# cohortMethodData <- getDbCohortMethodData(connectionDetails = connectionDetails, cdmDatabaseSchema = source_schema,oracleTempSchema = results_schema, targetId = 55, comparatorId = 6, outcomeIds = 1, studyStartDate = "", studyEndDate = "", exposureDatabaseSchema = results_schema, exposureTable = results_table, outcomeDatabaseSchema = results_schema, outcomeTable = results_table, cdmVersion = cdmVersion, firstExposureOnly = TRUE, removeDuplicateSubjects = TRUE, restrictToCommonPeriod = FALSE, washoutPeriod = 180, covariateSettings = covSettings) # 5.43 hours

# saveCohortMethodData(cohortMethodData, "./sepsisDrugCombos/sepsisWithCombo.zip")
#print("Loading cohort method data, with drug combo")
#cohortMethodData <- loadCohortMethodData("./sepsisDrugCombos/sepsisWithCombo.zip")
#
#print("Setting up studyPop and PS model")
#studyPop <- createStudyPopulation(cohortMethodData = cohortMethodData, outcomeId = 1,
#firstExposureOnly = FALSE, restrictToCommonPeriod = FALSE, washoutPeriod = 0, removeDuplicateSubjects = "remove all",
#removeSubjectsWithPriorOutcome = FALSE, minDaysAtRisk = 1,riskWindowStart = 0, startAnchor = "cohort start", riskWindowEnd = 30,endAnchor = "cohort end")
#ps <- createPs(cohortMethodData = cohortMethodData, population = studyPop,errorOnHighCorrelation=FALSE)
#psauc<-computePsAuc(ps)
#print("PS AUC")
#print(psauc)
#outcomeModel <- fitOutcomeModel(population = ps, modelType = "cox",inversePtWeighting = TRUE)
#print("Outcome model, IPW:")
#print(outcomeModel)
	#[1] "Outcome model, IPW:"
	#Model type: cox
	#Stratified: FALSE
	#Use covariates: FALSE
	#Use inverse probability of treatment weighting: TRUE
	#Status: OK 
	#
	#           Estimate lower .95 upper .95     logRr seLogRr
	#treatment  0.180711  0.054235  0.427320 -1.710856  0.5266
	#Using prior: None
	#Using 1 thread(s)
	#Fitting outcome model took 0.354 secs

#matchedPop <- matchOnPs(ps,caliper = 0.2, caliperScale = "standardized logit", maxRatio= 1)
#outcomeModel <- fitOutcomeModel(population = matchedPop, modelType = "cox",)
#print("Outcome model, Matched")
#print(outcomeModel)
	#[1] "Outcome model, Matched"
	#Model type: cox
	#Stratified: FALSE
	#Use covariates: FALSE
	#Use inverse probability of treatment weighting: FALSE
	#Status: OK
	#
	#          Estimate lower .95 upper .95    logRr seLogRr
	#treatment  0.86041   0.50671   1.42146 -0.15034  0.2631

## Repeat cohort method look without combo drugs
print("Starting cohort method")
cohortMethodData <- getDbCohortMethodData(connectionDetails = connectionDetails, cdmDatabaseSchema = source_schema,oracleTempSchema = results_schema, targetId = 2, comparatorId = 4, outcomeIds = 1, studyStartDate = "", studyEndDate = "", exposureDatabaseSchema = results_schema, exposureTable = results_table, outcomeDatabaseSchema = results_schema, outcomeTable = results_table, cdmVersion = cdmVersion, firstExposureOnly = TRUE, removeDuplicateSubjects = TRUE, restrictToCommonPeriod = FALSE, washoutPeriod = 180, covariateSettings = covSettings)

saveCohortMethodData(cohortMethodData, "./sepsisDrugCombos/sepsisNoCombo.zip")
print("Loading cohort method data, without drug combo")
cohortMethodData <- loadCohortMethodData("./sepsisDrugCombos/sepsisNoCombo.zip")

print("Setting up studyPop and PS model")
studyPop <- createStudyPopulation(cohortMethodData = cohortMethodData, outcomeId = 1,
firstExposureOnly = FALSE, restrictToCommonPeriod = FALSE, washoutPeriod = 0, removeDuplicateSubjects = "remove all",
removeSubjectsWithPriorOutcome = FALSE, minDaysAtRisk = 1,riskWindowStart = 0, startAnchor = "cohort start", riskWindowEnd = 30,endAnchor = "cohort end")
ps <- createPs(cohortMethodData = cohortMethodData, population = studyPop,errorOnHighCorrelation=FALSE)
psauc<-computePsAuc(ps)
print("PS AUC")
print(psauc)
outcomeModel <- fitOutcomeModel(population = ps, modelType = "cox",inversePtWeighting = TRUE)
print("Outcome model, IPW:")
print(outcomeModel)

matchedPop <- matchOnPs(ps,caliper = 0.2, caliperScale = "standardized logit", maxRatio= 1)
outcomeModel <- fitOutcomeModel(population = matchedPop, modelType = "cox",)
print("Outcome model, Matched")
print(outcomeModel)

#pdf("OptumNeuropath_PS_overlap_72920.pdf") 
#plotPs(ps, scale = "preference", showCountsLabel = TRUE, showAucLabel = TRUE, showEquiposeLabel = TRUE)
#dev.off()
#pdf("OptumNeuropath_FollwUpDist_72920.pdf") 
#plotFollowUpDistribution(population = matchedPop)
#dev.off()
#
