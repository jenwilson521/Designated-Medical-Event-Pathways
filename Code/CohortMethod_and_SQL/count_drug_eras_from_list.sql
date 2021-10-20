INSERT INTO @results_schema.@results_table (
        cohort_definition_id,
        cohort_start_date,
        cohort_end_date,
        subject_id)
SELECT DISTINCT @cohort_id, -- Look at condition occurrence for neuropathy outcomes
	c.drug_era_start_date,
	c.drug_era_end_date,
	c.person_id
	FROM @vocabulary_database_schema.DRUG_ERA c
	WHERE c.drug_concept_id IN
	(SELECT descendant_concept_id
                FROM @vocabulary_database_schema.concept_ancestor
                WHERE ancestor_concept_id in (@drug_list_str)); -- relevant clinical finding
