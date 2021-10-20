INSERT INTO @results_schema.@results_table (
        cohort_definition_id,
        cohort_start_date,
        cohort_end_date,
        subject_id)
SELECT DISTINCT @cohort_id, -- Look at condition occurrence for neuropathy outcomes
	c.condition_start_date,
	c.condition_end_date,
	c.person_id
	FROM @vocabulary_database_schema.CONDITION_OCCURRENCE c
	WHERE c.condition_concept_id IN
	(SELECT descendant_concept_id
                FROM @vocabulary_database_schema.concept_ancestor
                WHERE ancestor_concept_id = @rel_concept_id); -- relevant clinical finding
