INSERT INTO @results_schema.@results_table (
        cohort_definition_id,
        cohort_start_date,
        cohort_end_date,
        subject_id)
SELECT DISTINCT @cohort_id, -- Overlapping treatment, first reference then predicted combination
        f.cohort_start_date,
        f.cohort_end_date,
        f.subject_id
        FROM @results_schema.@results_table f -- sinlge or double drug exposure 
        INNER JOIN @results_schema.@results_table s -- DME outcomes
        ON s.cohort_start_date > f.cohort_start_date
        AND s.subject_id = f.subject_id
        AND s.cohort_definition_id = @dme_id  -- dme outcomes 
        And f.cohort_definition_id = @drug_id; -- single or double drug exposure 
