INSERT INTO @results_schema.@results_table (
        cohort_definition_id,
        cohort_start_date,
        cohort_end_date,
        subject_id)
SELECT DISTINCT @cohort_id, -- Overlapping treatment, first reference then predicted combination
        f.cohort_start_date,
        f.cohort_end_date,
        f.subject_id
        FROM @results_schema.@results_table f -- first treatment
        INNER JOIN @results_schema.@results_table s -- second treatment
        ON s.cohort_start_date > f.cohort_start_date
        AND s.cohort_start_date < f.cohort_end_date
        AND s.subject_id = f.subject_id
        AND s.cohort_definition_id = @second_cohort  -- predicted cotherapy
        And f.cohort_definition_id = @first_cohort; -- drugs with or without network mechanisms
