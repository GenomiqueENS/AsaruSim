process TRUNCATION_ESTIMATION {
    input:
    path alignment
    
    output:
    path "truncation_model.csv",  emit: truncation_model

    script:
    """
    python3.11 $projectDir/bin/AsaruSim.py truncation_estimator --paf $alignment 
    """
}