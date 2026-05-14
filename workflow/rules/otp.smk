# Snakemake rules to preprocess otp data

rule otp_drug:
    input:
        trait_metadata = "resources/metadata/trait_metadata_curated.xlsx",
        otp_drug_dir = "resources/nci/otp_output/{otp_ver}/evidence_clinical_precedence"
    output:
        drug = "results/otp/{otp_ver}/otp_drug.gz.parquet",
        drug_summary = "results/otp/{otp_ver}/otp_drug_summary.gz.parquet"
    params:
        output_dir = lambda x, output: Path(output[0]).parent.as_posix()
    conda: "renv"
    script: "snakescripts/otp/otp_drug.R"

rule otp_assoc:
    input:
        trait_metadata = "resources/metadata/trait_metadata_curated.xlsx",
        assoc_datatype_direct = "resources/nci/otp_output/{otp_ver}/association_by_datatype_direct",
        assoc_datatype_indirect = "resources/nci/otp_output/{otp_ver}/association_by_datatype_indirect",
        assoc_datasource_direct = "resources/nci/otp_output/{otp_ver}/association_by_datasource_direct",
        assoc_datasource_indirect = "resources/nci/otp_output/{otp_ver}/association_by_datasource_indirect",
        assoc_overall_direct = "resources/nci/otp_output/{otp_ver}/association_overall_direct",
        assoc_overall_indirect = "resources/nci/otp_output/{otp_ver}/association_overall_indirect"
    output:
        assoc_overall = "results/otp/{otp_ver}/otp_assoc_overall.gz.parquet",
        assoc_datasource = "results/otp/{otp_ver}/otp_assoc_datasource.gz.parquet",
        assoc_datatype = "results/otp/{otp_ver}/otp_assoc_datatype.gz.parquet"
    params:
        output_dir = lambda x, output: Path(output[0]).parent.as_posix()
    conda: "renv"
    script: "snakescripts/otp/otp_assoc.R"

rule otp_drug_summary:
    """
    get info for each drug in the otp drug data
    """
    input:
        trait_metadata = "resources/metadata/trait_metadata_curated.xlsx",
        otp_drug_dir = "resources/nci/otp_output/{otp_ver}/known_drug"
    output:
        drug_summary = "results/otp/{otp_ver}/all_drug_summary.tsv"
    params:
        output_dir = lambda x, output: Path(output[0]).parent.as_posix()
    conda: "renv"
    script: "snakescripts/otp/drug_summary.R"

rule otp_drug_mechanism_chembl:
    """
    get info for mechanism of drug in the otp chembl
    """
    input:
        otp_drug_dir = "resources/nci/otp_output/{otp_ver}/drug_molecule",
        otp_drug_mechanism_dir = "resources/nci/otp_output/{otp_ver}/drug_mechanism_of_action",
        otp_evidence_chembl_dir = "resources/nci/otp_output/{otp_ver}/evidence_clinical_precedence"
    output:
        drug_mechanism = "results/otp/{otp_ver}/all_drug_mechanism_chembl.tsv",
        drug_evidence = "results/otp/{otp_ver}/all_drug_evidence_chembl.tsv",
        drug_molecule = "results/otp/{otp_ver}/all_drug_molecule_chembl.tsv"
    params:
        output_dir = lambda x, output: Path(output[0]).parent.as_posix()
    conda: "renv"
    script: "snakescripts/otp/drug_mechanism_chembl.R"


rule otp_disease_summary:
    """
    get info for each disease in the otp disease data
    """
    input:
        trait_metadata = "resources/metadata/trait_metadata_curated.xlsx",
        otp_disease_dir = "resources/nci/otp_output/{otp_ver}/disease"
    output: "results/otp/{otp_ver}/all_disease.tsv"
    params:
        output_dir = lambda x, output: Path(output[0]).parent.as_posix()
    conda: "renv"
    script: "snakescripts/otp/disease_summary.R"


rule otp_all:
    input:
        drug = rules.otp_drug.output.drug,
        drug_summary = rules.otp_drug.output.drug_summary,
        assoc_overall = rules.otp_assoc.output.assoc_overall,
        assoc_datasource = rules.otp_assoc.output.assoc_datasource,
        assoc_datatype = rules.otp_assoc.output.assoc_datatype,
        drug_mechanism = rules.otp_drug_mechanism_chembl.output.drug_mechanism,
        drug_evidence = rules.otp_drug_mechanism_chembl.output.drug_evidence,
        drug_molecule = rules.otp_drug_mechanism_chembl.output,
        disease_summary = rules.otp_disease_summary.output
    output:
        touch("results/otp/{otp_ver}/.done_all")
    shell:
        "touch {output}"