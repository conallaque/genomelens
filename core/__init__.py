"""Input, output and provenance — everything that is not analysis.

    genome_input   Chip and VCF parsing, format sniffing, build detection
    imputation     Optional Beagle imputation with r2 gating
    qc             Call rate, sex inference, chip identification
    snp_registry   The unified SNP record store
    provenance     Per-variant source tagging (chip vs imputed) and build
    references     Citation collection and evidence levels
    build_stamp    Which commit produced an artefact — see its docstring
    fhir_export, pdf_export, emergency_card   Output formats
"""
