def dtypes():
    cell_counts = {
        "cluster_id": "int",
        "cell_id": "str",
        "read_count": "int"
    }
    library = {
        "prediction_id": "int",
        "num_reads": "int",
        "num_unique_reads": "int",
        "library": "str",
        "is_normal": "bool",
        "patient_id": "float"
    }
    breakpoints = {
        "prediction_id": "int",
        "chromosome_1": "str",
        "strand_1": "str",
        "position_1": "int",
        "chromosome_2": "str",
        "strand_2": "str",
        "position_2": "int",
        "homology": "int",
        "num_split": "int",
        "inserted": "str",
        "mate_score": "float",
        "template_length_1": "int",
        "log_cdf": "float",
        "template_length_2": "int",
        "log_likelihood": "float",
        "template_length_min": "int",
        "num_reads": "int",
        "num_unique_reads": "int",
        "type": "str",
        "num_inserted": "int",
        "sequence": "str",
        "gene_id_1": "str",
        "gene_name_1": "str",
        "gene_location_1": "str",
        "gene_id_2": "str",
        "gene_name_2": "str",
        "gene_location_2": "str",
        "dgv_ids": "float",
        "is_germline": "bool",
        "is_dgv": "bool",
        "num_patients": "int",
        "is_filtered": "bool",
        "dist_filtered": "float",
        "balanced": "bool",
        "rearrangement_type": "str"
    }

    dtypes = locals()

    return dtypes
