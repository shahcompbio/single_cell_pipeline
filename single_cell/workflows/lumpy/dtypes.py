def dtypes():
    evidence = {
        "breakpoint_id": "int",
        "cell_id": "str",
        "count": "int"
    }

    breakpoint = {
        "breakpoint_id": "int",
        "chrom1": "str",
        "start1": "int",
        "end1": "int",
        "strand1": "str",
        "max_chr1": "str",
        "max_pos1": "int",
        "confidence_interval_chr1": "str",
        "confidence_interval_start1": "int",
        "confidence_interval_end1": "int",
        "chrom2": "str",
        "start2": "int",
        "end2": "int",
        "strand2": "str",
        "max_chr2": "str",
        "max_pos2": "int",
        "confidence_interval_chr2": "str",
        "confidence_interval_start2": "int",
        "confidence_interval_end2": "int",
        "type": "str",
        "score": "float",
        "strands": "str",
        "normal_PE": "float",
        "tumour_PE": "float",
        "tumour_SR": "float",
        "normal_SR": "float",
    }
    dtypes = locals()

    return dtypes
