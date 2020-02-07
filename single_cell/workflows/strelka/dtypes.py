def dtypes():
    snv_strelka = {
        "chrom": "str",
        "coord": "int",
        "ref": "str",
        "alt": "str",
        "score": "int"
    }

    dtypes = locals()

    return dtypes
