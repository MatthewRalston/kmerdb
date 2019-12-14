import jsonschema

VERSION="0.0.1"

metadata_schema = {
            "type": "object",
            "properties": {
                "neighbors": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "id": {"type": "number"},
                            "count_delta": {"type": "number"}
                        }
                    } # int k-mer ids
                },
                "tags": {
                    "type": "array",
                    "items": {"type": "string"} # k-mer tags
                },
                "bools": {
                    "type": "array", 
                    "items": { # key-value boolean name, value
                        "type": "object",
                        "properties": {
                            "name": {"type": "string"},
                            "value": {"type": "boolean"}
                        }
                    }
                },
                "floats": {
                    "type": "array",
                    "items": { # numerical annotations, floating-point name, value
                        "type": "object",
                        "properties": {
                            "name": {"type": "string"},
                            "value": {"type": "number"}
                        }
                    }
                },
                "ints": {
                    "type": "array",
                    "items": { # numerical annotations, integer name, value (e.g. degree
                        "type": "object",
                        "properties": {
                            "name": {"type": "string"},
                            "value": {"type": "number"}
                        }
                    }
                }
            }
    }

node_schema = {
    "type": "object",
    "properties": {
        "id": {"type": "number"},
        "count": {"type": "number"},
        "metadata": metadata_schema
    }
}

header_schema = {
    "type": "object",
    "properties": {
        "kdb_ver": {"type": "string"},
        "k": {"type": "number"},
        "metadata": {"type": "boolean"},
        "tags": {
            "type": "array",
            "items": {"type": "string"}
        },
        "files": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "filename": {"type": "string"},
                    "sha256": {
                        "type": "string",
                        "minLength": 64,
                        "maxLength": 64
                    },
                    "md5": {
                        "type": "string",
                        "minLength": 32,
                        "maxLength": 32
                    },
                    "total_reads": {"type": "number"},
                    "total_kmers": {"type": "number"},
                    "unique_kmers": {"type": "number"},
                    "nullomers": {"type": "number"}
                }
            }
        },
        "comments": {
            "type": "array",
            "items": {"type": "string"}
        }
    }
}
