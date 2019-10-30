import jsonschema

VERSION="0.0.1"

header_schema = {
    "type": "object",
    "properties": {
        "kdb_ver": {"type": "string"},
        "k": {"type": "number"},
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
                    "unique_kmers": {"type": "number"}
                }
            }
        },
        "comments": {
            "type": "array",
            "items": {"type": "string"}
        }
    }
}
