VERSION="0.0.2"

header_schema = {
    "type": "object",
    "properties": {
        "version": {"type": "string"},
        "metadata_blocks": {"type": "number"},
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
                    "nullomers": {"type": "number"},
                    "mononucleotides": {
                        "type": "object",
                        "properties": {
                            "A": {"type": "number"},
                            "C": {"type": "number"},
                            "G": {"type": "number"},
                            "T": {"type": "number"}
                        }
                    }
                },
                "required": ["filename", "sha256", "md5", "total_reads", "total_kmers", "unique_kmers", "nullomers"]
            }
        },
        "comments": {
            "type": "array",
            "items": {"type": "string"}
        }
    },
    "required": ["version", "metadata_blocks", "k", "tags", "files"]
}
