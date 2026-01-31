# Why identifiers are messy

- **Version suffixes**: Ensembl IDs can appear as `ENSG... .7` and often must be stripped.
- **One-to-many mappings**: a stable ID may map to multiple symbols (aliases/synonyms).
- **Annotation drift**: mappings depend on genome build and Ensembl/HGNC releases.

Best practice: keep the stable ID as the key and *add* symbols as annotation columns.

