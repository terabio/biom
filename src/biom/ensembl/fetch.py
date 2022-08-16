import sys
from pathlib import Path
from subprocess import check_call

ENSEMBL_RELEASES = {
    107: "http://jul2022.archive.ensembl.org"
}


def transcripts(organism: str, version: int, saveto: Path, force: bool):
    if saveto.exists() and not force:
        return
    saveto.parent.mkdir(parents=True, exist_ok=True)
    query = f"""
        <?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
            <Dataset name = "{organism}_gene_ensembl" interface = "default" >
                <Attribute name = "ensembl_gene_id" />
                <Attribute name = "external_transcript_name" />
                <Attribute name = "gene_biotype" />
                <Attribute name = "transcript_biotype" />
                <Attribute name = "external_synonym" />
                <Attribute name = "transcript_gencode_basic" />
                <Attribute name = "transcript_tsl" />
                <Attribute name = "chromosome_name" />
                <Attribute name = "ensembl_gene_id_version" />
                <Attribute name = "ensembl_transcript_id" />
                <Attribute name = "ensembl_transcript_id_version" />
                <Attribute name = "transcript_start" />
                <Attribute name = "transcript_end" />
                <Attribute name = "strand" />
                <Attribute name = "transcript_is_canonical" />
                <Attribute name = "external_gene_name" />
            </Dataset>
        </Query>
    """
    query = [x.strip() for x in query.split('\n')]
    query = ''.join(query)

    if version not in ENSEMBL_RELEASES:
        raise FileNotFoundError(f"Too new/old Ensembl release: {version}.")
    url = f'{ENSEMBL_RELEASES[version]}/biomart/martservice?query={query}'
    print(url)

    assert saveto.name.endswith(".gz")
    # Trim gz suffix
    ungzipped = saveto.with_suffix('')
    try:
        check_call(['wget', '-O', ungzipped, url], stdout=sys.stdout)
        check_call(['gzip', ungzipped])
    except Exception as e:
        if ungzipped.exists():
            ungzipped.unlink()
        raise e
    assert saveto.exists()
