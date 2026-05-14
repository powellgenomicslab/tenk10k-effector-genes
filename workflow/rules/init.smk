# Initial file mapping
def get_file_map(host, type):
    config= f"config/path/{host}.yaml"
    with open(config, 'r') as f:
        filemap = yaml.full_load(f)

    d = {}
    
    if "directories" in filemap:
        if type == "target":
            d_dir = {k: directory(v[type]) for k,v in filemap['directories'].items()}
        else:
            d_dir = {k: v[type] for k,v in filemap['directories'].items()}
        d = {**d, **d_dir}
    
    if "files" in filemap:
        d_file = {k: v[type] for k,v in filemap['files'].items()}

        d = {**d, **d_file}
    
    return d
    
rule init_file_brenner:
    input: **get_file_map('brenner', 'source')
    output: **get_file_map('brenner', 'target')
    run:
        for i,o in zip(input, output):
            Path(str(o)).symlink_to(Path(i).resolve())

rule init_file_nci:
    input: **get_file_map('nci', 'source')
    output: **get_file_map('nci', 'target')
    run:
        for i,o in zip(input, output):
            Path(str(o)).symlink_to(Path(i).resolve())

# rule get_gene_universe:
#     input: lambda x: [str(f) for f in Path(f"resources/saige_eqtl/{x.study}").glob("**/common_raw.tsv")],
#     output: "resources/misc/{study}.genes.txt"
#     script: "snakescripts/get_gene_universe.R"

rule mk_gene_annot:
  """Make gene type annotation from gencode gtf"""
  input: "resources/misc/gencode.v44.basic.annotation.gtf"
  output: "resources/misc/gencode.{version}.gene_type.tsv"
  conda: "renv"
  script: "snakescripts/mk_gene_annot.R"

# Initialise trait metadata (update from Anne's and additional manual entries)
# rule init_trait_metadata:
#     output: "resources/metadata/trait_metadata_n.tsv"
#     conda: "renv"
#     script: "snakescripts/prep_trait_metadata.R"
