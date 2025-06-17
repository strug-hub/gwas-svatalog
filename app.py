import dash_bootstrap_components as dbc
import dash_daq as daq
import heapq
import numpy as np
import glob
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import polars as pl
import pyarrow

from dash import Dash, dcc, html, dash_table, Input, Output, State, ctx
from jupyter_dash import JupyterDash
from gtfparse import read_gtf # gtfparse==1.3.0


app = JupyterDash(__name__, 
           external_stylesheets = [dbc.themes.LUX, dbc.icons.BOOTSTRAP],
           title = "GWAS SVatalog")

app._favicon = "assets/favicon.ico"

server = app.server


DF_GWAS_FULL = pl.read_csv("data/gwas_catalog_v1.0-associations_e108.tsv",
                           separator = "\t",
                           quote_char = '',
                           dtypes = {"PUBMEDID": pl.Utf8, 
                                     "LINK": pl.Utf8,
                                     "P-VALUE": pl.Float64},
                           infer_schema_length = 1_000_000)


df_gwas = DF_GWAS_FULL.select([
    pl.col("CHR_ID"),
    pl.col("CHR_POS").cast(pl.Int32),
    pl.col("SNPS"),
    pl.col("STRONGEST SNP-RISK ALLELE"),
    pl.col("RISK ALLELE FREQUENCY"),
    pl.col("DISEASE/TRAIT"),
    pl.col("STUDY"),
    pl.col("PUBMEDID"),
    pl.col("LINK"),
    pl.col("P-VALUE")
]).drop_nulls(["CHR_ID"]).with_columns([
    (pl.lit("chr") + pl.col("CHR_ID")).alias("CHR_ID"),
    pl.col("SNPS").str.split(";").arr.first(),
    pl.col("STRONGEST SNP-RISK ALLELE").str.split("-").arr.last()
]).filter(pl.col("P-VALUE") <= 1)

# Rename columns in the output
df_gwas = df_gwas.rename({
    "CHR_ID": "Chromosome",
    "CHR_POS" : "SNP_Position",
    "SNPS": "SNP_Name_GWAS",
    "STRONGEST SNP-RISK ALLELE": "Risk_Allele",
    "RISK ALLELE FREQUENCY": "Risk_Allele_Frequency",
    "DISEASE/TRAIT": "Phenotype",
    "STUDY": "Study",
    "PUBMEDID": "Pubmed_ID",
    "LINK": "Link",
    "P-VALUE": "P-Value"
})


## SV Annotation data - created by Dr. Zhuozhi Wang on March 13, 2023
DF_ANNO_FULL = pl.read_csv("data/sv_annotations.tsv",
                           separator = "\t",
                           quote_char = '',
                           infer_schema_length = 1_000_000)

# Select and rename specific columns
df_anno = DF_ANNO_FULL.select([
    pl.col("Chromosome"),
    pl.col("Start").alias("SV_Start"),
    pl.col("End").cast(pl.Int64).alias("SV_End"),
    pl.col("ID").alias("SV_Name"),
    pl.col("Type").alias("SV_Type"),
    pl.col("Length").alias("SV_Length")
])


LD_file_paths = glob.glob(os.path.join("data/LD",
                                       "*_ld_stats.txt"))
LD_tables = [pl.read_csv(file_path,
                         separator = "\t",
                         has_header = True,
                         new_columns = ["SNP_Name", 
                                        "SNP_Position",
                                        "SV_Name",
                                        "SV_Position",
                                        "R2",
                                        "D'"])
             for file_path in LD_file_paths]

DF_LD = pl.concat(LD_tables)


## SV/SNP allele annotations done by Thomas Nalpathamkalam on April 10, 2023
allele_file_paths = glob.glob(os.path.join("data/LD",
                                            "*_allele_freq.txt"))
allele_tables = [pl.read_csv(file_path,
                             separator = "\t",
                             has_header = True,
                             new_columns = ["Chromosome",
                                            "SNP_Position",
                                            "SV_SNP_Name",
                                            "Reference_Allele",
                                            "Alternate_Allele",
                                            "Sample_AF",
                                            "gnomAD_nfe_AF",
                                            "dbSNP"])
                 for file_path in allele_file_paths]

DF_ALLELE = pl.concat(allele_tables)


## Gene/exon data from MANE.GRCh38.v1.0.ensembl_genomic.gtf extracted January 5, 2023
DF_GTF = read_gtf("data/MANE.GRCh38.v1.0.ensembl_genomic.gtf")

df_gene = pl.DataFrame(DF_GTF).filter(
    (pl.col("feature") == "transcript") | (pl.col("feature") == "exon")
).filter(pl.col("tag").str.contains("MANE_Select")).select([
    pl.col("seqname").alias("chr"),
    pl.col("start").cast(pl.Int64),
    pl.col("end").cast(pl.Int64),
    pl.col("strand"),
    pl.col("feature"),
    pl.col("gene_name").alias("gene")
])


#### DATA WRANGLING ####
## Subset allele data for only SV information
df_sv_allele = (DF_ALLELE
                .filter(pl.col("SV_SNP_Name").str.starts_with("P"))
                .rename({"SV_SNP_Name": "SV_Name"})
                .drop(["Reference_Allele", "SNP_Position", "gnomAD_nfe_AF", "dbSNP"])
                .unique())


## Subset allele data for only SNP information
df_snp_allele = (DF_ALLELE.filter(~pl.col("SV_SNP_Name").str.starts_with("P"))
                 .rename({"SV_SNP_Name": "SNP_Name"})
                 .unique())


## Subset DF_ANNO_FULL with columns that can be displayed to the public. Add SV AF to this table as well.
df_sv_anno = (DF_ANNO_FULL
              .join(
                  df_sv_allele.select(["SV_Name",
                                       "Alternate_Allele", 
                                       "Sample_AF"]),
                                       left_on = "ID",
                                       right_on = "SV_Name",
                                       how = "inner")
              .drop("Alternate_Allele")
              .select(["ID",
                       *DF_ANNO_FULL.columns[1:6],
                       "Sample_AF",
                       *DF_ANNO_FULL.columns[6:]])
              .rename({"ID": "SV Name",
                       "Sample_AF": "SV Sample AF"}))


## Creation of SV table for easier access of data required - join the dataframes based on SV names
df_sv_join = (DF_LD
              .join(df_anno, 
                    on = "SV_Name", 
                    how = "left")
              .select(["Chromosome", 
                       "SV_Name", 
                       "SV_Start", 
                       "SV_End",
                       "SNP_Name", 
                       "SNP_Position", 
                       "R2", 
                       "D'"])
              .join(df_sv_allele, 
                    on = ["Chromosome", 
                          "SV_Name"], 
                    how="left")
              .rename({"Alternate_Allele": "SV_Type",
                       "Sample_AF": "SV_AF"})
              .with_columns([pl.when(pl.col("SV_Type") == "<DEL>").then("Deletion")
                               .when(pl.col("SV_Type") == "<INS>").then("Insertion")
                               .when(pl.col("SV_Type") == "<DUP>").then("Duplication")
                               .when(pl.col("SV_Type") == "<INV>").then("Inversion")
                               .otherwise(pl.col("SV_Type"))
                               .alias("SV_Type")])
              .select(["Chromosome", 
                       "SV_Name", 
                       "SV_Start", 
                       "SV_End",
                       "SV_Type", 
                       "SV_AF", 
                       "SNP_Name", 
                       "SNP_Position", 
                       "R2", 
                       "D'"]))


## Merging SNP allele information with the LD dataset
df_sv_snp_join = (
    df_sv_join
    .join(
        df_snp_allele, 
        on=["Chromosome", "SNP_Name", "SNP_Position"], 
        how="left"
    )
    .drop("SNP_Name")
    .with_columns(
        pl.col("SNP_Position").cast(pl.Int32))
    .select([
        "Chromosome", "SV_Name", "SV_Start", "SV_End", "SV_Type", "SV_AF",
        "dbSNP", "SNP_Position", 
        "Reference_Allele", "Alternate_Allele", "Sample_AF", "gnomAD_nfe_AF", "R2", "D'"
    ])  # Reorder columns
    .rename({"dbSNP": "SNP_Name_dbSNP"})  # Rename column
)


## Creation of a full dataframe to make data extraction much easier
df_full_join = (df_sv_snp_join
                .join(df_gwas, 
                      on = ["Chromosome", 
                            "SNP_Position"], 
                      how = "left")
                .filter(pl.col("SV_AF").is_not_null())
                .with_columns(pl.when(pl.col("P-Value") < 1e-50)
                                .then(1e-50)
                                .otherwise(pl.col("P-Value"))
                                .alias("P-Value"))
                .with_columns((-pl.col("P-Value")
                                  .log10())
                                  .alias("P-Value_log10")))


# Subset df_anno for the SVs in df_full_foin
sv_names = df_full_join.select("SV_Name").unique()
df_anno_subset = df_anno.filter(pl.col("SV_Name")
                                  .is_in(sv_names["SV_Name"]))


## Creation of phenotype dictionary to encompass SV for each gene - computationally faster when creating SV table
df_phenotype = (
    df_full_join
    .groupby("Phenotype")
    .agg(pl.col("SV_Name").unique())  # Ensure SV_Name is unique within each group
)

dict_pheno = {
    row[0]: set(row[1]) for row in df_phenotype.iter_rows()
}


## Creation of y-coordinates to prevent overlap in genes when plotted
def assign_y_coordinates(df: pl.DataFrame,
                         y_start: float = -0.3,
                         y_increment: float = 0.5,
                         col_name: str = "y_coord") -> pl.DataFrame:
    
    # Priority queue to keep track of the end positions of genes for each y-coordinate
    df_filtered = (df.filter(pl.col("feature") == "transcript")
                     .sort(["chr", "start"]))
    
    rows = df_filtered.iter_rows()
    
     # Dictionary to store priority queues for each chromosome
    chromosome_heaps = {}

    # List to store the y-coordinate for each gene annotation
    y_coordinates = []

    for row in rows:
        chrom, start, end = row[0], row[1], row[2]

        # Initialize heap for new chromosomes
        if chrom not in chromosome_heaps:
            chromosome_heaps[chrom] = []

        heap = chromosome_heaps[chrom]

        # Check for available y-coordinates (no overlap)
        if heap and heap[0][0] <= start:
            _, y_coord = heapq.heappop(heap)
            heapq.heappush(heap, (end, y_coord))
        else:
            # Create a new y-coordinate if none available
            y_coord = heap[0][1] - y_increment if heap else y_start
            heapq.heappush(heap, (end, y_coord))

        y_coordinates.append(y_coord)

    # Add the y-coordinates as a new column
    return df_filtered.with_columns(pl.Series(col_name, y_coordinates))


df_transcript = assign_y_coordinates(df = df_gene,
                                     y_start = -0.3,
                                     y_increment = 0.15,
                                     col_name = "y_coord_nopheno")

df_transcript = assign_y_coordinates(df = df_transcript,
                                     y_start = -6.25,
                                     y_increment = 4,
                                     col_name = "y_coord_pheno")


# Add icon for help popup
question_icon = html.I(className = 'bi bi-question-diamond',
                       style = {"display" : "inline-block",
                                "color" : "lightsalmon",
                                "font-size" : "70%"})
                    #    id = 'question-icon')


#### FILTERS ####

# Make dropdown box for chromosomes and textbox for genomic range
REGION_LABEL = html.P("Genomic Region")

chromosomes = ["Any", *range(1, 23), "X", "Y"]

CHR_DD = dcc.Dropdown(id = 'chromosome-dropdown',
                      placeholder = "Select a chromosome",
                      options = [{"label" : str(i),
                                  "value" : str(i)} for i in chromosomes])

EMPTY_SPACE_REGION = html.Div(id = 'empty-space-region')

RANGE_START_TB = dcc.Input(id = 'range-start-textbox',
                           type = "number",
                           placeholder = "Start bp",
                           min = 1,
                           debounce = True,
                           className = 'class-range')

RANGE_END_TB = dcc.Input(id = 'range-end-textbox',
                         type = "number",
                         placeholder = "End bp",
                         min = 2,
                         debounce = True,
                         className = 'class-range')

REGION_DIV = html.Div(id = 'chromosome-dropdown-div',
                      children = [REGION_LABEL,
                                  CHR_DD,
                                  EMPTY_SPACE_REGION,
                                  RANGE_START_TB,
                                  RANGE_END_TB],
                      className = 'class-filter')


PHENO_LABEL = html.Div([html.P(["Phenotype",
                                question_icon],
                               id = 'pheno-label'),
                        dbc.Tooltip(dcc.Markdown("Select phenotype of interest availible in GWAS Catalog.\n\nSVs are selected if they have an LD calculation with an associated SNP for the phenotype.",
                                                 style = {'white-space':'pre-wrap'}),
                                    target = 'pheno-label',
                                    placement = "auto",
                                    delay = {"show" : 500,
                                             "hide" : 50},
                                    id = 'pheno-label-tooltip')])


phenotypes = sorted(df_full_join.select("Phenotype")
                                .unique()
                                .to_series()
                                .to_list())
phenotypes.insert(0, "Any")

PHENO_DD = dcc.Dropdown(id = 'phenotype-dropdown',
                        placeholder = "Select a phenotype of interest",
                        options = [{"label" : str(i),
                                  "value" : str(i)} for i in phenotypes],
                        optionHeight = 75)

PHENO_DD_DIV = html.Div(id = 'phenotype-dropdown-div',
                        children = [# PHENO_LINE,
                                    PHENO_LABEL,
                                    PHENO_DD],
                        className = 'class-section-filter-div')


# Make dropdown box for genes
GENE_LABEL = html.P("Gene")

genes = sorted(list(set(df_gene["gene"])))

GENE_DD = dcc.Dropdown(id = 'gene-dropdown',
                       placeholder = "Select a gene of interest",
                       options = [{"label" : str(i),
                                   "value" : str(i)} for i in genes])

GENE_DD_DIV = html.Div(id = 'gene-dropdown-div',
                       children = [GENE_LABEL,
                                   GENE_DD],
                       className = 'class-filter')


# Create a tab section for filter between gene and genomic region
TABS_FILTER_DIV = html.Div([dbc.Tabs([dbc.Tab(children = [GENE_DD_DIV],
                                         label = "by Gene",
                                         tab_id = "tab-gene"),
                                     dbc.Tab(children = [REGION_DIV],
                                             label = "by Region",
                                             tab_id = "tab-region")],
                                     id = "tab-filters",
                                     active_tab = "tab-gene")],
                           className = 'class-section-filter-div')

FILTER_DIV = html.Div(id = 'filter-div',
                      children = [TABS_FILTER_DIV,
                                  PHENO_DD_DIV],
                      className = 'class-filter-div')


#### TABLE ####

# Make DataTable
def make_table(chrom = chromosomes,
               range_start = 0,
               range_end = 250000000,
               gene = genes,
               pheno = phenotypes):

    global tab_show

    # Filter by chromosome
    if isinstance(chrom, list):
        sv_list_chr = df_anno_subset["SV_Name"].to_list()
    else:
        chrom = f"chr{chrom}"
        sv_list_chr = df_anno_subset.filter(pl.col("Chromosome") == chrom)["SV_Name"].to_list()

    # Filter by range
    sv_list_range = df_anno_subset.filter((pl.col("SV_Start") >= range_start) & \
                                          (pl.col("SV_End") <= range_end))["SV_Name"].to_list()

    # Filter by phenotype
    if isinstance(pheno, list):
        sv_list_pheno = df_anno_subset["SV_Name"].to_list()
    else:
        sv_list_pheno = dict_pheno[pheno]

    # Filter by gene
    if isinstance(gene, list):
        sv_list_gene = list(df_anno_subset["SV_Name"])
    else:
        df_1gene = df_gene.filter((pl.col("feature") == "transcript") & \
                                  (pl.col("gene") == gene))
        start = df_1gene["start"].min()
        end = df_1gene["end"].max()
        chrom = df_1gene.select("chr").to_series()[0]
        
        sv_list_gene = df_sv_anno.filter((pl.col("Chromosome") == chrom) & \
                                         (pl.col("Start") >= (start - 100000)) & \
                                         (pl.col("End") <= (end + 100000)))["SV Name"].to_list()

    # Find overlaps sv's from all lists
    sv_list_complete = set(sv_list_pheno).intersection(sv_list_chr, sv_list_gene, sv_list_range)

    tab_show = df_sv_anno.filter(pl.col("SV Name").is_in(sv_list_complete)).select([
        pl.col("SV Name").alias("ID"),
        pl.col("Chromosome").alias("Chrom"),
        pl.col("Start"),
        pl.col("End"),
        pl.col("Type"),
        pl.col("Length").alias("Size (bp)"),
        pl.col("SV Sample AF").alias("AF")])

    return dash_table.DataTable(id ='strvar-table',
                                data = tab_show.to_dicts(),
                                columns = [{'id': c, 'name': c} for c in ["Chrom", "Start", "End", "Type", "Size (bp)", "AF"]],
                                page_size = 100,
                                fixed_rows = {'headers': True},
                                style_table = {'height': '400px',
                                               'overflowY': 'auto'},
                                style_cell = {'textAlign': 'center',
                                              'font-family': 'Nunito Sans',
                                              'backgroundColor': 'whitesmoke',
                                              'minWidth' : '100px'},
                                style_header = {'backgroundColor': 'turquoise',
                                                'fontWeight': 'bold'},
                                row_selectable = 'single',
                                selected_rows = [])


# Make SV annotation table
def make_annotation_table(sv = "None"):
    if sv == "None":
        df_anno_subset_sv = pl.DataFrame()
        return dash_table.DataTable(id ='anno-table',
                                    data = [],
                                    columns = [{'id': c, 'name': c} for c in df_anno_subset_sv.columns],
                                    page_size = 400,
                                    fixed_rows = {'headers': True},
                                    style_table = {'height': '400px',
                                                   'overflowY': 'auto'},
                                    style_cell = {'textAlign': 'center',
                                                  'font-family': 'Nunito Sans',
                                                  'backgroundColor': 'whitesmoke'},
                                    style_header = {'backgroundColor': 'turquoise',
                                                    'fontWeight': 'bold'})

    else:
        df_anno_subset_sv = df_sv_anno.filter(pl.col("SV Name") == sv)
        df_anno_subset_sv = (df_anno_subset_sv.melt()
                                              .rename({"variable": "Header",
                                                       "value": "Information"}))
        
        return dash_table.DataTable(id ='anno-table',
                                    data = df_anno_subset_sv.to_dicts(),
                                    columns = [{'id': c, 'name': c} for c in df_anno_subset_sv.columns],
                                    page_size = 400,
                                    fixed_rows = {'headers': True},
                                    style_table = {'height': '400px',
                                                   'overflowX': 'auto',
                                                   'overflowY': 'auto'},
                                    style_data = {'whiteSpace' : 'normal',
                                                  'height' : 'auto'},
                                    style_cell = {'textAlign': 'left',
                                                  'font-family': 'Nunito Sans',
                                                  'backgroundColor': 'whitesmoke'},
                                    style_cell_conditional = [{'if': {'column_id': "Header"},
                                                              'backgroundColor': 'turquoise',
                                                              'fontWeight': 'bold',
                                                              'maxWidth' : '200px',
                                                              'whiteSpace' : 'normal',
                                                              'height' : 'auto'}],
                                    style_header = {'display' : 'none',
                                                    'height' : '0px'})


DTABLE = html.Div(id = 'strvar-table-div',
                  children = [make_table()],
                  className = 'sv-table')

# Create a header for the SV table
DTABLE_HEADER = html.Div([html.H4(["SV Information",
                                   question_icon],
                                   id = 'dtable-header',
                                   style = {'textAlign' : 'center',
                                            'color' : '#2F2E2D',
                                            'height' : '30 px',
                                            'text-transform': 'none'}),
                          dbc.Tooltip(dcc.Markdown("List of SVs subsetted by the filters applied.\n\nSelect an SV to visualize LD with GWAS-associated SNPs.",
                                                 style = {'white-space':'pre-wrap'}),
                                      target = 'dtable-header',
                                      placement = "auto",
                                      delay = {"show" : 500,
                                               "hide" : 50},
                                      id = 'dtable-header-tooltip')])

DTABLE_TEXTBOX = html.Div(id = 'dtable-textbox')

DTABLE_DIV = html.Div(id = 'dtable-div',
                      children = [DTABLE_HEADER,
                                  DTABLE,
                                  DTABLE_TEXTBOX],
                      className = 'class-filter-div')


ANNOTABLE = html.Div(id = 'anno-table-div',
                     children = [make_annotation_table()],
                     className = 'sv-anno-table')


# Create the SV details collapse and table
SV_TABLE_BUTTON = dbc.Button("SV Annotations",
                             id = "sv-button-collapse",
                             color = "secondary",
                             outline = True,
                             size = "sm",
                             n_clicks = 0)

# Create the filter division
ALL_FILTER_DIV = html.Div(id = 'all-filter-div',
                          children = [FILTER_DIV,
                                      DTABLE_DIV],
                                      style = {'display' : 'flex',
                                               'flexWrap' : "wrap",
                                               'justifyContent' : "space-evenly",
                                               'alignContent' : 'space-around'},
                                               className = 'sections')

SV_TABLE_BUTTON_DIV = html.Div(id = "sv-button-div",
                               children = [SV_TABLE_BUTTON])


SV_TABLE_COLLAPSE = dbc.Collapse(children = [ANNOTABLE],
                                 id = "sv-table-collapse",
                                 is_open = False)

SV_TABLE = html.Div(id = 'sv-table-collapse-div',
                    children = [SV_TABLE_COLLAPSE],
                    style = {'display' : 'flex',
                             'flexWrap' : "wrap",
                             'justifyContent' : "space-evenly",
                             'alignContent' : 'space-around'})


SV_TABLE_DIV = html.Div(id = 'sv-table-div',
                        children = [SV_TABLE,
                                    SV_TABLE_BUTTON_DIV],
                        className = 'sections')


#### PLOT ####

# Make scatter plot to compare r2 and D' values as seen in DataTable
def make_plot(sv = "None",
              pheno = "None",
              toggle = False,
              toggle_pheno = False):

    if toggle == False:
        linkage_value = "D'"
        size_value = "R2"
        legend_title = "D'"
    else:
        linkage_value = "R2"
        size_value = "D'"
        legend_title = "r<sup>2</sup>"

    global tab, sv_start_plot, sv_end_plot, extra_points, extra_points_pheno

    if sv == "None":
        tab = pl.DataFrame({"SNP_Position": [],
                            "-log10(P-Value)": []})
        
        sc_plot = px.scatter(data_frame = tab.to_pandas(),
                             x = "SNP_Position",
                             y = "-log10(P-Value)",
                             labels = {"-log10(P-Value)" : "-log<sub>10</sub> (P-Value)",
                                       "SNP_Position" : "Position in hg38 (bp)"},
                             width = 850,
                             height = 650,
                             template = "ggplot2")

        sc_plot.update_layout(margin = dict(l = 20, r = 20, t = 50, b = 50),
                              autosize = False,
                              font_family = "Nunito Sans",
                              font_size = 20)
    else:
        tab = df_full_join.filter(pl.col("SV_Name") == sv)

        # Display +/- 1KB from SV position
        sv_start = df_anno_subset.filter(pl.col("SV_Name") == sv).select("SV_Start").to_numpy()[0, 0]
        sv_end = df_anno_subset.filter(pl.col("SV_Name") == sv).select("SV_End").to_numpy()[0, 0]

        sv_start_plot = sv_start - 1000
        sv_end_plot = sv_end + 1000


        if pheno == "None":
            tab_lab = tab.to_dicts()

            def tab_label_maker(row):

                return [row["Phenotype"],
                        row["SNP_Name_dbSNP"],
                        row["Chromosome"],
                        row["SNP_Position"],
                        row["P-Value"],
                        row["R2"],
                        row["D'"],
                        row["SV_Name"],
                        row["SV_Start"],
                        row["SV_End"],
                        row["SV_Type"]]

            sc_plot = px.scatter(data_frame = tab.to_pandas(),
                                 x = "SNP_Position",
                                 y = linkage_value,
                                 labels = {"SNP_Position" : "Position in hg38 (bp)",
                                           "R2" : "r<sup>2</sup>"},
                                 width = 850,
                                 height = 650,
                                #  marginal_x = "rug",
                                 template = "ggplot2")

            sc_plot.update_traces(marker = dict(color = "rgba(57,201,187,0.5)",
                                                line = dict(color = "rgba(57,201,187,1)"),
                                                size = 10),
                                  hovertemplate = '<br><b>Phenotype:</b> %{customdata[0]}<br>' + 
                                                  '<br>SNP Name: %{customdata[1]}' + 
                                                  '<br>Chromosome: %{customdata[2]}' + 
                                                  '<br>SNP Position: %{customdata[3]}' + 
                                                  '<br>P-Value: %{customdata[4]}' + 
                                                  '<br>r2: %{customdata[5]}' + 
                                                  "<br>D': %{customdata[6]}",
                                  customdata = [tab_label_maker(row) for row in tab_lab])

            sc_plot.update_xaxes(showspikes = True,
                                 spikemode = "across")
            
            y_max = float(tab[linkage_value].max())
            ticks = np.linspace(start = 0, 
                                stop = y_max, 
                                num = 5)
            ticks = [round(t, 2) for t in ticks if t >= 0]
            sc_plot.update_yaxes(tickvals = ticks)

            sc_plot.update_layout(margin = dict(l = 20, r = 20, t = 50, b = 50),
                                  autosize = False,
                                  font_family = "Nunito Sans",
                                  font_size = 20)
                                #   xaxis_range = [sv_start_plot, sv_end_plot],
                                #   yaxis_range = [-0.6, 1])

             # Extract min and max SNP positions in chosen SV
            min_SNP = tab.select("SNP_Position").min().item() - 100000
            max_SNP = tab.select("SNP_Position").max().item() + 100000

            if sv_start < min_SNP:
                min_SNP = sv_start - 100000

            elif sv_start > max_SNP:
                max_SNP = sv_end  + 100000

            else:
                min_SNP = tab.select("SNP_Position").min()[0, 0] - 100000
                max_SNP = tab.select("SNP_Position").max()[0, 0] + 100000


            chromo = tab.select("Chromosome").unique().item()
            svname = tab.select("SV_Name").unique().item()


            sc_plot.add_shape(type = "rect",
                              x0 = sv_start,
                              x1 = sv_end,
                              y0 = -0.1,
                              y1 = -0.15,
                              line = dict(color = "darkslategrey",
                                           width = 1),
                              fillcolor = "lightsalmon")

            sc_plot.add_annotation(x = ((sv_end - sv_start) / 2) + sv_start,
                                   y = -0.01,
                                   text = svname,
                                   font = dict(size = 12),
                                   showarrow = False)

            sc_plot.add_trace(go.Scatter(x = np.array([((sv_end - sv_start) / 2) + sv_start]),
                                         y = np.array([-0.07]),
                                         mode = "markers",
                                         marker = dict(symbol = "triangle-down",
                                                       size = 15,
                                                       color = "lightsalmon"),
                                         showlegend = False,
                                         name = svname,
                                         hovertemplate = '<br><b>SV Name:</b> %{customdata[7]}<br>' + 
                                                         '<br>Chromosome: %{customdata[2]}' + 
                                                         '<br>Start: %{customdata[8]}' + 
                                                         '<br>End: %{customdata[9]}' + 
                                                         '<br>Type: %{customdata[10]}',
                                         customdata = [tab_label_maker(tab_lab[0])]))

            # Add extra points from df_sv with R2 and D' values
            extra_points = df_sv_snp_join.filter(pl.col("SV_Name") == sv)

            extra_points = extra_points.filter(~pl.col("SNP_Name_dbSNP").is_in(tab.select("SNP_Name_dbSNP").to_series()))

             # Create a list of labels for hovertext
            df_lab = extra_points.to_dicts()

            def extra_label_maker(row):

                return [row["SNP_Name_dbSNP"],
                        row["Chromosome"],
                        row["SNP_Position"],
                        row["R2"],
                        row["D'"]]

            sc_plot.add_trace(go.Scattergl(x = extra_points.select("SNP_Position").to_series().to_list(),
                                           y = extra_points.select(linkage_value).to_series().to_list(),
                                           mode = "markers",
                                           marker = dict(size = 10,
                                                         color = "rgba(168,168,168,0.5)",
                                                         line = dict(color = "rgba(168,168,168,1)")),
                                           hovertemplate = '<br><b>SNP Name: %{customdata[0]}</b>' +
                                                           '<br>Chromosome: %{customdata[1]}' + 
                                                           '<br>SNP Position: %{customdata[2]}' + 
                                                           '<br>r2: %{customdata[3]}' + 
                                                           "<br>D': %{customdata[4]}",
                                           customdata = [extra_label_maker(row) for row in df_lab],
                                           showlegend = False,
                                           name = "Not in \nGWAS Catalog"))

            # Filter exons and group them for gene-level information
            exons = df_gene.filter((pl.col("chr") == chromo) &
                                   (pl.col("start") >= min_SNP) &
                                   (pl.col("end") <= max_SNP))

            # Group to get gene-level start and end positions
            exons_genes = (exons.groupby(["chr", "gene", "strand"])
                                .agg([pl.col("start").min().alias("start"),
                                      pl.col("end").max().alias("end")]))

            # Filter for exons only
            exons_exons = exons.filter(pl.col("feature") == "exon")

            # Join to add y-coordinates (gene_line_position) to both dataframes
            exons_genes = exons_genes.join(df_transcript.select(["gene", "y_coord_nopheno"]),
                                           on = "gene",
                                           how = "left").rename({"y_coord_nopheno": "gene_line_position"})

            exons_exons = exons_exons.join(df_transcript.select(["gene", "y_coord_nopheno"]),
                                           on = "gene",
                                           how = "left").rename({"y_coord_nopheno": "gene_line_position"})

            # Prepare shapes for exons (rectangles)
            shapes = []
            for row in exons_exons.iter_rows(named=True):
                shapes.append(
                    dict(type = "rect",
                         x0 = row["start"],
                         x1 = row["end"],
                         y0 = row["gene_line_position"] + 0.02,
                         y1 = row["gene_line_position"] - 0.02,
                         line = dict(color = "mediumseagreen", 
                                     width = 2),
                         fillcolor = "mediumseagreen"))

            # Prepare shapes and annotations for genes (lines and text)
            annotations = []
            for row in exons_genes.iter_rows(named = True):
                # Line shape for the gene
                shapes.append(dict(type = "line",
                                   x0 = row["start"],
                                   x1 = row["end"],
                                   y0 = row["gene_line_position"],
                                   y1 = row["gene_line_position"],
                                   line = dict(width = 3, 
                                               color = "darkslategrey")))

                # Compute annotation details
                text_anno_x = round((row["end"] - row["start"]) / 2) + row["start"]
                direction = "\u2192" if row["strand"] == "+" else "\u2190"
                hover_text = f"{row['gene']}<br><b>Gene direction:</b> {'forward' if row['strand'] == '+' else 'reverse'}"

                # Add annotation
                annotations.append(dict(x = text_anno_x,
                                        y = row["gene_line_position"] + 0.060,
                                        text = f"{row['gene']}<span style='font-size:15px'>{direction}</span>",
                                        font = dict(size = 10),
                                        hovertext = hover_text,showarrow = False))

            # Add all shapes and annotations to the plot in bulk - computational time severely reduced
            sc_plot.update_layout(shapes = shapes)
            sc_plot.update_layout(annotations = annotations)
 

        else:
            tab = tab.filter(pl.col("Phenotype") == pheno)

            tab_lab = tab.to_dicts()

            def tab_label_maker(row):

                return [row["Phenotype"],
                        row["SNP_Name_dbSNP"],
                        row["Chromosome"],
                        row["SNP_Position"],
                        row["P-Value"],
                        row["R2"],
                        row["D'"],
                        row["SV_Name"],
                        row["SV_Start"],
                        row["SV_End"],
                        row["SV_Type"]]

            sc_plot = px.scatter(data_frame = tab.to_pandas(),
                                 x = "SNP_Position",
                                 y = "P-Value_log10",
                                 labels = {"SNP_Position" : "Position in hg38 (bp)",
                                           "P-Value_log10" : "-log<sub>10</sub> (P-Value)",
                                           linkage_value : legend_title},
                                 color = linkage_value,
                                 color_continuous_scale = [(0, "#4643CE"),
                                                           (0.2, "#4643CE"),
                                                           (0.2, "#90CAEE"),
                                                           (0.4, "#90CAEE"),
                                                           (0.4, "#46E17C"),
                                                           (0.6, "#46E17C"),
                                                           (0.6, "#EEC06B"),
                                                           (0.8, "#EEC06B"),
                                                           (0.8, "#CE5C5C"),
                                                           (1, "#CE5C5C")],
                                 range_color = [0, 1],
                                 width = 850,
                                 height = 650,
                                 template = "ggplot2")

            sc_plot.update_traces(hovertemplate = '<br><b>Phenotype:</b> %{customdata[0]}<br>' + 
                                                  '<br>SNP Name: %{customdata[1]}' + 
                                                  '<br>Chromosome: %{customdata[2]}' + 
                                                  '<br>SNP Position: %{customdata[3]}' + 
                                                  '<br>P-Value: %{customdata[4]}' + 
                                                  '<br>r2: %{customdata[5]}' + 
                                                  "<br>D': %{customdata[6]}",
                                  customdata = [tab_label_maker(row) for row in tab_lab],
                                  marker = dict(size = 10))

            sc_plot.update_xaxes(showspikes = True,
                                 spikemode = "across")
            
            y_max = float(tab["P-Value_log10"].max())
            ticks = np.linspace(start = 0, 
                                stop = y_max, 
                                num = 5)
            ticks = [round(t, 0) for t in ticks if t >= 0]
            sc_plot.update_yaxes(tickvals=ticks)

            sc_plot.update_layout(margin = dict(l = 20, r = 20, t = 50, b = 50),
                                 autosize = True,
                                 font_family = "Nunito Sans",
                                 font_size = 20,
                                 coloraxis_colorbar = dict(lenmode = "pixels",
                                                           len = 150,
                                                           yanchor = "top",
                                                           y = 1,
                                                           dtick = 0.2))

            # Extract min and max SNP positions in chosen SV
            min_SNP = tab.select("SNP_Position").min().item() - 100000
            max_SNP = tab.select("SNP_Position").max().item() + 100000

            if sv_start < min_SNP:
                min_SNP = sv_start - 100000

            elif sv_start > max_SNP:
                max_SNP = sv_end  + 100000

            else:
                min_SNP = tab.select("SNP_Position").min()[0, 0] - 100000
                max_SNP = tab.select("SNP_Position").max()[0, 0] + 100000


            chromo = tab.select("Chromosome").unique().item()
            svname = tab.select("SV_Name").unique().item()

            # Add points from gwas for pheno, without r2 and D' values
            extra_points = (df_gwas.filter((pl.col("Chromosome") == chromo) &
                                           (pl.col("SNP_Position") >= min_SNP) &
                                           (pl.col("SNP_Position") <= max_SNP) &
                                           (pl.col("Phenotype") == pheno))
                                   .filter(~pl.col("SNP_Name_GWAS").is_in(tab.select("SNP_Name_GWAS").to_series())))

            extra_points = (extra_points.with_columns(pl.when(pl.col("P-Value") < 1e-50)
                                                        .then(1e-50)
                                                        .otherwise(pl.col("P-Value"))
                                                        .alias("P-Value"),
                                                        (-pl.col("P-Value").log10()).alias("P-Value_log10")))

            # Create a list of labels for hovertext
            df_lab = extra_points.to_dicts()

            def label_maker(row):

                return [row["Phenotype"],
                        row["SNP_Name_GWAS"],
                        row["Chromosome"],
                        row["SNP_Position"],
                        row["P-Value"]]


            sc_plot.add_trace(go.Scattergl(x = extra_points.select("SNP_Position").to_series().to_list(),
                                           y = extra_points.select("P-Value_log10").to_series().to_list(),
                                           mode = "markers",
                                           marker = dict(size = 10,
                                                         color = "rgba(168,168,168,0.5)",
                                                         line = dict(color = "rgba(168,168,168,1)")),
                                           hovertemplate = '<br><b>Phenotype:</b> %{customdata[0]}<br>' + 
                                                           '<br>SNP Name: %{customdata[1]}' + 
                                                           '<br>Chromosome: %{customdata[2]}' + 
                                                           '<br>SNP Position: %{customdata[3]}' + 
                                                           '<br>P-Value: %{customdata[4]}',
                                           customdata = [label_maker(row) for row in df_lab],
                                           showlegend = False,
                                           name = "No LD Data"))

            # Add points from other phenos in LD with SV
            if toggle_pheno == True:
                extra_points_pheno = (df_full_join.filter((pl.col("Chromosome") == chromo) &
                                                          (pl.col("SNP_Position") >= min_SNP) &
                                                          (pl.col("SNP_Position") <= max_SNP) &
                                                          (pl.col("SV_Name") == sv))
                                                  .filter(~pl.col("SNP_Name_GWAS").is_in(tab.select("SNP_Name_GWAS").to_series())))


                # Create a list of labels for hovertext
                df_lab_2 = extra_points_pheno.to_dicts()

                def label_maker(row):

                    return [row["Phenotype"],
                            row["SNP_Name_GWAS"],
                            row["Chromosome"],
                            row["SNP_Position"],
                            row["P-Value"],
                            row["R2"],
                            row["D'"]]


                sc_plot.add_trace(go.Scattergl(x = extra_points_pheno.select("SNP_Position").to_series().to_list(),
                                               y = extra_points_pheno.select("P-Value_log10").to_series().to_list(),
                                               mode = "markers",
                                               marker = dict(size = 10,
                                                             color = "rgba(168,168,168,0.5)",
                                                             line = dict(color = "rgba(168,168,168,1)")),
                                               hovertemplate = '<br><b>Phenotype:</b> %{customdata[0]}<br>' + 
                                                               '<br>SNP Name: %{customdata[1]}' + 
                                                               '<br>Chromosome: %{customdata[2]}' + 
                                                               '<br>SNP Position: %{customdata[3]}' + 
                                                               '<br>P-Value: %{customdata[4]}' + 
                                                               '<br>r2: %{customdata[5]}' + 
                                                               "<br>D': %{customdata[6]}",
                                               customdata = [label_maker(row) for row in df_lab_2],
                                               showlegend = False,
                                               name = "Other Phenotypes"))
            else:
                pass


            sc_plot.add_shape(type = "rect",
                              x0 = sv_start,
                              x1 = sv_end,
                              y0 = -2,
                              y1 = -3,
                              line = dict(color = "darkslategrey",
                                          width = 1),
                              fillcolor = "lightsalmon")

            sc_plot.add_annotation(x = ((sv_end - sv_start) / 2) + sv_start,
                                   y = -0.5,
                                   text = svname,
                                   font = dict(size = 12),
                                   showarrow = False)

            sc_plot.add_trace(go.Scatter(x = np.array([((sv_end - sv_start) / 2) + sv_start]),
                                         y = np.array([-1.5]),
                                         mode = "markers",
                                         marker = dict(symbol = "triangle-down",
                                                       size = 15,
                                                       color = "lightsalmon"),
                                         showlegend = False,
                                         name = svname,
                                         hovertemplate = '<br><b>SV Name:</b> %{customdata[7]}<br>' + 
                                                         '<br>Chromosome: %{customdata[2]}' + 
                                                         '<br>Start: %{customdata[8]}' + 
                                                         '<br>End: %{customdata[9]}' + 
                                                         '<br>Type: %{customdata[10]}',
                                         customdata = [tab_label_maker(tab_lab[0])]))
            
            # Filter exons and group them for gene-level information
            exons = df_gene.filter((pl.col("chr") == chromo) &
                                   (pl.col("start") >= min_SNP) &
                                   (pl.col("end") <= max_SNP))

            # Group to get gene-level start and end positions
            exons_genes = (exons.groupby(["chr", "gene", "strand"])
                                .agg([pl.col("start").min().alias("start"),
                                      pl.col("end").max().alias("end")]))

            # Filter for exons only
            exons_exons = exons.filter(pl.col("feature") == "exon")

            # Join to add y-coordinates (gene_line_position) to both dataframes
            exons_genes = exons_genes.join(df_transcript.select(["gene", "y_coord_pheno"]),
                                           on = "gene",
                                           how = "left").rename({"y_coord_pheno": "gene_line_position"})

            exons_exons = exons_exons.join(df_transcript.select(["gene", "y_coord_pheno"]),
                                           on = "gene",
                                           how = "left").rename({"y_coord_pheno": "gene_line_position"})

            # Prepare shapes for exons (rectangles)
            shapes = []
            for row in exons_exons.iter_rows(named=True):
                shapes.append(
                    dict(type = "rect",
                         x0 = row["start"],
                         x1 = row["end"],
                         y0 = row["gene_line_position"] + 0.75,
                         y1 = row["gene_line_position"] - 0.75,
                         line = dict(color = "mediumseagreen", 
                                     width = 2),
                         fillcolor = "mediumseagreen"))

            # Prepare shapes and annotations for genes (lines and text)
            annotations = []
            for row in exons_genes.iter_rows(named = True):
                # Line shape for the gene
                shapes.append(dict(type = "line",
                                   x0 = row["start"],
                                   x1 = row["end"],
                                   y0 = row["gene_line_position"],
                                   y1 = row["gene_line_position"],
                                   line = dict(width = 3, 
                                               color = "darkslategrey")))

                # Compute annotation details
                text_anno_x = round((row["end"] - row["start"]) / 2) + row["start"]
                direction = "\u2192" if row["strand"] == "+" else "\u2190"
                hover_text = f"{row['gene']}<br><b>Gene direction:</b> {'forward' if row['strand'] == '+' else 'reverse'}"

                # Add annotation
                annotations.append(dict(x = text_anno_x,
                                        y = row["gene_line_position"] + 1.75,
                                        text = f"{row['gene']}<span style='font-size:15px'>{direction}</span>",
                                        font = dict(size = 10),
                                        hovertext = hover_text,showarrow = False))

            # Add all shapes and annotations to the plot in bulk - computational time severely reduced
            sc_plot.update_layout(shapes = shapes)
            sc_plot.update_layout(annotations = annotations)

    return sc_plot

SCATPLOT = dcc.Loading(id = "plot-loading",
                       children = [dcc.Graph(id = "scatter-plot",
                                             figure = make_plot(),
                                             config = {"modeBarButtonsToRemove": ["select2d",
                                                                                  "lasso2d",
                                                                                  "zoomIn2d",
                                                                                  "zoomOut2d"],
                                                        "toImageButtonOptions": {"format": "svg",
                                                                                 "filename": "gwas_svatalog_plot",
                                                                                 "height": 500,
                                                                                 "width": 700,
                                                                                 "scale": 1}})],
                       type = "dot",
                       color = "turquoise",
                       parent_className = "loading_wrapper")


TOGGLE_SWITCH = html.Div([dbc.Row([dbc.Col(["D'"],
                                           id = 'toggle-d'),
                                   dbc.Col(daq.ToggleSwitch(id = 'toggle-switch',
                                                            color = "lightsalmon")),
                                   dbc.Col(["r", html.Sup(2,
                                                          style = {'font-size' : '15px'})],
                                            id = 'toggle-r')],
                                   align = "center",
                                   id = 'toggle-region'),
                          dbc.Row([html.P(["Show Other Phenotypes",
                                           question_icon],
                                          id = 'toggle-pheno-header'),
                                   dbc.Tooltip(dcc.Markdown("Show GWAS-associated SNPs from other phenotypes in respect to selected SV.",
                                                            style = {'white-space':'pre-wrap'}),
                                               target = 'toggle-pheno-header',
                                               placement = "auto",
                                               delay = {"show" : 500,
                                                        "hide" : 50},
                                               id = 'toggle-pheno-header-tooltip')],
                                  id = 'toggle-pheno-header-row'),
                          dbc.Row([dbc.Col(["OFF"],
                                           id = 'toggle-off'),
                                   dbc.Col(daq.ToggleSwitch(id = 'other-pheno-switch',
                                                            color = "lightsalmon")),
                                   dbc.Col(["ON"],
                                            id = 'toggle-on')],
                                   align = "center",
                                   id = 'toggle-region-pheno'),
                          dbc.Row([dbc.Button("Export SNP Data to CSV",
                                               id = 'download-button',
                                               color = "secondary",
                                               outline = True,
                                               size = "sm",),
                                   dcc.Download(id = 'download-csv')],
                                   id = 'download-region')])


SCATPLOT_DIV = html.Div(id = 'scatter-plot-div',
                        children = [SCATPLOT,
                                    TOGGLE_SWITCH],
                        className = 'sections')


# Make SNP table from GWAS data based on clickData from plot
def make_SNP_table(snp = "None"):
    if snp == "None":
        df_snp = pl.DataFrame()
        return dash_table.DataTable(id ='snp-table',
                                    data = [],
                                    columns = [],
                                    page_size = 400,
                                    fixed_rows = {'headers': True},
                                    style_table = {'height': '500px', 'overflowY': 'auto'},
                                    style_cell = {'textAlign': 'center',
                                                  'font-family': 'Nunito Sans',
                                                  'backgroundColor': 'whitesmoke'},
                                    style_header = {'backgroundColor': 'turquoise',
                                                    'fontWeight': 'bold'})

    else:
        df_snp = (df_full_join.filter((pl.col("SNP_Name_dbSNP") == snp) | (pl.col("SNP_Name_GWAS") == snp))
                              .with_columns([
                pl.format("[{}](https://{})", pl.col("Pubmed_ID"), pl.col("Link")).alias("Link"),
                pl.format("[{}](https://www.ncbi.nlm.nih.gov/snp/{})", pl.col("SNP_Name_dbSNP"), pl.col("SNP_Name_dbSNP")).alias("SNP_Name_dbSNP"),
                pl.format("[{}](https://www.ebi.ac.uk/gwas/variants/{})", pl.col("SNP_Name_GWAS"), pl.col("SNP_Name_GWAS")).alias("SNP_Name_GWAS")])
                              .select([
                pl.col("Chromosome").alias("Chrom"),
                pl.col("SNP_Position").alias("SNP Position"),
                pl.col("SNP_Name_dbSNP").alias("SNP Name: dbSNP"),
                pl.col("SNP_Name_GWAS").alias("SNP Name: GWAS"),
                pl.col("Reference_Allele").alias("Reference Allele"),
                pl.col("Alternate_Allele").alias("Alternate Allele"),
                pl.col("Risk_Allele").alias("Risk Allele"),
                pl.col("Risk_Allele_Frequency").alias("Risk AF"),
                pl.col("Sample_AF").alias("Sample AF"),
                pl.col("gnomAD_nfe_AF").alias("gnomAD NFE AF"),
                pl.col("Phenotype"),
                pl.col("P-Value"),
                pl.col("Study"),
                pl.col("Link").alias("Pubmed Link")])
                               .unique())

        return dash_table.DataTable(id ='snp-table',
                                    data = df_snp.to_dicts(),
                                    columns = [{'id': c, 'name': c, 'presentation': 'markdown'} 
                                               if (c == "Pubmed Link" or c == "SNP Name: dbSNP" or c == "SNP Name: GWAS") 
                                               else ({'id': c, 'name': c, 'type':'numeric', 'format': {'specifier': '.3f'}} 
                                                     if (c == "Risk AF" or c == "Sample AF" or c == "gnomAD NFE AF") 
                                                     else {'id': c, 'name': c}) for c in df_snp.columns],
                                    markdown_options = {"html" : True,
                                                        "link_target": "_blank"},
                                    page_size = 400,
                                    fixed_rows = {'headers': True},
                                    style_table = {'height': '500px',
                                                   'minWidth': '100%',
                                                   'width' : '100%',
                                                   'maxWidth' : '100%',
                                                   'overflowX' : 'auto',
                                                   'overflowY': 'auto'},
                                    style_data = {'whiteSpace' : 'normal',
                                                  'height' : 'auto'},
                                    style_cell_conditional = [{'if': {'column_id': ["Chrom",
                                                                                    "SNP Position",
                                                                                    "SNP Name: dbSNP",
                                                                                    "SNP Name: GWAS",
                                                                                    "Reference Allele",
                                                                                    "Alternate Allele",
                                                                                    "Risk Allele",
                                                                                    "Risk AF",
                                                                                    "Sample AF",
                                                                                    "gnomAD NFE AF",
                                                                                    "P-Value",
                                                                                    "Pubmed Link"]},
                                                               'minWidth' : '100px',
                                                               'textAlign': 'left',
                                                               'font-family': 'Nunito Sans',
                                                               'backgroundColor': 'whitesmoke'},
                                                               {'if': {'column_id': ["Study",
                                                                                     "Phenotype"]},
                                                               'maxWidth' : '200px',
                                                               'textAlign': 'left',
                                                               'font-family': 'Nunito Sans',
                                                               'backgroundColor': 'whitesmoke',
                                                               'whiteSpace' : 'normal',
                                                               'height' : 'auto'}],
                                    style_header = {'backgroundColor': 'turquoise',
                                                    'fontWeight': 'bold',
                                                    'textAlign' : 'left',
                                                    'whiteSpace' : 'normal',
                                                    'height' : 'auto'})

SNP_TABLE = html.Div([html.Div(id = 'snp-table-div',
                               children = [make_SNP_table()],
                               className = 'sections'),
                      dbc.Modal([dbc.ModalHeader(),
                                 dbc.ModalBody("SNP table populated below.")],
                                 id = 'modal-snp-table',
                                 is_open = False)])


#### HEADERS ####

## Navigation Bar to match LocusFocus
NAVBAR = dbc.NavbarSimple(id = 'navbar-lf',
                          brand = "LocusFocus",
                          brand_href = "https://locusfocus.research.sickkids.ca/",
                          brand_external_link = "_blank",
                          color = "#40e0d0",
                          links_left = True,
                          expand = "lg",
                          children = [dbc.NavItem(dbc.NavLink("Colocalization",
                                                              href = "https://locusfocus.research.sickkids.ca/",
                                                              target = "_blank",
                                                              class_name = 'navbar-items')),
                                      dbc.NavItem(dbc.NavLink("Set-Based Test",
                                                              href = "https://locusfocus.research.sickkids.ca/setbasedtest",
                                                              target = "_blank",
                                                              class_name = 'navbar-items')),
                                      dbc.NavItem(dbc.NavLink("GWAS SVatalog",
                                                              href = "https://svatalog.research.sickkids.ca/",
                                                              class_name = 'navbar-items')),
                                      dbc.NavItem(dbc.NavLink("Documentation",
                                                              href = "https://gwas-svatalog-docs.readthedocs.io/en/latest/index.html",
                                                              target = "_blank",
                                                              class_name = 'navbar-items')),
                                      dbc.NavItem([dbc.Button("Contact Us",
                                                              id = 'button-modal-contact-us',
                                                              n_clicks = 0,
                                                              class_name = 'navbar-items'),
                                                   dbc.Modal([dbc.ModalHeader(dbc.ModalTitle("Contact Us")),
                                                              dbc.ModalBody(dcc.Markdown('''
                                                                                         <h5 children="GWAS SVatalog"  />
                                                                                         Shalvi Chirmade <span style="color: #2ba089" children="shalvi.chirmade@sickkids.ca" />

                                                                                         &nbsp;

                                                                                         <h5 children="LocusFocus" />
                                                                                         Mackenzie Frew <span style="color: #2ba089" children="mackenzie.frew@sickkids.ca" />
                                                                                        ''',
                                                                                        dangerously_allow_html = True,
                                                                                        id = 'contact-us-mrkdwn'))],
                                                              id = 'modal-contact-us',
                                                              is_open = False)]),
                                      dbc.NavItem(dbc.NavLink("Subscribe",
                                                              href = "https://mailchi.mp/752ab1c4d516/locusfocus",
                                                              target = "_blank",
                                                              class_name = 'navbar-items')),
                                      dbc.NavItem([dbc.Button("Citation",
                                                              id = 'button-modal-citation',
                                                              n_clicks = 0,
                                                              class_name = 'navbar-items'),
                                                   dbc.Modal([dbc.ModalHeader(dbc.ModalTitle("Citation")),
                                                              dbc.ModalBody(dcc.Markdown('''
                                                                                         *to be published*
                                                                                        '''))],
                                                              id = 'modal-citation',
                                                              is_open = False)])])


## Logo for the top of the webpage, disclaimer and info
GSV_LOGO = html.Img(src = "assets/gwas-svatalog-name.png",
                    id = 'gsv-logo-image')

TCAG_LOGO = html.Img(src = "assets/tcaglogo.png",
                     id = 'tcag-logo-image')

SK_LOGO = html.Img(src = "assets/SKlogo.svg",
                   id = 'sk-logo-image')

DISCLAIMER_BUTTON = dbc.Button("Disclaimer",
                               color = "secondary",
                               outline = True,
                               size = "sm",
                               id = 'disclaimer-button')

DISCLAIMER_COLLAPSE = dbc.Collapse(dbc.Card(dbc.CardBody("Database constructed from predominantly European population of 101 individuals with Cystic Fibrosis (CF). The alleles affected by CF aside, the remainder of the genome is comparable to a healthy population of European descent (citation TBD). Genomic location is referenced against GRCh38.")),
                                   id = 'disclaimer-collapse',
                                   is_open = False)

DISCLAIMER_DIV = html.Div(id = 'disclaimer-div',
                          children = [DISCLAIMER_BUTTON,
                                      DISCLAIMER_COLLAPSE])

GSV_LOGO_DIV = html.Div(id = 'gsv-logo-image-div',
                        children = [GSV_LOGO])

COMP_LOGO_DIV = html.Div(id = 'comp-logo-image-div',
                         children = [TCAG_LOGO,
                                     SK_LOGO,
                                     DISCLAIMER_DIV])

LOGO_DIV = html.Div(id = 'logo-image-div',
                    children = [GSV_LOGO_DIV,
                                COMP_LOGO_DIV])

HEADER_DIV = html.Div(id = 'header-line-div')

# Header for filters
FILTER_HEADER = html.H3("Search for Structural Variants",
                        id = 'filter-header',
                        style = {'textAlign' : 'center',
                                 'color' : '#2F2E2D',
                                 'height' : '70 px'})

FILTER_HEADER_DIV = html.Div(id = 'filter-header-div',
                             children = [FILTER_HEADER],
                             className = 'section-headers')



# Header for plot
PLOT_HEADER = html.Div([html.H3(["SVs and GWAS Hits Linkage Disequilibirum Plot",
                                 question_icon],
                                 id = 'plot-header',
                                 style = {'textAlign' : 'center',
                                          'color' : '#2F2E2D',
                                          'height' : '70 px'}),
                        dbc.Tooltip(dcc.Markdown("Interactive plot visualizing linkage disequilibrium (LD) between the selected SV and GWAS-associated SNPs from GWAS Catalog.",
                                                 style = {'white-space':'pre-wrap'}),
                                    target = 'plot-header',
                                    placement = "auto",
                                    delay = {"show" : 500,
                                             "hide" : 50},
                                    id = 'plot-header-tooltip')])

PLOT_HEADER_DIV = html.Div(id = 'plot-header-div',
                           children = [PLOT_HEADER],
                           className = 'section-headers')


#### CALLBACKS ####

# Incorporate together to create interactive webpage
app.layout = html.Div([NAVBAR,
                       LOGO_DIV,
                       HEADER_DIV,
                       FILTER_HEADER_DIV,
                       ALL_FILTER_DIV,
                       SV_TABLE_DIV,
                       PLOT_HEADER_DIV,
                       SCATPLOT_DIV,
                       SNP_TABLE])

# Contact Us modal
@app.callback(
    Output('modal-contact-us', 'is_open'),
    [Input('button-modal-contact-us', 'n_clicks')],
    [State('modal-contact-us', 'is_open')],
)
def toggle_modal(modal_contact, is_open):
    if modal_contact:
        return not is_open
    return is_open


# Citation modal
@app.callback(
    Output('modal-citation', 'is_open'),
    [Input('button-modal-citation', 'n_clicks')],
    [State('modal-citation', 'is_open')],
)
def toggle_modal(modal_citation, is_open):
    if modal_citation:
        return not is_open
    return is_open


# Disclaimer cutton opens collapsed text
@app.callback(
    Output('disclaimer-collapse', 'is_open'),
    [Input('disclaimer-button', 'n_clicks')],
    [State('disclaimer-collapse', 'is_open')]
)
def dislcaimer_collapse(clicks, is_open):
    if clicks:
        return not is_open
    return is_open


# Chromosome dropdown changes the genes and phenotypes displayed
@app.callback(
    Output('gene-dropdown', 'options'),
    Output('phenotype-dropdown', 'options'),
    Input('chromosome-dropdown', 'value')
)
def update_gene_dropdown(chrom):

    if chrom is None:
        new_genes = genes
        new_pheno = phenotypes
    elif chrom == "Any":
        new_genes = genes
        new_pheno = phenotypes
    else:
        chromosome = f"chr{chrom}"
        new_genes = df_gene.filter(pl.col("chr") == chromosome).select("gene").unique().to_series().to_list()
        new_genes = sorted(new_genes)
        new_pheno = df_full_join.filter(pl.col("Chromosome") == chromosome).select("Phenotype").unique().to_series().to_list()
        new_pheno = sorted(map(str, new_pheno))

    return [{'label': str(i), 'value': str(i)} for i in new_genes], [{'label': str(i), 'value': str(i)} for i in new_pheno]


# Filters alter the table displayed
@app.callback(
    Output('strvar-table-div', 'children'),
    Input('chromosome-dropdown', 'value'),
    Input('phenotype-dropdown', 'value'),
    Input('gene-dropdown', 'value'),
    Input('range-start-textbox', 'value'),
    Input('range-end-textbox', 'value')
)
def update_table(chrom, pheno, gene, range_start, range_end):

    chrom = chromosomes if chrom is None else chrom
    pheno = phenotypes if pheno is None else pheno
    gene = genes if gene is None else gene
    range_start = 0 if range_start is None else range_start
    range_end = 250000000 if range_end is None else range_end

    if chrom == "Any":
        chrom = chromosomes
    else:
        chrom = chrom
    
    if pheno == "Any":
        pheno = phenotypes
    else:
        pheno = pheno

    new_df = make_table(chrom = chrom,
                        range_start = range_start,
                        range_end = range_end,
                        pheno = pheno,
                        gene = gene)

    return new_df


@app.callback(
    Output('dtable-textbox', 'children'),
    Input('strvar-table', 'derived_virtual_data'),
    Input('strvar-table', 'selected_rows'),
    Input('range-start-textbox', 'value'),
    Input('range-end-textbox', 'value')
)
def update_sv_textbox(data, selected_rows, range_start, range_end):

    range_start = 0 if range_start is None else range_start
    range_end = 250000000 if range_end is None else range_end

    if range_end <= range_start:
        return html.P("Fix range values.")

    if data == []:
        return html.P("No SVs in selected filters.")
    
    elif selected_rows:
        return html.P("Plot generated below.")
    
    else:
        return html.P("")


@app.callback(
    Output('anno-table-div', 'children'),
    Input('strvar-table', 'selected_rows')
)
def update_anno_table(selected_rows):

    if selected_rows:
        sv_selected = tab_show[selected_rows[0], 0]
        new_anno_table = make_annotation_table(sv = sv_selected)

    else:
        new_anno_table = make_annotation_table(sv = "None")

    return new_anno_table


# Collapse SV annotation table
@app.callback(
    Output("sv-table-collapse", "is_open"),
    [Input("sv-button-collapse", "n_clicks")],
    [State("sv-table-collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


# Only row selection on table changes plot, SV div textbox and button enabling
@app.callback(
    Output('scatter-plot', 'figure'),
    Output('scatter-plot', 'clickData'),
    Output('download-button', 'disabled'),
    Output('other-pheno-switch', 'disabled'),
    Input('strvar-table', 'selected_rows'),
    Input('phenotype-dropdown', 'value'),
    Input('toggle-switch', 'value'),
    Input('other-pheno-switch', 'value')
)
def update_plot(selected_rows, pheno, toggle, toggle_pheno):

    pheno = "None" if pheno is None else pheno
    toggle = False if toggle is None else toggle
    toggle_pheno = False if toggle_pheno is None else toggle_pheno

    if selected_rows:
        sv_name = tab_show[selected_rows[0], 0]
        button = False

    else:
        sv_name = "None"
        button = True

    if pheno == "None":
        pheno_toggle_disable = True

    else:
        pheno_toggle_disable = False

    new_plot = make_plot(sv = sv_name,
                         pheno = pheno,
                         toggle = toggle,
                         toggle_pheno = toggle_pheno)

    return new_plot, {}, button, pheno_toggle_disable


# Update SNP table based on data point clicked on plot
@app.callback(
        Output('snp-table-div', 'children'),
        Output('modal-snp-table', 'is_open'),
        Input('scatter-plot', 'clickData'),
        Input('strvar-table', 'selected_rows'),
        State('modal-snp-table', 'is_open'),
    prevent_initial_call = True
)
def update_SNP_table(clickData, selected_rows, is_open):

    if clickData:
        if clickData == {}:
            rsID = "None"
            is_open = False
        else:
            rsID = clickData["points"][0]["customdata"][1]
            is_open = True

    else:
        rsID = "None"
        is_open = False

    if not selected_rows:
        rsID = "None"
        is_open = False

    new_snp_table = make_SNP_table(snp = rsID)

    return new_snp_table, is_open


# Download SNP data seen in plot
@app.callback(
    Output('download-csv', 'data'),
    Input('download-button', 'n_clicks'),
    Input('scatter-plot', 'relayoutData'),
    Input('toggle-switch', 'value'),
    Input('other-pheno-switch', 'value'),
    prevent_initial_call = True
)
def download_plot_data(n_clicks, relayoutData, toggle, toggle_pheno):
    # Example from https://towardsdatascience.com/building-a-dashboard-in-plotly-dash-c748588e2920
    global df_download

    if toggle == False:
        linkage = "D'"
    elif toggle == True:
        linkage = "R2"
    else:
        linkage = "D'"
        
    # Function to concat all dfs using polars
    def align_and_reorder_columns(dfs, reference_df):
    
        reference_columns = reference_df.columns
        reference_schema = reference_df.schema

        aligned_dfs = []
        for df in dfs:
            missing_columns = [col for col in reference_columns if col not in df.columns]
            df_with_all_columns = df.with_columns([pl.lit(None, dtype = reference_schema[col]).alias(col)
                                                   for col in missing_columns])
            
            df_with_all_columns = df_with_all_columns.select(['Chromosome',
                                                                'SV_Name',
                                                                'SV_Start',
                                                                'SV_End',
                                                                'SV_Type',
                                                                'SV_AF',
                                                                'SNP_Name_GWAS',
                                                                'SNP_Name_dbSNP',
                                                                'SNP_Position',
                                                                'Reference_Allele',
                                                                'Alternate_Allele',
                                                                'Sample_AF',
                                                                'gnomAD_nfe_AF',
                                                                'Risk_Allele',
                                                                'Risk_Allele_Frequency',
                                                                "D'",
                                                                'R2',
                                                                'Phenotype',
                                                                'P-Value',
                                                                'P-Value_log10',
                                                                'Pubmed_ID',
                                                                'Study',
                                                                'Link'])
            aligned_dfs.append(df_with_all_columns)

        concatenated_df = pl.concat(aligned_dfs)

        return concatenated_df

        
    # If button was triggered
    if ctx.triggered[0]['prop_id'] == 'download-button.n_clicks':
        # Default range (no zoom)
        if relayoutData is None or "xaxis.autorange" in relayoutData:
            df_download = tab
            df_extra = extra_points

            if toggle_pheno:
                df_extra_pheno = extra_points_pheno
                df_download = align_and_reorder_columns([df_download, df_extra, df_extra_pheno], 
                                                        reference_df = df_download)
            else:
                df_download = align_and_reorder_columns([df_download, df_extra], 
                                                        reference_df = df_download)

        # Filter by x-axis and y-axis ranges
        elif "xaxis.range[0]" in relayoutData:
            x_min, x_max = relayoutData["xaxis.range[0]"], relayoutData["xaxis.range[1]"]

            if "yaxis.range[0]" in relayoutData:
                y_min, y_max = relayoutData["yaxis.range[0]"], relayoutData["yaxis.range[1]"]

                if toggle_pheno:
                    df_download = tab.filter((pl.col("SNP_Position") >= x_min) &
                                             (pl.col("SNP_Position") <= x_max) &
                                             (pl.col("P-Value_log10") >= y_min) &
                                             (pl.col("P-Value_log10") <= y_max))
                    
                    df_extra = extra_points.filter((pl.col("SNP_Position") >= x_min) &
                                                   (pl.col("SNP_Position") <= x_max) &
                                                   (pl.col("P-Value_log10") >= y_min) &
                                                   (pl.col("P-Value_log10") <= y_max))
                    
                    df_extra_pheno = extra_points_pheno.filter((pl.col("SNP_Position") >= x_min) &
                                                               (pl.col("SNP_Position") <= x_max) &
                                                               (pl.col("P-Value_log10") >= y_min) &
                                                               (pl.col("P-Value_log10") <= y_max))
                else:
                    
                    if "D'" in extra_points.columns:
                        df_download = tab.filter((pl.col("SNP_Position") >= x_min) &
                                                (pl.col("SNP_Position") <= x_max) &
                                                (pl.col(linkage) >= y_min) &
                                                (pl.col(linkage) <= y_max))
                        
                        df_extra = extra_points.filter((pl.col("SNP_Position") >= x_min) &
                                                    (pl.col("SNP_Position") <= x_max) &
                                                    (pl.col(linkage) >= y_min) &
                                                    (pl.col(linkage) <= y_max))
                    
                    else:
                        df_download = tab.filter((pl.col("SNP_Position") >= x_min) &
                                                (pl.col("SNP_Position") <= x_max) &
                                                (pl.col("P-Value_log10") >= y_min) &
                                                (pl.col("P-Value_log10") <= y_max))
                        
                        df_extra = extra_points.filter((pl.col("SNP_Position") >= x_min) &
                                                    (pl.col("SNP_Position") <= x_max) &
                                                    (pl.col("P-Value_log10") >= y_min) &
                                                    (pl.col("P-Value_log10") <= y_max))
                    
            else:  # Only x-axis range
                df_download = tab.filter((pl.col("SNP_Position") >= x_min) &
                                         (pl.col("SNP_Position") <= x_max))
                
                df_extra = extra_points.filter((pl.col("SNP_Position") >= x_min) &
                                               (pl.col("SNP_Position") <= x_max))
                
                if toggle_pheno:
                    df_extra_pheno = extra_points_pheno.filter((pl.col("SNP_Position") >= x_min) &
                                                               (pl.col("SNP_Position") <= x_max))
                    
                    df_download = align_and_reorder_columns([df_download, df_extra, df_extra_pheno], 
                                                             reference_df = df_download)
                else:
                    df_download = align_and_reorder_columns([df_download, df_extra], 
                                                             reference_df = df_download)
        else:  # Only y-axis range
            y_min, y_max = relayoutData["yaxis.range[0]"], relayoutData["yaxis.range[1]"]
            
            if toggle_pheno:
                df_download = tab.filter((pl.col("P-Value_log10") >= y_min) &
                                         (pl.col("P-Value_log10") <= y_max))
            
                df_extra = extra_points.filter((pl.col("P-Value_log10") >= y_min) &
                                               (pl.col("P-Value_log10") <= y_max))
                
                df_extra_pheno = extra_points_pheno.filter((pl.col("P-Value_log10") >= y_min) &
                                                           (pl.col("P-Value_log10") <= y_max))
                
                df_download = align_and_reorder_columns([df_download, df_extra, df_extra_pheno], 
                                                         reference_df = df_download)
            else:
                
                if "D'" in extra_points.columns:
                    df_download = tab.filter((pl.col(linkage) >= y_min) &
                                                (pl.col(linkage) <= y_max))
                
                    df_extra = extra_points.filter((pl.col(linkage) >= y_min) &
                                                    (pl.col(linkage) <= y_max))
                else:
                    df_download = tab.filter((pl.col("P-Value_log10") >= y_min) &
                                             (pl.col("P-Value_log10") <= y_max))
                
                    df_extra = extra_points.filter((pl.col("P-Value_log10") >= y_min) &
                                                   (pl.col("P-Value_log10") <= y_max))
                
                df_download = align_and_reorder_columns([df_download, df_extra], 
                                                         reference_df = df_download)
            #raise PreventUpdate

        return dcc.send_data_frame(df_download.to_pandas().to_csv, 
                                   "gwas_svatalog_snp_data.csv")


if __name__ == '__main__':
    app.run_server(debug = True,
                   host = "127.0.0.1",
                   port = 4321)
