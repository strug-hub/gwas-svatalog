import dash_bootstrap_components as dbc
import dash_daq as daq
import heapq
import numpy as np
import glob
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from dash import Dash, dcc, html, dash_table, Input, Output, State, ctx
from gtfparse import read_gtf # gtfparse==1.3.0


app = Dash(__name__, 
           external_stylesheets = [dbc.themes.LUX, dbc.icons.BOOTSTRAP],
           title = "GWAS SVatalog")

app._favicon = "assets/favicon.ico"


#### READ IN REQUIRED DATA ####

## GWAS Catalog data v1.0 - downloaded on January 14, 2023 and edited by Dr. Zhuozhi Wang.
    # Expanded SNP haplotypes to match the multiple chromosome locations
    # Used dbSNP151 to correspond the correct rsID to the chromosome location when multiple rsIDs were provided
DF_GWAS_FULL = pd.read_csv("data/gwas_catalog_v1.0-associations_e108.tsv",
                           sep = "\t",
                           header = 0,
                           dtype = {"PUBMEDID" : "object",
                                    "LINK" : "object",
                                    "P-VALUE" : "float64"})

df_gwas = DF_GWAS_FULL[["CHR_ID",
                        "CHR_POS",
                        "SNPS",
                        "STRONGEST SNP-RISK ALLELE",
                        "RISK ALLELE FREQUENCY",
                        "DISEASE/TRAIT",
                        "STUDY",
                        "PUBMEDID",
                        "LINK",
                        "P-VALUE"]]


# Remove the rows with empty chromosome location
df_gwas = df_gwas.dropna(subset = "CHR_ID")

# Use the first rsID only
df_gwas["SNPS"] = df_gwas["SNPS"].str.split(";").str[0]

# Change dtype of chromosome position - prevent decimals. Could not do before as there were ; in rows
df_gwas["CHR_POS"] = df_gwas["CHR_POS"].astype(int)

df_gwas.columns = ["Chromosome",
                   "SNP_Position",
                   "SNP_Name_GWAS",
                   "Risk_Allele",
                   "Risk_Allele_Frequency",
                   "Phenotype",
                   "Study",
                   "Pubmed_ID",
                   "Link",
                   "P-Value"]


df_gwas["Chromosome"] = "chr" + df_gwas["Chromosome"].astype(str)

df_gwas["Risk_Allele"] = df_gwas["Risk_Allele"].str.split("-").str[1]

df_gwas = df_gwas[df_gwas["P-Value"] <= 1]


## SV Annotation data - created by Dr. Zhuozhi Wang on March 13, 2023
DF_ANNO_FULL = pd.read_csv("data/sv_annotations.tsv",
                           sep = "\t",
                           header = 0,
                           dtype = {"startAnn" : "Int64",
                                    "endAnn" : "Int64"})

df_anno = DF_ANNO_FULL[["chrAnn",
                        "startAnn",
                        "endAnn",
                        "ID",
                        "variantTypeAnn",
                        "sizeAnn"]]

df_anno.columns =["Chromosome",
                  "SV_Start",
                  "SV_End",
                  "SV_Name",
                  "SV_Type",
                  "SV_Length"]


## SV-SNP LD data created by Dr. Zhuozhi Wang on April 6, 2023
LD_file_paths = glob.glob(os.path.join("data/LD",
                                       "*_ld_stats.txt"))
LD_tables = [pd.read_table(file_path,
                           header = 0,
                           names = ["SNP_Name", 
                                    "SNP_Position",
                                    "SV_Name",
                                    "SV_Position",
                                    "R2",
                                    "D'"])
             for file_path in LD_file_paths]

DF_LD = pd.concat(LD_tables)


## SV/SNP allele annotations done by Thomas Nalpathamkalam on April 10, 2023
allele_file_paths = glob.glob(os.path.join("data/LD",
                                            "*_allele_freq.txt"))
allele_tables = [pd.read_table(file_path,
                               header = 0,
                               names = ["Chromosome",
                                        "SNP_Position",
                                        "SV_SNP_Name",
                                        "Reference_Allele",
                                        "Alternate_Allele",
                                        "Sample_AF",
                                        "gnomAD_nfe_AF",
                                        "dbSNP"])
                 for file_path in allele_file_paths]

DF_ALLELE = pd.concat(allele_tables)


## Gene/exon data from MANE.GRCh38.v1.0.ensembl_genomic.gtf extracted January 5, 2023

DF_GTF = read_gtf("data/MANE.GRCh38.v1.0.ensembl_genomic.gtf")
DF_GTF = pd.DataFrame(data = DF_GTF)

df_gene = DF_GTF[(DF_GTF["feature"] == "transcript") | \
                 (DF_GTF["feature"] == "exon")]

df_gene = df_gene[df_gene["tag"].str.contains("MANE_Select")]
df_gene = df_gene[["seqname",
                    "start", 
                    "end", 
                    "strand", 
                    "feature", 
                    "gene_name"]]

df_gene.columns = ["chr",
                   "start",
                   "end",
                   "strand",
                   "feature",
                   "gene"]


#### DATA WRANGLING ####
## Subset allele data for only SV information
df_sv_allele = DF_ALLELE[DF_ALLELE["SV_SNP_Name"].str.startswith("P")]

df_sv_allele = df_sv_allele.rename(columns = {"SV_SNP_Name" : "SV_Name"})

df_sv_allele = df_sv_allele.drop(["Reference_Allele",
                                  "SNP_Position",
                                  "gnomAD_nfe_AF",
                                  "dbSNP"],
                                 axis = 1)

df_sv_allele = df_sv_allele.drop_duplicates()


## Subset allele data for only SNP information
df_snp_allele = DF_ALLELE[~DF_ALLELE["SV_SNP_Name"].str.startswith("P")]

df_snp_allele = df_snp_allele.rename(columns = {"SV_SNP_Name" : "SNP_Name"})

df_snp_allele = df_snp_allele.drop_duplicates()



## Subset DF_ANNO_FULL with columns that can be displayed to the public. Add SV AF to this table as well.
df_sv_anno = DF_ANNO_FULL

df_sv_anno = df_sv_anno.rename(columns = {"chrAnn" : "Chromosome",
                                          "startAnn" : "SV_Start",
                                          "endAnn" : "SV_End",
                                          "variantTypeAnn" : "SV_Type",
                                          "sizeAnn" : "SV_Length"})

df_sv_anno = df_sv_anno.merge(df_sv_allele.iloc[:,[1,2,3]],
                              left_on = "ID",
                              right_on = "SV_Name")

df_sv_anno = df_sv_anno.drop(columns = "Alternate_Allele")

columns_range = list(range(1,33))
columns_range.insert(0, 33)
columns_range.insert(6, 34)

df_sv_anno = df_sv_anno.iloc[:,columns_range]



## Creation of SV table for easier access of data required - join the dataframes based on SV names
df_sv_join = DF_LD.merge(df_anno,
                              on = "SV_Name",
                              how = "left")

df_sv_join = df_sv_join[["Chromosome",
                         "SV_Name",
                         "SV_Start",
                         "SV_End",
                         "SNP_Name",
                         "SNP_Position",
                         "R2",
                         "D'"]]


## Merging SV allele information with the LD dataset
df_sv_join = df_sv_join.merge(df_sv_allele,
                              on = ["Chromosome",
                                    "SV_Name"],
                              how = "left")

df_sv_join = df_sv_join.rename(columns = {"Alternate_Allele" : "SV_Type",
                                          "Sample_AF" : "SV_AF"})

df_sv_join.loc[df_sv_join["SV_Type"] == "<DEL>", "SV_Type"] = "Deletion"
df_sv_join.loc[df_sv_join["SV_Type"] == "<INS>", "SV_Type"] = "Insertion"
df_sv_join.loc[df_sv_join["SV_Type"] == "<DUP>", "SV_Type"] = "Duplication"
df_sv_join.loc[df_sv_join["SV_Type"] == "<INV>", "SV_Type"] = "Inversion"

df_sv_join = df_sv_join.iloc[:, [0,1,2,3,8,9,4,5,6,7]]


## Merging SNP allele information with the LD dataset
df_sv_snp_join = df_sv_join.merge(df_snp_allele,
                                  on = ["Chromosome",
                                        "SNP_Name",
                                        "SNP_Position"],
                                        how = "left")

df_sv_snp_join = df_sv_snp_join.drop("SNP_Name", axis = 1)

df_sv_snp_join = df_sv_snp_join.iloc[:, [0,1,2,3,4,5,13,6,9,10,11,12,7,8]]

df_sv_snp_join = df_sv_snp_join.rename(columns = {"dbSNP" : "SNP_Name_dbSNP"})



## Creation of a full dataframe to make data extraction much easier
df_full_join = df_sv_snp_join.merge(df_gwas,
                                    how = "left",
                                    on = ["Chromosome",
                                          "SNP_Position"])

# Remove rows where SV AF is NA - these are SVs with AF < 0.1
df_full_join = df_full_join[df_full_join["SV_AF"].notna()]

# Change p-values to log10 scale
    # Smallest pval was 1-e323. After looking at the distribution, all pvalues under 1e-50 were changed to 1e-50.
df_full_join.loc[df_full_join["P-Value"] < 1e-50, "P-Value"] = 1e-50

df_full_join["P-Value_log10"] = np.log10(df_full_join["P-Value"])
df_full_join["P-Value_log10"] = -df_full_join["P-Value_log10"]


# Subset df_anno for the SVs in df_full_foin
sv_names = df_full_join["SV_Name"].unique()
df_anno_subset = df_anno[df_anno["SV_Name"].isin(sv_names)]


## Creation of phenotype dictionary to encompass SV for each gene - computationally faster when creating SV table
dict_pheno = {}

df_dict = df_full_join[["SV_Name", "SV_Start", "SV_End", "Phenotype"]]
df_dict = df_dict.to_dict("records")

for row in df_dict:
    if row["Phenotype"] not in dict_pheno:
        dict_pheno[row["Phenotype"]] = []
        dict_pheno[row["Phenotype"]].append(row["SV_Name"])
    else:
        dict_pheno[row["Phenotype"]].append(row["SV_Name"])

for k, v in dict_pheno.items():
    dict_pheno[k] = set(v)


## Creation of y-coordinates to prevent overlap in genes when plotted
def assign_y_coordinates(df = df_gene,
                         y_start = -0.3,
                         y_increment = 0.5,
                         col_name = "y_coord"):
    # Priority queue to keep track of the end positions of genes for each y-coordinate
    df = df[df["feature"] == "transcript"]
    df = df.sort_values(by = ["chr", "start"])
    df = df.reset_index(drop = True) #reset index 
     # Dictionary to store priority queues for each chromosome
    chromosome_heaps = {}

    # List to store the y-coordinate for each gene annotation
    y_coordinates = [y_start] * len(df)

    for index, row in df.iterrows():
        chrom = row['chr']

        # If the chromosome is not in the dictionary, initialize its heap
        if chrom not in chromosome_heaps:
            chromosome_heaps[chrom] = []

        # Work with the heap for the current chromosome
        heap = chromosome_heaps[chrom]

        # If there's an available y-coordinate whose gene ends before the current gene starts,
        # reuse that y-coordinate and update its gene end position
        if heap and heap[0][0] <= row['start']:
            previous_gene_end, y_coordinate = heapq.heappop(heap)
            heapq.heappush(heap, (row['end'], y_coordinate))
        else:
            # If no y-coordinate is available, create a new one
            if heap:
                y_coordinate = min(y_coordinate for _, y_coordinate in heap) - y_increment
            else:
                y_coordinate = y_start
            heapq.heappush(heap, (row['end'], y_coordinate))

        # Assign the y-coordinate to the gene annotation
        y_coordinates[index] = y_coordinate

    # Add the y-coordinates to the dataframe
    df[col_name] = y_coordinates

    return df

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


# Make dropdown box for phenotypes
# PHENO_LINE = html.Div(id = 'pheno-line-div')

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

phenotypes = sorted(list(df_full_join["Phenotype"].unique()))
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
        sv_list_chr = list(df_anno_subset["SV_Name"])
    else:
        chrom = str("chr" + chrom)
        sv_list_chr = df_anno_subset[(df_anno_subset["Chromosome"] == chrom)]
        sv_list_chr = list(sv_list_chr["SV_Name"])


    # Filter by range
    sv_list_range = df_anno_subset[(df_anno_subset["SV_Start"] >= range_start) & \
                                   (df_anno_subset["SV_End"] <= range_end)]

    sv_list_range = list(sv_list_range["SV_Name"])


    # Filter by phenotype
    if isinstance(pheno, list):
        sv_list_pheno = list(df_anno_subset["SV_Name"])
    else:
        sv_list_pheno = dict_pheno[pheno]


    # Filter by gene
    if isinstance(gene, list):
        sv_list_gene = list(df_anno_subset["SV_Name"])
    else:
        df_1gene = df_gene[(df_gene["feature"] == "transcript") & \
                           (df_gene["gene"] == gene)]
        start = df_1gene["start"].min()
        end = df_1gene["end"].max()
        chrom = df_1gene.iat[0,0]
        strand = df_1gene.iat[0,3]

        sv_list_gene = df_sv_anno[(df_sv_anno["Chromosome"] == chrom) & \
                                  (df_sv_anno["SV_Start"] >= (start - 100000)) & \
                                  (df_sv_anno["SV_End"] <= (end + 100000))]

        sv_list_gene = list(sv_list_gene["SV_Name"])


    # Find overlaps sv's from all lists
    sv_list_complete = set(sv_list_pheno).intersection(sv_list_chr, sv_list_gene, sv_list_range)

    tab_show = df_sv_anno[df_sv_anno["SV_Name"].isin(sv_list_complete)]
    tab_show = tab_show[["SV_Name", "Chromosome", "SV_Start", "SV_End", "SV_Type", "SV_Length", "Sample_AF"]]
    tab_show = tab_show.rename(columns = {"SV_Name" : "ID",
                                          "Chromosome" : "Chrom",
                                          "SV_Start" : "Start",
                                          "SV_End" : "End",
                                          "SV_Type" : "Type",
                                          "SV_Length" : "Size (bp)",
                                          "Sample_AF" : "AF"})

    return dash_table.DataTable(id ='strvar-table',
                                data = tab_show.to_dict("records"),
                                columns = [{'id': c, 'name': c} for c in tab_show.loc[:,[ "Chrom", "Start", "End", "Type", "Size (bp)", "AF"]]],
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
        df_anno_subset_sv = pd.DataFrame()
        return dash_table.DataTable(id ='anno-table',
                                    data = df_anno_subset_sv.to_dict("records"),
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
        df_anno_subset_sv = df_sv_anno[df_sv_anno["SV_Name"] == sv]
        df_anno_subset_sv = df_anno_subset_sv.transpose()
        df_anno_subset_sv = df_anno_subset_sv.rename_axis(" ").reset_index()
        df_anno_subset_sv.columns = ["Header", "Information"]
        return dash_table.DataTable(id ='anno-table',
                                    data = df_anno_subset_sv.to_dict("records"),
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
        tab = pd.DataFrame(columns = ["SNP Position",
                                      "-log10(P-Value)"])
        sc_plot = px.scatter(data_frame = tab,
                             x = "SNP Position",
                             y = "-log10(P-Value)",
                             labels = {"-log10(P-Value)" : "-log<sub>10</sub> (P-Value)",
                                       "SNP Position" : "Position in hg38 (bp)"},
                             width = 850,
                             height = 650,
                            #  marginal_x = "rug",
                             template = "ggplot2")

        sc_plot.update_layout(margin = dict(l = 20, r = 20, t = 50, b = 50),
                              autosize = False,
                              font_family = "Nunito Sans",
                              font_size = 20)
    else:
        tab = df_full_join[df_full_join["SV_Name"] == sv]

        # Display +/- 1KB from SV position
        sv_start = df_anno_subset[df_anno_subset["SV_Name"] == sv]["SV_Start"].values[0]
        sv_end = df_anno_subset[df_anno_subset["SV_Name"] == sv]["SV_End"].values[0]

        sv_start_plot = sv_start - 1000
        sv_end_plot = sv_end + 1000


        if pheno == "None":
            tab_lab = tab.to_dict("records")

            def tab_label_maker(row):

                label_pheno = row["Phenotype"]
                label_snp_name = row["SNP_Name_dbSNP"]
                label_chr = row["Chromosome"]
                label_snp_pos = row["SNP_Position"]
                label_pval = row["P-Value"]
                label_r2 = row["R2"]
                label_d = row["D'"]
                label_svname = row["SV_Name"]
                label_svstart = row["SV_Start"]
                label_svend = row["SV_End"]
                label_svtype = row["SV_Type"]

                label = [label_pheno, label_snp_name, label_chr, label_snp_pos, label_pval, label_r2, label_d, label_svname, label_svstart, label_svend, label_svtype]

                return label

            sc_plot = px.scatter(data_frame = tab,
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
                                  hovertemplate = '<br><b>Phenotype:</b> %{customdata[0]}<br>' + '<br>SNP Name: %{customdata[1]}' + '<br>Chromosome: %{customdata[2]}' + '<br>SNP Position: %{customdata[3]}' + '<br>P-Value: %{customdata[4]}' + '<br>R2: %{customdata[5]}' + "<br>D': %{customdata[6]}",
                                  customdata = [tab_label_maker(row) for row in tab_lab])

            sc_plot.update_xaxes(showspikes = True,
                                 spikemode = "across")

            sc_plot.update_layout(margin = dict(l = 20, r = 20, t = 50, b = 50),
                                  autosize = False,
                                  font_family = "Nunito Sans",
                                  font_size = 20)
                                #   xaxis_range = [sv_start_plot, sv_end_plot],
                                #   yaxis_range = [-0.6, 1])

             # Extract min and max SNP positions in chosen SV
            min_SNP = min(tab["SNP_Position"]) - 100000
            max_SNP = max(tab["SNP_Position"]) + 100000

            if sv_start < min_SNP:
                min_SNP = sv_start - 100000

            elif sv_start > max_SNP:
                max_SNP = sv_end  + 100000

            else:
                min_SNP = min(tab["SNP_Position"]) - 100000
                max_SNP = max(tab["SNP_Position"]) + 100000


            chromo = list(set(tab["Chromosome"]))[0]
            svname = list(set(tab["SV_Name"]))[0]


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
                                         hovertemplate = '<br><b>SV Name:</b> %{customdata[7]}<br>' + '<br>Chromosome: %{customdata[2]}' + '<br>Start: %{customdata[8]}' + '<br>End: %{customdata[9]}' + '<br>Type: %{customdata[10]}',
                                         customdata = [tab_label_maker(tab_lab[0])]))

            # Add extra points from df_sv with R2 and D' values
            extra_points = df_sv_snp_join[df_sv_snp_join["SV_Name"] == sv]

            extra_points = extra_points[~extra_points["SNP_Name_dbSNP"].isin(tab["SNP_Name_dbSNP"])]

             # Create a list of labels for hovertext
            df_lab = extra_points.to_dict("records")

            def extra_label_maker(row):

                label_snp_name = row["SNP_Name_dbSNP"]
                label_chr = row["Chromosome"]
                label_snp_pos = row["SNP_Position"]
                label_r2 = row["R2"]
                label_d = row["D'"]

                label = [label_snp_name, label_chr, label_snp_pos, label_r2, label_d]

                return label

            sc_plot.add_trace(go.Scattergl(x = extra_points["SNP_Position"],
                                           y = extra_points[linkage_value],
                                           mode = "markers",
                                           marker = dict(size = 10,
                                                         color = "rgba(168,168,168,0.5)",
                                                         line = dict(color = "rgba(168,168,168,1)")),
                                           hovertemplate = '<br><b>SNP Name: %{customdata[0]}</b>' +'<br>Chromosome: %{customdata[1]}' + '<br>SNP Position: %{customdata[2]}' + '<br>R2: %{customdata[3]}' + "<br>D': %{customdata[4]}",
                                           customdata = [extra_label_maker(row) for row in df_lab],
                                           showlegend = False,
                                           name = "Not in \nGWAS Catalog"))

            # Use min and max SNP positions to extract required rows from df_exons to obtain gene information
            exons = df_gene[(df_gene["chr"] == chromo) & \
                            (df_gene["start"] >= min_SNP) & \
                            (df_gene["end"] <= max_SNP)]

            exons_genes = exons.groupby(["chr", "gene", "strand"])["start"].min()
            exons_genes = exons_genes.to_frame().reset_index()

            exons_end = exons.groupby(["chr", "gene", "strand"])["end"].max()
            exons_end = exons_end.to_frame().reset_index()

            exons_genes["end"] = exons_end["end"]

            # gene_line_position = -0.3

            for i in range(0,len(exons_genes)):
                row = exons_genes.iloc[i]
                gene_line_position = float(df_transcript["y_coord_nopheno"][df_transcript["gene"] == row["gene"]].iloc[0]) 

                sc_plot.add_shape(type = "line",
                                x0 = (row["start"]),
                                x1 = (row["end"]),
                                y0 = gene_line_position,
                                y1 = gene_line_position,
                                line_width = 3,
                                line_color = "darkslategrey")

                text_anno_x = round((row["end"] - row["start"]) / 2) + row["start"]

                if row["strand"] == "+":
                    sc_plot.add_annotation(x = text_anno_x,
                                           y = gene_line_position + 0.060,
                                           text = row["gene"] + '<span style="font-size:15px">\u2192</span>',
                                           font = dict(size = 10),
                                           hovertext = row["gene"] + '<br><b>Gene direction:</b> forward',
                                           name = "Gene",
                                           showarrow = False)

                else:
                    sc_plot.add_annotation(x = text_anno_x,
                                           y = gene_line_position + 0.060,
                                           text = row["gene"] + '<span style="font-size:15px">\u2190</span>',
                                           font = dict(size = 10),
                                           hovertext = row["gene"] + '<br><b>Gene direction:</b> reverse',
                                           name = "Gene",
                                           showarrow = False)


            exons_exons = exons[exons["feature"] == "exon"]

            for i in range(0,len(exons_exons)):
                row = exons_exons.iloc[i]
                gene_line_position = float(df_transcript["y_coord_nopheno"][df_transcript["gene"] == row["gene"]].iloc[0])
                sc_plot.add_shape(type = "rect",
                                x0 = (row["start"]),
                                x1 = (row["end"]),
                                y0 = gene_line_position + 0.020,
                                y1 = gene_line_position - 0.020,
                                line = dict(color = "mediumseagreen",
                                            width = 2),
                                fillcolor = "mediumseagreen")


        else:
            tab = tab[tab["Phenotype"] == pheno]

            tab_lab = tab.to_dict("records")

            def tab_label_maker(row):

                label_pheno = row["Phenotype"]
                label_snp_name = row["SNP_Name_dbSNP"]
                label_chr = row["Chromosome"]
                label_snp_pos = row["SNP_Position"]
                label_pval = row["P-Value"]
                label_r2 = row["R2"]
                label_d = row["D'"]
                label_svname = row["SV_Name"]
                label_svstart = row["SV_Start"]
                label_svend = row["SV_End"]
                label_svtype = row["SV_Type"]

                label = [label_pheno, label_snp_name, label_chr, label_snp_pos, label_pval, label_r2, label_d, label_svname, label_svstart, label_svend, label_svtype]

                return label

            sc_plot = px.scatter(data_frame = tab,
                                 x = "SNP_Position",
                                 y = "P-Value_log10",
                                 labels = {"SNP_Position" : "Position in hg38 (bp)",
                                           "P-Value_log10" : "-log<sub>10</sub> (P-Value)",
                                           linkage_value : legend_title},
                                 color = linkage_value,
                                #  color_continuous_scale = px.colors.sequential.Viridis_r,
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
                                #  size = size_value,
                                 width = 850,
                                 height = 650,
                                #  marginal_x = "rug",
                                 template = "ggplot2")

            sc_plot.update_traces(hovertemplate = '<br><b>Phenotype:</b> %{customdata[0]}<br>' + '<br>SNP Name: %{customdata[1]}' + '<br>Chromosome: %{customdata[2]}' + '<br>SNP Position: %{customdata[3]}' + '<br>P-Value: %{customdata[4]}' + '<br>r2: %{customdata[5]}' + "<br>D': %{customdata[6]}",
                                  customdata = [tab_label_maker(row) for row in tab_lab],
                                  marker = dict(size = 10))

            sc_plot.update_xaxes(showspikes = True,
                                 spikemode = "across")

            sc_plot.update_layout(margin = dict(l = 20, r = 20, t = 50, b = 50),
                                 autosize = True,
                                 font_family = "Nunito Sans",
                                 font_size = 20,
                                #  xaxis_range = [sv_start_plot, sv_end_plot],
                                #  yaxis_range = [-20, 50],
                                 coloraxis_colorbar = dict(lenmode = "pixels",
                                                           len = 150,
                                                           yanchor = "top",
                                                           y = 1,
                                                           dtick = 0.2))

             # Extract min and max SNP positions in chosen SV
            min_SNP = min(tab["SNP_Position"]) - 100000
            max_SNP = max(tab["SNP_Position"]) + 100000

            if sv_start < min_SNP:
                min_SNP = sv_start - 100000

            elif sv_start > max_SNP:
                max_SNP = sv_end  + 100000

            else:
                min_SNP = min(tab["SNP_Position"]) - 100000
                max_SNP = max(tab["SNP_Position"]) + 100000


            chromo = list(set(tab["Chromosome"]))[0]
            svname = list(set(tab["SV_Name"]))[0]

            # Add points from gwas for pheno, without r2 and D' values
            extra_points = df_gwas[(df_gwas["Chromosome"] == chromo) & \
                                   (df_gwas["SNP_Position"] >= min_SNP) & \
                                   (df_gwas["SNP_Position"] <= max_SNP) & \
                                   (df_gwas["Phenotype"] == pheno)]

            extra_points = extra_points[~extra_points["SNP_Name_GWAS"].isin(tab["SNP_Name_GWAS"])]

            extra_points.loc[extra_points["P-Value"] < 1e-50, "P-Value"] = 1e-50

            extra_points["P-Value_log10"] = np.log10(extra_points["P-Value"])
            extra_points["P-Value_log10"] = -extra_points["P-Value_log10"]

            # Create a list of labels for hovertext
            df_lab = extra_points.to_dict("records")

            def label_maker(row):

                label_pheno = row["Phenotype"]
                label_snp_name = row["SNP_Name_GWAS"]
                label_chr = row["Chromosome"]
                label_snp_pos = row["SNP_Position"]
                label_pval = row["P-Value"]

                label = [label_pheno, label_snp_name, label_chr, label_snp_pos, label_pval]

                return label


            sc_plot.add_trace(go.Scattergl(x = extra_points["SNP_Position"],
                                           y = extra_points["P-Value_log10"],
                                           mode = "markers",
                                           marker = dict(size = 10,
                                                         color = "rgba(168,168,168,0.5)",
                                                         line = dict(color = "rgba(168,168,168,1)")),
                                           hovertemplate = '<br><b>Phenotype:</b> %{customdata[0]}<br>' + '<br>SNP Name: %{customdata[1]}' + '<br>Chromosome: %{customdata[2]}' + '<br>SNP Position: %{customdata[3]}' + '<br>P-Value: %{customdata[4]}',
                                           customdata = [label_maker(row) for row in df_lab],
                                           showlegend = False,
                                           name = "No LD Data"))

            # Add points from other phenos in LD with SV
            if toggle_pheno == True:
                extra_points_pheno = df_full_join[(df_full_join["Chromosome"] == chromo) & \
                                                  (df_full_join["SNP_Position"] >= min_SNP) & \
                                                  (df_full_join["SNP_Position"] <= max_SNP) & \
                                                  (df_full_join["SV_Name"] == sv)]

                extra_points_pheno = extra_points_pheno[~extra_points_pheno["SNP_Name_GWAS"].isin(tab["SNP_Name_GWAS"])]


                # Create a list of labels for hovertext
                df_lab_2 = extra_points_pheno.to_dict("records")

                def label_maker(row):

                    label_pheno = row["Phenotype"]
                    label_snp_name = row["SNP_Name_GWAS"]
                    label_chr = row["Chromosome"]
                    label_snp_pos = row["SNP_Position"]
                    label_pval = row["P-Value"]
                    label_r2 = row["R2"]
                    label_d = row["D'"]

                    label = [label_pheno, label_snp_name, label_chr, label_snp_pos, label_pval, label_r2, label_d]

                    return label


                sc_plot.add_trace(go.Scattergl(x = extra_points_pheno["SNP_Position"],
                                               y = extra_points_pheno["P-Value_log10"],
                                               mode = "markers",
                                               marker = dict(size = 10,
                                                             color = "rgba(168,168,168,0.5)",
                                                             line = dict(color = "rgba(168,168,168,1)")),
                                               hovertemplate = '<br><b>Phenotype:</b> %{customdata[0]}<br>' + '<br>SNP Name: %{customdata[1]}' + '<br>Chromosome: %{customdata[2]}' + '<br>SNP Position: %{customdata[3]}' + '<br>P-Value: %{customdata[4]}' + '<br>r2: %{customdata[5]}' + "<br>D': %{customdata[6]}",
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
                                         hovertemplate = '<br><b>SV Name:</b> %{customdata[7]}<br>' + '<br>Chromosome: %{customdata[2]}' + '<br>Start: %{customdata[8]}' + '<br>End: %{customdata[9]}' + '<br>Type: %{customdata[10]}',
                                         customdata = [tab_label_maker(tab_lab[0])]))

            # Use min and max SNP positions to extract required rows from df_exons to obtain gene information
            exons = df_gene[(df_gene["chr"] == chromo) & \
                            (df_gene["start"] >= min_SNP) & \
                            (df_gene["end"] <= max_SNP)]

            exons_genes = exons.groupby(["chr", "gene", "strand"])["start"].min()
            exons_genes = exons_genes.to_frame().reset_index()

            exons_end = exons.groupby(["chr", "gene", "strand"])["end"].max()
            exons_end = exons_end.to_frame().reset_index()

            exons_genes["end"] = exons_end["end"]

            # gene_line_position = -6.25

            for i in range(0,len(exons_genes)):
                row = exons_genes.iloc[i]
                gene_line_position = float(df_transcript["y_coord_pheno"][df_transcript["gene"] == row["gene"]].iloc[0]) 
                sc_plot.add_shape(type = "line",
                                x0 = (row["start"]),
                                x1 = (row["end"]),
                                y0 = gene_line_position,
                                y1 = gene_line_position,
                                line_width = 3,
                                line_color = "darkslategrey")

                text_anno_x = round((row["end"] - row["start"]) / 2) + row["start"]

                if row["strand"] == "+":
                    # sc_plot.add_trace(go.Scatter(x = np.array([text_anno_x]),
                    #                             y = np.array([gene_line_position - 2.25]),
                    #                             mode = "markers+lines",
                    #                             marker = dict(symbol = "triangle-right-open",
                    #                                         size = 12,
                    #                                         color = "darkslategrey"),
                    #                             showlegend = False,
                    #                             hovertemplate = row["gene"] + '<br><b>Gene direction:</b> forward',
                    #                             name = "gene direction"))

                    sc_plot.add_annotation(x = text_anno_x,
                                           y = gene_line_position + 1.75,
                                           text = row["gene"] + '<span style="font-size:15px">\u2192</span>',
                                           font = dict(size = 10),
                                           hovertext = row["gene"] + '<br><b>Gene direction:</b> forward',
                                           name = "Gene",
                                           showarrow = False)

                else:
                    # sc_plot.add_trace(go.Scatter(x = np.array([text_anno_x]),
                    #                             y = np.array([gene_line_position - 1.25]),
                    #                             mode = "markers+lines",
                    #                             marker = dict(symbol = "triangle-left-open",
                    #                                         size = 12,
                    #                                         color = "darkslategrey"),
                    #                             showlegend = False,
                    #                             hovertemplate = row["gene"] + '<br><b>Gene direction:</b> reverse',
                    #                             name = "gene direction"))

                    sc_plot.add_annotation(x = text_anno_x,
                                           y = gene_line_position + 1.75,
                                           text = row["gene"] + '<span style="font-size:15px">\u2190</span>',
                                           font = dict(size = 10),
                                           hovertext = row["gene"] + '<br><b>Gene direction:</b> reverse',
                                           name = "Gene",
                                           showarrow = False)


            exons_exons = exons[exons["feature"] == "exon"]

            for i in range(0,len(exons_exons)):
                row = exons_exons.iloc[i]
                gene_line_position = float(df_transcript["y_coord_pheno"][df_transcript["gene"] == row["gene"]].iloc[0]) 
                sc_plot.add_shape(type = "rect",
                                x0 = (row["start"]),
                                x1 = (row["end"]),
                                y0 = gene_line_position + 0.75,
                                y1 = gene_line_position - 0.75,
                                line = dict(color = "mediumseagreen",
                                            width = 2),
                                fillcolor = "mediumseagreen")


    return sc_plot


# SCATPLOT = dcc.Graph(id = "scatter-plot",
#                      figure = make_plot())

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
        df_snp = pd.DataFrame()
        return dash_table.DataTable(id ='snp-table',
                                    data = df_snp.to_dict("records"),
                                    columns = [{'id': c, 'name': c} for c in df_snp.columns],
                                    page_size = 400,
                                    fixed_rows = {'headers': True},
                                    style_table = {'height': '500px', 'overflowY': 'auto'},
                                    style_cell = {'textAlign': 'center',
                                                  'font-family': 'Nunito Sans',
                                                  'backgroundColor': 'whitesmoke'},
                                    style_header = {'backgroundColor': 'turquoise',
                                                    'fontWeight': 'bold'})

    else:
        df_snp = df_full_join[(df_full_join["SNP_Name_dbSNP"] == snp) | (df_full_join["SNP_Name_GWAS"] == snp)]

        df_snp["Link"] = "[" + df_snp["Pubmed_ID"] + "]" + "(https://" + df_snp["Link"] + ")"

        df_snp["SNP_Name_dbSNP"] = "[" + df_snp["SNP_Name_dbSNP"] + "]" + "(https://www.ncbi.nlm.nih.gov/snp/" + df_snp["SNP_Name_dbSNP"] + ")"

        df_snp["SNP_Name_GWAS"] = "[" + df_snp["SNP_Name_GWAS"] + "]" + "(https://www.ebi.ac.uk/gwas/variants/" + df_snp["SNP_Name_GWAS"] + ")"

        df_snp = df_snp[["Chromosome",
                         "SNP_Position",
                         "SNP_Name_dbSNP",
                         "SNP_Name_GWAS",
                         "Reference_Allele",
                         "Alternate_Allele",
                         "Risk_Allele",
                         "Risk_Allele_Frequency",
                         "Sample_AF",
                         "gnomAD_nfe_AF",
                         "Phenotype",
                         "P-Value",
                         "Study",
                         "Link"]]
        df_snp = df_snp.rename(columns = {"Link" : "Pubmed Link",
                                          "Chromosome" : "Chrom",
                                          "SNP_Position" : "SNP Position",
                                          "SNP_Name_dbSNP" : "SNP Name: dbSNP",
                                          "SNP_Name_GWAS" : "SNP Name: GWAS",
                                          "Risk_Allele" : "Risk Allele",
                                          "Risk_Allele_Frequency" : "Risk AF",
                                          "Reference_Allele" : "Reference Allele",
                                          "Alternate_Allele" : "Alternate Allele",
                                          "Sample_AF" : "Sample AF",
                                          "gnomAD_nfe_AF" : "gnomAD NFE AF"})

        df_snp = df_snp.drop_duplicates()

        return dash_table.DataTable(id ='snp-table',
                                    data = df_snp.to_dict("records"),
                                    columns = [{'id': c, 'name': c, 'presentation': 'markdown'} if (c == "Pubmed Link" or c == "SNP Name: dbSNP" or c == "SNP Name: GWAS") else ({'id': c, 'name': c, 'type':'numeric', 'format': {'specifier': '.3f'}} if (c == "Risk AF" or c == "Sample AF" or c == "gnomAD NFE AF") else {'id': c, 'name': c}) for c in df_snp.columns],
                                    markdown_options = {"html" : True,
                                                        "link_target": "_blank"},
                                    page_size = 400,
                                    fixed_rows = {'headers': True},
                                    # fixed_columns = {'headers' : True,
                                    #                  'data' : 4},
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
                                                              href = "https://locusfocus.research.sickkids.ca/colocalization",
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

# DOCUMENTATION = dcc.Markdown('''
# Documentation for GWAS SVatalog can be found [here](https://gwas-svatalog-docs.readthedocs.io/en/latest/index.html).
# Colocalization testing of GWAS loci across various datasets can be conducted at [LocusFocus] (https://locusfocus.research.sickkids.ca/).
# ''',
#                              id = 'doc-markdown')

# DISC_DOC_DIV = html.Div(id = 'disc-doc-div',
#                         children = [DOCUMENTATION,
#                                     DISCLAIMER_DIV])

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
        chromosome = str("chr" + str(chrom))
        new_genes = df_gene[df_gene["chr"] == chromosome]
        new_genes = sorted(list(new_genes["gene"]))
        new_pheno = df_full_join[df_full_join["Chromosome"] == chromosome]
        new_pheno = sorted(list(new_pheno["Phenotype"].unique().astype(str)))

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
        sv_selected = tab_show.iloc[selected_rows[0],[0]][0]
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
        sv_name = tab_show.iloc[selected_rows[0],[0]][0]
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
        State('modal-snp-table', 'is_open')
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


    # If button was triggered
    if ctx.triggered[0]['prop_id'] == 'download-button.n_clicks':

        if (relayoutData == None or "xaxis.autorange" in relayoutData):
            # This means the plot is in the default starting position
            df_download = tab
            df_extra = extra_points

            if toggle_pheno == True:
                df_extra_pheno = extra_points_pheno
                df_download = pd.concat([df_download, df_extra, df_extra_pheno])
            elif toggle_pheno == False:
                df_download = pd.concat([df_download, df_extra])
            else:
                df_download = pd.concat([df_download, df_extra])


        elif "xaxis.range[0]" in relayoutData:

            if "yaxis.range[0]" in relayoutData:

                if "P-Value" in extra_points.columns:
                    df_download = tab[(tab["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                      (tab["SNP_Position"] <= relayoutData["xaxis.range[1]"]) & \
                                      (tab["P-Value_log10"] >= relayoutData["yaxis.range[0]"]) & \
                                      (tab["P-Value_log10"] <= relayoutData["yaxis.range[1]"])]

                    df_extra = extra_points[(extra_points["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                            (extra_points["SNP_Position"] <= relayoutData["xaxis.range[1]"]) &
                                            (extra_points["P-Value_log10"] >= relayoutData["yaxis.range[0]"]) & \
                                            (extra_points["P-Value_log10"] <= relayoutData["yaxis.range[1]"])]

                    df_extra_pheno = extra_points_pheno[(extra_points_pheno["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                                        (extra_points_pheno["SNP_Position"] <= relayoutData["xaxis.range[1]"]) &
                                                        (extra_points_pheno["P-Value_log10"] >= relayoutData["yaxis.range[0]"]) & \
                                                        (extra_points_pheno["P-Value_log10"] <= relayoutData["yaxis.range[1]"])]

                else:
                    df_download = tab[(tab["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                      (tab["SNP_Position"] <= relayoutData["xaxis.range[1]"]) & \
                                      (tab[linkage] >= relayoutData["yaxis.range[0]"]) & \
                                      (tab[linkage] <= relayoutData["yaxis.range[1]"])]

                    df_extra = extra_points[(extra_points["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                            (extra_points["SNP_Position"] <= relayoutData["xaxis.range[1]"]) &
                                            (extra_points[linkage] >= relayoutData["yaxis.range[0]"]) & \
                                            (extra_points[linkage] <= relayoutData["yaxis.range[1]"])]

            else:

                if toggle_pheno == True:

                    df_download = tab[(tab["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                      (tab["SNP_Position"] <= relayoutData["xaxis.range[1]"])]

                    df_extra = extra_points[(extra_points["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                            (extra_points["SNP_Position"] <= relayoutData["xaxis.range[1]"])]

                    df_extra_pheno = extra_points_pheno[(extra_points_pheno["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                                        (extra_points_pheno["SNP_Position"] <= relayoutData["xaxis.range[1]"])]

                    df_download = pd.concat([df_download, df_extra, df_extra_pheno])

                elif toggle_pheno == False:

                    df_download = tab[(tab["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                      (tab["SNP_Position"] <= relayoutData["xaxis.range[1]"])]

                    df_extra = extra_points[(extra_points["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                            (extra_points["SNP_Position"] <= relayoutData["xaxis.range[1]"])]

                    df_download = pd.concat([df_download, df_extra])

                else:

                    df_download = tab[(tab["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                      (tab["SNP_Position"] <= relayoutData["xaxis.range[1]"])]

                    df_extra = extra_points[(extra_points["SNP_Position"] >= relayoutData["xaxis.range[0]"]) & \
                                            (extra_points["SNP_Position"] <= relayoutData["xaxis.range[1]"])]

                    df_download = pd.concat([df_download, df_extra])



        else: # yaxis.range[0] only

            if toggle_pheno == True:

                    df_download = tab[(tab["P-Value_log10"] >= relayoutData["yaxis.range[0]"]) & \
                                      (tab["P-Value_log10"] <= relayoutData["yaxis.range[1]"])]

                    df_extra = extra_points[(extra_points["P-Value_log10"] >= relayoutData["yaxis.range[0]"]) & \
                                            (extra_points["P-Value_log10"] <= relayoutData["yaxis.range[1]"])]

                    df_extra_pheno = extra_points_pheno[(extra_points_pheno["P-Value_log10"] >= relayoutData["yaxis.range[0]"]) & \
                                                        (extra_points_pheno["P-Value_log10"] <= relayoutData["yaxis.range[1]"])]

                    df_download = pd.concat([df_download, df_extra, df_extra_pheno])

            elif toggle_pheno == False:

                    df_download = tab[(tab[linkage] >= relayoutData["yaxis.range[0]"]) & \
                                      (tab[linkage] <= relayoutData["yaxis.range[1]"])]

                    df_extra = extra_points[(extra_points[linkage] >= relayoutData["yaxis.range[0]"]) & \
                                            (extra_points[linkage] <= relayoutData["yaxis.range[1]"])]

                    df_download = pd.concat([df_download, df_extra])

            else:

                    df_download = tab[(tab[linkage] >= relayoutData["yaxis.range[0]"]) & \
                                      (tab[linkage] <= relayoutData["yaxis.range[1]"])]

                    df_extra = extra_points[(extra_points[linkage] >= relayoutData["yaxis.range[0]"]) & \
                                            (extra_points[linkage] <= relayoutData["yaxis.range[1]"])]

                    df_download = pd.concat([df_download, df_extra])
            #raise PreventUpdate

        return dcc.send_data_frame(df_download.to_csv, "gwas_svatalog_snp_data.csv")



if __name__ == '__main__':
    app.run_server(debug = True,
                   host = "127.0.0.1",
                   port = 4321)
