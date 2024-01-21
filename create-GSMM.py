from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import gzip
import cobra

# Define your input files
fasta_file = "C:/Users/ash24/Downloads/sce.fasta.gz"
gff_file = "C:/Users/ash24/Downloads/genomic_annotation.gff"
output_file = "C:/Users/ash24/Downloads/annotated_output.fasta"

# Create a dictionary to store gene coordinates
gene_coordinates = {}

# Parse GFF file and store gene coordinates
with open(gff_file, "r") as gff:
    for line in gff:
        if not line.startswith("#"):
            fields = line.strip().split("\t")
            feature_type = fields[2]
            if feature_type == "gene":
                gene_name = fields[8].split(";")[0].split("=")[1]
                start = int(fields[3]) - 1  # Convert to 0-based indexing
                end = int(fields[4])
                gene_coordinates[gene_name] = (start, end)

# Update FASTA file with gene names
records = []
with gzip.open(fasta_file, "rt") as fasta_handle:
    for record in SeqIO.parse(fasta_handle, "fasta"):
        for gene_name, (start, end) in gene_coordinates.items():
            if record.id == gene_name:
                # Add gene name to the description
                record.description = f"gene={gene_name} {record.description}"

                # Add a SeqFeature to store the gene location
                feature_location = FeatureLocation(start, end)
                gene_feature = SeqFeature(location=feature_location, type="gene", id=gene_name)
                record.features.append(gene_feature)

        records.append(record)

# Write the annotated records to a new FASTA file
with open(output_file, "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")

################################### Developing GSMM ##################################################
fasta_file="C:/Users/ash24/Downloads/annotated_output.fasta"

from Bio import SeqIO
def extract_gene_names(fasta_file):
    gene_ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract the 'ID' directly from the record
        gene_id = record.name
        gene_ids.append(gene_id)

    return gene_ids

#list of genes in the annotated fasta file
gene_list = extract_gene_names(fasta_file)

# extract reactions in gene from gene list

##################### Extract information from KEGG XML ############
import xml.etree.ElementTree as ET
xml_file = "C:/Users/ash24/Downloads/sce-kegg.xml"
gene_names = []
# Parse the XML file
tree = ET.parse(xml_file)
root = tree.getroot()

reaction_info_list = []

# Iterate through reactions in the KGML file
for reaction in root.findall(".//reaction"):
    reaction_id = reaction.get("id")
    name = reaction.get("name")

    substrates = [substrate.get("name") for substrate in reaction.findall(".//substrate")]
    products = [product.get("name") for product in reaction.findall(".//product")]

    # Store reaction information in a dictionary
    reaction_info = {
        "reaction_id": reaction_id,
        "name": name,
        "substrates": substrates,
        "products": products
    }

    reaction_info_list.append(reaction_info)










#########################################################################################################
############# Check the model inforamtion ###############
from cobra.io import read_sbml_model

# Replace 'your_model.xml' with the actual path to your SBML file
sbml_file_path = "C:/Users/ash24/Downloads/new_output_gsmm.xml"

# Read the SBML model
model = read_sbml_model(sbml_file_path)

# Get the number of genes, reactions, and metabolites
num_genes = len(model.genes)
num_reactions = len(model.reactions)
num_metabolites = len(model.metabolites)

# Print the results
print(f"Number of Genes: {num_genes}")
print(f"Number of Reactions: {num_reactions}")
print(f"Number of Metabolites: {num_metabolites}")
















######### Stoichiometric Matrix #####################
# Load your GSMM (replace 'your_model.xml' with the actual file path)
model = cobra.io.read_sbml_model("C:/Users/ash24/Downloads/iMM904.xml")

# Print reaction IDs and names
for reaction in model.reactions:
    print(f"Reaction ID: {reaction.id}, Name: {reaction.name}")

# Set the objective function (for example, biomass production)
model.objective = model.reactions.get_by_id('THRA')


# Create an empty matrix with rows representing metabolites and columns representing reactions
stoichiometric_matrix = []

# Create a list to keep track of metabolite indices
metabolite_indices = {}

# Iterate through metabolites and reactions to populate the matrix
for i, metabolite in enumerate(model.metabolites):
    metabolite_indices[metabolite.id] = i  # Keep track of metabolite index
    row = [0] * len(model.reactions)  # Initialize a row with zeros
    for reaction in metabolite.reactions:
        coefficient = model.reactions.get_by_id(reaction.id).metabolites[metabolite]
        row[metabolite_indices[metabolite.id]] = coefficient
    stoichiometric_matrix.append(row)

# Convert the list of lists to a numpy array for easier manipulation (optional)
import numpy as np
stoichiometric_matrix = np.array(stoichiometric_matrix)

# Now, 'stoichiometric_matrix' represents the stoichiometry of each metabolite in each reaction
print(stoichiometric_matrix)






