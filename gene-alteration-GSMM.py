import cobra

# Load your GSMM (replace 'your_model.xml' with the actual file path)
model = cobra.io.read_sbml_model("C:/Users/ash24/Downloads/iMM904.xml")

############ find gene information and check how many genes are involved in more than one pathway ###########

# Create lists to store gene-related information
gene_names = []
num_reactions = []
associated_reactions_list = []

# Iterate through genes in the model
for gene in model.genes:
    # Get reactions associated with the gene
    associated_reactions = [reaction.id for reaction in gene.reactions]

    # Append information to lists
    gene_names.append(gene.id)
    num_reactions.append(len(associated_reactions))
    associated_reactions_list.append(associated_reactions)

# Create a DataFrame from the lists
gene_reaction_df = pd.DataFrame({
    'Gene Name': gene_names,
    'Number of Reactions': num_reactions,
    'Associated Reactions': associated_reactions_list
})

# Print the DataFrame
gene_reaction_df=gene_reaction_df.sort_values(by='Number of Reactions', ascending=False)
print(gene_reaction_df)

#chossing the gene with maximum nmber of reactions YKR009C

# Specify the gene ID to knockout
gene_to_knockout = 'YKR009C'

# Find the corresponding reactions associated with the gene
reactions_to_knockout = [reaction.id for reaction in model.genes.get_by_id(gene_to_knockout).reactions]

# Knock out the gene by setting the flux of associated reactions to zero
for reaction_id in reactions_to_knockout:
    model.reactions.get_by_id(reaction_id).lower_bound = 0
    model.reactions.get_by_id(reaction_id).upper_bound = 0

# Set the objective function to maximize biomass production
model.objective = 'BIOMASS_SC5_notrace'

# Perform Flux Balance Analysis (FBA)
solution = model.optimize()

# Get the biomass production rate (objective value)
biomass_production_rate = solution.objective_value

# no difference in the biomass production

###### Knocking out multiple genes
# Create a dictionary to store genes and their associated reactions
genes_to_knockout = list(gene_reaction_df[gene_reaction_df['Number of Reactions']>10]['Gene Name'])

# Iterate through the list of genes
for gene_id in genes_to_knockout:
    # Find the corresponding reactions associated with the gene
    reactions_to_knockout = [reaction.id for reaction in model.genes.get_by_id(gene_id).reactions]

    # Create a copy of the DataFrame row and append it to the DataFrame
    new_row = pd.DataFrame({'Gene Name': [gene_id], 'Number of Reactions': [len(reactions_to_knockout)], 'Associated Reactions': [reactions_to_knockout]})
    gene_reaction_df = pd.concat([gene_reaction_df, new_row], ignore_index=True, sort=False)

    # Knock out the gene by setting the flux of associated reactions to zero
    for reaction_id in reactions_to_knockout:
        model.reactions.get_by_id(reaction_id).lower_bound = 0
        model.reactions.get_by_id(reaction_id).upper_bound = 0

# Set the objective function to maximize biomass production
model.objective = 'BIOMASS_SC5_notrace'

# Perform Flux Balance Analysis (FBA)
solution = model.optimize()

# Get the biomass production rate (objective value)
biomass_production_rate = solution.objective_value

#saw the difference in the biomass production

