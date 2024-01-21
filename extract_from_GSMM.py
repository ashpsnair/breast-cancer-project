import cobra
# Load your GSMM (replace 'your_model.xml' with the actual file path)
model = cobra.io.read_sbml_model("C:/Users/ash24/Downloads/iMM904.xml")

# Search for the biomass reaction by checking reaction names and IDs
biomass_reaction = None
for reaction in model.reactions:
    if 'biomass' in reaction.id.lower() or 'biomass' in reaction.name.lower():
        biomass_reaction = reaction
        break

# Print the biomass reaction information
if biomass_reaction:
    print(f"Biomass Reaction ID: {biomass_reaction.id}")
    print(f"Biomass Reaction Name: {biomass_reaction.name}")
    print(f"Biomass Reaction Equation: {biomass_reaction.reaction}")
else:
    print("Biomass reaction not found in the model.")


# Set the objective function (for example, biomass production)
model.objective = model.reactions.get_by_id('BIOMASS_SC5_notrace')

# Optimize the model
solution = model.optimize()
# Print the results
print(solution)

# Get the flux distribution for all reactions
flux_distribution = solution.fluxes
print("Flux Distribution:")
print(flux_distribution)

############# Stoichiometric Matrix ####################
import pandas as pd
# Create an empty matrix with rows representing metabolites and columns representing reactions
stoichiometric_matrix = []

# Create an empty matrix with rows representing metabolites and columns representing reactions
stoichiometric_matrix = pd.DataFrame(index=[metabolite.id for metabolite in model.metabolites],
                                     columns=[reaction.id for reaction in model.reactions],
                                     dtype=float)

# Populate the matrix with stoichiometric coefficients
for reaction in model.reactions:
    for metabolite, coefficient in reaction.metabolites.items():
        stoichiometric_matrix.at[metabolite.id, reaction.id] = coefficient

# Print the stoichiometric matrix with reaction and metabolite names
print(stoichiometric_matrix)

