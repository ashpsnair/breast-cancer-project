"""
**1.How many genes, metabolites, and reactions does this genome-scale metabolic model have?**
"""
import cobra

model = cobra.io.read_sbml_model("C:/Users/ash24/Downloads/iMM904.xml")
# Get the number of genes, metabolites and reactions
num_genes = len(model.genes)
num_metabolites = len(model.metabolites)
num_reactions = len(model.reactions)

# Print the results
print(f"Number of genes: {num_genes}")
print(f"Number of metabolites: {num_metabolites}")
print(f"Number of reactions: {num_reactions}")

"""
**2.What is the definition of biomass in this genome-scale metabolic model?**
"""
# Check for the existence of an objective function
if model.objective:
    # Get the objective function (biomass equation)
    biomass_equation = model.objective.expression

    # Print the biomass equation
    print("Biomass Equation:")
    print(biomass_equation)
else:
    print("No objective function found in the model.")

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


#answer:
# Maximize: 1.0*added_biomass_sink - 1.0*added_biomass_sink_reverse_63abe
# forward flux of the "added_biomass_sink" reaction is not considered in the optimization
# term 2 has -1 coefficient which means the reverse flux of the "added_biomass_sink" reaction (indicated by "added_biomass_sink_reverse_63abe") is included in the objective.
#The optimization aims to minimize or reduce the reverse flux of this reaction.
# "biomass sink" often refers to a reaction that represents the overall production of biomass or cellular growth.
# This type of reaction is typically included in the model to simulate the process of cellular growth and reproduction.

# Check for biomass-related reaction IDs
biomass_related_ids = [reaction.id for reaction in model.reactions if "biomass" in reaction.id.lower()]
print("Biomass-Related Reaction IDs:", biomass_related_ids)

# Identify the biomass sink reaction (replace 'biomass_reaction' with the actual reaction ID)
biomass_reaction = model.reactions.get_by_id('added_biomass_sink')
# Get the stoichiometry of the biomass sink reaction
biomass_components = biomass_reaction.metabolites
# Formulate the biomass sink equation
biomass_equation = ' + '.join([f"{coefficient} {metabolite.id}" for metabolite, coefficient in biomass_components.items()])

# Print the biomass equation
print("Biomass Sink Equation:")
print(biomass_equation)

for metabolite, coefficient in biomass_components.items():
    print(f"{coefficient} {metabolite.id}")


"""
**3.How many exchange metabolites does this genome-scale metabolic model have?**
"""
import pandas as pd
# Get the stoichiometric matrix
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

# Find exchange metabolites
num_exchange_metabolites = 0
exchange_metabolites = []

for reaction in model.reactions:
    # Check if the reaction is an exchange reaction
    if len(reaction.reactants) == 1 and len(reaction.products) == 0:
        num_exchange_metabolites += 1
        exchange_metabolites.append(reaction.reactants[0].id)
    elif len(reaction.reactants) == 0 and len(reaction.products) == 1:
        num_exchange_metabolites += 1
        exchange_metabolites.append(reaction.products[0].id)

# Print the number of exchange metabolites
print("No. of Exchange Metabolites:", num_exchange_metabolites)

# Print the list of exchange metabolite IDs
print("Exchange Metabolite IDs:", exchange_metabolites)

"""
** 4. Which exchange metabolites are important for biomass production? Can you rank them? **
"""
#set the objective
biomass_reaction_id = "BIOMASS_SC5_notrace"
biomass_reaction = model.reactions.get_by_id(biomass_reaction_id)
model.objective = model.reactions.get_by_id('BIOMASS_SC5_notrace')

# Perform sensitivity analysis
exchange_metabolite_sensitivities = {}

for exchange_reaction in model.exchanges:
    with model:
        # Knock out the exchange reaction
        exchange_reaction.knock_out()

        # Optimize the biomass reaction
        solution = model.optimize()

        # Record the change in biomass production
        biomass_change = solution.objective_value - model.optimize().objective_value

        # Store the sensitivity for the exchange metabolite
        exchange_metabolite_sensitivities[exchange_reaction.id] = biomass_change

# Rank the exchange metabolites based on sensitivity
ranked_exchange_metabolites = sorted(exchange_metabolite_sensitivities.items(), key=lambda x: x[1], reverse=True)

# Print the ranked exchange metabolites
print("Ranked Exchange Metabolites based on Sensitivity:")
for exchange_id, sensitivity in ranked_exchange_metabolites:
    print(f"{exchange_id}: {sensitivity}")

#large absolute flux value suggests that the uptake of metabolite is essential for biomass production.
#Both negative and positive flux values are important.
# Negative flux values indicate uptake of a metabolite, while positive flux values indicate secretion.

"""
** 6. Can you describe what a possible optimal culture medium consists of?  **
"""
# Optimize the model to find the optimal nutrient composition
solution = model.optimize()

# Display the optimized nutrient composition
print("Optimal Nutrient Composition:")
for reaction in model.exchanges:
    flux = solution.fluxes.get(reaction.id, 0.0)
    print(f"{reaction.id}: {flux}")

# Alternatively, you can print only the nutrients that are taken up (negative fluxes)
print("\nNutrients Uptake in Optimal Composition:")
for reaction in model.exchanges:
    flux = solution.fluxes.get(reaction.id, 0.0)
    if flux < 0:
        print(f"{reaction.id}: {flux}")


#EX_glc__D_e: -10.0
#The negative flux value indicates glucose uptake.
# The model predicts that the microorganism should take up 10 units of glucose per unit of time to maximize biomass production.

"""
** 6. Suggest a possible strategy to improve this genome-scale metabolic model.   **
"""

