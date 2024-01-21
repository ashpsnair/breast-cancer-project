import cobra

# Load your GSMM (replace 'your_model.xml' with the actual file path)
model = cobra.io.read_sbml_model("C:/Users/ash24/Downloads/iMM904.xml")

# Set the objective function to maximize biomass production
model.objective = 'BIOMASS_SC5_notrace'

# Perform Flux Balance Analysis (FBA)
solution = model.optimize()

# Get the biomass production rate (objective value)
biomass_production_rate = solution.objective_value

# Print the result
print("Biomass Production Rate:", biomass_production_rate)

############ Sensitivity analysis ############

import pandas as pd

# Get the baseline biomass production rate
baseline_biomass_production_rate = solution.objective_value

# Create an empty DataFrame to store sensitivity results
sensitivity_df = pd.DataFrame(columns=['Metabolite', 'Sensitivity'])

# Create an empty DataFrame to store sensitivity results
sensitivity_df = pd.DataFrame(columns=['Metabolite', 'Sensitivity'])

# Perform sensitivity analysis for each metabolite
for metabolite in model.metabolites:
    # Perturb the flux of each reaction associated with the metabolite
    for reaction in metabolite.reactions:
        original_bounds = reaction.bounds
        reaction.bounds = (0, 0)  # Knock out the reaction

        # Perform FBA after perturbing the flux
        perturbed_solution = model.optimize()
        perturbed_biomass_production_rate = perturbed_solution.objective_value

        # Calculate sensitivity as the change in biomass production rate
        sensitivity = perturbed_biomass_production_rate - baseline_biomass_production_rate

        # Store the sensitivity result in the DataFrame
        sensitivity_df = pd.concat([sensitivity_df, pd.DataFrame([[metabolite.id, sensitivity]], columns=['Metabolite', 'Sensitivity'])])

        # Reset the reaction bounds to the original values
        reaction.bounds = original_bounds

# Sort the DataFrame by sensitivity in descending order
sensitivity_df = sensitivity_df.sort_values(by='Sensitivity', ascending=False)

# Print the top contributors to biomass production
print("Top Contributors to Biomass Production:")
print(sensitivity_df.head())