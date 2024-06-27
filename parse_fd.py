import pandas as pd

# Load the CSV file into a DataFrame
df = pd.read_csv('function_data.csv')

# Compute the Tcycle/Count for each function
df['Tcycle_per_Count'] = df['Tcycle'] / df['Count']

# Select the relevant columns
result_df = df[['Function Name', 'Tcycle', 'Count', 'Tcycle_per_Count']]

# Save the resulting DataFrame to a new CSV file
result_df.to_csv('processed_function_data.csv', index=False)
