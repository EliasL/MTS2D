from tabulate import tabulate
import pandas as pd

# Example DataFrame
data = {
    "Nr_Unused_Cores": ["24/64", "16/64", "7/64", "38/64", "128/128", "8/128", "80/128", "8/128", "8/128"]
}
index = ["galois", "pascal", "schwartz", "lagrange", "condorcet", "dalembert", "poincare", "fourier", "descartes"]
df = pd.DataFrame(data, index=index)

# Function to add color based on the value
def colorize(value):
    numerator, denominator = map(int, value.split('/'))
    percentage = numerator / denominator
    # Using ANSI escape codes for coloring: Red for low, Green for high
    if percentage > 0.5:
        return f"\033[92m{value}\033[0m"  # Green
    else:
        return f"\033[91m{value}\033[0m"  # Red

# Applying colorize function to each value in the DataFrame
df['Nr_Unused_Cores'] = df['Nr_Unused_Cores'].apply(colorize)

# Convert DataFrame to a format that tabulate can use
data_for_tabulate = [df.index.tolist()] + df.T.values.tolist()

# Printing the table using tabulate
print(tabulate(data_for_tabulate, headers='firstrow', tablefmt='grid'))
