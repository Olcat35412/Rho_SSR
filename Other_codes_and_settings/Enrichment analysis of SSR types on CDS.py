import pandas as pd
import scipy.stats as stats

# put your input file here
df = pd.read_csv('', sep='\t')
species_names = df['Species_name'].unique()
ssr_types = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6']
genomic_data = df[df['type'] == 'genomic'].set_index(['Species_name'])
cds_data = df[df['type'] == 'CDS'].set_index(['Species_name'])
results = []
for species in species_names:
    for ssr_type in ssr_types:
        try:
            genomic_count = genomic_data.loc[species, ssr_type]
            cds_count = cds_data.loc[species, ssr_type]
            total_genomic = genomic_data.loc[species, ssr_types].sum()
            total_cds = cds_data.loc[species, ssr_types].sum()
            contingency_table = [[cds_count, total_cds - cds_count],
                                 [genomic_count, total_genomic - genomic_count]]
            print(contingency_table)
            odds_ratio, p_value = stats.fisher_exact(contingency_table)
            results.append([species, ssr_type, odds_ratio, p_value])

        except KeyError:
            print(f"Error occurred for {species} and SSR type {ssr_type}. Skipping...")

results_df = pd.DataFrame(results, columns=['Species_name', 'SSR_type', 'Odds_Ratio', 'P_value'])

#change your output_file here
output_file = ''
results_df.to_csv(output_file, sep='\t', index=False)
print(results_df)
