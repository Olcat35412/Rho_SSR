import pandas as pd
import scipy.stats as stats

#put your inputfile here
df = pd.read_csv('', sep='\t')
genomic_data = df[df['type'] == 'CDS']
print(df)
def chi_square_test(row):
    observed_counts = row[['P1', 'P2', 'P3', 'P4', 'P5', 'P6']]
    total_ssr = observed_counts.sum()
    expected_counts = [total_ssr / len(observed_counts)] * len(observed_counts)
    
    chi2, p = stats.chisquare(observed_counts, f_exp=expected_counts)
    return pd.Series({'Chi-squared': chi2, 'P-value': p})
results = genomic_data.apply(chi_square_test, axis=1)
genomic_data = genomic_data.join(results)
genomic_data['Species_name'] = genomic_data['Species_name']
#put your outputfile here
output_file = ''
genomic_data[['Species_name', 'Chi-squared', 'P-value']].to_csv(output_file, sep='\t', index=False)
print(f"Results saved to {output_file}")
print(genomic_data[['Species_name', 'Chi-squared', 'P-value']])
