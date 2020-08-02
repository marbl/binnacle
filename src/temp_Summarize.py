import pandas as pd
from os import listdir

filepath = '/Users/harihara/Research-Activities/Data/Binnacle-Op/Scaffold_Coverage_After_Delinking/'
files = listdir(filepath)

for f in files:
	if f[0] != '.' and 'Coverages' in f:
		print(f)
		df = pd.read_csv(filepath+f)
		df_mean = df.groupby(['Connected_Component']).mean()[['Coverage']]
		df_mean = df_mean.rename(columns = {'Coverage':'Mean'})
		df_std = df.groupby(['Connected_Component']).std()[['Coverage']]
		df_std = df_std.rename(columns = {'Coverage':'Dev'})
		df_cov = df_mean.join(df_std)

		print(df_cov.head())

		df_cov.to_csv(filepath+f.replace("_Coverage.csv", "")+'_Summary.csv')

