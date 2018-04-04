from sys import argv
import pandas
import numpy

script, output_ped, output_map = argv

# baseline = pandas.read_table(input_baseline, header=None, sep='\t')

eradicated = range(0,122)
persistent = range(0,35)


series_list = []
for i in persistent:
	for j in eradicated:
		absence_pers = 34-i
		absence_erad = 121-j

		col1 = ["T"]*i + ["G"]*absence_pers + ["T"]*j + ["G"]*absence_erad
		col2 = ["T"]*i + ["G"]*absence_pers + ["T"]*j + ["G"]*absence_erad

		series_col1 = pandas.Series(col1, name=str(i)+'.'+str(j)+'.1')
		series_col2 = pandas.Series(col2, name=str(i)+'.'+str(j)+'.2')

		series_list.append(series_col1)
		series_list.append(series_col2)

output_df = pandas.concat(series_list, axis=1)

isolate_list = range(1,156)
zeroes = ["0"]*155
ones = ["1"]*155
phenotypes = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

output_df.insert(0, column="phenotype", value=phenotypes)
output_df.insert(0, column="sex", value=zeroes)
output_df.insert(0, column="motherid", value=zeroes)
output_df.insert(0, column="fatherid", value=zeroes)
output_df.insert(0, column="iid", value=isolate_list)
output_df.insert(0, column="familyid", value=isolate_list)

# Write output_df to file

output_df.to_csv(path_or_buf=output_ped, sep='\t', header=False, index=False)


map_df = pandas.concat([pandas.Series(["1"]*4270), pandas.Series(range(0,4270)), pandas.Series(["0"]*4270), pandas.Series(range(0,4270))], axis=1)

# Write map_df to file
map_df.to_csv(path_or_buf=output_map, sep='\t', header=False, index=False)