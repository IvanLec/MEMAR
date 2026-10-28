import random as rd
import pandas as pd
def readfile(FILE,dic_sites):
	header = FILE.columns
	for l in FILE.itertuples(index=False):  
		name = l[3]
		pos = l[1]
		if name not in dic_sites:
			dic_sites[name] = {}
		if pos not in dic_sites[name]:
			dic_sites[name][pos] = []
		dic_sites[name][pos].append(l)
	return dic_sites, header


def split_TVT(dic_sites):
	l_train = 0
	l_test = 0
	l_validate = 0
	gene = [g for g in list(dic_sites.keys())]
	rd.shuffle(gene)
	set_train = [gene[0]]
	l_train += len(dic_sites[gene[0]])
	set_test = []
	set_validate = []

	for g in gene[1:]:
		if l_test < 0.15 * (l_train + l_test + l_validate):
			set_test.append(g)
			l_test += len(dic_sites[g])
		elif l_validate < 0.15 * (l_train + l_test + l_validate):
			set_validate.append(g)
			l_validate += len(dic_sites[g])
		else:
			set_train.append(g)
			l_train += len(dic_sites[g])

	s = l_train+l_test+l_validate
	return set_train,set_validate, set_test 

def write_set(sett, dic_sites, header):
	df = pd.DataFrame(columns=header)
	for g in sett:
		for p in dic_sites[g]:
			if len(dic_sites[g][p]) == 1:
				line = dic_sites[g][p][0]
				df.loc[len(df)] = [line[0], line[1], line[2], line[3], float(line[4]), int(line[5]), int(line[6]), float(line[7]), float(line[8]), float(line[9]), float(line[10])]
			else:
				line = dic_sites[g][p][0]
				chromosome, sites, strand, gene, cr1, agcov1, acov1, genecov1, ratio1, pvalue1, p_adjust1 = line[0], line[1], line[2], line[3], float(line[4]), int(line[5]), int(line[6]), float(line[7]), float(line[8]), float(line[9]), float(line[10])
				line = dic_sites[g][p][1]
				chromosome, sites, strand, gene, cr2, agcov2, acov2, genecov2, ratio2, pvalue2, p_adjust2 = line[0], line[1], line[2], line[3], float(line[4]), int(line[5]), int(line[6]), float(line[7]), float(line[8]), float(line[9]), float(line[10])
				r1 = agcov1
				r2 = agcov2
				s = r1 + r2
				cr, agcov, acov, genecov, ratio, pvalue, p_adjust = str((r1*cr1 + r2*cr2)/s), str(agcov1 + agcov2), str(acov1 + acov2), str(genecov1 + genecov2), str((r1*ratio1 + r2*ratio2)/s), str((r1*pvalue1 + r2*pvalue2)/s), str((r1*p_adjust1 + r2*p_adjust2)/s)
				df.loc[len(df)] = [chromosome, sites, strand,gene,cr,agcov,acov,genecov,ratio,pvalue,p_adjust]
	return df



