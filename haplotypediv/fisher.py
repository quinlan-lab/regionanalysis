import scipy.stats as stats
import sys
odds_ratio, pvalue = stats.fisher_exact([[int(sys.argv[1]), int(sys.argv[2])], [int(sys.argv[3]), int(sys.argv[4])]],alternative='greater')
print odds_ratio
print pvalue
