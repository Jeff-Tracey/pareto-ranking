source("paretoRankFunctions.R")

criteria <- matrix(c(0.1, 0.1, 0.1,
                     0.2, 0.2, 0.3,
                     0.3, 0.2, 0.2,
                     0.3, 0.2, 0.2,
                     0.3, 0.2, 0.4,
                     0.3, 0.3, 0.4),
                   nrow=6, ncol=3, byrow=T)
pareto.ranks <- cppFasterParetoRank(criteria)
print(cbind(criteria, pareto.ranks))

# Result:
# 0.1 0.1 0.1            1
# 0.2 0.2 0.3            2
# 0.3 0.2 0.2            2
# 0.3 0.2 0.2            2
# 0.3 0.2 0.4            3
# 0.3 0.3 0.4            4

