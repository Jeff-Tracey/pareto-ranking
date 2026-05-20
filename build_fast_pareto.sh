export PKG_CXXFLAGS=`Rscript -e "Rcpp:::CxxFlags()"`
export PKG_LIBS=`Rscript -e "Rcpp:::LdFlags()"`

R CMD SHLIB fast_pareto_rank.cpp
