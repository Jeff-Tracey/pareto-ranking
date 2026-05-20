#include <Rcpp.h>
#include <vector>
#include <algorithm> // for iter_swap
#include <stack>
#include <queue>
#include <cstdlib>

using namespace Rcpp;

const int NO_RANK = -1;
const int NO_INDEX = -1;

/******************************************************************************
 * Testing functions
 *****************************************************************************/
void printArray(Rcpp::NumericMatrix &array, std::vector<long> &index, 
                std::vector<long> &rank) {
  long nr = array.nrow();
  long nc = array.ncol();
	for(long i = 0; i < nr; i++) {
    Rcpp::Rcout << index[i] << "[" << rank[index[i]] << "]:\t(";
    for (long j = 0; j < nc; j++) {
       Rcpp::Rcout << array(index[i], j);
       if (j < (nc - 1)) {
         Rcpp::Rcout << ", ";
       }
    }
    Rcpp::Rcout << ")" << std::endl;
	}
}

void printIndex(std::vector<long> &index, long first, long last) {
  long n = index.size();
  for (long i = 0; i < n; i++) {
    if ((i >= first) && (i <= last)) {
      Rcpp::Rcout << index[i];
    } else {
      Rcpp::Rcout << "*";
    }
    if (i == (n - 1)) {
      Rcpp::Rcout << std::endl;
    } else {
      Rcpp::Rcout << "\t";
    }
  }
}

/******************************************************************************
 * Functions to compare vectors of criteria for sorting
 *****************************************************************************/

/*
 * Comparison function to test whether b weakly dominates a
 * If b dominates a, return true
 * Otherwise, return false
 * "Weak Pareto dominace"
 *   Condition 1: b <= a for all criteria
 *   Condition 2: b < a for at least one criterion
 */
inline bool isDominatedBy(Rcpp::NumericMatrix &array, long a, long b) {
  long ncol = array.ncol();
  // initial condition implies the elements in the two arrays are equal
  bool allLessEq = true;  // condition 1 is true
  bool oneBetter = false; // condition 2 is false
  for (int i = 0; i < ncol; i++) {
    if (array(b, i) > array(a, i)) { // condition 1 false
      allLessEq = false;
      break; // b does not dominate a, no point in going further
    }
    if (array(b, i) < array(a, i)) { // condition 2 true
      oneBetter = true;
    }
  }
  if (allLessEq && oneBetter) { // both conditions true
    return true;
  } else {
    return false;
  }
}

inline bool dominates(Rcpp::NumericMatrix &array, long a, long b) {
  bool oneBetter = false; // condition 2 is false
  for (int i = 0; i < array.ncol(); i++) {
    if (array(a, i) > array(b, i)) return false; // condition 1 false
    if (array(a, i) < array(b, i)) oneBetter = true; // condition 2 true
  }
  return oneBetter;
}

// mostly for testing with R
RcppExport SEXP isDominated(SEXP obj, SEXP I, SEXP J) {
  Rcpp::NumericMatrix objectives(obj);
  long i = Rcpp::as<long>(I);
  i--; // shift to C++ index
  long j = Rcpp::as<long>(J);
  j--; // shift to C++ index
  bool res = isDominatedBy(objectives, i, j);
  return Rcpp::wrap(res);
}

RcppExport SEXP dominates(SEXP obj, SEXP I, SEXP J) {
  Rcpp::NumericMatrix objectives(obj);
  long i = Rcpp::as<long>(I);
  i--; // shift to C++ index
  long j = Rcpp::as<long>(J);
  j--; // shift to C++ index
  bool res = dominates(objectives, i, j);
  return Rcpp::wrap(res);
}

/******************************************************************************
 * Quicksort functions
 * Note: quicksort doesn't work correctly for Pareto ranking
 *****************************************************************************/

long partition(Rcpp::NumericMatrix &array, std::vector<long> &index, 
              long first, long last) {
  long pivot = first + (last - first)/2;
  std::iter_swap(index.begin()+pivot, index.begin()+last);
  //printIndex(index, first, last); // FOR TESTING
  long storeIndex = first;
  for (int i = first; i < last; i++) { // from first to (last - 1)
    if (isDominatedBy(array, index[last], index[i])) {
      std::iter_swap(index.begin()+i, index.begin()+storeIndex);
      storeIndex++;
      //printIndex(index, first, last); // FOR TESTING
    }
  }
  std::iter_swap(index.begin()+storeIndex, index.begin()+last);
  //printIndex(index, first, last); // FOR TESTING
  return storeIndex;
}

void qsort(Rcpp::NumericMatrix &array, std::vector<long> &index, 
          long first, long last) {
  if (first < last) {
    long p = partition(array, index, first, last);
    //Rcpp::Rcout << "partion = " << p << std::endl; // FOR TESTING
    qsort(array, index, first, p-1);
    qsort(array, index, p+1, last);
  }
}

std::vector<long> quicksortPareto(Rcpp::NumericMatrix &array) {
  int n = array.nrow();
  std::vector<long> index(n);
  // initialize in current order
  for (long i = 0; i < n; i++) {
    index[i] = i;
  }
  // start quicksort
  qsort(array, index, 0, n-1);
  //printIndex(index, 0, n-1); // FOR TESTING
  return index;
}

/******************************************************************************
 * The final algorithms
 *****************************************************************************/
/*
 * NOTE: AS AN OPTION, WE COULD HAVE A MAXIMUM RANK OR A MINIMUM NUMBER (OR
 * AREA) OF UNITS RANKED, AFTER WHICH WE STOP.
 */

RcppExport SEXP paretoRank(SEXP obj) {
    Rcpp::NumericMatrix objectives(obj);
    int nSol = objectives.nrow();
    int nObj = objectives.ncol();
    std::vector<int> ranks(nSol);
    std::vector<int> dup(nSol);
    int remaining = nSol;
    bool isDominated = false;
    int currRank = 1;
    int d1 = 0;
    Rcpp::Rcout << "pareto_rank: Scanning for duplicate solutions..." << 
      std::endl;
    // set all ranks to invalid value
    for (int i = 0; i < nSol; i++) {
        ranks[i] = NO_RANK;
        dup[i] = NO_INDEX;
    }
    // scan for duplicates (-1 implies reference solution)
    for (int i = 0; i < nSol; i++) { // reference
        if (dup[i] == NO_INDEX) {
            for (int j = 0; j < nSol; j++) {
                if ((j != i) && (dup[j] == NO_INDEX)) {
                    d1 = 0;
                    for (int k = 0; k < nObj; k++) {
                        if (objectives(i, k) == objectives(j, k)) {
                            d1++;
                        }
                    }
                    if (d1 == nObj) {
                        dup[j] = i;
                        remaining--; // don't have to rank this one
                        Rcpp::Rcout << "solution " << j <<
                            " matches solution " << i << std::endl;
                    }
                }
            }
        }
    }
    // rank solutions
    Rcpp::Rcout << "pareto_rank: Ranking solutions..." << 
      std::endl;
    while (remaining > 0) {
        Rcpp::Rcout << "Current Rank: " << currRank << " (" <<
            remaining << " remaining)" << std::endl;
        for (int i = 0; i < nSol; i++) {
            if ((ranks[i] == NO_RANK) && (dup[i] == NO_INDEX)) {
                // --- test to see if solution i is dominated
                isDominated = false;
                for (int j = 0; j < nSol; j++) {
                    if ((j != i) && (dup[j] == NO_INDEX) &&
                            ((ranks[j] == NO_RANK) || (ranks[j] == currRank))) {
                        d1 = 0;
                        for (int k = 0; k < nObj; k++) {
                            if (objectives(i, k) < objectives(j, k)) {
                                d1++;
                            }
                        }
                        if (d1 == 0) {
                            isDominated = true;
                            break; // has to be dominated by only one other solution
                        }
                    }
                } // --- end test
                // update if appropriate
                if (!isDominated) {
                    ranks[i] = currRank;
                    remaining--;
                }
            }
        }
        currRank++;
    }
    // assign ranks to duplicates
    Rcpp::Rcout << "pareto_rank: Assigning ranks to duplicate solutions..." << 
      std::endl;
    for (int i = 1; i < nSol; i++) {
        if (dup[i] != NO_INDEX) {
            ranks[i] = ranks[dup[i]];
        }
    }
    return Rcpp::wrap(ranks);
}

RcppExport SEXP fastParetoRank(SEXP obj) {
    Rcpp::NumericMatrix objectives(obj);
    long nSol = objectives.nrow();
    long nObj = objectives.ncol();
    
    // 0. Compare elements
    std::vector<long> ndom(nSol);
    std::vector< std::stack<long> > domstack(nSol);
    std::queue<long> ranktree;
    for (long i = 0; i < nSol; i++) {
      if (i%2000 == 0) {
	Rcpp::Rcout << "Comparing element " << i << " (of " << nSol << ")" << std::endl;
      }
      //
      ndom[i] = 0;
      for (long j = 0; j < nSol; j++) {
	if (dominates(objectives, i, j)) {
	  domstack[i].push(j);
	} else if (dominates(objectives, j, i)) {
	  ndom[i]++;
	} // note: both will be false if i == j
      }
      if (ndom[i] == 0) ranktree.push(i);
    }
    // 1. Assign Ranks
    std::vector<long> rank(nSol);
    long currentRank = 1;
    while (!ranktree.empty()) {
      Rcpp::Rcout << "Assigning rank " << currentRank;
      Rcpp::Rcout << " to " << ranktree.size() << " elements..." << std::endl;
      std::queue<long> nexttree;
      while (!ranktree.empty()) {
	long currentIndex = ranktree.front(); // get ref to first element
	ranktree.pop(); // remove first element
	rank[currentIndex] = currentRank;
	
	while(!domstack[currentIndex].empty()) {
	  long childIndex = domstack[currentIndex].top();
	  domstack[currentIndex].pop();
	  ndom[childIndex]--;
	  if (ndom[childIndex] == 0) {
	    nexttree.push(childIndex);
	  }
	}
      }
      ranktree = nexttree;
      currentRank++;
    }
    return Rcpp::wrap(rank);
}

RcppExport SEXP fasterParetoRank(SEXP obj) {
    Rcpp::NumericMatrix objectives(obj);
    long nSol = objectives.nrow();
    long nObj = objectives.ncol();
    std::vector<long> rank(nSol);
    std::vector<long> remaining(nSol);
    for (long i = 0; i < nSol; i++) {
      remaining[i] = i;
    }
    long insertPos = 0, rankPos = 0, currentRank = 1;
    bool oneBetter = false, noneWorse = true, iDom = false;
    while (insertPos < remaining.size()) {
      Rcpp::Rcout << "Current Rank: " << currentRank << " (" <<
	(nSol-rankPos) << " remaining)" << std::endl;
      for (long i = rankPos; i < nSol; i++) {
	iDom = false;
	for (long j = rankPos; j < nSol; j++) {
	  // does j weakly dominate i?
	  // note: the following is faster than calling dominates()
	  oneBetter = false;
	  noneWorse = true;
	  for (int k = 0; k < nObj; k++) {
	    if (objectives(remaining[j], k) > objectives(remaining[i], k)) {
	      noneWorse = false;
	      break;
	    }
	    if (objectives(remaining[j], k) < objectives(remaining[i], k)) oneBetter = true;
	  }
	  if (oneBetter && noneWorse) {
	    iDom = true;
	    break; // no need to check others
	  }
	  //if (dominates(objectives, remaining[j], remaining[i])) {iDom = true; break;}
	}
	if (!iDom) {
	  rank[remaining[i]] = currentRank;
	  std::iter_swap(remaining.begin()+insertPos, remaining.begin()+i);
	  insertPos++;
	}
      }
      currentRank++;
      rankPos = insertPos;
    }
    return Rcpp::wrap(rank);
}


