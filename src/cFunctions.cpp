#include <Rcpp.h>
using namespace Rcpp;

//#define NANCHECK(val) if (ISNAN(val)) {stop("Illegal NaN.");} // THIS HAS SIGNIFCANT COMPUTATIONAL BURDEN
#define STARTTIME 0.


//' @title Compute stage frequencies from simulated stage durations.
//' @param sampleTimes Vector of times at which the cohort was sampled.
//' @param sDur A matrix of stage durations where each row corresponds to the stage duration of an individual.
//'             If there is no exit from the last stage, the stage duration of that stage can be set to Inf.
//' @param mortRate A vector of mortality rates for each stage. Should be of length equal to the number of 
//'                 columns of sDur.
//' @return A matrix of stage counts. Each row gives the number of individuals in each stage for a given sampling time.             
//' @export
// [[Rcpp::export]]
NumericMatrix stageFreqMC(const NumericVector sampleTimes, const NumericMatrix sDur, const NumericVector mortRate) {
  int nStage = sDur.ncol();
  int N = sDur.nrow();
  int nSample = sampleTimes.length();
  double state, lSurvProb, prevT;
  NumericMatrix stageFreq(nSample, nStage);
  for (int ind = 0; ind < N; ind++) {
    int stage = 0;
    int samp = 0;
    lSurvProb = 0.;
    prevT = STARTTIME;
    state = STARTTIME + sDur(ind, stage);
    while (samp < nSample) {
      while (stage < (nStage - 1) && state < sampleTimes[samp]) {
        lSurvProb -= mortRate[stage] * (state - prevT);
        prevT = state;
        stage ++;
        state += sDur(ind, stage); 
      }
    lSurvProb -= mortRate[stage]*(sampleTimes[samp] - prevT);
    prevT = sampleTimes[samp];
    stageFreq(samp, stage) += exp(lSurvProb);
    samp++; 
    }
  } 
  return stageFreq;
}
