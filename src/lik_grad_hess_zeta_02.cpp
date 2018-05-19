// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Function to calculate the negative loglikelihood and gradient
// [[Rcpp::export]]
List neglikgrad(const arma::colvec& Omega, // parameters
                const arma::colvec& D, // disease status (binary)
                const arma::mat& G, // genetic variables
                const arma::mat& E, // environmental variables
                const double pi1 // disease rate in the population
) {
  //declarations
  const unsigned int p = Omega.n_elem; // number of parameters
  const unsigned int n = D.n_elem; // number of patients
  const unsigned int n1 = accu(D); // number of cases
  const unsigned int n0 = n-n1; // number of controls
  const unsigned int nG = G.n_cols; // number of genetic variables
  const unsigned int nE = E.n_cols; // number of environmental variables
  const double pi0 = 1 - pi1; // proportion of controls in the population
  const double cst = (pi1*n0)/(pi0*n1); // the constant from the denominator of S
  const double pibyn [2] = {pi0/n0, pi1/n1}; // (pi0/ncontrol, pi1/ncase)
  unsigned int counter; // counter to create G*E interaction terms
  arma::rowvec GEcombo(p); // vector to hold (1, G, E, G*E)
  double m; // exp(dot product of Omega and GEcombo)
  double cm; // cst*m
  double Denom; // denominator of S
  double Denompibyn; // Denom * pibyn(D(j))
  double RAccumScalar = 0; // accumulates n values to get R
  arma::rowvec ROmegaAccumVec(p); // accumulates n values to get R_Omega
  double LogLikS = 0; // loglikelihood contribution from S
  double LogLikR = 0; // loglikelihood contribution from R
  arma::rowvec GradientLHS(p); // sums (S_Omega/S) over n patients
  arma::rowvec GradientRHS(p); // sums (R_Omega/R) over n patients
  Rcpp::NumericVector GradientVec(p); // to convert from a 1xp matrix to a p-vector

  ROmegaAccumVec.zeros();
  GradientLHS.zeros();
  GradientRHS.zeros();
  GEcombo.at(0)=1;


  for (unsigned int i=0; i<n; i++){ // begin i loop  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    // Create a row-vector with the model matrix for subject i
    for (unsigned int k=0; k<nG; k++){GEcombo.at(1+k) = G.at(i,k);}
    for (unsigned int k=0; k<nE; k++){GEcombo.at(1+nG+k) = E.at(i,k);}
    counter=0;
    for (unsigned int k=0; k<nE; k++){
      for (unsigned int l=0; l<nG; l++){
        GEcombo.at(1+nG+nE+counter) = G.at(i,l)*E.at(i,k);
        counter++;
      } // end l loop
    } // end k loop

    m = exp(dot(GEcombo,Omega));
    cm = cst*m;
    Denom = 1/(1+cm);

    // Calculate the left-hand side of the loglikelihood
    if(D.at(i)==0){
      LogLikS +=log(Denom);
      GradientLHS -= GEcombo*cm*Denom;
    } else {
      LogLikS +=log(m*Denom);
      GradientLHS += GEcombo*Denom;
    }

    // Calculate the right-hand side of the loglikelihood
    RAccumScalar = 0;
    ROmegaAccumVec.zeros();

    for (unsigned int j=0; j<n; j++){ // begin j loop  jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
      for (unsigned int k=0; k<nG; k++){GEcombo.at(1+k) = G.at(j,k);}
      // for (unsigned int k=0; k<nE; k++){GEcombo(1+nG+k) = E(i,k);}  // can comment this out because these values are never reset
      counter=0;
      for (unsigned int k=0; k<nE; k++){
        for (unsigned int l=0; l<nG; l++){
          GEcombo.at(1+nG+nE+counter) = G.at(j,l)*E.at(i,k);
          counter++;
        }
      }
      m = exp(dot(GEcombo,Omega));
      cm = cst*m;
      Denom = 1/(1+cm);
      if(D.at(j)==0){
        Denompibyn = Denom*pibyn[0];
      } else {
        Denompibyn = Denom*pibyn[1];
      }

      RAccumScalar += (1+m)*Denompibyn;
      ROmegaAccumVec += GEcombo*(m-cm)*Denom*Denompibyn; // R_hat_Omega
    } // end of the j loop  jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj

    LogLikR += log(RAccumScalar);
    GradientRHS += ROmegaAccumVec/RAccumScalar;
  } // end of the i loop  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  GradientVec = -1 * (GradientLHS-GradientRHS);
  GradientVec.attr("dim") = R_NilValue;
  return List::create(Named("objective") = -1*(LogLikS-LogLikR),
                      Named("gradient") = GradientVec);
}

// Function to calculate the negative hessian
// [[Rcpp::export]]
arma::mat neghess(const arma::colvec& Omega, // parameters
                  const arma::colvec& D, // disease status (binary)
                  const arma::mat& G, // genetic variables
                  const arma::mat& E, // environmental variables
                  const double pi1 // disease rate in the population
) {
  //declarations
  const unsigned int p = Omega.n_elem; // number of parameters
  const unsigned int n = D.n_elem; // number of patients
  const unsigned int n1 = accu(D); // number of cases
  const unsigned int n0 = n-n1; // number of controls
  const unsigned int nG = G.n_cols; // number of genetic variables
  const unsigned int nE = E.n_cols; // number of environmental variables
  const double pi0 = 1 - pi1; // proportion of controls in the population
  const double cst = (pi1*n0)/(pi0*n1); // the constant from the denominator of S
  const double pibyn [2] = {pi0/n0, pi1/n1}; // (pi0/ncontrol, pi1/ncase)
  unsigned int counter; // counter to create G*E interaction terms
  arma::rowvec GEcombo(p); // vector to hold (1, G, E, G*E)
  double m; // exp(dot product of Omega and GEcombo)
  double cm; // cst*m
  double Denom; // denominator of S
  double DenomSquared; // Denom squared
  double Denompibyn; // Denom * pibyn(D(j))
  double dcm; // Denom*(1-cm)
  double RAccumScalar = 0; // accumulates n values to get R
  double RSquared; // R squared
  arma::rowvec ROmegaAccumTemp; // stores values used to calculate R_Omega
  arma::rowvec ROmegaAccumVec(p); // accumulates n values to get R_Omega
  arma::mat ROmegaOmegaAccum(p,p);  // accumulates over n iterations to get the Hessian of R wrt Omega
  arma::mat LeftHessAccum(p,p);  // accumulates over n iterations to get the LHS of the Hessian of the loglikelihood
  arma::mat RightHessAccum(p,p);  // accumulates over n iterations to get the RHS of the Hessian of the loglikelihood
  arma::mat Hessian(p,p);  // hessian to be reported

  ROmegaAccumVec.zeros();
  LeftHessAccum.zeros();
  RightHessAccum.zeros();
  GEcombo.at(0)=1;


  for (unsigned int i=0; i<n; i++){ // begin i loop  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    // Create a row-vector with the model matrix for subject i
    for (unsigned int k=0; k<nG; k++){GEcombo.at(1+k) = G.at(i,k);}
    for (unsigned int k=0; k<nE; k++){GEcombo.at(1+nG+k) = E.at(i,k);}
    counter=0;
    for (unsigned int k=0; k<nE; k++){
      for (unsigned int l=0; l<nG; l++){
        GEcombo.at(1+nG+nE+counter) = G.at(i,l)*E.at(i,k);
        counter++;
      } // end l loop
    } // end k loop

    m = exp(dot(GEcombo,Omega));
    cm = cst*m;
    Denom = 1/(1+cm);
    DenomSquared = Denom*Denom;
    // LeftHessAccum -= GEcombo.t() * GEcombo * cm*Denom*Denom;  // commented out because only the lower triangle is needed
    for (unsigned int k=p; k>0; k--){  // lower triangle of LeftHessAccum
      for (unsigned int l=0; l<k; l++){
        LeftHessAccum.at(k-1,l) -= GEcombo.at(k-1) * GEcombo.at(l) * cm * DenomSquared;
      } // end l loop
    } // end k loop

    // Calculate the right-hand side of the loglikelihood
    RAccumScalar = 0;
    ROmegaAccumVec.zeros();
    ROmegaOmegaAccum.zeros();

    for (unsigned int j=0; j<n; j++){ // begin j loop  jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
      for (unsigned int k=0; k<nG; k++){GEcombo.at(1+k) = G.at(j,k);}
      // for (unsigned int k=0; k<nE; k++){GEcombo(1+nG+k) = E(i,k);}  // can comment this out because these values are never reset
      counter=0;
      for (unsigned int k=0; k<nE; k++){
        for (unsigned int l=0; l<nG; l++){
          GEcombo.at(1+nG+nE+counter) = G.at(j,l)*E.at(i,k);
          counter++;
        }
      }
      m = exp(dot(GEcombo,Omega));
      cm = cst*m;
      Denom = 1/(1+cm);
      if(D.at(j)==0){
        Denompibyn = Denom*pibyn[0];
      } else {
        Denompibyn = Denom*pibyn[1];
      }
      dcm = Denom*(1-cm);

      RAccumScalar += (1+m)*Denompibyn;
      ROmegaAccumTemp = GEcombo*(m-cm)*Denom*Denompibyn;
      for (unsigned int k=0; k<p; k++){ROmegaAccumVec.at(k) += ROmegaAccumTemp.at(k);}

      for (unsigned int k=p; k>0; k--){  // lower triangle of ROmegaOmegaAccum
        for (unsigned int l=0; l<k; l++){
          ROmegaOmegaAccum.at(k-1,l) += GEcombo.at(k-1) * ROmegaAccumTemp.at(l)*dcm;
        } // end l loop
      } // end k loop
    } // end of the j loop  jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj

    RSquared = RAccumScalar*RAccumScalar;

    // RightHessAccum += ROmegaOmegaAccum/RAccumScalar - (ROmegaAccumVec.t()*ROmegaAccumVec)/RSquared;  // commented out because only the lower triangle is needed
    for (unsigned int k=p; k>0; k--){
      for (unsigned int l=0; l<k; l++){
        RightHessAccum.at(k-1,l) += ROmegaOmegaAccum.at(k-1,l)/RAccumScalar-(ROmegaAccumVec.at(k-1)*ROmegaAccumVec.at(l))/RSquared;
      } // end l loop
    } // end k loop

  } // end of the i loop  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  Hessian = -1 * (LeftHessAccum-RightHessAccum);
  for (unsigned int k=0; k<p-1; k++){  // Upper triangle of Hessian
    for (unsigned int l=k+1; l<p; l++){
      Hessian.at(k,l) = Hessian.at(l,k);
    } // end l loop
  } // end k loop
  return Hessian;
}


// Function to calculate the negative Hessian, zeta0, and zeta1
// [[Rcpp::export]]
List hesszeta(const arma::colvec& Omega, // parameters
              const arma::colvec& D, // disease status (binary)
              const arma::mat& G, // genetic variables
              const arma::mat& E, // environmental variables
              const double pi1 // disease rate in the population
) {
  //declarations
  const unsigned int p = Omega.n_elem; // number of parameters
  const unsigned int n = D.n_elem; // number of patients
  const unsigned int n1 = accu(D); // number of cases
  const unsigned int n0 = n-n1; // number of controls
  const unsigned int nG = G.n_cols; // number of genetic variables
  const unsigned int nE = E.n_cols; // number of environmental variables
  const double pi0 = 1 - pi1; // proportion of controls in the population
  const double cst = (pi1*n0)/(pi0*n1); // the constant from the denominator of S
  const double pibyn [2] = {pi0/n0, pi1/n1}; // (pi0/ncontrol, pi1/ncase)
  unsigned int counter; // counter to create G*E interaction terms
  unsigned int countcase = 0; // counter to create zeta1
  unsigned int countcontrol = 0;  // counter to create zeta0
  arma::rowvec GEcombo(p); // vector to hold (1, G, E, G*E)
  // arma::mat GEmat(p,p); // matrix of GEcombo times its transpose
  double m; // exp(dot product of Omega and GEcombo)
  double cm; // cst*m
  double Denom; // denominator of S
  double DenomSquared; // Denom squared
  double dcm;
  double mcmdenomdenom;
  arma::mat SZeroOneMat(n, n);  // matrix of S(d=0, ...) + S(d=1, ...) for all combos of G & E
  arma::rowvec SOmegaByS(p);  // (S_Omega/S)
  arma::cube SOmegaZeroOneCube(p, n, n);  // matrix of S_Omega(d=0, ...) + S_Omega(d=1, ...) for all combos of G & E
  arma::mat RMat(n, n); // Matrix of the components of R for all combinations of G & E
  arma::rowvec RAccum(n); // accumulates the column sums of RMat to get a vector of R values
  arma::rowvec RSquared(n); // R squared
  arma::cube ROmegaCube(p, n, n); // Matrix of the components of R_Omega for all combinations of G & E
  arma::mat ROmegaAccum(p, n); // accumulates the slices of ROmegaCube to get a matrix of R_Omega values
  arma::mat ROmegaOmegaAccum(p,p);  // accumulates over n iterations to get the Hessian of R wrt Omega
  arma::mat LeftHessAccum(p,p);  // accumulates over n iterations to get the LHS of the Hessian of the loglikelihood
  arma::mat RightHessAccum(p,p);  // accumulates over n iterations to get the RHS of the Hessian of the loglikelihood
  arma::mat Hessian(p,p);  // hessian to be reported
  arma::rowvec RightZetaAccum(p);  // accumulates over n iterations to get the RHS of zeta
  arma::mat zeta0(n0, p);  // zeta0 matrix
  arma::mat zeta1(n1, p);  // zeta1 matrix

  RAccum.zeros();
  ROmegaAccum.zeros();
  LeftHessAccum.zeros();
  RightHessAccum.zeros();
  GEcombo.at(0)=1;

  for (unsigned int i=0; i<n; i++){ // begin i loop iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ROmegaOmegaAccum.zeros();
    for (unsigned int j=0; j<n; j++){ // begin j loop jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
      for (unsigned int k=0; k<nG; k++){GEcombo.at(1+k) = G.at(j,k);}
      for (unsigned int k=0; k<nE; k++){GEcombo.at(1+nG+k) = E.at(i,k);}
      counter=0;
      for (unsigned int k=0; k<nE; k++){
        for (unsigned int l=0; l<nG; l++){
          GEcombo.at(1+nG+nE+counter) = G.at(j,l)*E.at(i,k);
          counter++;
        } // end l loop
      } // end k loop
      m = exp(dot(GEcombo,Omega));
      cm = cst*m;
      Denom = 1/(1+cm);
      dcm = Denom*(1-cm);
      mcmdenomdenom = (m-cm)*Denom*Denom;
      SZeroOneMat(i,j) = (1+m) * Denom;
      for (unsigned int k=0; k<p; k++){SOmegaZeroOneCube.at(k,i,j) = GEcombo.at(k)*mcmdenomdenom;}
      // RMat(i,j) = SZeroOneMat(i,j)*pibyn(D(j));
      // for (unsigned int k=0; k<p; k++){ROmegaCube(k,i,j) = SOmegaZeroOneCube(k,i,j) * pibyn(D(j));}
      if(D.at(j)==0){
        RMat.at(i,j) = SZeroOneMat.at(i,j)*pibyn[0];
        for (unsigned int k=0; k<p; k++){ROmegaCube.at(k,i,j) = SOmegaZeroOneCube.at(k,i,j) * pibyn[0];}
      } else {
        RMat.at(i,j) = SZeroOneMat.at(i,j)*pibyn[1];
        for (unsigned int k=0; k<p; k++){ROmegaCube.at(k,i,j) = SOmegaZeroOneCube.at(k,i,j) * pibyn[1];}
      }

      for (unsigned int k=p; k>0; k--){  // lower triangle of ROmegaOmegaAccum
        for (unsigned int l=0; l<k; l++){
          ROmegaOmegaAccum.at(k-1,l) += GEcombo.at(k-1) * ROmegaCube.at(l,i,j)*dcm;
        } // end l loop
      } // end k loop
    } // end j loop jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj

    for (unsigned int j=0; j<n; j++){RAccum.at(i) += RMat.at(i,j);}
    for (unsigned int k=0; k<p; k++){
      for (unsigned int j=0; j<n; j++){
        ROmegaAccum.at(k,i) += ROmegaCube.at(k,i,j);
      } // end j loop
    } // end k loop
    RSquared.at(i) = RAccum.at(i)*RAccum.at(i);

    for (unsigned int k=p; k>0; k--){
      for (unsigned int l=0; l<k; l++){
        RightHessAccum.at(k-1,l) += ROmegaOmegaAccum.at(k-1,l)/RAccum.at(i)-(ROmegaAccum.at(k-1,i)*ROmegaAccum.at(l,i))/RSquared.at(i);
      } // end l loop
    } // end k loop
  } // end i loop iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii


  for (unsigned int i=0; i<n; i++){ // begin i loop
    RightZetaAccum.zeros();
    for (unsigned int k=0; k<p; k++){
      for (unsigned int j=0; j<n; j++){
        RightZetaAccum.at(k) += (SOmegaZeroOneCube.at(k,j,i)/RAccum.at(j)) - (ROmegaAccum.at(k,j)*SZeroOneMat.at(j,i)/RSquared.at(j));
      } // end j loop
    } // end k loop

    for (unsigned int k=0; k<nG; k++){GEcombo.at(1+k) = G.at(i,k);}
    for (unsigned int k=0; k<nE; k++){GEcombo.at(1+nG+k) = E.at(i,k);}
    counter=0;
    for (unsigned int k=0; k<nE; k++){
      for (unsigned int l=0; l<nG; l++){
        GEcombo.at(1+nG+nE+counter) = G.at(i,l)*E.at(i,k);
        counter++;
      } // end l loop
    } // end k loop
    m = exp(dot(GEcombo,Omega));
    cm = cst*m;
    Denom = 1/(1+cm);
    // Calculate S-terms
    if(D[i]==0){
      SOmegaByS = -GEcombo*cm*Denom;
      countcontrol++;
    } else {
      SOmegaByS = GEcombo*Denom;
      countcase++;
    }

    DenomSquared = Denom*Denom;
    // LeftHessAccum -= GEcombo.t() * GEcombo * cm*Denom*Denom;
    for (unsigned int k=p; k>0; k--){  // lower triangle of LeftHessAccum
      for (unsigned int l=0; l<k; l++){
        LeftHessAccum.at(k-1,l) -= GEcombo.at(k-1) * GEcombo.at(l) * cm * DenomSquared;
      } // end l loop
    } // end k loop

    if(D[i]==0){
      for (unsigned int k=0; k<p; k++){
        zeta0.at(countcontrol-1,k) = SOmegaByS.at(k) - (ROmegaAccum.at(k,i)/RAccum.at(i)) - pibyn[0]*RightZetaAccum.at(k);
      } // end k loop
    } else {
      for (unsigned int k=0; k<p; k++){
        zeta1.at(countcase-1,k) = SOmegaByS.at(k) - (ROmegaAccum.at(k,i)/RAccum.at(i)) - pibyn[1]*RightZetaAccum.at(k);
      } // end k loop
    }
  } // end i loop
  Hessian = -1 * (LeftHessAccum-RightHessAccum);
  for (unsigned int k=0; k<p-1; k++){  // Upper triangle of Hessian
    for (unsigned int l=k+1; l<p; l++){
      Hessian.at(k,l) = Hessian.at(l,k);
    } // end l loop
  } // end k loop
  return List::create(Named("hessian")=Hessian,
                      Named("zeta0")=zeta0,
                      Named("zeta1")=zeta1);
}
