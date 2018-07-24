// Functions for the bivariate Gaussian copula

#include <Rcpp.h>
using namespace Rcpp;

// Global variables -------------------------
extern double TOL; // defined in correst_c.cpp


// This is wrong. Need to multiply by derivatives of inverse cdf's of normals.
// Density function of the bivariate normal copula.
//
// @param x The value of the first element. Must be between 0 and 1.
// @param y The value of the second element. Must be between 0 and 1.
// @param rho The correlation. Must be between 0 and 1.
// @param log_p A logical. Should we return the log-density (\code{TRUE}) or not (\code{FALSE})?
//
// @author David Gerard
//
// double dnormcop(double x, double y, double rho, bool log_p = true) {
//   // check input -------------------------------------------------
//   if ((x > 1.0 + TOL) | (x < -TOL)) {
//     Rcpp::stop("dnormcop: x should be between 0 and 1.");
//   }
//   if ((y > 1.0 + TOL) | (y < -TOL)) {
//     Rcpp::stop("dnormcop: y should be between 0 and 1.");
//   }
//   if ((rho > 1.0) | (rho < 0.0)) {
//     Rcpp::stop("dnormcop: rho should be between 0 and 1.");
//   }
//
//   double ldense;
//   if ((x < TOL) | (y < TOL) | ((1.0 - y) < TOL) | ((1.0 - x) < TOL)) {
//     ldense = R_NegInf;
//   } else {
//     // Invert standard normal ---------------------------------------
//     double phi_inv_x = R::qnorm5(x, 0.0, 1.0, 1, 0); // Phi^{-1}(x)
//     double phi_inv_y = R::qnorm5(y, 0.0, 1.0, 1, 0); // Phi^{-1}(y)
//
//     // Bivariate normal density -------------------------------------
//     ldense = -1.0 * std::log(2.0 * M_PI) -
//       std::log(1.0 - std::pow(rho, 2.0)) / 2.0 -
//       (std::pow(phi_inv_x, 2.0) + std::pow(phi_inv_y, 2.0) - 2.0 * rho * phi_inv_x * phi_inv_y) / (2.0 * (1.0 - std::pow(rho, 2.0)));
//   }
//
//   // Exp if desired -----------------------------------------------
//   double return_val = 0.0;
//   if (log_p) {
//     return_val = ldense;
//   } else {
//     return_val = std::exp(ldense);
//   }
//
//   return return_val;
// }

//' A C++ implementation of bivariate normal probabilities, derived from the MATLAB code of Alan Genz.
//'
//' Calculates the probability that x > dh and y > dk.
//'
//' Original Matlab code written by:
//' Alan Genz, Department of Mathematics,
//' Washington State University,
//' Pullman, Wa 99164-3113, Email: alangenz@wsu.edu,
//'
//' @param dh 1st lower integration limit
//' @param dk 2nd lower integration limit
//' @param r Correlation coefficient
//'
//'
//' @references Drezner, Z., & Wesolowsky, G. O. (1990). On the computation of the bivariate normal integral. Journal of Statistical Computation and Simulation, 35(1-2), 101-107.
//'
//' @author David Gerard
// [[Rcpp::export]]
double bvnu(double dh, double dk, double r) {
  double p;
  if (dh == R_PosInf | dk == R_PosInf) {
    p = 0;
  } else if (dh == R_NegInf) {
    if (dk == R_NegInf) {
      p = 1;
    } else {
      R::pnorm5(-dk, 0.0, 1.0, 1, 0);
    }
  } else if (dk == R_NegInf) {
    p = R::pnorm5(-dh, 0.0, 1.0, 1, 0);
  } else if (std::abs(r) < TOL) {
    p = R::pnorm5(-dh, 0.0, 1.0, 1, 0) * R::pnorm5(-dk, 0.0, 1.0, 1, 0);
  } else {
    double tp = 2.0 * M_PI;
    double h = dh;
    double k = dk;
    double hk = h * k;
    double bvn = 0.0;
    double hs;
    double asr;
    NumericVector sn;
    NumericVector x;
    NumericVector w;
    if (std::abs(r) < 0.3) {
      // Gauss Legendre points and weights, n = 6
      w = NumericVector::create(0.1713244923791705, 0.3607615730481384, 0.4679139345726904,
                                0.1713244923791705, 0.3607615730481384, 0.4679139345726904);
      x = NumericVector::create(1.0 - 0.9324695142031522, 1.0 - 0.6612093864662647, 1.0 - 0.2386191860831970,
                                1.0 + 0.9324695142031522, 1.0 + 0.6612093864662647, 1.0 + 0.2386191860831970);
    } else if (std::abs(r) < 0.75) {
      // Gauss Legendre points and weights, n = 12
      w = NumericVector::create(0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                                0.2031674267230659, 0.2334925365383547, 0.2491470458134029,
                                0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                                0.2031674267230659, 0.2334925365383547, 0.2491470458134029);
      x = NumericVector::create(1.0 - 0.9815606342467191, 1.0 - 0.9041172563704750, 1.0 - 0.7699026741943050,
                                1.0 - 0.5873179542866171, 1.0 - 0.3678314989981802, 1.0 - 0.1252334085114692,
                                1.0 + 0.9815606342467191, 1.0 + 0.9041172563704750, 1.0 + 0.7699026741943050,
                                1.0 + 0.5873179542866171, 1.0 + 0.3678314989981802, 1.0 + 0.1252334085114692);
    } else {
      // Gauss Legendre points and weights, n = 20
      w = NumericVector::create(0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
                                0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
                                0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
                                0.1527533871307259,
                                0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
                                0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
                                0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
                                0.1527533871307259);
      x = NumericVector::create(1.0 - 0.9931285991850949, 1.0 - 0.9639719272779138, 1.0 - 0.9122344282513259,
                                1.0 - 0.8391169718222188, 1.0 - 0.7463319064601508, 1.0 - 0.6360536807265150,
                                1.0 - 0.5108670019508271, 1.0 - 0.3737060887154196, 1.0 - 0.2277858511416451,
                                1.0 - 0.07652652113349733,
                                1.0 + 0.9931285991850949, 1.0 + 0.9639719272779138, 1.0 + 0.9122344282513259,
                                1.0 + 0.8391169718222188, 1.0 + 0.7463319064601508, 1.0 + 0.6360536807265150,
                                1.0 + 0.5108670019508271, 1.0 + 0.3737060887154196, 1.0 + 0.2277858511416451,
                                1.0 + 0.07652652113349733);
    }
    if (std::abs(r) < 0.925) {
      hs = (h * h + k * k) / 2.0;
      asr = std::asin(r) / 2.0;
      sn = Rcpp::sin(asr * x);
      bvn = Rcpp::sum(Rcpp::exp((sn * hk - hs) / (1 - Rcpp::pow(sn, 2.0))) * w);
      bvn = bvn * asr / tp + R::pnorm5(-h, 0.0, 1.0, 1, 0) * R::pnorm5(-k, 0.0, 1.0, 1, 0);
    } else {
      NumericVector edave;
      NumericVector wsub;
      NumericVector rs;
      NumericVector ep;
      NumericVector asr_vec;
      NumericVector sp_vec;
      NumericVector xs;
      NumericVector xsub;
      LogicalVector ix;
      double as;
      double a;
      double bs;
      double b;
      double c;
      double d;
      double sp;
      if (r < 0.0) {
        k = -k;
        hk = -hk;
      }
      if (std::abs(r) < 1.0) {
        as = 1.0 - r * r;
        a = std::sqrt(as);
        bs = std::pow(h - k, 2.0);
        asr = -(bs / as + hk) / 2.0;
        c = (4.0 - hk) / 8.0;
        d = (12.0 - hk) / 80.0;
        if (asr > -100.0) {
          bvn = a * std::exp(asr) * (1.0 - c * (bs - as) * (1.0 - d * bs) / 3.0 + c * d * std::pow(as, 2.0));
        }
        if (hk > -100.0) {
          b = std::sqrt(bs);
          sp = std::sqrt(tp) * R::pnorm5(-b / a, 0.0, 1.0, 1, 0);
          bvn = bvn - std::exp(-hk / 2) * sp * b * (1.0 - c * bs * (1.0 - d * bs) / 3);
        }
        a = a / 2.0;
        xs = Rcpp::pow(a * x, 2.0);
        asr_vec = -(bs / xs + hk) / 2.0;
        ix =  asr_vec > -100.0;
        xsub = xs[ix];
        sp_vec = (1.0 + c * xsub * (1.0 + 5.0 * d * xsub));
        rs = Rcpp::sqrt(1.0 - xsub);
        ep = Rcpp::exp(-(hk / 2.0) * xsub / Rcpp::pow(1.0 + rs, 2.0)) / rs;
        edave = Rcpp::exp(asr_vec[ix]) * (sp_vec - ep);
        wsub = w[ix];
        bvn = (a * Rcpp::sum(edave * wsub) - bvn) / tp;
      }
      if (r > 0.0) {
        bvn =  bvn + R::pnorm5(-std::max(h, k), 0.0, 1.0, 1, 0);
      } else if (h >= k) {
        bvn = -bvn;
      } else {
        double L;
        if (h < 0.0) {
          L = R::pnorm5(k, 0.0, 1.0, 1, 0) - R::pnorm5(h, 0.0, 1.0, 1, 0);
        } else {
          L = R::pnorm5(-h, 0.0, 1.0, 1, 0) - R::pnorm5(-k, 0.0, 1.0, 1, 0);
        }
        bvn = L - bvn;
      }
    }
    p = std::max(0.0, std::min(1.0, bvn));
  }
  return p;
}

//' A C++ implementation of the bivariate normal CDF, derived from the MATLAB code of Alan Genz
//'
//' Calculates the probability that x < dh and y < dk
//'
//' Original Matlab code written by: Alan Genz,
//' Department of Mathematics, Washington State University,
//' Pullman, Wa 99164-3113, Email: alangenz@wsu.edu,
//'
//' @param dh 1st upper integration limit
//' @param dk 2nd upper integration limit
//' @param r Correlation coefficient
//'
//'
//' @references Drezner, Z., & Wesolowsky, G. O. (1990). On the computation of the bivariate normal integral. Journal of Statistical Computation and Simulation, 35(1-2), 101-107.
//'
//'
//' @author David Gerard
// [[Rcpp::export]]
double bvnl(double dh, double dk, double r) {
  double p = bvnu(-dh, -dk, r);
  return p;
}


//' Distribution function of the bivariate Gaussian copula.
//'
//' @param x The value of the first element. Must be between 0 and 1.
//' @param y The value of the second element. Must be between 0 and 1.
//' @param rho The correlation. Must be between 0 and 1.
//'
//' @author David Gerard
//'
//' @seealso \code{\link{bvnl}} For the cdf of the bivariate normal.
//'
//'
// [[Rcpp::export]]
double pnormcop(double x, double y, double rho) {
  // check input -------------------------------------------------
  if ((x > 1.0 + TOL) | (x < -TOL)) {
    Rcpp::stop("dnormcop: x should be between 0 and 1.");
  }
  if ((y > 1.0 + TOL) | (y < -TOL)) {
    Rcpp::stop("dnormcop: y should be between 0 and 1.");
  }
  if ((rho > 1.0) | (rho < -1.0)) {
    Rcpp::stop("dnormcop: rho should be between 0 and 1.");
  }

  double phi_inv_x;
  double phi_inv_y;

  double lower_prob;
  if ((x < TOL) | (y < TOL)) {
    lower_prob = 0.0;
  } else if (((1.0 - x) < TOL) & ((1.0 - y) < TOL)) {
    lower_prob = 1.0;
  } else if ((1.0 - x) < TOL) {
    phi_inv_y = R::qnorm5(y, 0.0, 1.0, 1, 0); // Phi^{-1}(y)
    lower_prob = R::pnorm5(phi_inv_y, 0.0, 1.0, 1, 0);
  } else if ((1.0 - y) < TOL) {
    phi_inv_x = R::qnorm5(x, 0.0, 1.0, 1, 0); // Phi^{-1}(x)
    lower_prob = R::pnorm5(phi_inv_x, 0.0, 1.0, 1, 0);
  } else {
    phi_inv_x = R::qnorm5(x, 0.0, 1.0, 1, 0); // Phi^{-1}(x)
    phi_inv_y = R::qnorm5(y, 0.0, 1.0, 1, 0); // Phi^{-1}(y)
    lower_prob = bvnl(phi_inv_x, phi_inv_y, rho);
  }
  return(lower_prob);
}
