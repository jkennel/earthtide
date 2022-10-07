// This function is faster than the R based ones though.

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(BH)]]

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <RcppEigen.h>
#include <Eigen/StdVector>

#include <RcppThread.h>

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::ArrayXd;
using Eigen::ArrayXi;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using Eigen::Vector3d;
using Eigen::Vector2d;


//==============================================================================
// [[Rcpp::export]]
Eigen::MatrixXd time_mat(const Eigen::ArrayXd& time) {

  size_t n = time.size();
  MatrixXd t_mat = MatrixXd::Ones(n,5);

  t_mat.col(1) = time;
  t_mat.col(2) = t_mat.col(1).array() * time;
  t_mat.col(3) = t_mat.col(2).array() * time;
  t_mat.col(4) = t_mat.col(3).array() * time;

  return(t_mat);
}

// [[Rcpp::export]]
Eigen::MatrixXd time_der_mat(const Eigen::ArrayXd& time) {

  size_t n = time.size();
  MatrixXd t_mat = MatrixXd::Zero(n,5);

  t_mat.col(1).setOnes();
  t_mat.col(2) = time * 2.0;
  t_mat.col(3) = t_mat.col(2).array() * time * 3.0; // time*time*3.0 need to check
  t_mat.col(4) = t_mat.col(3).array() * time * 4.0;

  return(t_mat);

}

// [[Rcpp::export]]
Eigen::MatrixXd astro(const Eigen::ArrayXd& t_astro,
                      const Eigen::MatrixXd simon,
                      double longitude,
                      Eigen::RowVectorXd hours,
                      Eigen::RowVectorXd ddt) {

  MatrixXd at_mat = simon * time_mat(t_astro).transpose();


  // First row needs correction
  at_mat.row(0) = at_mat.row(2).array() - at_mat.row(1).array() +
    (longitude + hours.array() * 15.0) -
    (0.0027 * ddt.array() * 15.0 / 3600.0);

  // modulus
  at_mat = at_mat.array() - (at_mat.array() / 360).floor() * 360;

  return(at_mat);
}



// [[Rcpp::export]]
Eigen::MatrixXd astro_der(const Eigen::ArrayXd& t_astro,
                          const Eigen::MatrixXd simon) {

  double time_scale = 1.0 / (365250.0 * 24.0);
  MatrixXd at_mat = simon * time_der_mat(t_astro).transpose() * time_scale;

  // First row needs correction
  at_mat.row(0) = at_mat.row(2).array() - at_mat.row(1).array() + 15.0;

  return(at_mat);
}



// [[Rcpp::export]]
double legendre_bh(int l, int m, double x, int csphase = -1) {

  return(boost::math::legendre_p(l, m, x) * std::pow(csphase, m));

}

// [[Rcpp::export]]
double legendre_deriv_bh(int l, int m, double x) {

  double pm1 = legendre_bh(l, m-1,  x);
  double pp1 = legendre_bh(l, m+1,  x);

  return(0.5 * ((l+m) * (l-m+1) * pm1 - pp1));

}

// [[Rcpp::export]]
double scale_legendre_bh(int l, int m) {

  double k;

  if(m == 0) {
    k = 1.0;
  } else {
    k = 2.0;
  }

  double num   = boost::math::factorial<double>(l - m);
  double denom = boost::math::factorial<double>(l + m);

  return(std::sqrt((double) (k * (2.0 * l + 1.0) * num /  denom)));

}

// [[Rcpp::export]]
Eigen::MatrixXd legendre(int l_max, double x) {

  double scale;

  size_t n = VectorXi::LinSpaced(l_max - 1, 3, l_max + 1).sum();
  MatrixXd out(n, 4);


  int i = 0;

  for (int l=2; l <= l_max; l++){
    for(int m=0; m <= l; m++){
      scale = scale_legendre_bh(l, m);
      out(i, 0) = l;
      out(i, 1) = m;
      out(i, 2) = legendre_bh(l, m, x) * scale;
      out(i, 3) = legendre_deriv_bh(l, m, x) * scale;
      i += 1;
    }
  }

  return(out);
}

// [[Rcpp::export]]
Eigen::MatrixXi get_catalog_indices(Eigen::VectorXi index, size_t ng) {

  size_t nw = index.size();
  size_t counter = 1;

  MatrixXi inds(ng, 2);

  inds(0, 0) = 0;
  inds(ng - 1, 1) = nw - 1;

  for (size_t i = 1; i < nw; ++i) {
    if(index[i] != index[i-1]) {
      inds(counter, 0) = i;
      inds(counter - 1, 1) = i - 1;
      counter = counter + 1;
    }
  }

  return(inds);
}

// [[Rcpp::export]]
Eigen::VectorXi subset_2_eigen(const Eigen::VectorXi& input)
{
  size_t n = input.size();
  size_t counter = 0;
  Eigen::VectorXi out(n);

  for (size_t i=0; i < n; ++i)
  {
    if(input[i] + 1 == 2) {
      out[counter] = i;
      counter += 1;
    }
  }

  out.conservativeResize(counter);

  return out;
}

// [[Rcpp::export]]
Eigen::ArrayXd subset_eigen(const Eigen::ArrayXd& input,
                            const Eigen::VectorXi& subs)
{
  size_t n = subs.size();
  Eigen::VectorXd out(n);


  for (size_t i=0; i < n; ++i)
  {
    out[i] = input[subs[i]];
  }

  return out;
}


//[[Rcpp::export]]
Eigen::VectorXi unique_eigen(Eigen::VectorXi index) {

  std::vector<int> v(index.data(), index.data() + index.size());

  std::sort(v.begin(), v.end());
  std::vector<int>::iterator it;
  it = std::unique(v.begin(), v.end());
  v.resize(std::distance(v.begin(), it));

  Eigen::VectorXi index_unique = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(v.data(), v.size());

  return(index_unique);

}

//[[Rcpp::export]]
Eigen::ArrayXd calc_dc2(const Eigen::MatrixXd& k_mat,
                        const Eigen::VectorXd& astro,
                        const Eigen::ArrayXd& pk) {

  double to_rad = M_PI / 180.0;

  // is there a way to vectorize this?  matrix size issue
  ArrayXd dc2 = (k_mat * astro).array() + pk + 360.0;
  dc2 = (dc2 - (dc2 / 360).floor() * 360) * to_rad;

  return(dc2);
}

// [[Rcpp::export]]
Eigen::VectorXd set_fac(const Eigen::ArrayXd& body,
                        const Eigen::ArrayXi& body_inds,
                        const Eigen::MatrixXd& k_mat,
                        const Eigen::VectorXd& astro_der,
                        double delta,
                        double deltar,
                        double o1,
                        double resonance,
                        size_t max_amp
)
{

  // size_t n = body.size();
  size_t n = body_inds.size();
  Eigen::ArrayXd out = body;

  // size_t k;
  double dc3;

  for(size_t i = 0; i < n; ++i) {
    dc3 = k_mat.row(body_inds[i]) * astro_der;
    out(body_inds[i]) = delta + deltar * (dc3 - o1) / (resonance - dc3);
  }

  out = out / out[max_amp];

  return(out);
}

// [[Rcpp::export]]
Eigen::MatrixXd et_analyze_one(const Eigen::VectorXd& astro,
                               const Eigen::VectorXd& astro_der,
                               const Eigen::MatrixXd& k_mat,
                               const Eigen::ArrayXd& pk,
                               const Eigen::ArrayXd& body,
                               const Eigen::ArrayXi& body_inds,
                               double delta,
                               double deltar,
                               const Eigen::MatrixXd& x,
                               const Eigen::MatrixXd& y,
                               double j2000,
                               double o1,
                               double resonance,
                               size_t max_amp,
                               double update_coef,
                               bool scale) {


  MatrixXd output(1, 2);


  // is there a way to vectorize this?  matrix size issue
  ArrayXd dc2 = calc_dc2(k_mat, astro, pk);

  const ArrayXd fac = set_fac(body,
                              body_inds,
                              k_mat,
                              astro_der,
                              delta,
                              deltar,
                              o1,
                              resonance,
                              max_amp);

  const Eigen::Vector3d v(1.0, j2000, j2000 * j2000);
  ArrayXd fac_x = fac * (x * v).array();
  ArrayXd fac_y = fac * (y * v).array();


  RowVectorXd dtham = (fac_x * fac_x + fac_y * fac_y).sqrt();
  ArrayXd dthph = dc2 -  fac_y.binaryExpr(fac_x, [] (double a, double b) { return std::atan2(a, b);} );

  // determine phase correction
  VectorXd cos_dc2 = dthph.cos(); //dc0
  VectorXd sin_dc2 = dthph.sin(); //ds0

  double cc = dtham * cos_dc2;
  double ss = dtham * sin_dc2;

  if (scale) {
    cc = cc / dtham.maxCoeff();
    ss = ss / dtham.maxCoeff();
  }

  output(0, 0) = cc;
  output(0, 1) = ss;


  // output = (fac * (((x0 + x1 * j2000 + x2 * j2000_sq) * cos_dc2) +
  //                  ((y0 + y1 * j2000 + y2 * j2000_sq) * sin_dc2))).sum();


  return(output);

}

// [[Rcpp::export]]
double et_predict_one(const Eigen::VectorXd& astro,
                      const Eigen::VectorXd& astro_der,
                      const Eigen::MatrixXd& k_mat,
                      const Eigen::ArrayXd& pk,
                      const Eigen::ArrayXd& body,
                      const Eigen::ArrayXi& body_inds,
                      double delta,
                      double deltar,
                      const Eigen::MatrixXd& x,
                      const Eigen::MatrixXd& y,
                      double j2000,
                      double o1,
                      double resonance,
                      size_t max_amp,
                      double update_coef) {


  double output;


  // is there a way to vectorize this?  matrix size issue
  ArrayXd dc2 = calc_dc2(k_mat, astro, pk);

  const ArrayXd fac = set_fac(body,
                              body_inds,
                              k_mat,
                              astro_der,
                              delta,
                              deltar,
                              o1,
                              resonance,
                              max_amp);

  const Eigen::Vector3d v(1.0, j2000, j2000 * j2000);

  output = (fac * ((x * v).array().colwise() * dc2.cos() +
    (y * v).array().colwise() * dc2.sin())).sum();


  // output = (fac * (((x0 + x1 * j2000 + x2 * j2000_sq) * cos_dc2) +
  //                  ((y0 + y1 * j2000 + y2 * j2000_sq) * sin_dc2))).sum();


  return(output);

}


//[[Rcpp::export]]
Eigen::MatrixXd et_calculate(const Eigen::MatrixXd& astro,
                             const Eigen::MatrixXd& astro_der,
                             const Eigen::MatrixXd& k_mat,
                             const Eigen::ArrayXd& phases,
                             const Eigen::ArrayXd& delta,
                             double deltar,
                             const Eigen::MatrixXd& cc,
                             const Eigen::MatrixXd& ss,
                             const Eigen::ArrayXd& dgk,
                             const Eigen::VectorXi& jcof,
                             const Eigen::ArrayXd& j2000,
                             double o1,
                             double resonance,
                             const Eigen::VectorXi& index,
                             size_t astro_update,
                             double update_coef,
                             const Eigen::ArrayXd& multiplier,
                             bool predict,
                             bool scale) {



  Eigen::Index max_elem;
  size_t i_max;

  // number of times
  size_t nt = astro.cols();

  // number of wave groups
  const VectorXi un = unique_eigen(index);
  size_t ng = un.size();
  size_t start_seg, n_seg;

  // sin and cos terms
  const ArrayXd dgk_sub = jcof.unaryExpr(dgk);
  const MatrixXd x = cc.array().colwise() * dgk_sub * 1.e-10;
  const MatrixXd y = ss.array().colwise() * dgk_sub * 1.e-10;

  // get equilibrium wave amplitude
  const VectorXd amplitude = (pow(x.col(0).array(), 2) +
                              pow(y.col(0).array(), 2)).sqrt();

  // get the subsets for each wave group
  const MatrixXi sub = get_catalog_indices(index, ng);

  MatrixXd output;

  if(predict) {
    output = Eigen::MatrixXd::Zero(nt,1);
  } else {
    output = Eigen::MatrixXd::Zero(nt,ng*2);
  }

  // subset for each wave group
  for(std::size_t j = 0; j < ng; ++j) {

    start_seg = sub(j, 0);
    n_seg = sub(j, 1) - sub(j, 0) + 1;

    const ArrayXd pk   = subset_eigen(phases, jcof.segment(start_seg, n_seg));
    const ArrayXd body = subset_eigen(delta, jcof.segment(start_seg, n_seg));

    amplitude.segment(start_seg, n_seg).maxCoeff(&max_elem);
    i_max = max_elem;
    const ArrayXi body_inds = subset_2_eigen(jcof.segment(start_seg, n_seg));

    const MatrixXd k_mat_sub = k_mat.middleRows(start_seg, n_seg);
    const MatrixXd x_sub = x.middleRows(start_seg, n_seg);
    const MatrixXd y_sub = y.middleRows(start_seg, n_seg);
    double mult = multiplier[j];

    if (predict) {
      // subset for each time
      RcppThread::parallelFor(0, nt, [&] (size_t k) {

        output(k,0) += mult * et_predict_one(
          astro.col(k),
          astro_der.col(k),
          k_mat_sub,
          pk,
          body,
          body_inds,
          delta(1),
          deltar,
          x_sub,
          y_sub,
          j2000(k),
          o1,
          resonance,
          i_max,
          update_coef);
      });
    } else {
      // subset for each time
      RcppThread::parallelFor(0, nt, [&] (size_t k) {

        output.block(k, j * 2, 1, 2) = mult * et_analyze_one(
          astro.col(k),
          astro_der.col(k),
          k_mat_sub,
          pk,
          body,
          body_inds,
          delta(1),
          deltar,
          x_sub,
          y_sub,
          j2000(k),
          o1,
          resonance,
          i_max,
          update_coef,
          scale);
      });
    }
  }

  return(output);
}


// ***********Just for testing*********
// //[[Rcpp::export]]
// Eigen::ArrayXd fac_xy(Eigen::MatrixXd x, Eigen::ArrayXd fac, double j2000 ){
//   const Eigen::Vector3d v(1.0, j2000, j2000 * j2000);
//   ArrayXd fac_x = fac * (x * v).array();
//   return(fac_x);
// }
//
// //[[Rcpp::export]]
// Eigen::RowVectorXd fac_2(Eigen::ArrayXd fac_x, Eigen::ArrayXd fac_y){
//   RowVectorXd dtham = (fac_x * fac_x + fac_y * fac_y).sqrt();
//   return(dtham);
// }
//
// //[[Rcpp::export]]
// Eigen::VectorXd phase(Eigen::ArrayXd dc2, Eigen::ArrayXd fac_x, Eigen::ArrayXd fac_y){
//
//   VectorXd dthph = dc2 -  fac_y.binaryExpr(fac_x, [] (double a, double b) { return std::atan2(a, b);} );
//
//   return(dthph.cos());
// }

/*** R

set.seed(10)
dc2 <- earthtide:::calc_dc2(matrix(rnorm(1100), ncol = 11), matrix(rnorm(11), ncol = 1), 1:100)
a   <- earthtide:::fac_xy(matrix(rnorm(300), ncol = 3), 1:100, 0.5)
b   <- earthtide:::fac_xy(matrix(rnorm(300), ncol = 3), 1:100, 0.5)
amp <- earthtide:::fac_2(a,b);
earthtide:::phase(dc2, a,b);



library(earthtide)
tms <- seq.POSIXt(as.POSIXct('1995-01-01', tz = 'UTC'), as.POSIXct('1995-01-01 00:00:30', tz = 'UTC'),1)
wave_groups = data.frame(start = c(0, 1), end = c(1,2))
et <- Earthtide$new(utc = tms,
                    latitude = 49.00937,
                    longitude = 8.40444,
                    elevation = 120,
                    gravity = 9.8127,
                    cutoff = 1.0e-10,
                    wave_groups = wave_groups)
x <- 1:1e6
astro_args <- as.matrix(earthtide:::simon_coef_1994[, -1])
longitude <- 8.40444
tms_2 <- earthtide:::.prepare_datetime(tms)
# recipes6:::time_mat(tms)


bench::mark(
  a <- t(recipes6:::time_mat(tms)),
  b <- earthtide:::time_mat(tms),
  check = TRUE
)

bench::mark(
  a <- t(recipes6:::time_der_mat(tms)),
  b <- earthtide:::time_der_mat(tms),
  check = TRUE
)


bench::mark(
  a <-  recipes6:::astro(tms_2$t_astro, astro_args, longitude, tms_2$hours, tms_2$ddt),
  b <- earthtide:::astro(tms_2$t_astro, astro_args, longitude, tms_2$hours, tms_2$ddt),
  check = TRUE
)

bench::mark(
  a <-  recipes6:::astro_der(tms_2$t_astro, astro_args),
  b <- earthtide:::astro_der(tms_2$t_astro, astro_args),
  check = TRUE
)


bench::mark(
  a <- recipes6:::legendre(5, 0.1),
  b <- earthtide:::legendre(5, 0.1),
  check = TRUE
)




library(data.table)
library(recipes6)
library(earthtide)

wave_groups = na.omit(eterna_wavegroups[eterna_wavegroups$time == '1 month',])

tms <- seq.POSIXt(as.POSIXct('1995-01-01', tz = 'UTC'), as.POSIXct('1995-01-01 02:00:00', tz = 'UTC'), 1)

et <- Earthtide$new(utc = tms,
                    latitude = 52.3868,
                    longitude = 9.7144,
                    elevation = 110,
                    gravity = 9.8127,
                    cutoff = 1.0e-10,
                    catalog = 'hw95s',
                    wave_groups = wave_groups)
et$analyze(method = 'tidal_potential',
           scale = FALSE)
bench::mark(
et$analyze(method = 'tidal_potential')
)
et$tides

et$delta[1:12] <- et$love_params$dklat
et$deltar <- et$love_params$dkr

x <- cbind(et$catalog$c0,
           et$catalog$c1,
           et$catalog$c2)
y <- cbind(et$catalog$s0,
           et$catalog$s1,
           et$catalog$s2)

a <- bench::mark(
  cc <- recipes6:::et_calculate(et$astro$astro,
                                et$astro$astro_der,
                                et$catalog$k,
                                et$pk,
                                et$delta,
                                et$deltar,
                                x,
                                y,
                                et$station$dgk,
                                et$catalog$jcof - 1,
                                et$datetime$j2000,
                                et$love_params$dom0,
                                et$love_params$domr,
                                et$catalog$id,
                                1,
                                et$update_coef,
                                et$catalog$wave_groups$multiplier,
                                TRUE,
                                TRUE),
  dd <- as.numeric(et$calculate()),
  check = TRUE
)
setDT(a)
a
head(dd)
head(cc)
plot(tidal_potential~datetime, et$tides, type = 'l')



# const Eigen::ArrayXd& fac,
# const Eigen::MatrixXd& x,
# const Eigen::MatrixXd& y,
# double j2000,
# const Eigen::ArrayXd& cos_dc2,
# const Eigen::ArrayXd& sin_dc2
n <- 10000
fac <- rnorm(n)
x <- matrix(rnorm(3*n), ncol = 3)
y <- matrix(rnorm(3*n), ncol = 3)
cos_dc2 <-  rnorm(n)
tmp <- matrix(cos_dc2, nrow = 1)
sin_dc2 <-  rnorm(n)
tmp2 <- matrix(sin_dc2, nrow = 1)

sum((fac*tmp) %*% x * c(1,22,22*22) + (fac*tmp2) %*% y * c(1,22,22*22))

a <- bench::mark(
  b <- recipes6:::test_2(fac,
                         x,
                         y,
                         22.0,
                         cos_dc2,
                         sin_dc2),
  d <- recipes6:::test_3(fac,
                         x,
                         y,
                         22.0,
                         tmp,
                         tmp2),
  check = FALSE
)
data.table::setDT(a)
a


c <- bench::mark(
  d <- recipes6:::test_1(fac,
                         x[,1],
                         y[,1],
                         x[,2],
                         y[,2],
                         x[,3],
                         y[,3],
                         22.0,
                         cos_dc2,
                         sin_dc2)
)
setDT(c)
c

*/

