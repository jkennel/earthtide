// [[Rcpp::depends(RcppEigen)]]
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
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::VectorXcd;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::ArrayXd;
using Eigen::ArrayXi;


//==============================================================================
// [[Rcpp::export]]
Eigen::MatrixXd time_mat(const Eigen::ArrayXd& time) {

  size_t n = time.size();
  MatrixXd t_mat = MatrixXd::Ones(n, 5);

  t_mat.col(1) = time;
  t_mat.col(2) = t_mat.col(1).array() * time;
  t_mat.col(3) = t_mat.col(2).array() * time;
  t_mat.col(4) = t_mat.col(3).array() * time;

  return(t_mat);
}

// [[Rcpp::export]]
Eigen::MatrixXd time_der_mat(const Eigen::ArrayXd& time) {

  size_t n = time.size();
  MatrixXd t_mat = MatrixXd::Zero(n, 5);

  t_mat.col(1).setOnes();
  t_mat.col(2) = time * 2.0;
  t_mat.col(3) = t_mat.col(2).array() * time * 3.0; // time * time * 3.0 need to check
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

  double pm1 = legendre_bh(l, m - 1,  x);
  double pp1 = legendre_bh(l, m + 1,  x);

  return(0.5 * ((l + m) * (l - m + 1) * pm1 - pp1));

}

// [[Rcpp::export]]
double scale_legendre_bh(int l, int m) {

  double k;

  if (m == 0) {
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

  for (int l = 2; l <= l_max; ++l){
    for (int m = 0; m <= l; ++m){
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
    if (index[i] != index[i - 1]) {
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
  VectorXi out(n);

  for (size_t i = 0; i < n; ++i)
  {
    if (input[i] + 1 == 2) {
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
  VectorXd out(n);


  for (size_t i = 0; i < n; ++i)
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

  VectorXi index_unique = Eigen::Map<VectorXi, Eigen::Unaligned>(v.data(), v.size());

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

  size_t n = body_inds.size();
  ArrayXd out = body;

  double dc3;

  for (size_t i = 0; i < n; ++i) {
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
                               bool scale) {


  RowVector2d output;


  // is there a way to vectorize this?  matrix size issue
  const ArrayXd dc2 = calc_dc2(k_mat, astro, pk);

  const ArrayXd fac = set_fac(body,
                              body_inds,
                              k_mat,
                              astro_der,
                              delta,
                              deltar,
                              o1,
                              resonance,
                              max_amp);

  const Vector3d v(1.0, j2000, j2000 * j2000);
  const ArrayXd fac_x = fac * (x * v).array();
  const ArrayXd fac_y = fac * (y * v).array();

  const RowVectorXd dtham = (fac_x * fac_x + fac_y * fac_y).sqrt();
  const ArrayXd dthph = dc2 - fac_y.binaryExpr(fac_x, [] (double a, double b) { return std::atan2(a, b);} ).array();

  // determine phase correction
  const VectorXd cos_dc2 = Eigen::cos(dthph); //dc0
  const VectorXd sin_dc2 = Eigen::sin(dthph); //ds0

  double cc = dtham * cos_dc2;
  double ss = dtham * sin_dc2;

  if (scale) {
    cc = cc / dtham.maxCoeff();
    ss = ss / dtham.maxCoeff();
  }

  output(0) = cc;
  output(1) = ss;


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
                      size_t max_amp) {


  double output;


  // is there a way to vectorize this?  matrix size issue
  const ArrayXd dc2 = calc_dc2(k_mat, astro, pk);

  const ArrayXd fac = set_fac(body,
                              body_inds,
                              k_mat,
                              astro_der,
                              delta,
                              deltar,
                              o1,
                              resonance,
                              max_amp);

  const Vector3d v(1.0, j2000, j2000 * j2000);

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
  const ArrayXd amplitude = (x.col(0).array() * x.col(0).array() +
                             y.col(0).array() * y.col(0).array()).sqrt();

  // get the subsets for each wave group
  const MatrixXi sub = get_catalog_indices(index, ng);

  MatrixXd output;

  if (predict) {
    output = MatrixXd::Zero(nt, 1);
  } else {
    output = MatrixXd::Zero(nt, ng * 2);
  }

  // subset for each wave group
  for (std::size_t j = 0; j < ng; ++j) {

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
      // single curve - subset for each time
      RcppThread::parallelFor(0, nt, [&] (size_t k) {

        output(k) += mult * et_predict_one(
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
          i_max);
      });
    } else {
      // separate curves - subset for each time
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
          scale);
      });
    }
  }

  return(output);
}

/*** R

*/
