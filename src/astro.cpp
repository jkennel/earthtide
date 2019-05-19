
#define BOOST_DISABLE_ASSERTS
#define ARMA_DONT_PRINT_ERRORS

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]


#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace RcppParallel;
using namespace Rcpp;


/**
 * Possible scenarios
 *   - single sum of entire earth tides
 *   - sum for each earthtide group
 *   - sin, cos for each earthtide group
 */



//==============================================================================
/**
 * Extend division reminder to vectors
 * https://stackoverflow.com/questions/27686319/c-armadillo-modulus-function
 *
 * @param   a       Dividend
 * @param   n       Divisor
 */
template<typename T>
T mod(T a, int n)
{
  return a - floor(a / n) * n;
}



//==============================================================================
// @title
// time_mat
//
// @description
// Matrix to multiply by astrological.  Code adapted from ETERNA.
//
// @param time the times to calculate astrological parameters
//
// @return matrix of times
//
// @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
//
// @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// processing package ETERNA 3.3. Bulletin d'Informations
// Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
// 
// @keywords internal
//
// [[Rcpp::export]]
arma::mat time_mat(const arma::rowvec time) {
  
  unsigned len = time.n_elem;
  arma::mat t_mat(5, len);
  
  t_mat.ones();
  t_mat.row(1) = time;
  t_mat.row(2) = t_mat.row(1) % time;
  t_mat.row(3) = t_mat.row(2) % time;
  t_mat.row(4) = t_mat.row(3) % time;
  
  return(t_mat);
}


//==============================================================================
// @title
// time_der_mat
//
// @description
// Matrix to multiply by astrological.  Code adapted from ETERNA.
//
// @param time the times to calculate the derivatives of the astrological parameters
//
// @return matrix of time derivatives
//
// @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
//
// @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// processing package ETERNA 3.3. Bulletin d'Informations
// Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
// 
// @keywords internal
//
// [[Rcpp::export]]
arma::mat time_der_mat(const arma::rowvec time) {
  
  unsigned len = time.n_elem;
  arma::mat td_mat(5, len);
  
  td_mat.zeros();
  td_mat.row(1).ones();
  td_mat.row(2) = time * 2.0;
  td_mat.row(3) = td_mat.row(1) % time * 3.0;
  td_mat.row(4) = td_mat.row(2) % time * 4.0;
  
  return(td_mat);
}


//==============================================================================
// @title
// astro
//
// @description
// Calculate astronomical parameters.  Code adapted from ETERNA.
//
// @param t_astro astronomical time
// @param simon coefficient matrix
// @param longitude the longitude
// @param hours hour or measurement
// @param ddt time in ddt
//
// @return astonomical parameters
//
// @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
//
// @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// processing package ETERNA 3.3. Bulletin d'Informations
// Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
// 
// @keywords internal
//
// [[Rcpp::export]]
arma::mat astro(const arma::rowvec t_astro,
                const arma::mat simon,
                double longitude,
                const arma::rowvec hours,
                const arma::rowvec ddt) {
  
  arma::mat at_mat = simon * time_mat(t_astro);
  
  at_mat.row(0) = at_mat.row(2) - at_mat.row(1) + longitude + hours * 15.0 -
    0.0027 * ddt * 15.0 / 3600.0;
  
  at_mat = mod(at_mat, 360);
  
  return(at_mat);
}





//==============================================================================
// @title
// astro_der
//
// @description
// Calculate derivatives of astronomical parameters.  Code adapted from ETERNA.
//
// @param t_astro astronomical time
// @param simon coefficient matrix
//
// @return derivative astonomical parameters
//
// @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
//
// @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// processing package ETERNA 3.3. Bulletin d'Informations
// Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
//
//
// [[Rcpp::export]]
arma::mat astro_der(const arma::rowvec t_astro,
                    const arma::mat simon) {
  
  double time_scale = 1.0 / (365250.0 * 24.0);
  
  arma::mat at_mat = simon * time_der_mat(t_astro) * time_scale;
  
  at_mat.row(0) = at_mat.row(2) - at_mat.row(1) + 15.0;
  
  return(at_mat);
  
}


// helpful content https://www.mat.univie.ac.at/~westra/associatedlegendrefunctions.pdf
//==============================================================================
// @title
// legendre_bh
//
// @description
// legendre function
//
// @param l maximum degree
// @param m the m value
// @param x the x value
// @param csphase the csphase value
//
// @return legendre scale
// 
// @keywords internal
//
// [[Rcpp::export]]
double legendre_bh(int l, int m, double x, int csphase = -1) {
  
  
  return(boost::math::legendre_p(l, m, x) * std::pow(csphase, m));
  
}


// helpful content https://www.mat.univie.ac.at/~westra/associatedlegendrefunctions.pdf
//==============================================================================
// @title
// legendre_deriv_bh
//
// @description
// Scaling for legendre function
//
// @param l maximum degree
// @param m the x value
// @param x the x value
//
// @return legendre scale
//
// @keywords internal
//
// [[Rcpp::export]]
double legendre_deriv_bh(int l, int m, double x) {
  
  double pm1 = legendre_bh(l, m-1,  x);
  double pp1 = legendre_bh(l, m+1,  x);
  
  return(0.5 * ((l+m) * (l-m+1) * pm1 - pp1));
  
}


//==============================================================================
// @title
// scale_legendre_bh
//
// @description
// Scaling for legendre function
//
// @param l maximum degree
// @param m the x value
//
// @return legendre scale
//
// @keywords internal
//
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


//==============================================================================
// @title
// legendre
//
// @description
// Calculate legendre
//
// @param l_max maximum degree
// @param x the x value
//
// @return legendre polynomials
// 
// @keywords internal
//
// [[Rcpp::export]]
Rcpp::NumericMatrix legendre(int l_max, double x) {  
  
  double scale;
  
  int n = Rcpp::sum(Rcpp::seq(3, l_max + 1));
  Rcpp::NumericMatrix out(n, 4);
  
  
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



//==============================================================================
// @title
// et_analyze
//
// @description
// Calculate tidal potential for a single time and multiple waves.  Code adapted from ETERNA.
//
// @param astro vector astronomical parameters
// @param astro_der vector derivative of astronomical parameters
// @param k_mat matrix tidal catalog values
// @param pk vector of phases for each group
// @param body vector of body tide for each group
// @param body_inds indices of body max for each group
// @param delta vector body
// @param deltar double gravimentric factor
// @param x0 vector tidal catalog values
// @param y0 vector tidal catalog values
// @param x1 vector tidal catalog values
// @param y1 vector tidal catalog values
// @param x2 vector tidal catalog values
// @param y2 vector tidal catalog values
// @param j2000 double Julian date
// @param o1 double frequency
// @param resonance double frequency
// @param max_amp int index of the wave with maximum amplitude in wave group
// @param update_coef double constant for phase updating
//
// @return synthetic gravity
//
// @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
// 
// @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// processing package ETERNA 3.3. Bulletin d'Informations
// Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
// 
// @keywords internal
//
// [[Rcpp::export]]
arma::mat et_analyze(const arma::mat astro,
                     const arma::mat astro_der,
                     const arma::mat k_mat,
                     const arma::vec pk,
                     const arma::vec body,
                     const arma::uvec body_inds,
                     double delta,
                     double deltar,
                     const arma::vec x0, 
                     const arma::vec y0, 
                     const arma::vec x1, 
                     const arma::vec y1,
                     const arma::vec x2, 
                     const arma::vec y2,
                     const arma::vec j2000,
                     double o1,
                     double resonance,
                     const arma::uword max_amp,
                     double update_coef,
                     bool scale) {
  
  
  int nr = k_mat.n_rows;  // number of constituents
  int nt = astro.n_cols;  // number of times
  
  double j2000_sq;
  double to_rad = M_PI / 180.0; 
  //  double coef =  to_rad * 1.0 / 3600.0;
  
  arma::vec dc2(nr), dc3(nr);
  arma::vec cos_dc2(nr), sin_dc2(nr);
  arma::vec dummy(nr), cos_c(nr), sin_c(nr);
  arma::vec dtham(nr), dthph(nr);
  arma::vec fac = body;
  arma::vec fac_x(nr), fac_y(nr);
  arma::mat output(nt, 2);
  output.zeros();
  
  
  dc2 = arma::vectorise(k_mat * astro.col(0) + pk + 360.0);
  dc2 = mod(dc2, 360) * to_rad;
  dc3 = arma::vectorise(k_mat * astro_der.col(0));
  
  
  // normalize to max of group
  fac(body_inds) = delta + deltar * (dc3(body_inds) - o1) / (resonance - dc3(body_inds));
  fac = fac / fac(max_amp);
  
  j2000_sq = j2000[0] * j2000[0];
  
  fac_x = fac % (x0 + x1 * j2000[0] + x2 * j2000_sq);
  fac_y = fac % (y0 + y1 * j2000[0] + y2 * j2000_sq);
  
  dtham = arma::sqrt(fac_x % fac_x + fac_y % fac_y);
  dthph = dc2 - arma::atan2(fac_y, fac_x);
  
  // determine phase correction
  cos_dc2 = cos(dthph); //dc0
  sin_dc2 = sin(dthph); //ds0
  
  
  if (nt == 1) {
    
    if (scale) {
      output(0, 0) = arma::dot(dtham, cos_dc2) / arma::max(dtham);
      output(0, 1) = arma::dot(dtham, sin_dc2) / arma::max(dtham); 
    } else {
      output(0, 0) = arma::dot(dtham, cos_dc2);
      output(0, 1) = arma::dot(dtham, sin_dc2); 
    }
    
  } else {
    
    dc3 = dc3 * update_coef; // tidal frequencies in radian per hour. dthfr
    
    // speed enhancement but sacrifices precision
    cos_c = arma::cos(dc3); // ddc
    sin_c = arma::sin(dc3); // dds
    
    // loop through each time group
    for (int k = 0; k < nt; k++) {
      
      if (scale) {
        output(k, 0) = arma::dot(dtham, cos_dc2) / arma::max(dtham);
        output(k, 1) = arma::dot(dtham, sin_dc2) / arma::max(dtham); 
      } else {
        output(k, 0) = arma::dot(dtham, cos_dc2);
        output(k, 1) = arma::dot(dtham, sin_dc2); 
      }
      
      // update
      dummy   = cos_dc2 % cos_c - sin_dc2 % sin_c;
      sin_dc2 = sin_dc2 % cos_c + cos_dc2 % sin_c;  
      cos_dc2 = dummy;  
      
    }
  }
  
  
  return(output);
  
}



//==============================================================================
// @title
// et_predict
//
// @description
// Calculate tidal potential for a single time and multiple waves.  Code adapted from ETERNA.
//
// @param astro vector astronomical parameters
// @param astro_der vector derivative of astronomical parameters
// @param k_mat matrix tidal catalog values
// @param pk vector of phases for each group
// @param body vector of body tide for each group
// @param body_inds indices of body max for each group
// @param delta vector body
// @param deltar double gravimentric factor
// @param x0 vector tidal catalog values
// @param y0 vector tidal catalog values
// @param x1 vector tidal catalog values
// @param y1 vector tidal catalog values
// @param x2 vector tidal catalog values
// @param y2 vector tidal catalog values
// @param j2000 double Julian date
// @param o1 double frequency
// @param resonance double frequency
// @param max_amp int index of the wave with maximum amplitude in wave group
// @param update_coef double constant for phase updating
//
// @return synthetic gravity
//
// @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
// 
// @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// processing package ETERNA 3.3. Bulletin d'Informations
// Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
// 
// 
// @keywords internal
//
// [[Rcpp::export]]
arma::mat et_predict(const arma::mat astro,
                     const arma::mat astro_der,
                     const arma::mat k_mat,
                     const arma::vec pk,
                     const arma::vec body,
                     const arma::uvec body_inds,
                     double delta,
                     double deltar,
                     const arma::vec x0, 
                     const arma::vec y0, 
                     const arma::vec x1, 
                     const arma::vec y1,
                     const arma::vec x2, 
                     const arma::vec y2,
                     const arma::vec j2000,
                     double o1,
                     double resonance,
                     const arma::uword max_amp,
                     double update_coef) {
  
  
  int nr = k_mat.n_rows;  // number of constituents
  int nt = astro.n_cols;  // number of times
  
  double j2000_sq;
  double to_rad = M_PI / 180.0; 
  
  arma::vec dc2(nr), dc3(nr);
  arma::vec cos_dc2(nr), sin_dc2(nr);
  arma::vec dummy(nr), cos_c(nr), sin_c(nr);
  arma::vec fac = body;
  arma::mat output(nt, 1);
  
  
  dc2 = arma::vectorise(k_mat * astro.col(0) + pk + 360.0);
  dc2 = mod(dc2, 360) * to_rad;
  dc3 = arma::vectorise(k_mat * astro_der.col(0));
  
  
  // normalize to max of group
  fac(body_inds) = delta + deltar * (dc3(body_inds) - o1) / (resonance - dc3(body_inds));
  fac = fac / fac(max_amp);
  
  
  // determine phase correction
  cos_dc2 = arma::cos(dc2);
  sin_dc2 = arma::sin(dc2);
  
  
  
  if (nt == 1) {
    
    j2000_sq = j2000[0] * j2000[0];
    
    output(0, 0) = arma::dot(fac % (x0 + x1 * j2000[0] + x2 * j2000_sq), cos_dc2)+
      arma::dot(fac % (y0 + y1 * j2000[0] + y2 * j2000_sq), sin_dc2);
    
  } else {
    
    dc3 = dc3 * update_coef;
    
    // speed enhancement but sacrifices precision
    cos_c = arma::cos(dc3);
    sin_c = arma::sin(dc3);
    
    // loop through each time group
    for (int k = 0; k < nt; k++) {
      
      j2000_sq = j2000[k] * j2000[k];
      
      output(k, 0) = arma::dot(fac % (x0 + x1 * j2000[0] + x2 * j2000_sq), cos_dc2) +
        arma::dot(fac % (y0 + y1 * j2000[0] + y2 * j2000_sq), sin_dc2);
      
      // update
      dummy   = cos_dc2 % cos_c - sin_dc2 % sin_c;
      sin_dc2 = sin_dc2 % cos_c + cos_dc2 % sin_c;  
      cos_dc2 = dummy;  
      
    }
  }
  
  
  return(output);
  
}



//==============================================================================
// Main parallel structure for time and wave groups
struct earthtide_worker : public Worker
{
  const arma::mat astro;
  const arma::mat astro_der;
  const arma::mat k_mat;
  const arma::field<arma::vec> pk;
  const arma::field<arma::vec> body;
  const arma::field<arma::uvec> body_inds;
  double delta;
  double deltar;
  const arma::vec x0;
  const arma::vec y0;
  const arma::vec x1;
  const arma::vec y1;
  const arma::vec x2;
  const arma::vec y2;
  const arma::vec j2000;
  double o1;
  double resonance;
  const arma::field<arma::uvec> inds;
  const arma::uvec i_max;
  int astro_update;
  double update_coef;
  const arma::vec multiplier;
  bool predict;
  bool scale;
  arma::mat& output;
  
  earthtide_worker(const arma::mat astro,
                   const arma::mat astro_der,
                   const arma::mat k_mat,
                   const arma::field<arma::vec> pk,
                   const arma::field<arma::vec> body,
                   const arma::field<arma::uvec> body_inds,
                   double delta,
                   double deltar,
                   const arma::vec x0,
                   const arma::vec y0,
                   const arma::vec x1,
                   const arma::vec y1,
                   const arma::vec x2,
                   const arma::vec y2,
                   const arma::vec j2000,
                   double o1,
                   double resonance,
                   const arma::field<arma::uvec> inds,
                   const arma::uvec i_max,
                   int astro_update,
                   double update_coef,
                   const arma::vec multiplier,
                   bool predict,
                   bool scale,
                   arma::mat& output)
    : astro(astro), astro_der(astro_der), k_mat(k_mat), pk(pk), body(body),body_inds(body_inds), delta(delta), deltar(deltar),
      x0(x0),  y0(y0), x1(x1), y1(y1), x2(x2), y2(y2),
      j2000(j2000), o1(o1), resonance(resonance), inds(inds), i_max(i_max),
      astro_update(astro_update), update_coef(update_coef), multiplier(multiplier),
      predict(predict), scale(scale), output(output) {}
  
  
  void operator()(std::size_t begin_row, std::size_t end_row) {
    
    int start;
    int end;
    arma::uword r_start;
    arma::uword r_end;
    
    if (predict) {
      // loop through each time chunk
      for (std::size_t k = begin_row; k < end_row; k += astro_update) {
        
        start = k;
        end   = std::min(k + astro_update - 1, end_row-1);
        
        // loop through each time wave group
        for(std::size_t j = 0; j < inds.n_elem; j++) {
          
          // r_start = (j * 2);
          // r_end = (j * 2 + 1);
          
          output(arma::span(start, end), 0) += multiplier[j] * et_predict(
            astro.cols(arma::span(start, end)),
            astro_der.cols(arma::span(start, end)),
            k_mat.rows(inds(j)),
            pk(j),
            body(j),
            body_inds(j),
            delta,
            deltar,
            x0(inds(j)),
            y0(inds(j)),
            x1(inds(j)),
            y1(inds(j)),
            x2(inds(j)),
            y2(inds(j)),
            j2000(arma::span(start, end)),
            o1,
            resonance,
            i_max(j),
            update_coef);
        }
      }
    } else {
      
      // loop through each time chunk
      for (std::size_t k = begin_row; k < end_row; k += astro_update) {
        
        start = k;
        end   = std::min(k + astro_update - 1, end_row-1);
        
        // loop through each time wave group
        for(std::size_t j = 0; j < inds.n_elem; j++) {
          
          r_start = (j * 2);
          r_end = (j * 2 + 1);
          
          output(arma::span(start, end), arma::span(r_start, r_end)) = multiplier[j] * et_analyze(
            astro.cols(arma::span(start, end)),
            astro_der.cols(arma::span(start, end)),
            k_mat.rows(inds(j)),
            pk(j),
            body(j),
            body_inds(j),
            delta,
            deltar,
            x0(inds(j)),
            y0(inds(j)),
            x1(inds(j)),
            y1(inds(j)),
            x2(inds(j)),
            y2(inds(j)),
            j2000(arma::span(start, end)),
            o1,
            resonance,
            i_max(j),
            update_coef,
            scale);
        }
      }
    }
  }
};


//==============================================================================
// @title
// et_calculate
//
// @description
// Parallel calculation of tidal potential for a single time and multiple waves.  Code adapted from ETERNA.
//
// @param astro matrix astronomical parameters
// @param astro_der matrix derivative of astronomical parameters
// @param k_mat matrix tidal catalog values
// @param phases vector phases
// @param delta vector body
// @param deltar double gravimentric factor
// @param c0 vector tidal catalog values
// @param s0 vector tidal catalog values
// @param c1 vector tidal catalog values
// @param s1 vector tidal catalog values
// @param c2 vector tidal catalog values
// @param s2 vector tidal catalog values
// @param dgk vector geodetic coefficients
// @param jcof vector wave index
// @param j2000 double Julian date
// @param o1 double frequency
// @param resonance double frequency
// @param index vector wave group index
// @param astro_update how often to recalculate astronomical parameters
// @param update_coef time for approx
// @param magnifier scale the wave group value
// @param predict predict or analyze
//
// @return synthetic gravity
//
// @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
// 
// @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// processing package ETERNA 3.3. Bulletin d'Informations
// Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
// 
// @keywords internal
//
// [[Rcpp::export]]
arma::mat et_calculate(const arma::mat astro,
                       const arma::mat astro_der,
                       const arma::mat k_mat,
                       const arma::vec phases,
                       const arma::vec delta,
                       double deltar,
                       const arma::vec c0, 
                       const arma::vec s0, 
                       const arma::vec c1, 
                       const arma::vec s1,
                       const arma::vec c2, 
                       const arma::vec s2,
                       const arma::vec dgk,
                       const arma::uvec jcof,
                       const arma::vec j2000,
                       double o1,
                       double resonance,
                       const arma::ivec index,
                       int astro_update,
                       double update_coef,
                       const arma::vec magnifier,
                       bool predict,
                       bool scale) {
  
  
  // number of times
  int nt = astro.n_cols;
  
  
  // number of wave groups
  arma::ivec un = arma::unique(index);
  int ng = un.n_elem;
  
  
  // two columns for each wave group and one row for each time
  arma::mat output;
  
  if(predict) {
    output = arma::zeros(nt,1);
  } else {
    output = arma::zeros(nt, 2 * ng);
  };
  
  
  arma::uvec i_max(ng);                   // index of wave group max
  arma::field<arma::uvec> inds(ng);       // indices for each wave group
  arma::field<arma::uvec> body_inds(ng);  // indices for each wave group
  arma::field<arma::vec> body(ng);        // indices for each wave group
  arma::field<arma::vec> pk(ng);          // indices for each wave group
  
  arma::vec x0 = c0 % dgk(jcof) * 1.e-10;
  arma::vec y0 = s0 % dgk(jcof) * 1.e-10;
  arma::vec x1 = c1 % dgk(jcof) * 1.e-10;
  arma::vec y1 = s1 % dgk(jcof) * 1.e-10;
  arma::vec x2 = c2 % dgk(jcof) * 1.e-10;
  arma::vec y2 = s2 % dgk(jcof) * 1.e-10;
  
  
  // get the subsets for each wave group
  // get equilibrium wave amplitude
  arma::vec amplitude = arma::sqrt(pow(x0, 2) + pow(y0, 2));
  
  for(int j = 0; j < ng; j++) {
    arma::uvec sub = arma::find(index == j+1);
    inds(j)        = sub;
    pk(j)          = phases.elem(jcof.elem(sub));
    body(j)        = delta.elem(jcof.elem(sub));
    i_max(j)       = amplitude.elem(sub).index_max();
    body_inds(j)   = find((jcof.elem(sub) + 1) == 2);
  }
  
  
  earthtide_worker ew(astro,
                      astro_der,
                      k_mat,
                      pk,
                      body,
                      body_inds,
                      delta(1),
                      deltar,
                      x0,
                      y0,
                      x1,
                      y1,
                      x2,
                      y2,
                      j2000,
                      o1,
                      resonance,
                      inds,
                      i_max,
                      astro_update,
                      update_coef,
                      magnifier,
                      predict,
                      scale,
                      output);
  
  RcppParallel::parallelFor(0, nt, ew);
  
  
  return(output);
}



// //==============================================================================
// // The below prediction method typically 10-50% faster for cutoff > 1e-5
// //==============================================================================

// //==============================================================================
// //' @title
// //' et_predict_one
// //'
// //' @description
// //' Calculate tidal potential for a single time and multiple waves.  Code adapted from ETERNA.
// //'
// //' @param astro vector astronomical parameters
// //' @param astro_der vector derivative of astronomical parameters
// //' @param k_mat matrix tidal catalog values
// //' @param phases vector phases
// //' @param delta vector body
// //' @param deltar double gravimentric factor
// //' @param x0 vector tidal catalog values
// //' @param y0 vector tidal catalog values
// //' @param x1 vector tidal catalog values
// //' @param y1 vector tidal catalog values
// //' @param x2 vector tidal catalog values
// //' @param y2 vector tidal catalog values
// //' @param j2000 double Julian date
// //' @param o1 double frequency
// //' @param resonance double frequency
// //' @param max_amp int index of the wave with maximum amplitude in wave group
// //'
// //' @return synthetic gravity
// //'
// //' @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
// //' 
// //' @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// //' processing package ETERNA 3.3. Bulletin d'Informations
// //' Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
// //' 
// //' @export
// //'
// // [[Rcpp::export]]
// arma::vec et_predict_one(const arma::mat astro,
//                          const arma::mat astro_der,
//                          const arma::mat k_mat,
//                          const arma::vec pk,
//                          const arma::vec body,
//                          double delta,
//                          double deltar,
//                          const arma::vec x0, 
//                          const arma::vec y0, 
//                          const arma::vec x1, 
//                          const arma::vec y1,
//                          const arma::vec x2, 
//                          const arma::vec y2,
//                          const arma::vec j2000,
//                          double o1,
//                          double resonance,
//                          const arma::uvec indices,
//                          const arma::uvec max_amp,
//                          double update_coef) {
//   
//   
//   int nr = k_mat.n_rows;  // number of constituents
//   int nt = astro.n_cols;  // number of times
//   
//   double j2000_sq;
//   double to_rad = M_PI / 180.0; 
//   // double coef =  to_rad * 1.0 / 3600.0;
//   
//   arma::vec dc2(nr), dc3(nr);
//   arma::vec cos_dc2(nr), sin_dc2(nr);
//   arma::vec dummy(nr), cos_c(nr), sin_c(nr);
//   arma::vec fac = body;
//   arma::vec output(nt);
//   
//   
//   
//   dc2 = arma::vectorise(k_mat * astro.col(0) + pk + 360.0);
//   dc2 = mod(dc2, 360) * to_rad;
//   dc3 = arma::vectorise(k_mat * astro_der.col(0));
//   
//   
//   // normalize to max of group
//   fac(indices) = delta + deltar * (dc3(indices) - o1) / (resonance - dc3(indices));
//   fac = fac / fac(max_amp);
//   
//   
//   // determine phase correction
//   cos_dc2 = arma::cos(dc2);
//   sin_dc2 = arma::sin(dc2);
//   
//   
//   if (nt == 1) {
//     
//     j2000_sq = j2000[0] * j2000[0];
//     
//     output(0) = arma::dot(fac % (x0 + x1 * j2000[0] + x2 * j2000_sq), cos_dc2) +
//                 arma::dot(fac % (y0 + y1 * j2000[0] + y2 * j2000_sq), sin_dc2);
//     
//   } else {
//     
//     dc3 = dc3 * update_coef;
//     
//     // speed enhancement but sacrifices precision
//     cos_c = arma::cos(dc3);
//     sin_c = arma::sin(dc3);
//     
//     // loop through each time group
//     for (int k = 0; k < nt; k++) {
//       
//       j2000_sq = j2000[k] * j2000[k];
//       
//       output(k) = arma::dot(fac % (x0 + x1 * j2000[k] + x2 * j2000_sq), cos_dc2) +
//         arma::dot(fac % (y0 + y1 * j2000[k] + y2 * j2000_sq), sin_dc2);
//       
//       // update
//       dummy   = cos_dc2 % cos_c - sin_dc2 % sin_c;
//       sin_dc2 = sin_dc2 % cos_c + cos_dc2 % sin_c;  
//       cos_dc2 = dummy;  
//       
//     }
//   }
//   
//   
//   return(output);
//   
// }
// 
// 
// 
// 
// //==============================================================================
// // Main parallel structure for time and wave groups
// struct earthtide_predict_worker : public Worker
// {
//   const arma::mat astro;
//   const arma::mat astro_der;
//   const arma::mat k_mat;
//   const arma::vec pk;
//   const arma::vec body;
//   double delta;
//   double deltar;
//   const arma::vec x0;
//   const arma::vec y0;
//   const arma::vec x1;
//   const arma::vec y1;
//   const arma::vec x2;
//   const arma::vec y2;
//   const arma::vec j2000;
//   double o1;
//   double resonance;
//   const arma::uvec indices;
//   const arma::uvec i_max;
//   int astro_update;
//   double update_coef;
//   arma::vec& output;
//   
//   earthtide_predict_worker(const arma::mat astro,
//                         const arma::mat astro_der,
//                         const arma::mat k_mat,
//                         const arma::vec pk,
//                         const arma::vec body,
//                         double delta,
//                         double deltar,
//                         const arma::vec x0,
//                         const arma::vec y0,
//                         const arma::vec x1,
//                         const arma::vec y1,
//                         const arma::vec x2,
//                         const arma::vec y2,
//                         const arma::vec j2000,
//                         double o1,
//                         double resonance,
//                         const arma::uvec indices,
//                         const arma::uvec i_max,
//                         int astro_update,
//                         double update_coef,
//                         arma::vec& output)
//     : astro(astro), astro_der(astro_der), k_mat(k_mat), pk(pk), body(body), delta(delta), deltar(deltar),
//       x0(x0),  y0(y0), x1(x1), y1(y1), x2(x2), y2(y2),
//       j2000(j2000), o1(o1), resonance(resonance), indices(indices), i_max(i_max),
//       astro_update(astro_update),update_coef(update_coef), output(output) {}
//   
//   
//   void operator()(std::size_t begin_row, std::size_t end_row) {
//     
//     int start;
//     int end;
//     
//     // loop through each time chunk
//     for (std::size_t k = begin_row; k < end_row; k += astro_update) {
//       
//       start = k;
//       end   = std::min(k + astro_update - 1, end_row-1);
//       
//       output(arma::span(start, end)) = et_predict_one(
//         astro.cols(arma::span(start, end)),
//         astro_der.cols(arma::span(start, end)),
//         k_mat,
//         pk,
//         body,
//         delta,
//         deltar,
//         x0,
//         y0,
//         x1,
//         y1,
//         x2,
//         y2,
//         j2000(arma::span(start, end)),
//         o1,
//         resonance,
//         indices,
//         i_max,
//         update_coef);
//     }
//     
//   }
// };
// 
// 
// //==============================================================================
// //' @title
// //' et_predict
// //'
// //' @description
// //' Parallel calculation of tidal potential for a single time and multiple waves.  Code adapted from ETERNA.
// //'
// //' @param astro matrix astronomical parameters
// //' @param astro_der matrix derivative of astronomical parameters
// //' @param k_mat matrix tidal catalog values
// //' @param pk vector phases
// //' @param delta vector body
// //' @param deltar double gravimentric factor
// //' @param c0 vector tidal catalog values
// //' @param s0 vector tidal catalog values
// //' @param c1 vector tidal catalog values
// //' @param s1 vector tidal catalog values
// //' @param c2 vector tidal catalog values
// //' @param s2 vector tidal catalog values
// //' @param dgk vector geodetic coefficients
// //' @param jcof vector wave index
// //' @param j2000 double Julian date
// //' @param o1 double frequency
// //' @param resonance double frequency
// //' @param index vector wave group index
// //' @param astro_update how often to recalculate astronomical parameters
// //'
// //' @return synthetic gravity
// //'
// //' @author Jonathan Kennel, \email{jkennel@uoguelph.ca}
// //' 
// //' @references Wenzel, H.-G. (1996): The nanogal software: Earth tide data
// //' processing package ETERNA 3.3. Bulletin d'Informations
// //' Marees Terrestres vol. 124, 9425-9439, Bruxelles 1996.
// //' 
// //' @export
// //'
// // [[Rcpp::export]]
// arma::vec et_predict(const arma::mat astro,
//                      const arma::mat astro_der,
//                      const arma::mat k_mat,
//                      const arma::vec phases,
//                      const arma::vec delta,
//                      double deltar,
//                      const arma::vec c0, 
//                      const arma::vec s0, 
//                      const arma::vec c1, 
//                      const arma::vec s1,
//                      const arma::vec c2, 
//                      const arma::vec s2,
//                      const arma::vec dgk,
//                      const arma::uvec jcof,
//                      const arma::vec j2000,
//                      double o1,
//                      double resonance,
//                      const arma::ivec index,
//                      int astro_update,
//                      double update_coef) {
//   
//   
//   // number of times
//   int nt = astro.n_cols;
//   
//   
//   // number of wave groups
//   arma::ivec un = arma::unique(index);
//   int ng = un.n_elem;
//   
//   
//   // two columns for each wave group and one row for each time
//   arma::vec output(nt);
//   output.fill(0.0);
//   
//   
//   arma::uvec i_max(index.size());    // index of wave group max
//   
//   arma::vec x0 = c0 % dgk(jcof) * 1.e-10;
//   arma::vec y0 = s0 % dgk(jcof) * 1.e-10;
//   arma::vec x1 = c1 % dgk(jcof) * 1.e-10;
//   arma::vec y1 = s1 % dgk(jcof) * 1.e-10;
//   arma::vec x2 = c2 % dgk(jcof) * 1.e-10;
//   arma::vec y2 = s2 % dgk(jcof) * 1.e-10;
//   
//   
//   // get the subsets for each wave group
//   // get equilibrium wave amplitude
//   arma::vec amplitude = arma::sqrt(pow(x0, 2) + pow(y0, 2));
//   
//   int n_count = 0;
//   for(int j = 0; j < ng; j++) {
//     arma::uvec sub = arma::find(index == j+1);
//     int n_sub = sub.size();
//     arma::uvec max_ind(n_sub);
//     max_ind.fill(n_count + amplitude.elem(sub).index_max());
//     i_max.elem(sub) = max_ind; // need to remove empty subsets
//     n_count += n_sub;
//   }
//   
//   
//   arma::vec pk = phases(jcof);
//   arma::vec body = delta(jcof);
//   arma::uvec indices = find((jcof + 1) == 2);
//   
//   earthtide_predict_worker ew(astro,
//                            astro_der,
//                            k_mat,
//                            pk,
//                            body,
//                            delta(1),
//                            deltar,
//                            x0,
//                            y0,
//                            x1,
//                            y1,
//                            x2,
//                            y2,
//                            j2000,
//                            o1,
//                            resonance,
//                            indices,
//                            i_max,
//                            astro_update,
//                            update_coef,
//                            output);
//   
//   RcppParallel::parallelFor(0, nt, ew);
//   
//   
//   return(output);
// }



//==============================================================================
/*** R
*/
