
#include "r2d.h"

#include <stdint.h>
#include <stdbool.h>


/** \struct r2d_parabola
 *  \brief A parabola: n ⋅ (x - x0) + (kappa0/2) * (t ⋅ (x - x0))^2 = 0.
 */
typedef struct {
	r2d_rvec2 n; /*!< Unit-length normal vector. */
	r2d_real kappa0;
} r2d_parabola;

/**
 * \brief Clip a polygon against a parabola and compute the zeroth and first moment.
 *
 * \param [in, out] poly
 * The polygon to be clipped. Upon return the poly represents the intersected polygon, except that parabolic
 * faces are approximated by a line
 *
 * \param [in] parabola
 * An array of planes against which to clip this polygon.
 * 
 * \param [in, out] grad_s
 * The derivative of the plane constant w.r.t. the normal angle phi
 * 	-	if not given (equal to NaN) then it is computed such as to enforce volume conservation (hence derivative[0] = 0.)
 * 	- if given (not equal to NaN) then it is used to compute the derivative of the zeroth moment
 *
 * \param[out] moments
 * The zeroth and first moment of the clipped polygon (the moments are exact, in contrary to the
 * approximate output poly)
 *
 * \param[out] derivative
 * The derivative of the zeroth and first moment (w.r.t. the normal angle) of the clipped polygon (derivative[0:2])
 * and the derivative w.r.t. the curvature (derivative[3])
 */
void r2d_clip_parabola_cmpMoments(r2d_poly* poly, r2d_parabola *parabola, r2d_real *grad_s, r2d_real *moments, r2d_real *derivative, bool *compute_derivative);

r2d_int parabola_line_intersection(r2d_parabola *parabola, r2d_rvec2 pos1, r2d_rvec2 pos2, r2d_real *roots);

void adjust_moments_for_parabola(r2d_parabola *parabola, r2d_rvec2 pos1, r2d_rvec2 pos2, r2d_real *moments);

void compute_moment_derivatives(r2d_parabola *parabola, r2d_vertex *vertbuffer, r2d_int *nverts, bool *is_on_parabola, r2d_real *grad_s, r2d_real *derivative, bool compute_grad_s);

void real_roots(r2d_real *coeff, r2d_real *roots);

void r2d_split_ptr(r2d_poly* inpolys, r2d_int npolys, r2d_plane *plane, r2d_poly* out_pos, r2d_poly* out_neg);