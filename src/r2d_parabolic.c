#include "r2d_parabolic.h"

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdbool.h>

// useful macros
#define NR_MONOMIALS 6
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)
#define dot_relative(va, vb, vr) ((va.x - vr.x)*vb.x + (va.y - vr.y)*vb.y)
#define dot_relative_rotate(va, vb, vr) ((va.y - vr.y)*vb.x - (va.x - vr.x)*vb.y)
#define wav(va, wa, vb, wb, vr) {			\
	vr.x = ((wa)*va.x + (wb)*vb.x)/(wa + wb);	\
	vr.y = ((wa)*va.y + (wb)*vb.y)/(wa + wb);	\
}

// Clip poly with parabola, and compute zeroth and first moment
void r2d_clip_parabola_moments_01(r2d_poly* poly, r2d_parabola *parabola, r2d_real *grad_s, r2d_real *moments, r2d_real *derivative, bool *compute_derivative) {

	bool compute_grad_s;

	moments[0] = 0.0;
	moments[1] = 0.0;
	moments[2] = 0.0;

	// The derivative of the zeroth and first moment w.r.t. the normal angle
	derivative[0] = 0.0;
	derivative[1] = 0.0;
	derivative[2] = 0.0;

	compute_grad_s = isnan(grad_s[0]);
	if (compute_grad_s) {
		grad_s[0] = 0.;
		grad_s[1] = 0.;
	}

	if (parabola->kappa0 > 0.0) {
		printf("Kappa0 can not be positive\n");
		return;
	}

	// variable declarations
	r2d_int v, p, np, old_nverts, vstart, vcur, vnext, new_nverts;

	// direct access to vertex buffer
	r2d_vertex* vertbuffer = poly->verts;
	r2d_int* nverts = &poly->nverts;
	if(*nverts <= 0) return;

	// signed distance to the clipping plane
	r2d_real sdist;

	// for marking clipped vertices
	r2d_int clipped[R2D_MAX_VERTS];

	// for marking parabola vertices
	bool is_on_parabola[R2D_MAX_VERTS];
	r2d_real roots[2], root, tau;
	r2d_int nr_roots, sgn;

	// calculate signed distances to the clip plane
	old_nverts = *nverts;
	memset(&clipped, 0, sizeof(clipped));
	for(v = 0; v < old_nverts; ++v) {
		// NB this does no longer represent a distance (unless kapp0 = 0), but the
		// sign still indicates on which side of the parabola we are on
		tau = dot_relative_rotate(vertbuffer[v].pos, parabola->n, parabola->x0);
		sdist = dot_relative(vertbuffer[v].pos, parabola->n, parabola->x0) + (parabola->kappa0 / 2.0) * tau * tau;
		if(sdist >= 0.0) clipped[v] = 1;
	}

	memset(&is_on_parabola, false, sizeof(is_on_parabola));
	// check all edges and insert new vertices on the edges
	for(vcur = 0; vcur < old_nverts; ++vcur) {
		if(clipped[vcur]) continue;
		for(np = 0; np < 2; ++np) {
			vnext = vertbuffer[vcur].pnbrs[np];
			if(!clipped[vnext]) continue;
			// Case 1) an edge has one node clipped: bisected edge
			nr_roots = parabola_line_intersection(parabola, vertbuffer[vcur].pos, vertbuffer[vnext].pos, roots);
			if (nr_roots==0) continue;

			// link new vertex with unclipped vertex (negative value of clipped vertex may be used for finding neighbors)
			root = roots[0];
			vertbuffer[*nverts].pnbrs[1-np] = vcur;
			vertbuffer[*nverts].pnbrs[np] = -1;
			vertbuffer[vcur].pnbrs[np] = *nverts;
			vertbuffer[vnext].pnbrs[1-np] = *nverts;
			wav(vertbuffer[vcur].pos, 1.0 - root,
				  vertbuffer[vnext].pos, root,
				  vertbuffer[*nverts].pos);
			is_on_parabola[*nverts] = true;
			// printf("Case 1: pos1 = (%e, %e)\n", vertbuffer[vcur].pos.x, vertbuffer[vcur].pos.y);
			// printf("        pos2 = (%e, %e)\n", vertbuffer[vnext].pos.x, vertbuffer[vnext].pos.y);
			// printf("        root, new pos = %e, (%e, %e)\n", root, vertbuffer[*nverts].pos.x, vertbuffer[*nverts].pos.y);

			(*nverts)++;

		}
	}

	for(vcur = 0; vcur < old_nverts; ++vcur) {
		vnext = vertbuffer[vcur].pnbrs[0];
		if (vnext >= old_nverts) continue;
		if (!clipped[vcur] && !clipped[vnext] && parabola->kappa0 < 0.0) {
			// Case 3) an edge has neither vertex clipped: possibly trisected edge (kappa0 < 0)
			nr_roots = parabola_line_intersection(parabola, vertbuffer[vcur].pos, vertbuffer[vnext].pos, roots);
			if (nr_roots == 2) {
				// link new vertices with unclipped vertices (negative value of clipped vertex is not used for finding neighbors)
				vertbuffer[*nverts].pnbrs[1] = vcur;
				vertbuffer[*nverts].pnbrs[0] = -1;
				vertbuffer[vcur].pnbrs[0] = *nverts;
				wav(vertbuffer[vcur].pos, 1.0 - roots[0],
					  vertbuffer[vnext].pos, roots[0],
					  vertbuffer[*nverts].pos);
				is_on_parabola[*nverts] = true;
				// printf("Case 3: pos1 = (%e, %e)\n", vertbuffer[vcur].pos.x, vertbuffer[vcur].pos.y);
				// printf("        pos2 = (%e, %e)\n", vertbuffer[vnext].pos.x, vertbuffer[vnext].pos.y);
				// printf("        root, new pos = %e, (%e, %e)\n", roots[0], vertbuffer[*nverts].pos.x, vertbuffer[*nverts].pos.y);

				vertbuffer[(*nverts)+1].pnbrs[0] = vnext;
				vertbuffer[(*nverts)+1].pnbrs[1] = -1;
				vertbuffer[vnext].pnbrs[1] = (*nverts)+1;
				wav(vertbuffer[vcur].pos, 1.0 - roots[1],
					  vertbuffer[vnext].pos, roots[1],
					  vertbuffer[(*nverts)+1].pos);
				is_on_parabola[(*nverts)+1] = true;
				// printf("        root, new pos = %e, (%e, %e)\n", roots[1], vertbuffer[(*nverts)+1].pos.x, vertbuffer[(*nverts)+1].pos.y);

				(*nverts) += 2;
			}
		}
	}

	// for each new vertex, search around the poly for its new neighbors
	// and doubly-link everything
	for(vstart = old_nverts; vstart < *nverts; ++vstart) {
		if(vertbuffer[vstart].pnbrs[1] >= 0) continue;
		vcur = vertbuffer[vstart].pnbrs[0];
		do {
			// printf("vcur = %d\n", vcur);
			vcur = vertbuffer[vcur].pnbrs[0];
		} while(vcur < old_nverts);
		// printf("doubly link (%d,%d)\n", vcur, vstart);
		vertbuffer[vstart].pnbrs[1] = vcur;
		vertbuffer[vcur].pnbrs[0] = vstart;
	}

	// go through and compress the vertex list, removing clipped verts
	// and re-indexing accordingly (reusing `clipped` to re-index everything)
	new_nverts = 0;
	for(v = 0; v < *nverts; ++v) {
		if(!clipped[v]) {
			vertbuffer[new_nverts] = vertbuffer[v];
			is_on_parabola[new_nverts] = is_on_parabola[v];
			clipped[v] = new_nverts++;
		}
	}
	*nverts = new_nverts;
	for(v = 0; v < *nverts; ++v) {
		vertbuffer[v].pnbrs[0] = clipped[vertbuffer[v].pnbrs[0]];
		vertbuffer[v].pnbrs[1] = clipped[vertbuffer[v].pnbrs[1]];
	}


	// for(r2d_int v = 0; v < poly->nverts; ++v) {
	// 	printf("  vertex %d: pos = ( %.10e , %.10e ), nbrs = (%d, %d), is_on_parabola = %d\n",
	// 			v, poly->verts[v].pos.x, poly->verts[v].pos.y, poly->verts[v].pnbrs[0], poly->verts[v].pnbrs[1], is_on_parabola[v]);
	// }

	// first we compute the moments of the polygon
	r2d_reduce(poly, moments, 1);

	// then we adjust the moments with the parabolic contributions
	if (parabola->kappa0 < 0.0){
		for(vcur = 0; vcur < *nverts; ++vcur) {
			vnext = vertbuffer[vcur].pnbrs[0];
			if (is_on_parabola[vcur] && is_on_parabola[vnext]) {
				adjust_moments_for_parabola(parabola, vertbuffer[vcur].pos, vertbuffer[vnext].pos, moments);
			}
		}
	}

	if (*compute_derivative) compute_moment_derivatives(parabola, vertbuffer, nverts, is_on_parabola, grad_s, derivative, compute_grad_s);
}

void compute_moment_derivatives(r2d_parabola *parabola, r2d_vertex *vertbuffer, r2d_int *nverts, bool *is_on_parabola, r2d_real *grad_s, r2d_real *derivative, bool compute_grad_s){
	r2d_real tau[2], kappa0, shift;

	r2d_real tauPower[2], monomial, monomials_sum[NR_MONOMIALS], der1_eta, der1_tau;

	r2d_int vcur, vnext;

	if (*nverts==0) return;

	// TODO for efficiency we could compute the monomials_sum once and re-use it for
	// computation of the moments

	kappa0 = parabola->kappa0;
	shift = dot(parabola->x0, parabola->n);

	for (int i = 0; i < NR_MONOMIALS; ++i){
		monomials_sum[i] = 0.0;
	}

	// The derivative is composed of integrals over the interface
	for(vcur = 0; vcur < *nverts; ++vcur) {
		vnext = vertbuffer[vcur].pnbrs[0];

		if (is_on_parabola[vcur] && is_on_parabola[vnext]) {
			tau[0] = dot_relative_rotate(vertbuffer[vcur].pos, parabola->n, parabola->x0);
			tau[1] = dot_relative_rotate(vertbuffer[vnext].pos, parabola->n, parabola->x0);

			// each of the intergrals is over the total interface (which may be split into parts)
			// hence we compute the integrals
			// 	int τ^i dτ,
			// for i = 0, 5
			tauPower[0] = 1.0;
			tauPower[1] = 1.0;
			for (int i = 0; i < NR_MONOMIALS; ++i){
				tauPower[0] *= tau[0];
				tauPower[1] *= tau[1];
				monomial = (tauPower[1] - tauPower[0]) / (i + 1);
				monomials_sum[i] = monomials_sum[i] + monomial;
			}
		}
	}

	if (monomials_sum[0] == 0.0) return;

	// grad_s is the gradient of the shift (w.r.t. the normal angle and curvature) and is given by
	// 	-int [(s - κ/2 τ^2) κ τ - τ] dτ / int 1 dτ
	// 	int [τ^2/2] dτ / int 1 dτ
	// respectively
	if (compute_grad_s){
		grad_s[0] = -(monomials_sum[1] * shift * kappa0 - monomials_sum[3] * kappa0 * kappa0 / 2 - monomials_sum[1]) / monomials_sum[0];
		grad_s[1] = (monomials_sum[2]/2) / monomials_sum[0];

	}
	else
	{
		derivative[0] = monomials_sum[0] * grad_s[0] - monomials_sum[1] + monomials_sum[1] * shift * kappa0 - monomials_sum[3] * kappa0 * kappa0 / 2;
		derivative[3] = monomials_sum[0] * grad_s[1] - monomials_sum[2]/2;
	}

	// we write the derivative in terms of its normal and tangential components
	// the normal component is given by
	// 	int [grad_s - τ + κ τ(s - κ/2 τ^2)][s - κ/2 τ^2] dτ
	der1_eta = shift * grad_s[0] * monomials_sum[0]  + monomials_sum[3] * kappa0 / 2 - monomials_sum[2] * grad_s[0] * kappa0 / 2 - shift * monomials_sum[1] + kappa0 * (monomials_sum[1] * shift * shift + monomials_sum[5] * kappa0 * kappa0 / 4 - shift * kappa0 * monomials_sum[3]);

	// and the tangential component
	// 	int [grad_s - τ + κ τ(s - κ/2 τ^2)]τ dτ
	der1_tau = grad_s[0] * monomials_sum[1] - monomials_sum[2] + kappa0 * shift * monomials_sum[2] - monomials_sum[4] * kappa0 * kappa0 / 2;
	
	derivative[1] = parabola->n.x * der1_eta - parabola->n.y * der1_tau;
	derivative[2] = parabola->n.y * der1_eta + parabola->n.x * der1_tau;
}

void adjust_moments_for_parabola(r2d_parabola *parabola, r2d_rvec2 pos1, r2d_rvec2 pos2, r2d_real *moments){
	// pos in coordinates of the parabola: tangential, normal
	r2d_real eta[2], tau[2], kappa0, c[2], volume_correction, first_moment_correction[2], tauPower[2];

	r2d_real monomials[5];

	kappa0 = parabola->kappa0;
	tau[0] = dot_relative_rotate(pos1, parabola->n, parabola->x0);
	tau[1] = dot_relative_rotate(pos2, parabola->n, parabola->x0);

	if (tau[0] == tau[1]) return;

	eta[0] = dot_relative(pos1, parabola->n, parabola->x0);
	eta[1] = dot_relative(pos2, parabola->n, parabola->x0);

	// in the local coordinates the polygon face is given by
	// x_η = c_0 * x_τ + c_1
	c[0] = (eta[1] - eta[0]) / (tau[1] - tau[0]);
	c[1] = (eta[1] + eta[0]) / 2 - c[0] * (tau[0] + tau[1]) / 2;

	// and the parabola is given by
	// x_η = -(κ_0/2) x_τ^2

	tauPower[0] = 1.0;
	tauPower[1] = 1.0;
	for (int i = 0; i < 5; ++i){
		tauPower[0] *= tau[0];
		tauPower[1] *= tau[1];
		monomials[i] = (tauPower[1] - tauPower[0]) / (i + 1);
	}

	volume_correction = -(kappa0/2) * monomials[2] - (c[0] * monomials[1] + c[1] * monomials[0]);
	moments[0] += volume_correction;
	// printf("Parabola volume correction = %e for tau1, tau2 = (%e, %e)\n", volume_correction, tau[0], tau[1]);

	first_moment_correction[0] = -(kappa0/2) * monomials[3] - (c[0] * monomials[2] + c[1] * monomials[1]);
	first_moment_correction[1] = (kappa0*kappa0/8) * monomials[4] - (c[0] * c[0] * monomials[2] + 2 * c[0] * c[1] * monomials[1] + c[1] * c[1] * monomials[0]) / 2;
	// printf("Parabola first moment correction = (%e, %e)\n", first_moment_correction[0], first_moment_correction[1]);

	// transform the first moment correction back to global coordinates, and increment first moment
	moments[1] += parabola->n.x * first_moment_correction[1] - parabola->n.y * first_moment_correction[0] + parabola->x0.x * volume_correction;
	moments[2] += parabola->n.y * first_moment_correction[1] + parabola->n.x * first_moment_correction[0] + parabola->x0.y * volume_correction;
}

// Given a parabola and a line connecting the points pos1, pos2; find the intersection
r2d_int parabola_line_intersection(r2d_parabola *parabola, r2d_rvec2 pos1, r2d_rvec2 pos2, r2d_real *roots) {
	r2d_real coeff[3];
	r2d_int nr_roots, root_is_good[2];

	// We parametrize the line as: l(t) = pos1 + t * (pos2 - pos1),
	// and solve for t
	coeff[0] = (parabola->kappa0/2) * pow(dot_relative_rotate(pos2, parabola->n, pos1), 2);
	coeff[1] = parabola->kappa0 * dot_relative_rotate(pos1, parabola->n, parabola->x0) * dot_relative_rotate(pos2, parabola->n, pos1) + dot_relative(pos2, parabola->n, pos1);
	coeff[2] = dot_relative(pos1, parabola->n, parabola->x0) + (parabola->kappa0/2) * pow(dot_relative_rotate(pos1, parabola->n, parabola->x0), 2);

	real_roots(coeff, roots);

	nr_roots = 0;
	for (int i = 0; i < 2; ++i){
		root_is_good[i] = !isnan(roots[i]) && roots[i] >= 0.0 && roots[i] <= 1.0;
		if (root_is_good[i]) ++nr_roots;
	}

	if (root_is_good[1] && !root_is_good[0]) roots[0] = roots[1];

	return nr_roots;
}

// Finds real solution of c[0] * t^2 + c[1] * t + c[2] = 0
void real_roots(r2d_real *coeff, r2d_real *roots){
	r2d_real det, scaled[2], maxSqrt, tmp;

	if (coeff[0] == 0.0){
		if (coeff[1] == 0.0){
			roots[0] = NAN;
			roots[1] = NAN;
		}
		else{
			roots[0] = -coeff[2] / coeff[1];
			roots[1] = NAN;
		}
	}
	else if (coeff[2] == 0.0){
		roots[0] = -coeff[1] / coeff[0];
		roots[1] = 0;
	}
	else {
		// Citarduaq's formula is used to allow for small c[0] in a numerically stable way
		scaled[0] = coeff[1] / coeff[0];
		scaled[1] = coeff[2] / coeff[0];

		maxSqrt = sqrt(DBL_MAX);
		if (scaled[0] > maxSqrt || scaled[0] < -maxSqrt) {
			// scaled[0]^2 would overflow, so we let √D ≈ |scaled[0]| instead
			roots[0] = -scaled[1] / scaled[0];
			roots[1] = scaled[1] / roots[0];
		}
		else {
			det = scaled[0] * scaled[0] - 4 * scaled[1];
		}

		if (det < 0.0) {
			roots[0] = NAN;
			roots[1] = NAN;
		}
		else if (det == 0.0) {
			roots[0] = -scaled[1] / 2;
			roots[1] = NAN;
		}
		else {
			roots[0] = -2 * scaled[1] / (scaled[0] + (scaled[0] > 0 ? sqrt(det) : -sqrt(det)));
			roots[1] = scaled[1] / roots[0];
		}
	}

	if (!isnan(roots[0]) && !isnan(roots[1])){
		if (roots[0] > roots[1]) {
			tmp = roots[0];
			roots[0] = roots[1];
			roots[1] = tmp;
		}
	}
}