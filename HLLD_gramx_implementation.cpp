#ifndef GRHYDROX_HLLE_HXX
#define GRHYDROX_HLLE_HXX

#include "GRHydroX_var_groups.hxx"

#include <cctk.h>

#include <algorithm>
#include <cmath>

namespace GRHydroX
{

	// Build up multiple element MIN/MAX from pairwise recursive
	// comparison. These are ~10 times faster than the MIN(MIN(MIN(...
	// approach, presumably since they allow for out-of-order evalutation
	// of intermediate results.
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T min3(T a, T b, T c)
	{
		using std::min;
		return min(min(a, b), c);
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T min4(T a, T b, T c, T d)
	{
		using std::min;
		return min(min(a, b), min(c, d));
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T min8(T a, T b, T c, T d, T e, T f, T g, T h)
	{
		using std::min;
		return min(min4(a, b, c, d), min4(e, f, g, h));
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T min10(T a, T b, T c, T d, T e, T f, T g, T h, T i,
																   T j)
	{
		using std::min;
		return min(min8(a, b, c, d, e, f, g, h), min(i, j));
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T min11(T a, T b, T c, T d, T e, T f, T g, T h, T i, T j,
																   T k)
	{
		using std::min;
		return min(min8(a, b, c, d, e, f, g, h), min3(i, j, k));
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T min10_0(const T (&restrict a)[10])
	{
		return min11(T(0), a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8],
					 a[9]);
	}

	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T max3(T a, T b, T c)
	{
		using std::max;
		return max(max(a, b), c);
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T max4(T a, T b, T c, T d)
	{
		using std::max;
		return max(max(a, b), max(c, d));
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T max8(T a, T b, T c, T d, T e, T f, T g, T h)
	{
		using std::max;
		return max(max4(a, b, c, d), max4(e, f, g, h));
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T max10(T a, T b, T c, T d, T e, T f, T g, T h, T i,
																   T j)
	{
		using std::max;
		return max(max8(a, b, c, d, e, f, g, h), max(i, j));
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T max11(T a, T b, T c, T d, T e, T f, T g, T h, T i, T j,
																   T k)
	{
		using std::max;
		return max(max8(a, b, c, d, e, f, g, h), max3(i, j, k));
	}
	template <typename T>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T max10_0(const T (&restrict a)[10])
	{
		return max11(T(0), a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8],
					 a[9]);
	}

	// Used when we divide by a quantity such that we do not get any 0 divison errors
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
	prevent_zero_division(CCTK_REAL value)
	{
		const CCTK_REAL min_value = 1e-10;
		if (value >= 0)
		{
			return (value < min_value) ? min_value : value;
		}
		else
		{
			return (value > -min_value) ? -min_value : value;
		}
	}

	// calculation of additional values needed by flux_x() function
	// Total flops = 119
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	vlow_blow(const prim_point &P, const metric_point &metric, CCTK_REAL &ab0,
			  CCTK_REAL &b2, CCTK_REAL &w, CCTK_REAL &bxlow, CCTK_REAL &bylow,
			  CCTK_REAL &bzlow)
	{

		// Calculates v_i (see Anton Eq. 5) and B_i (Bvecxlow)- undensitized!
		// calculates B^i v_i [Anton eq. 44] and b^2 [LHS of Anton eq. 46]
		// Calculates w (Lorentz factor) as (1-v^i v_i)^{-1/2}
		// Calculates b_i (bxlow)

		// vel_i  = g_ij v^j
		// B_i = g_ij B^i
		const CCTK_REAL gxx = metric.gxx;
		const CCTK_REAL gxy = metric.gxy;
		const CCTK_REAL gxz = metric.gxz;
		const CCTK_REAL gyy = metric.gyy;
		const CCTK_REAL gyz = metric.gyz;
		const CCTK_REAL gzz = metric.gzz;
		const CCTK_REAL velx = P.velx;
		const CCTK_REAL vely = P.vely;
		const CCTK_REAL velz = P.velz;
		const CCTK_REAL Bvecx = P.Bvecx;
		const CCTK_REAL Bvecy = P.Bvecy;
		const CCTK_REAL Bvecz = P.Bvecz;

		const CCTK_REAL velxlow = gxx * velx + gxy * vely + gxz * velz;
		const CCTK_REAL velylow = gxy * velx + gyy * vely + gyz * velz;
		const CCTK_REAL velzlow = gxz * velx + gyz * vely + gzz * velz;
		const CCTK_REAL Bvecxlow = gxx * Bvecx + gxy * Bvecy + gxz * Bvecz;
		const CCTK_REAL Bvecylow = gxy * Bvecx + gyy * Bvecy + gyz * Bvecz;
		const CCTK_REAL Bveczlow =
			gxz * Bvecx + gyz * Bvecy + gzz * Bvecz; // flops = 30

		// B^i v_i (= b^0/u^0)
		const CCTK_REAL Bdotv = velxlow * Bvecx + velylow * Bvecy + velzlow * Bvecz;

		// v^2 = v_i v^i; w=(1-v^2)^{-1/2}
		const CCTK_REAL v2 = velxlow * velx + velylow * vely + velzlow * velz;
		w = 1. / sqrt(1. - v2); // flops = 31

		// b^2 = B^i B_i / w^2 + (b^0/u^0)^2
		b2 = (Bvecx * Bvecxlow + Bvecy * Bvecylow + Bvecz * Bveczlow) / SQR(w) +
			 SQR(Bdotv); // flops = 18

		// b_i = B_i/w +w*(B dot v)*v_i
		bxlow = Bvecxlow / w + w * Bdotv * velxlow;
		bylow = Bvecylow / w + w * Bdotv * velylow;
		bzlow = Bveczlow / w + w * Bdotv * velzlow;

		// alpha b_0 = w_lorentz B^i vel_i
		ab0 = w * Bdotv; // flops = 40
	}

	/***************************************************************************************************
	*************************************** X-DIRECTION FLUXES  ****************************************
	***************************************************************************************************/

	// flux calculation for left and right states (along x direction)
	// Total flops = 108
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE cons_point
	flux_x(const prim_point &P, const cons_point &U,
		   const metric_point &metric)
	{
		CCTK_REAL ab0, b2, w, bsubx, bsuby, bsubz;
		vlow_blow(P, metric, ab0, b2, w, bsubx, bsuby, bsubz);

		const CCTK_REAL velx = P.velx;
		const CCTK_REAL vely = P.vely;
		const CCTK_REAL velz = P.velz;
		const CCTK_REAL press = P.press;

		const CCTK_REAL dens = U.dens;
		const CCTK_REAL sx = U.sx;
		const CCTK_REAL sy = U.sy;
		const CCTK_REAL sz = U.sz;
		const CCTK_REAL tau = U.tau;
		const CCTK_REAL Bconsx = U.Bconsx;
		const CCTK_REAL Bconsy = U.Bconsy;
		const CCTK_REAL Bconsz = U.Bconsz;

		const CCTK_REAL betax = metric.betax;
		const CCTK_REAL betay = metric.betay;
		const CCTK_REAL betaz = metric.betaz;
		const CCTK_REAL alp = metric.alp;
		const CCTK_REAL sdet = sqrt(calculate_detg(metric)); // flops = 10

		const CCTK_REAL velmbetainvalpx =
			velx - betax / alp;								  // alp and beta come from metric
		const CCTK_REAL velmbetainvalpy = vely - betay / alp; // vel comes from prim
		const CCTK_REAL velmbetainvalpz = velz - betaz / alp;
		const CCTK_REAL pressstar = press + 0.5 * b2; // b2~Bvec^2
		const CCTK_REAL sqrtdetpressstar =
			sdet * pressstar; // sdet comes from metric //flops = 36

		cons_point F; // flux
		F.dens = dens * velmbetainvalpx;
		F.sx = sx * velmbetainvalpx + sqrtdetpressstar - bsubx * Bconsx / w;
		F.sy =
			sy * velmbetainvalpx - bsuby * Bconsx / w;		 // calculate w inside function
		F.sz = sz * velmbetainvalpx - bsubz * Bconsx / w; // TODO: bsub ~ Bvec_low
		F.tau = tau * velmbetainvalpx + sqrtdetpressstar * velx - ab0 * Bconsx / w;
		F.Bconsx = 0.0;
		F.Bconsy =
			Bconsy * velmbetainvalpx - Bconsx * velmbetainvalpy; // By=Bconsy
		F.Bconsz =
			Bconsz * velmbetainvalpx - Bconsx * velmbetainvalpz; // flops = 62
		return F;
	}

	// Calculate the eigenvalues along X-direction
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	eigenvalues_x(const prim_point &P, const metric_point &metric, CCTK_REAL cs2,
				  CCTK_REAL *const restrict lam)
	{

		// calculate w and u (uppermetric) here itself
		const CCTK_REAL gxx = metric.gxx;
		const CCTK_REAL gxy = metric.gxy;
		const CCTK_REAL gxz = metric.gxz;
		const CCTK_REAL gyy = metric.gyy;
		const CCTK_REAL gyz = metric.gyz;
		const CCTK_REAL gzz = metric.gzz;
		const CCTK_REAL velx = P.velx;
		const CCTK_REAL vely = P.vely;
		const CCTK_REAL velz = P.velz;
		const CCTK_REAL press = P.press;
		const CCTK_REAL rho = P.rho;
		const CCTK_REAL eps = P.eps;
		const CCTK_REAL Bvecx = P.Bvecx;
		const CCTK_REAL Bvecy = P.Bvecy;
		const CCTK_REAL Bvecz = P.Bvecz;
		const CCTK_REAL betax = metric.betax;
		const CCTK_REAL alp = metric.alp;
		const CCTK_REAL detg = calculate_detg(metric);

		// calculate uppermetric for x-direction: uxx
		const CCTK_REAL invdetg = 1.0 / detg;
		const CCTK_REAL uxx = (-gyz * gyz + gyy * gzz) * invdetg;
		const CCTK_REAL u = uxx; // TODO:u=uyy for y-direction and uzz for z-direction

		const CCTK_REAL vlowx = gxx * velx + gxy * vely + gxz * velz;
		const CCTK_REAL vlowy = gxy * velx + gyy * vely + gyz * velz;
		const CCTK_REAL vlowz = gxz * velx + gyz * vely + gzz * velz;
		const CCTK_REAL v2 = vlowx * velx + vlowy * vely + vlowz * velz;
		const CCTK_REAL w = 1. / sqrt(1. - v2);

		const CCTK_REAL boa = betax / alp; // TODO:will be betay for y-direction and betaz for z-direction
		const CCTK_REAL Bvecxlow = gxx * Bvecx + gxy * Bvecy + gxz * Bvecz;
		const CCTK_REAL Bvecylow = gxy * Bvecx + gyy * Bvecy + gyz * Bvecz;
		const CCTK_REAL Bveczlow = gxz * Bvecx + gyz * Bvecy + gzz * Bvecz;
		const CCTK_REAL B2 = Bvecxlow * Bvecx + Bvecylow * Bvecy + Bveczlow * Bvecz;
		const CCTK_REAL Bdotv = Bvecxlow * velx + Bvecylow * vely + Bveczlow * velz;
		const CCTK_REAL Bdotv2 = Bdotv * Bdotv;
		const CCTK_REAL w2 = w * w;
		const CCTK_REAL b2 = B2 / w2 + Bdotv2;
		const CCTK_REAL rhos = rho * (1.0 + eps) + press + b2;
		const CCTK_REAL va2 = b2 / rhos;
		const CCTK_REAL u2 = va2 + cs2 * (1.0 - va2);

		lam[1] = velx - boa;
		lam[2] = lam[1];
		lam[3] = lam[1];

		const CCTK_REAL lam_tmp1 =
			(velx * (1.0 - u2) -
			 sqrt(u2 * (1.0 - v2) *
				  (u * (1.0 - v2 * u2) - velx * velx * (1.0 - u2)))) /
			(1.0 - v2 * u2);
		const CCTK_REAL lam_tmp2 =
			(velx * (1.0 - u2) +
			 sqrt(u2 * (1.0 - v2) *
				  (u * (1.0 - v2 * u2) - velx * velx * (1.0 - u2)))) /
			(1.0 - v2 * u2);

		lam[0] = lam_tmp1 - boa;
		lam[4] = lam_tmp2 - boa;

		return;
	}

	/***************************************************************************************************
	*************************************** Y-DIRECTION FLUXES  ****************************************
	***************************************************************************************************/
	// flux calculation for left and right states (along y direction)
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE cons_point
	flux_y(const prim_point &P, const cons_point &U,
		   const metric_point &metric)
	{
		CCTK_REAL ab0, b2, w, bsubx, bsuby, bsubz;
		vlow_blow(P, metric, ab0, b2, w, bsubx, bsuby, bsubz);

		const CCTK_REAL velx = P.velx;
		const CCTK_REAL vely = P.vely;
		const CCTK_REAL velz = P.velz;
		const CCTK_REAL press = P.press;

		const CCTK_REAL dens = U.dens;
		const CCTK_REAL sx = U.sx;
		const CCTK_REAL sy = U.sy;
		const CCTK_REAL sz = U.sz;
		const CCTK_REAL tau = U.tau;
		const CCTK_REAL Bconsx = U.Bconsx;
		const CCTK_REAL Bconsy = U.Bconsy;
		const CCTK_REAL Bconsz = U.Bconsz;

		const CCTK_REAL betax = metric.betax;
		const CCTK_REAL betay = metric.betay;
		const CCTK_REAL betaz = metric.betaz;
		const CCTK_REAL alp = metric.alp;
		const CCTK_REAL sdet = sqrt(calculate_detg(metric));

		const CCTK_REAL velmbetainvalpx = velx - betax / alp;
		const CCTK_REAL velmbetainvalpy = vely - betay / alp;
		const CCTK_REAL velmbetainvalpz = velz - betaz / alp;
		const CCTK_REAL pressstar = press + 0.5 * b2; // b2~Bvec^2
		const CCTK_REAL sqrtdetpressstar = sdet * pressstar;

		cons_point F; // flux
		F.dens = dens * velmbetainvalpy;
		F.sx = sx * velmbetainvalpy - bsubx * Bconsy / w;
		F.sy = sy * velmbetainvalpy + sqrtdetpressstar - bsuby * Bconsy / w;
		F.sz = sz * velmbetainvalpy - bsubz * Bconsy / w; // bsub ~ Bvec_low
		F.tau = tau * velmbetainvalpy + sqrtdetpressstar * vely - ab0 * Bconsy / w;
		F.Bconsx = Bconsx * velmbetainvalpy - Bconsy * velmbetainvalpx;
		F.Bconsy = 0.0;
		F.Bconsz = Bconsz * velmbetainvalpy - Bconsy * velmbetainvalpz;
		return F;
	}

	// Calculate the eigenvalues along Y-direction
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	eigenvalues_y(const prim_point &P, const metric_point &metric, CCTK_REAL cs2,
				  CCTK_REAL *const restrict lam)
	{

		// calculate w and u (uppermetric) here itself
		const CCTK_REAL gxx = metric.gxx;
		const CCTK_REAL gxy = metric.gxy;
		const CCTK_REAL gxz = metric.gxz;
		const CCTK_REAL gyy = metric.gyy;
		const CCTK_REAL gyz = metric.gyz;
		const CCTK_REAL gzz = metric.gzz;
		const CCTK_REAL velx = P.velx;
		const CCTK_REAL vely = P.vely;
		const CCTK_REAL velz = P.velz;
		const CCTK_REAL press = P.press;
		const CCTK_REAL rho = P.rho;
		const CCTK_REAL eps = P.eps;
		const CCTK_REAL Bvecx = P.Bvecx;
		const CCTK_REAL Bvecy = P.Bvecy;
		const CCTK_REAL Bvecz = P.Bvecz;
		const CCTK_REAL betay = metric.betay;
		const CCTK_REAL alp = metric.alp;
		const CCTK_REAL detg = calculate_detg(metric);

		// calculate uppermetric for y-direction: uyy
		const CCTK_REAL invdetg = 1.0 / detg;
		const CCTK_REAL uyy = (-gxz * gxz + gxx * gzz) * invdetg;
		const CCTK_REAL u = uyy;

		const CCTK_REAL vlowx = gxx * velx + gxy * vely + gxz * velz;
		const CCTK_REAL vlowy = gxy * velx + gyy * vely + gyz * velz;
		const CCTK_REAL vlowz = gxz * velx + gyz * vely + gzz * velz;
		const CCTK_REAL v2 = vlowx * velx + vlowy * vely + vlowz * velz;
		const CCTK_REAL w = 1. / sqrt(1. - v2);

		const CCTK_REAL boa = betay / alp;
		const CCTK_REAL Bvecxlow = gxx * Bvecx + gxy * Bvecy + gxz * Bvecz;
		const CCTK_REAL Bvecylow = gxy * Bvecx + gyy * Bvecy + gyz * Bvecz;
		const CCTK_REAL Bveczlow = gxz * Bvecx + gyz * Bvecy + gzz * Bvecz;
		const CCTK_REAL B2 = Bvecxlow * Bvecx + Bvecylow * Bvecy + Bveczlow * Bvecz;
		const CCTK_REAL Bdotv = Bvecxlow * velx + Bvecylow * vely + Bveczlow * velz;
		const CCTK_REAL Bdotv2 = Bdotv * Bdotv;
		const CCTK_REAL w2 = w * w;
		const CCTK_REAL b2 = B2 / w2 + Bdotv2;
		const CCTK_REAL rhos = rho * (1.0 + eps) + press + b2;
		const CCTK_REAL va2 = b2 / rhos;
		const CCTK_REAL u2 = va2 + cs2 * (1.0 - va2);

		lam[1] = vely - boa;
		lam[2] = lam[1];
		lam[3] = lam[1];

		const CCTK_REAL lam_tmp1 =
			(vely * (1.0 - u2) -
			 sqrt(u2 * (1.0 - v2) *
				  (u * (1.0 - v2 * u2) - vely * vely * (1.0 - u2)))) /
			(1.0 - v2 * u2);
		const CCTK_REAL lam_tmp2 =
			(vely * (1.0 - u2) +
			 sqrt(u2 * (1.0 - v2) *
				  (u * (1.0 - v2 * u2) - vely * vely * (1.0 - u2)))) /
			(1.0 - v2 * u2);

		lam[0] = lam_tmp1 - boa;
		lam[4] = lam_tmp2 - boa;

		return;
	}

	/***************************************************************************************************
	*************************************** Z-DIRECTION FLUXES  ****************************************
	***************************************************************************************************/
	// flux calculation for left and right states (along z direction)
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE cons_point
	flux_z(const prim_point &P, const cons_point &U,
		   const metric_point &metric)
	{
		CCTK_REAL ab0, b2, w, bsubx, bsuby, bsubz;
		vlow_blow(P, metric, ab0, b2, w, bsubx, bsuby, bsubz);

		const CCTK_REAL velx = P.velx;
		const CCTK_REAL vely = P.vely;
		const CCTK_REAL velz = P.velz;
		const CCTK_REAL press = P.press;

		const CCTK_REAL dens = U.dens;
		const CCTK_REAL sx = U.sx;
		const CCTK_REAL sy = U.sy;
		const CCTK_REAL sz = U.sz;
		const CCTK_REAL tau = U.tau;
		const CCTK_REAL Bconsx = U.Bconsx;
		const CCTK_REAL Bconsy = U.Bconsy;
		const CCTK_REAL Bconsz = U.Bconsz;

		const CCTK_REAL betax = metric.betax;
		const CCTK_REAL betay = metric.betay;
		const CCTK_REAL betaz = metric.betaz;
		const CCTK_REAL alp = metric.alp;
		const CCTK_REAL sdet = sqrt(calculate_detg(metric));

		const CCTK_REAL velmbetainvalpx = velx - betax / alp;
		const CCTK_REAL velmbetainvalpy = vely - betay / alp;
		const CCTK_REAL velmbetainvalpz = velz - betaz / alp;
		const CCTK_REAL pressstar = press + 0.5 * b2; // b2~Bvec^2
		const CCTK_REAL sqrtdetpressstar = sdet * pressstar;

		cons_point F; // flux
		F.dens = dens * velmbetainvalpz;
		F.sx = sx * velmbetainvalpz - bsubx * Bconsz / w;
		F.sy = sy * velmbetainvalpz - bsuby * Bconsz / w;
		F.sz = sz * velmbetainvalpz + sqrtdetpressstar - bsubz * Bconsz / w; // bsub ~ Bvec_low
		F.tau = tau * velmbetainvalpz + sqrtdetpressstar * velz - ab0 * Bconsz / w;
		F.Bconsx = Bconsx * velmbetainvalpz - Bconsz * velmbetainvalpx;
		F.Bconsy = Bconsy * velmbetainvalpz - Bconsz * velmbetainvalpy;
		F.Bconsz = 0.0;
		return F;
	}

	// Calculate the eigenvalues along Z-direction
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	eigenvalues_z(const prim_point &P, const metric_point &metric, CCTK_REAL cs2,
				  CCTK_REAL *const restrict lam)
	{

		// calculate w and u (uppermetric) here itself
		const CCTK_REAL gxx = metric.gxx;
		const CCTK_REAL gxy = metric.gxy;
		const CCTK_REAL gxz = metric.gxz;
		const CCTK_REAL gyy = metric.gyy;
		const CCTK_REAL gyz = metric.gyz;
		const CCTK_REAL gzz = metric.gzz;
		const CCTK_REAL velx = P.velx;
		const CCTK_REAL vely = P.vely;
		const CCTK_REAL velz = P.velz;
		const CCTK_REAL press = P.press;
		const CCTK_REAL rho = P.rho;
		const CCTK_REAL eps = P.eps;
		const CCTK_REAL Bvecx = P.Bvecx;
		const CCTK_REAL Bvecy = P.Bvecy;
		const CCTK_REAL Bvecz = P.Bvecz;
		const CCTK_REAL betaz = metric.betaz;
		const CCTK_REAL alp = metric.alp;
		const CCTK_REAL detg = calculate_detg(metric);

		// calculate uppermetric for z-direction: uzz
		const CCTK_REAL invdetg = 1.0 / detg;
		const CCTK_REAL uzz = (-gxy * gxy + gxx * gyy) * invdetg;
		const CCTK_REAL u = uzz;

		const CCTK_REAL vlowx = gxx * velx + gxy * vely + gxz * velz;
		const CCTK_REAL vlowy = gxy * velx + gyy * vely + gyz * velz;
		const CCTK_REAL vlowz = gxz * velx + gyz * vely + gzz * velz;
		const CCTK_REAL v2 = vlowx * velx + vlowy * vely + vlowz * velz;
		const CCTK_REAL w = 1. / sqrt(1. - v2);

		const CCTK_REAL boa = betaz / alp;
		const CCTK_REAL Bvecxlow = gxx * Bvecx + gxy * Bvecy + gxz * Bvecz;
		const CCTK_REAL Bvecylow = gxy * Bvecx + gyy * Bvecy + gyz * Bvecz;
		const CCTK_REAL Bveczlow = gxz * Bvecx + gyz * Bvecy + gzz * Bvecz;
		const CCTK_REAL B2 = Bvecxlow * Bvecx + Bvecylow * Bvecy + Bveczlow * Bvecz;
		const CCTK_REAL Bdotv = Bvecxlow * velx + Bvecylow * vely + Bveczlow * velz;
		const CCTK_REAL Bdotv2 = Bdotv * Bdotv;
		const CCTK_REAL w2 = w * w;
		const CCTK_REAL b2 = B2 / w2 + Bdotv2;
		const CCTK_REAL rhos = rho * (1.0 + eps) + press + b2;
		const CCTK_REAL va2 = b2 / rhos;
		const CCTK_REAL u2 = va2 + cs2 * (1.0 - va2);

		lam[1] = velz - boa;
		lam[2] = lam[1];
		lam[3] = lam[1];

		const CCTK_REAL lam_tmp1 =
			(velz * (1.0 - u2) -
			 sqrt(u2 * (1.0 - v2) *
				  (u * (1.0 - v2 * u2) - velz * velz * (1.0 - u2)))) /
			(1.0 - v2 * u2);
		const CCTK_REAL lam_tmp2 =
			(velz * (1.0 - u2) +
			 sqrt(u2 * (1.0 - v2) *
				  (u * (1.0 - v2 * u2) - velz * velz * (1.0 - u2)))) /
			(1.0 - v2 * u2);

		lam[0] = lam_tmp1 - boa;
		lam[4] = lam_tmp2 - boa;

		return;
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	calculate_fluxes(const prim_point &P_L, const prim_point &P_R, const cons_point &U_L, const cons_point &U_R,
					 const metric_point &metric, CCTK_REAL c_s_squared_L, CCTK_REAL c_s_squared_R,
					 CCTK_REAL *const restrict lambda_left, CCTK_REAL *const restrict lambda_right,
					 cons_point &F_L, cons_point &F_R, const int flux_direction)
	{
		switch (flux_direction)
		{
		case 0:
			F_L = flux_x(P_L, U_L, metric);
			eigenvalues_x(P_L, metric, c_s_squared_L, lambda_left);

			F_R = flux_x(P_R, U_R, metric);
			eigenvalues_x(P_R, metric, c_s_squared_R, lambda_right);
			break;

		case 1:
			F_L = flux_y(P_L, U_L, metric);
			eigenvalues_y(P_L, metric, c_s_squared_L, lambda_left);

			F_R = flux_y(P_R, U_R, metric);
			eigenvalues_y(P_R, metric, c_s_squared_R, lambda_right);
			break;

		case 2:
			F_L = flux_z(P_L, U_L, metric);
			eigenvalues_z(P_L, metric, c_s_squared_L, lambda_left);

			F_R = flux_z(P_R, U_R, metric);
			eigenvalues_z(P_R, metric, c_s_squared_R, lambda_right);
			break;

		default:
			assert(0);
		}
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE cons_point
	calculate_R(const cons_point &U, const cons_point &F, const CCTK_REAL lamda)
	{
		cons_point R;
		R.dens = lamda * U.dens - F.dens;
		R.sx = lamda * U.sx - F.sx;
		R.sy = lamda * U.sy - F.sy;
		R.sz = lamda * U.sz - F.sz;
		R.tau = lamda * U.tau - F.tau;
		R.Bconsx = lamda * U.Bconsx - F.Bconsx;
		R.Bconsy = lamda * U.Bconsy - F.Bconsy;
		R.Bconsz = lamda * U.Bconsz - F.Bconsz;
		return R;
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
	calculate_p_hll(const cons_point &F_L, const cons_point &F_R,
					const cons_point &U_L, const cons_point &U_R, const CCTK_REAL lamda_L, const CCTK_REAL lamda_R)
	{
		CCTK_REAL numerator = lamda_R * F_L.sx - lamda_L * F_R.sx + lamda_L * lamda_R * (U_R.sx - U_L.sx);
		CCTK_REAL denominator = lamda_R - lamda_L;
		denominator = prevent_zero_division(denominator);

		CCTK_REAL p_hll = numerator / denominator;
		return p_hll;
	}

	// R.dens, R.sx, R.sy, R.sz, R.tau, R.Bconsx, R.Bconsy, R.Bconsz
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	calculate_v(CCTK_REAL v[3], const cons_point &R, const prim_point &P, const CCTK_REAL lamda, const CCTK_REAL p_guess)
	{
		const CCTK_REAL A = R.sx - lamda * R.tau + p_guess * (1 - SQR(lamda));
		const CCTK_REAL G = SQR(R.Bconsy) + SQR(R.Bconsz);
		const CCTK_REAL C = R.sy * R.Bconsy + R.sz * R.Bconsz;
		const CCTK_REAL Q = -A - G + SQR(P.Bvecx) * (1 - SQR(lamda));
		CCTK_REAL X = P.Bvecx * (A * lamda * P.Bvecx + C) - (A + G) * (lamda * p_guess + R.tau);

		// Making sure we do not get any 0 division errors
		X = prevent_zero_division(X);
		v[0] = (P.Bvecx * (A * P.Bvecx + C * lamda) - (A + G) * (p_guess + R.sx)) / (X);
		v[1] = (Q * R.sy + R.Bconsy * (C + P.Bvecx * (lamda * R.sx - R.tau))) / (X);
		v[2] = (Q * R.sz + R.Bconsz * (C + P.Bvecx * (lamda * R.sx - R.tau))) / (X);
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	calculate_B(CCTK_REAL B[3], const cons_point &R, const prim_point &P, const CCTK_REAL v[3], const CCTK_REAL lamda)
	{
		CCTK_REAL denominator = lamda - v[0];
		denominator = prevent_zero_division(denominator);

		B[0] = P.Bvecx;									  // B_x is constant across the discontinuity
		B[1] = (R.Bconsy - P.Bvecx * v[1]) / denominator; // B_y = R_b2 - B_x * v_y / (lamda - v_x)
		B[2] = (R.Bconsz - P.Bvecx * v[2]) / denominator; // B_z = R_b3 - B_x * v_z / (lamda - v_x)
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
	calculate_w(const cons_point &R, const prim_point &P, const CCTK_REAL v[3], const CCTK_REAL lamda, const CCTK_REAL p_guess)
	{
		CCTK_REAL denominator = lamda - v[0];
		denominator = prevent_zero_division(denominator);

		CCTK_REAL w = p_guess + (R.tau - (v[0] * R.sx + v[1] * R.sy + v[2] * R.sz)) / denominator; // w = w + R_4 - (v_x * R_1 + v_y * R_2 + v_z * R_3) / (lamda - v_x)
		return w;
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
	calculate_eta(const prim_point &P, const CCTK_REAL w, const CCTK_REAL sign)
	{
		const CCTK_REAL sign_Bx = (P.Bvecx >= 0.0) ? 1.0 : -1.0;
		const CCTK_REAL eta = sign * sign_Bx * sqrt(std::max(CCTK_REAL(1e-10), w));
		return eta;
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	calculate_K(const cons_point &R, const CCTK_REAL B[3], const CCTK_REAL lamda, const CCTK_REAL p_guess, const CCTK_REAL eta, CCTK_REAL K[3])
	{
		// Denominator for Equation 43: lamda * p + R_E + Bx * eta
		// R[4] corresponds to R_E (Energy component)
		CCTK_REAL denominator = lamda * p_guess + R.tau + B[0] * eta;
		denominator = prevent_zero_division(denominator);

		// R[1] corresponds to R_mx
		K[0] = (R.sx + p_guess + R.Bconsx * eta) / (denominator);

		// R[2] corresponds to R_my, R[6] corresponds to R_By
		K[1] = (R.sy + R.Bconsy * eta) / (denominator);

		// R[3] corresponds to R_mz, R[7] corresponds to R_Bz
		K[2] = (R.sz + R.Bconsz * eta) / (denominator);
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	calculate_B_c(CCTK_REAL B_c[3], const CCTK_REAL v_L[3], const CCTK_REAL v_R[3], const CCTK_REAL B_L[3], const CCTK_REAL B_R[3],
				  const CCTK_REAL K_L_x, const CCTK_REAL K_R_x, const CCTK_REAL P_L_B_x)
	{
		B_c[0] = P_L_B_x; // B_x is constant across the discontinuity

		CCTK_REAL denominator = K_R_x - K_L_x;
		denominator = prevent_zero_division(denominator);

		// B_c = (B_left * (lamda_R - v_left) + B_right * (v_right - lamda_L)) / (lamda_R - lamda_L)
		for (int i = 1; i < 3; i++)
		{
			B_c[i] = ((B_R[i] * (K_R_x - v_R[0]) + B_R[0] * v_R[i]) - (B_L[i] * (K_L_x - v_L[0]) + B_L[0] * v_L[i])) / (denominator);
		}
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
	calculate_Y_R(const CCTK_REAL K_L[3], const CCTK_REAL K_R[3], const CCTK_REAL B_c[3], const CCTK_REAL eta)
	{
		const CCTK_REAL top = 1 - (K_R[0] * K_R[0] + K_R[1] * K_R[1] + K_R[2] * K_R[2]);

		// for the actual value you need to divide by delta_K_x as well, but we dont do this here because we only use Y_L/Y_R to calculate f(p) (and there it cancels)
		// And dividing it by this quantity will lead in numerical instability
		CCTK_REAL denominator = eta * (K_R[0] * B_c[0] + K_R[1] * B_c[1] + K_R[2] * B_c[2]);
		denominator = prevent_zero_division(denominator);

		const CCTK_REAL Y_R = top / denominator;
		return Y_R;
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
	calculate_Y_L(const CCTK_REAL K_L[3], const CCTK_REAL K_R[3], const CCTK_REAL B_c[3], const CCTK_REAL eta)
	{
		const CCTK_REAL top = 1 - (K_L[0] * K_L[0] + K_L[1] * K_L[1] + K_L[2] * K_L[2]);

		// for the actual value you need to divide by delta_K_x as well, but we dont do this here because we only use Y_L/Y_R to calculate f(p) (and there it cancels)
		// And dividing it by this quantity will lead in numerical instability
		CCTK_REAL denominator = eta - (K_L[0] * B_c[0] + K_L[1] * B_c[1] + K_L[2] * B_c[2]);
		denominator = prevent_zero_division(denominator);
		const CCTK_REAL Y_L = top / denominator;
		return Y_L;
	}

	// TODO: Add the metric to the solver
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_REAL
	calculate_f_of_p(const prim_point &P_L, const prim_point &P_R,
					 const cons_point &R_L, const cons_point &R_R, CCTK_REAL lamda_L, CCTK_REAL lamda_R, CCTK_REAL p_guess)
	{

		CCTK_REAL v_L[3];
		calculate_v(v_L, R_L, P_L, lamda_L, p_guess);

		CCTK_REAL v_R[3];
		calculate_v(v_R, R_R, P_R, lamda_R, p_guess);

		CCTK_REAL B_L[3];
		calculate_B(B_L, R_L, P_L, v_L, lamda_L);
		CCTK_REAL B_R[3];
		calculate_B(B_R, R_R, P_R, v_R, lamda_R);

		const CCTK_REAL w_L = calculate_w(R_L, P_L, v_L, lamda_L, p_guess);
		const CCTK_REAL w_R = calculate_w(R_R, P_R, v_R, lamda_R, p_guess);

		// step 3: Calculate K and the B_c field in the intermediate state
		const CCTK_REAL eta_L = calculate_eta(P_L, w_L, -1.0);	  // eta is -1 for the left state
		const CCTK_REAL eta_R = calculate_eta(P_R, w_R, 1.0); // eta is 1 for the right state

		CCTK_REAL K_L[3];
		calculate_K(R_L, B_L, lamda_L, p_guess, eta_L, K_L);
		CCTK_REAL K_R[3];
		calculate_K(R_R, B_R, lamda_R, p_guess, eta_R, K_R);

		CCTK_REAL B_c[3];
		calculate_B_c(B_c, v_L, v_R, B_L, B_R, K_L[0], K_R[0], P_L.Bvecx);

		const CCTK_REAL Y_L = calculate_Y_L(K_L, K_R, B_c, eta_L);	  // eta is -1 for the left state
		const CCTK_REAL Y_R = calculate_Y_R(K_L, K_R, B_c, eta_R); // eta is 1 for the right state

		const CCTK_REAL delta_K_x = K_R[0] - K_L[0];
		const CCTK_REAL B_x = P_L.Bvecx;								 // B_x is constant across the discontinuity
		const CCTK_REAL f_of_p = (delta_K_x - B_x * (Y_R - Y_L)); // As mentioned in Y_L/Y_R, we dont do the delta_K_x since it gets divided out
		return f_of_p;
	}
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
	calculate_v_c(CCTK_REAL v_c[3], const CCTK_REAL K[3], const CCTK_REAL B[3], const CCTK_REAL eta)
	{

		const CCTK_REAL K_squared = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];
		const CCTK_REAL K_times_B = K[0] * B[0] + K[1] * B[1] + K[2] * B[2];

		CCTK_REAL denominator = eta - K_times_B;
		denominator = prevent_zero_division(denominator);
		v_c[0] = K[0] - B[0] * (1 - K_squared) / (denominator);
		v_c[1] = K[1] - B[1] * (1 - K_squared) / (denominator);
		v_c[2] = K[2] - B[2] * (1 - K_squared) / (denominator);
	}

	// Conditions from the mignone HLLD paper, if they match we can use
	// the solution from the HLLD solver, else not.
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE bool
	HLLD_conditions(const CCTK_REAL p, const CCTK_REAL w_L, const CCTK_REAL w_R, const CCTK_REAL v_x_aL, const CCTK_REAL v_x_aR,
					const CCTK_REAL lamda_L, const CCTK_REAL lamda_R, const CCTK_REAL v_x_cL, const CCTK_REAL v_x_cR, const CCTK_REAL lamda_a_L, const CCTK_REAL lamda_a_R)
	{
		if (w_L < p || w_R < p)
		{
			return false;
		}
		else if (v_x_aL < lamda_L || v_x_aR > lamda_R)
		{
			return false;
		}
		else if (v_x_cL < lamda_a_L || v_x_cR > lamda_a_R)
		{
			return false;
		}

		return true; // conditions all matched
	}

	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE cons_point
	calculate_HLL_flux(const CCTK_REAL lamda_R, const CCTK_REAL lamda_L, const cons_point F_L, const cons_point F_R,
					   const cons_point U_L, const cons_point U_R, cons_point &hll_flux)
	{
		CCTK_REAL denominator = lamda_R - lamda_L;
		denominator = prevent_zero_division(denominator);

		cons_point hll_flux;

		hll_flux.dens = (lamda_R * F_L.dens - lamda_L * F_R.dens + lamda_L * lamda_R * (U_R.dens - U_L.dens)) / denominator;
		hll_flux.sx = (lamda_R * F_L.sx - lamda_L * F_R.sx + lamda_L * lamda_R * (U_R.sx - U_L.sx)) / denominator;
		hll_flux.sy = (lamda_R * F_L.sy - lamda_L * F_R.sy + lamda_L * lamda_R * (U_R.sy - U_L.sy)) / denominator;
		hll_flux.sz = (lamda_R * F_L.sz - lamda_L * F_R.sz + lamda_L * lamda_R * (U_R.sz - U_L.sz)) / denominator;
		hll_flux.tau = (lamda_R * F_L.tau - lamda_L * F_R.tau + lamda_L * lamda_R * (U_R.tau - U_L.tau)) / denominator;
		hll_flux.Bconsx = (lamda_R * F_L.Bconsx - lamda_L * F_R.Bconsx + lamda_L * lamda_R * (U_R.Bconsx - U_L.Bconsx)) / denominator;
		hll_flux.Bconsy = (lamda_R * F_L.Bconsy - lamda_L * F_R.Bconsy + lamda_L * lamda_R * (U_R.Bconsy - U_L.Bconsy)) / denominator;
		hll_flux.Bconsz = (lamda_R * F_L.Bconsz - lamda_L * F_R.Bconsz + lamda_L * lamda_R * (U_R.Bconsz - U_L.Bconsz)) / denominator;

		return hll_flux;
	}

	/***************************************************************************************************
	***************************  MAIN FLUX CALCULATION ROUTINE  ****************************************
	***************************************************************************************************/
	// Flux calculation on a face
	template <int flux_direction>
	static inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE cons_point
	GRHydroX_HLLE_flux(const prim_point &P_L, const cons_point &U_L,
					   const CCTK_REAL c_s_squared_L, const prim_point &P_R,
					   const cons_point &U_R, const CCTK_REAL c_s_squared_R,
					   const metric_point &metric)
	{

		const CCTK_REAL min_value = 1e-10; // Small value to prevent division by zero and negative pressures
		CCTK_REAL lambda_leftright[10];
		CCTK_REAL *const lambda_temp_L = lambda_leftright;
		CCTK_REAL *const lambda_temp_R = lambda_leftright + 5;

		cons_point F_L;
		cons_point F_R;

		calculate_fluxes(P_L, P_R, U_L, U_R, metric, c_s_squared_L, c_s_squared_R, lambda_temp_L, lambda_temp_R, F_L, F_R, flux_direction);

		// This finds the lamda_R and lamda_L wavespeeds
		const CCTK_REAL lamda_R = max10_0(lambda_leftright);
		const CCTK_REAL lamda_L = min10_0(lambda_leftright);

		// Now we calculate the R_left and R_right according to mignone
		cons_point R_L = calculate_R(U_L, F_L, lamda_L);
		cons_point R_R = calculate_R(U_R, F_R, lamda_R);

		const CCTK_REAL p_hll = calculate_p_hll(F_L, F_R, U_L, U_R, lamda_L, lamda_R);
		CCTK_REAL p_used;
		const CCTK_REAL B_x_squared = SQR(U_L.Bconsx); // We can use either U_L or U_R since Bconsx is the same for both states

		// Mignone tells us to fall back to p_0 for weak magnetic fields proportionally to p_hll
		// This code essentially calculates the initial guess for the secant solver
		if ((B_x_squared / prevent_zero_division(p_hll)) < 0.1)
		{
			CCTK_REAL denominator = lamda_R - lamda_L;
			denominator = prevent_zero_division(denominator);

			// Here we the HLL states
			const CCTK_REAL E_hll = (lamda_R * U_R.tau - lamda_L * U_L.tau + F_L.tau - F_R.tau) / denominator;
			const CCTK_REAL m_x_hll = (lamda_R * U_R.sx - lamda_L * U_L.sx + F_L.sx - F_R.sx) / denominator;

			const CCTK_REAL F_E_hll = (lamda_R * F_L.tau - lamda_L * F_R.tau + lamda_L * lamda_R * (U_L.tau - U_R.tau)) / denominator;
			const CCTK_REAL F_mx_hll = (lamda_R * F_L.sx - lamda_L * F_R.sx + lamda_L * lamda_R * (U_L.sx - U_R.sx)) / denominator;

			const CCTK_REAL A = 1.0;
			const CCTK_REAL B = E_hll - F_mx_hll;
			const CCTK_REAL C = m_x_hll * F_E_hll - F_mx_hll * E_hll;

			// ABC formula
			CCTK_REAL discriminant = SQR(B) - 4.0 * A * C;
			if (discriminant < 0.0)
			{
				discriminant = 0.0;
			} // To prevent negative roots due to numerical errors
			CCTK_REAL p_0 = (-B + sqrt(discriminant)) / (2.0 * A); // Only plus root because you do not want negative pressure

			p_used = std::max(min_value, p_0); // clamping pressure, because you want the simulation to spit out actual values
		}
		else
		{
			p_used = std::max(min_value, p_hll); // clamping pressure, because you want the simulation to spit out actual values
		}

		CCTK_REAL p_old = 0.99 * p_used;
		CCTK_REAL p_new = 1.01 * p_used;

		CCTK_REAL f_of_p_old = calculate_f_of_p(P_left, P_right, R_left, R_right, lamda_left, lamda_right, p_old);

		// This is a simple convergence criterion, we will replace this with the actual convergence criterion later on
		while (std::abs(p_old - p_new) > 1e-6)
		{
			// Here we use the secant method with f(p) and update p till our convergence criterion meets
			const CCTK_REAL f_of_p_new = calculate_f_of_p(P_L, P_R, R_L, R_R, lamda_L, lamda_R, p_new);

			CCTK_REAL denominator = f_of_p_new - f_of_p_old;
			denominator = prevent_zero_division(denominator);
			CCTK_REAL p_next = p_new - f_of_p_new * (p_new - p_old) / (denominator); // This is a placeholder for the actual update of p, we will replace this with the actual calculation later on
			p_next = std::max(min_value, p_next);									 // prevents negative pressure
			p_old = p_new;
			f_of_p_old = f_of_p_new;

			p_new = p_next;
		}

		CCTK_REAL v_a_L[3];
		calculate_v(v_a_L, R_L, P_L, lamda_L, p_guess);
		CCTK_REAL v_a_R[3];
		calculate_v(v_a_R, R_R, P_R, lamda_R, p_guess);

		CCTK_REAL B_a_L[3];
		calculate_B(B_a_L, R_L, P_L, v_a_L, lamda_L);
		CCTK_REAL B_a_R[3];
		calculate_B(B_a_R, R_R, P_R, v_a_R, lamda_R);

		const CCTK_REAL w_a_L = calculate_w(R_L, P_L, v_a_L, lamda_L, p_guess);
		const CCTK_REAL w_a_R = calculate_w(R_R, P_R, v_a_R, lamda_R, p_guess);

		// step 3: Calculate K and the B_c field in the intermediate state
		const CCTK_REAL eta_a_L = calculate_eta(P_L, w_a_L, -1.0);  // eta is -1 for the left state
		const CCTK_REAL eta_a_R = calculate_eta(P_R, w_a_R, 1.0); // eta is 1 for the right state

		CCTK_REAL K_a_L[3];
		calculate_K(R_L, B_a_L, lamda_L, p_guess, eta_a_L, K_a_L);
		CCTK_REAL K_a_R[3];
		calculate_K(R_R, B_a_R, lamda_R, p_guess, eta_a_R, K_a_R);

		CCTK_REAL B_c[3];
		calculate_B_c(B_c, v_a_L, v_a_R, B_a_L, B_a_R, K_a_L[0], K_a_R[0], P_L.Bvecx);

		const CCTK_REAL Y_a_L = calculate_Y_L(K_a_L, K_a_R, B_c, eta_a_L);  // eta is -1 for the left state
		const CCTK_REAL Y_a_R = calculate_Y_R(K_a_L, K_a_R, B_c, eta_a_R); // eta is 1 for the right state

		CCTK_REAL v_c[3];
		calculate_v_c(v_c, K_a_L, B_c, eta_a_L); // We use the left values but we can also use the right values since at convergence they should be the same

		const CCTK_REAL lamda_c = v_c[0];
		const CCTK_REAL lamda_a_L = K_a_L[0];
		const CCTK_REAL lamda_a_R = K_a_R[0];

		cons_point final_flux;

		if (lamda_L > 0)
		{
			final_flux = F_L; // flux_chooser = F_left
		}
		else if (lamda_L < 0 && lamda_a_L > 0)
		{
			CCTK_REAL denominator = lamda_L - v_a_L[0];
			denominator = prevent_zero_division(denominator);
			const CCTK_REAL D_a_L = R_L.dens / (denominator); // D = R_D / (lamda - v_x), where R_D is the first component of R

			// TEMPORARY: Need to add the metric to this
			const prim_point P_a_L = calculate_final_P(D_a_L, v_a_L, B_a_L, w_a_L, p_found);
			final_flux = calculate_F(P_a_L); // F_left
		}
		else if (lamda_a_L < 0 && lamda_c > 0)
		{
			CCTK_REAL denominator = lamda_L - v_a_L[0];
			denominator = prevent_zero_division(denominator);
			CCTK_REAL D_a_L = R_L.dens / (denominator); // D = R_D / (lamda - v_x), where R_D is the first component of R

			const prim_point P_a_L = calculate_final_P(D_a_L, v_a_L, B_a_L, w_a_L, p_found);
			const cons_point U_a_L = calculate_U(P_a_L); // U_left
			const cons_point F_a_L = calculate_F(P_a_L); // F_left

			denominator = prevent_zero_division(denominator);

			const CCTK_REAL D_c_L = D_a_L * (lamda_a_L - v_a_L[0]) / (denominator);																						 // D_c = D_a * (lamda_a - v_a) / (lamda_c - v_a)
			const CCTK_REAL E_a_L = U_a_L.tau;																																	 // E_a = U[4]
			const CCTK_REAL m_a_x_L = F_a_L.tau;																																 // m_a_x = F[4]
			const CCTK_REAL E_c_L = (lamda_a_L * E_a_L - m_a_x_L + p_found * v_c[0] - (v_c[0] * B_c[0] + v_c[1] * B_c[1] + v_c[2] * B_c[2]) * B_a_L[0]) / (denominator); // E_c = lamda_a * E_a - m_a_x + p * v_c

			CCTK_REAL m_c_L[3];
			const CCTK_REAL v_c_dot_B_c = v_c[0] * B_c[0] + v_c[1] * B_c[1] + v_c[2] * B_c[2];
			m_c_L[0] = (E_c_L + p_found) * v_c[0] - v_c_dot_B_c * B_c[0];
			m_c_L[1] = (E_c_L + p_found) * v_c[1] - v_c_dot_B_c * B_c[1];
			m_c_L[2] = (E_c_L + p_found) * v_c[2] - v_c_dot_B_c * B_c[2];

			const cons_point U_c_L = calculate_U_intermediate_region(D_c_L, E_c_L, B_c, m_c_L);
			final_flux = calculate_flux_intermediate_region(U_a_L, U_c_L, F_a_L, lamda_c);
		}
		else if (lamda_c < 0 && lamda_a_R > 0)
		{
			CCTK_REAL denominator = lamda_R - v_a_R[0];
			denominator = prevent_zero_division(denominator);
			const CCTK_REAL D_a_R = R_R[0] / (denominator); // D = R_D / (lamda - v_x), where R_D is the first component of R

			const prim_point P_a_R = calculate_final_P(D_a_R, v_a_R, B_a_R, w_a_R, p_found);
			const cons_point U_a_R = calculate_U(P_a_R); // U_right
			const cons_point F_a_R = calculate_F(P_a_R); // F_right

			denominator = lamda_a_R - v_c[0];
			denominator = prevent_zero_division(denominator);

			const CCTK_REAL D_c_R = D_a_R * (lamda_a_R - v_a_R[0]) / (denominator);																							 // D_c = D_a * (lamda_a - v_a) / (lamda_c - v_a)
			const CCTK_REAL E_a_R = U_a_R.tau;																																		 // E_a = U[4]
			const CCTK_REAL m_a_x_R = F_a_R.tau;																																		 // m_a_x = F[4]
			const CCTK_REAL E_c_R = (lamda_a_R * E_a_R - m_a_x_R + p_found * v_c[0] - (v_c[0] * B_c[0] + v_c[1] * B_c[1] + v_c[2] * B_c[2]) * B_a_R[0]) / (denominator); // E_c = lamda_a * E_a - m_a_x + p * v_c

			CCTK_REAL m_c_R[3];
			const CCTK_REAL v_c_dot_B_c = v_c[0] * B_c[0] + v_c[1] * B_c[1] + v_c[2] * B_c[2];
			m_c_R[0] = (E_c_R + p_found) * v_c[0] - v_c_dot_B_c * B_c[0];
			m_c_R[1] = (E_c_R + p_found) * v_c[1] - v_c_dot_B_c * B_c[1];
			m_c_R[2] = (E_c_R + p_found) * v_c[2] - v_c_dot_B_c * B_c[2];

			const cons_point U_c_R = calculate_U_intermediate_region(D_c_R, E_c_R, B_c, m_c_R);
			final_flux = calculate_flux_intermediate_region(U_a_R, U_c_R, F_a_R, lamda_c);
		}
		else if (lamda_a_R < 0 && lamda_R > 0)
		{
			const CCTK_REAL D_a_R = R_right[0] / (lamda_R - v_a_R[0]); // D = R_D / (lamda - v_x), where R_D is the first component of R
			const prim_point P_a_R = calculate_final_P(D_a_R, v_a_R, B_a_R, w_a_R, p_found);
			const cons_point F_a_R = calculate_F(P_a_R); // F_right
			final_flux = F_a_R;								 // flux_chooser = F_a_R
		}
		else if (lamda_R < 0)
		{
			final_flux = F_R; // flux_chooser = F_right
		}

		// This checks whether the conditions matched or not
		bool HLLD_conditioned_matched = HLLD_conditions(p_found, w_a_L, w_a_R, v_a_left[0], v_a_R[0],
														lamda_L, lamda_R, v_c[0], v_c[0], lamda_a_L, lamda_a_R);
		// mignone tells us that if the conditions are too extreme, we should fall back to the HLL solver
		if (HLLD_conditioned_matched)
		{
			return final_flux;
		}
		else
		{
			const cons_point hll_flux = calculate_HLL_flux(lamda_R, lamda_L, F_left, F_R, U_L, U_R);
			return hll_flux;
		}
	}

} // namespace GRHydroX

#endif // #ifndef GRHYDROX_HLLD_HXX