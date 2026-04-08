#include <iostream>
#include <vector>
#include <cmath>

const double GAMMA_EOS = 5.0 / 3.0;
const double EPSILON = 1e-10;

struct conserved_variables
{
	double D;
	double v_x;
	double v_y;
	double v_z;
	double B_x;
	double B_y;
	double B_z;
	double w;
	double p;
};

std::vector<double> HLLD_solver(std::vector<double> P_left, std::vector<double> P_right);

int main()
{
	// our starting values for the left and right states
	double rho_left;
	std::vector<double> v_left;
	std::vector<double> B_left;
	double p_g_left;
	double rho_right;
	std::vector<double> v_right;
	std::vector<double> B_right;
	double p_g_right;

	// This should be replaced with the actual values for the left and right states, for now we will just use some placeholder values
	conserved_variables P_left = calculate_P(rho_left, v_left, B_left, p_g_left);
	conserved_variables P_right = calculate_P(rho_right, v_right, B_right, p_g_right);

	// Here we call the HLLD solver with the left and right states and get the result
	std::vector<double> result = HLLD_solver(P_left, P_right, p_g_left, p_g_right);

	return 0;
}

double calculate_gamma_factor(double v_x, double v_y, double v_z)
{
	double denominator = std::sqrt(1 - (v_x * v_x + v_y * v_y + v_z * v_z));
	denominator = prevent_zero_division(denominator);
	double gamma_factor = 1 / denominator;
	return gamma_factor;
}

double calculate_b0(double B_x, double B_y, double B_z, double v_x, double v_y, double v_z, double gamma_factor)
{
	double b0 = (B_x * v_x + B_y * v_y + B_z * v_z) * gamma_factor;
	return b0;
}

double calculate_eta(conserved_variables P, double w, double sign)
{	
	double sign_Bx = (P.B_x >= 0.0) ? 1.0 : -1.0;
	double eta = sign * sign_Bx * std::sqrt(w);
	return eta;
}

std::vector<double> calculate_bk(double b_0, double B_x, double B_y, double B_z, double v_x, double v_y, double v_z, double gamma_factor)
{
	std::vector<double> bk(3);
	bk[0] = B_x / gamma_factor + b_0 * v_x;
	bk[1] = B_y / gamma_factor + b_0 * v_y;
	bk[2] = B_z / gamma_factor + b_0 * v_z;
	return bk;
}
// P = {D, v_x, v_y, v_z, B_x, B_y, B_z, w, p}
std::vector<double> calculate_U(conserved_variables P)
{
	// First we compute the relevant quantities needed for the conservative variables
	double gamma_factor{calculate_gamma_factor(P.v_x, P.v_y, P.v_z)};
	double b_0{calculate_b0(P.B_x, P.B_y, P.B_z, P.v_x, P.v_y, P.v_z, gamma_factor)};
	std::vector<double> bk{calculate_bk(b_0, P.B_x, P.B_y, P.B_z, P.v_x, P.v_y, P.v_z, gamma_factor)};

	// now we define U and calculate the quantities needed for the conservative variables
	std::vector<double> U(8);
	U[0] = P.D;														// U[0] = D
	U[1] = P.w * gamma_factor * gamma_factor * P.v_x - b_0 * bk[0]; // U[1] = w * gamma_factor * gamma_factor * v_x - b_0 * bk[0]
	U[2] = P.w * gamma_factor * gamma_factor * P.v_y - b_0 * bk[1]; // U[2] = w * gamma_factor * gamma_factor * v_y - b_0 * bk[1]
	U[3] = P.w * gamma_factor * gamma_factor * P.v_z - b_0 * bk[2]; // U[3] = w * gamma_factor * gamma_factor * v_z - b_0 * bk[2]
	U[4] = P.w * gamma_factor * gamma_factor - P.p - b_0 * b_0;		// U[4] = w * gamma_factor * gamma_factor - p - b_0^2
	U[5] = P.B_x;													// U[5] = B_x
	U[6] = P.B_y;													// U[6] = B_y
	U[7] = P.B_z;													// U[7] = B_z

	return U;
}

std::vector<double> calculate_F(conserved_variables P)
{
	// First we compute the relevant quantities needed for the conservative variables
	double gamma_factor{calculate_gamma_factor(P.v_x, P.v_y, P.v_z)};
	double b_0{calculate_b0(P.B_x, P.B_y, P.B_z, P.v_x, P.v_y, P.v_z, gamma_factor)};
	std::vector<double> bk{calculate_bk(b_0, P.B_x, P.B_y, P.B_z, P.v_x, P.v_y, P.v_z, gamma_factor)};

	std::vector<double> F(8);
	F[0] = P.D * P.v_x; // F[0] = D * v_x
	F[1] = P.w * gamma_factor * gamma_factor * P.v_x * P.v_x - bk[0] * bk[0] + P.p;
	F[2] = P.w * gamma_factor * gamma_factor * P.v_y * P.v_x - bk[1] * bk[0];
	F[3] = P.w * gamma_factor * gamma_factor * P.v_z * P.v_x - bk[2] * bk[0];
	F[4] = P.w * gamma_factor * gamma_factor * P.v_x - b_0 * bk[0];
	F[5] = 0;
	F[6] = P.B_y * P.v_x - P.B_x * P.v_y;
	F[7] = P.B_z * P.v_x - P.B_x * P.v_z;
	return F;
}
// P = {D, v_x, v_y, v_z, B_x, B_y, B_z, w, p}
conserved_variables calculate_P(double rho, std::vector<double> v, std::vector<double> B, double p_g)
{
	double gamma = calculate_gamma_factor(v[0], v[1], v[2]);
	double D = rho * gamma; // D = rho * gamma
	double b_squared = (B[0] * B[0] + B[1] * B[1] + B[2] * B[2]) / (gamma * gamma) +
					   (B[0] * v[0] + B[1] * v[1] + B[2] * v[2]) * (B[0] * v[0] + B[1] * v[1] + B[2] * v[2]); // b^2 = (B^2 / gamma^2) + (B dot v)^2
	double w_g = rho + p_g * (GAMMA_EOS / (GAMMA_EOS - 1));													  // w = rho + p * (gamma_eos / (gamma_eos - 1))
	double w = w_g + b_squared;																				  // w = w_g + b^2
	double p = p_g + b_squared / 2;																			  // p = p_g + b^2 / 2

	conserved_variables P;
	P.D = D;
	P.v_x = v[0];
	P.v_y = v[1];
	P.v_z = v[2];
	P.B_x = B[0];
	P.B_y = B[1];
	P.B_z = B[2];
	P.w = w;
	P.p = p;

	return P;
}

// Used when we divide by a quantity such that we do not get any 0 divison errors
double prevent_zero_division(double value) {
	if (std::abs(value) < EPSILON) {
		if (value >= 0) {
			return EPSILON;
		} else {
			return -EPSILON;
		}
	} else {
		return value;
	}
}

std::vector<double> calculate_R(std::vector<double> U, std::vector<double> F, double lamda)
{
	std::vector<double> R(8);
	for (int i = 0; i < 8; i++)
	{
		R[i] = lamda * U[i] - F[i];
	}
	return R;
}

double quartic_function(double lamda, conserved_variables P, double p_g)
{	
	P.w = prevent_zero_division(P.w);
	double c_s_squared = GAMMA_EOS * p_g / P.w; // c_s^2 = gamma_eos * p_g / w
	double gamma = calculate_gamma_factor(P.v_x, P.v_y, P.v_z);
	double a = gamma * (lamda - P.v_x);
	double B_squared = P.B_x * P.B_x + P.B_y * P.B_y + P.B_z * P.B_z;
	double v_dot_B = P.v_x * P.B_x + P.v_y * P.B_y + P.v_z * P.B_z;

	// No 0 division check needed here as gamma always is greater than 1
	double b_x = gamma * (B_squared / (gamma * gamma) + P.v_x * (P.B_x * P.v_x + P.B_y * P.v_y + P.B_z * P.v_z));
	double b_0 = gamma * v_dot_B;
	double special_B_squared = (b_x - lamda * b_0) * (b_x - lamda * b_0);
	double abs_b_squared = std::abs(P.B_x * P.B_x + P.B_y * P.B_y + P.B_z * P.B_z) / (gamma * gamma) +
					   (P.B_x * P.v_x + P.B_y * P.v_y + P.B_z * P.v_z) * (P.B_x * P.v_x + P.B_y * P.v_y + P.B_z * P.v_z);

	double F_of_lamda = P.w * (1- c_s_squared) * std::pow(a, 4.0) - (1- lamda * lamda) * ((abs_b_squared + P.w * c_s_squared) * a * a - c_s_squared) * special_B_squared;
	return F_of_lamda;
}

double bisection_loop_lamda(conserved_variables P, double left_lamda, double right_lamda, double p_g, double v_x) {
	double error_margin = 1e-6;
	double midpoint_lamda = (left_lamda + right_lamda) / 2.0;

	while (std::abs(left_lamda - right_lamda) > error_margin)
	{
		double midpoint_lamda = (left_lamda + right_lamda) / 2.0;
		double f_midpoint = quartic_function(midpoint_lamda, P, p_g);
		double f_left = quartic_function(left_lamda, P, p_g);
		
		if (f_midpoint * f_left < 0) {
			right_lamda = midpoint_lamda;
		} else if (f_midpoint * f_left > 0) {
			left_lamda = midpoint_lamda;
		} else {
			return midpoint_lamda; // Found an exact root
		}

	}
	return midpoint_lamda;
}
// I decided to use the bisection method for now
// In the future, I might replace this with a faster algorithm
std::vector<double> calculate_lamdas(conserved_variables P_left, conserved_variables P_right, double p_g_left, double p_g_right)
{	
	// For the left side of the boundary
	double left_plus_lamda = P_left.v_x;
	double right_plus_lamda = 1;

	double left_plus_lamda_found = bisection_loop_lamda(P_left, left_plus_lamda, right_plus_lamda, p_g_left, P_left.v_x);

	double left_minus_lamda = -1;
	double right_minus_lamda = P_left.v_x;
	double left_minus_lamda_found = bisection_loop_lamda(P_left, left_minus_lamda, right_minus_lamda, p_g_left, P_left.v_x);

	// For the right side of the boundary
	double left_plus_lamda = P_right.v_x;
	double right_plus_lamda = 1;

	double left_plus_lamda_found = bisection_loop_lamda(P_right, left_plus_lamda, right_plus_lamda, p_g_right, P_right.v_x);

	double left_minus_lamda = -1;
	double right_minus_lamda = -P_right.v_x;
	double left_minus_lamda_found = bisection_loop_lamda(P_right, left_minus_lamda, right_minus_lamda, p_g_right, P_right.v_x);

	double lamda_left = std::min(left_minus_lamda, right_minus_lamda);
	double lamda_right = std::max(left_plus_lamda, right_plus_lamda);

	return {lamda_left, lamda_right};
}


// R_D = R_0
// R_mk = R_123
// R_E = R_4
// R_bk = R_5, R_6, R_7
std::vector<double> calculate_v(std::vector<double> R, conserved_variables P, double lamda, double p_guess)
{
	std::vector<double> v(3);

	double A{R[1] - lamda * R[4] + p_guess * (1 - lamda * lamda)};
	double G{R[6] * R[6] + R[7] * R[7]};
	double C{R[2] * R[6] + R[3] * R[7]};
	double Q = {-A - G + P.B_x * P.B_x * (1 - lamda * lamda)};
	double X = {P.B_x * (A * lamda * P.B_x + C) - (A + G) * (lamda * p_guess + R[4])};


	// Making sure we do not get any 0 division errors
	X = prevent_zero_division(X);
	v[0] = (P.B_x * (A * P.B_x + C * lamda) - (A + G) * (p_guess + R[1])) / (X);
	v[1] = (Q * R[2] + R[6] * (C + P.B_x* (lamda * R[1] - R[4]))) / (X);
	v[2] = (Q * R[3] + R[7] * (C + P.B_x * (lamda * R[1] - R[4]))) / (X);

	return v;
}

std::vector<double> calculate_B(std::vector<double> R, conserved_variables P, std::vector<double> v, double lamda)
{
	std::vector<double> B(3);

	double denominator = lamda - v[0];
	denominator = prevent_zero_division(denominator);
	

	B[0] = P.B_x; // B_x is constant across the discontinuity
	B[1] = (R[6] - P.B_x * v[1]) / denominator; // B_y = R_b2 - B_x * v_y / (lamda - v_x)
	B[2] = (R[7] - P.B_x * v[2]) / denominator; // B_z = R_b3 - B_x * v_z / (lamda - v_x)

	return B;
}

double calculate_w(std::vector<double> R, conserved_variables P, std::vector<double> v, double lamda, double p_guess)
{
	double denominator = lamda - v[0];
	denominator = prevent_zero_division(denominator);

	double w = p_guess + (R[4] - (v[0] * R[1] + v[1] * R[2] + v[2] * R[3])) / denominator; // w = w + R_4 - (v_x * R_1 + v_y * R_2 + v_z * R_3) / (lamda - v_x)
	return w;
}

std::vector<double> calculate_K(std::vector<double> R, std::vector<double> B, double lamda, double p_guess, double eta)
{
	std::vector<double> K(3);

	// Denominator for Equation 43: lamda * p + R_E + Bx * eta
    // R[4] corresponds to R_E (Energy component)
    double denominator = lamda * p_guess + R[4] + B[0] * eta;
	denominator = prevent_zero_division(denominator);

    // R[1] corresponds to R_mx
    K[0] = (R[1] + p_guess + R[5] * eta) / (denominator);

    // R[2] corresponds to R_my, R[6] corresponds to R_By
    K[1] = (R[2] + R[6] * eta) / (denominator);

    // R[3] corresponds to R_mz, R[7] corresponds to R_Bz
    K[2] = (R[3] + R[7] * eta) / (denominator);

	return K;
}

std::vector<double> calculate_B_c(std::vector<double> v_left, std::vector<double> v_right, std::vector<double> B_left, std::vector<double> B_right, 
	double lamda_L, double lamda_R, double K_left_x, double K_right_x, double P_left_B_x)
{
	std::vector<double> B_c(3);

	B_c[0] = P_left_B_x; // B_x is constant across the discontinuity

	double denominator = K_right_x - K_left_x;
	denominator = prevent_zero_division(denominator);

	for (int i = 1; i < 3; i++)
	{
		B_c[i] = ((B_right[i] * (lamda_R - v_right[0]) + B_right[0] * v_right[i]) - (B_left[i] * (lamda_L - v_left[0]) + B_left[0] * v_left[i])) / (denominator); // B_c = (B_left * (lamda_R - v_left) + B_right * (v_right - lamda_L)) / (lamda_R - lamda_L)
	}

	return B_c;
}

// PLACEHOLDER CODE FOR CALCULATING p_hll, we will replace this with the actual calculation later on
double calculate_p_hll(conserved_variables P)
{
	double p_hll = P.p;
	return p_hll;
}

double calculate_Y_R(std::vector<double> K_left, std::vector<double> K_right, std::vector<double> B_c, double eta)
{
	double Y_R;
	double delta_K_x = K_right[0] - K_left[0];
	double top = 1 - (K_right[0] * K_right[0] + K_right[1] * K_right[1] + K_right[2] * K_right[2]);
	double denominator = eta * delta_K_x - delta_K_x * (K_right[0] * B_c[0] + K_right[1] * B_c[1] + K_right[2] * B_c[2]);
	denominator = prevent_zero_division(denominator);

	Y_R = top / denominator;
	return Y_R;
}

double calculate_Y_L(std::vector<double> K_left, std::vector<double> K_right, std::vector<double> B_c, double eta)
{
	double Y_L;
	double delta_K_x = K_right[0] - K_left[0];
	double top = 1 - (K_left[0] * K_left[0] + K_left[1] * K_left[1] + K_left[2] * K_left[2]);
	double denominator = eta * delta_K_x - delta_K_x * (K_left[0] * B_c[0] + K_left[1] * B_c[1] + K_left[2] * B_c[2]);
	denominator = prevent_zero_division(denominator);
	Y_L = top / denominator;
	return Y_L;
}

conserved_variables calculate_final_P(double D, std::vector<double> v, std::vector<double> B, double w, double p_found)
{
	conserved_variables final_P;

	final_P.D = D;
	final_P.v_x = v[0];
	final_P.v_y = v[1];
	final_P.v_z = v[2];
	final_P.B_x = B[0]; // this value should be const across the whole calculation
	final_P.B_y = B[1];
	final_P.B_z = B[2];
	final_P.w = w;
	final_P.p = p_found;
	return final_P;
}

double calculate_f_of_p(conserved_variables P_left, conserved_variables P_right,
						std::vector<double> R_left, std::vector<double> R_right, double lamda_left, double lamda_right, double p_guess)
{

	std::vector<double> v_left = calculate_v(R_left, P_left, lamda_left, p_guess);	  
	std::vector<double> v_right = calculate_v(R_right, P_right, lamda_right, p_guess); 

	std::vector<double> B_left = calculate_B(R_left, P_left, v_left, lamda_left);	  
	std::vector<double> B_right = calculate_B(R_right, P_right, v_right, lamda_right);

	double w_left = calculate_w(R_left, P_left, v_left, lamda_left, p_guess);	  
	double w_right = calculate_w(R_right, P_right, v_right, lamda_right, p_guess);

	// step 3: Calculate K and the B_c field in the intermediate state
	double eta_left = calculate_eta(P_left, w_left, -1); // eta is -1 for the left state
	double eta_right = calculate_eta(P_right, w_right, 1); // eta is 1 for the right state

	std::vector<double> K_left = calculate_K(R_left, B_left, lamda_left, p_guess, eta_left); 
	std::vector<double> K_right = calculate_K(R_right, B_right, lamda_right, p_guess, eta_right);

	std::vector<double> B_c = calculate_B_c(v_left, v_right, B_left, B_right, lamda_left, lamda_right, K_left[0], K_right[0]);

	double Y_left = calculate_Y_L(K_left, K_right, B_c, eta_left); // eta is -1 for the left state
	double Y_right = calculate_Y_R(K_left, K_right, B_c, eta_right); // eta is 1 for the right state

	double delta_K_x = K_right[0] - K_left[0];
	double B_x = P_left.B_x; // B_x is constant across the discontinuity
	double f_of_p = (delta_K_x * (1 - B_x * (Y_right - Y_left))); // This is a placeholder for the actual calculation of f(p), we will replace this with the actual calculation later on
	return f_of_p;
}

std::vector<double> calculate_v_c(std::vector<double> K, std::vector<double> B, double eta)
{
	std::vector<double> v_c(3);

	double K_squared = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];
	double K_times_B = K[0] * B[0] + K[1] * B[1] + K[2] * B[2];

	double denominator = eta - K_times_B;
	denominator = prevent_zero_division(denominator);
	v_c[0] = K[0] - B[0] * (1 - K_squared) / (denominator);
	v_c[1] = K[1] - B[1] * (1 - K_squared) / (denominator);
	v_c[2] = K[2] - B[2] * (1 - K_squared) / (denominator);

	return v_c;
}

std::vector<double> calculate_flux_intermediate_region(std::vector<double> U_a, std::vector<double> U_c, std::vector<double> F_a, double lamda_c)
{
	std::vector<double> flux(8);
	for (int i = 0; i < 8; i++)
	{
		flux[i] = F_a[i] + lamda_c * (U_c[i] - U_a[i]);
	}
	return flux;
}

std::vector<double> calculate_U_intermediate_region(double D_c, double E_c, std::vector<double> B_c, std::vector<double> m_c) 
{
	std::vector<double> U_c(8);
	U_c[0] = D_c;
	U_c[1] = m_c[0];
	U_c[2] = m_c[1];
	U_c[3] = m_c[2];
	U_c[4] = E_c;
	U_c[5] = B_c[0];
	U_c[6] = B_c[1];
	U_c[7] = B_c[2];

	return U_c;
}

std::vector<double> HLLD_solver(conserved_variables P_left, conserved_variables P_right, double p_g_left, double p_g_right)
{

	// Here we calculate the conservative variables for the left and right states
	std::vector<double> U_left = calculate_U(P_left);
	std::vector<double> U_right = calculate_U(P_right);

	std::vector<double> F_left = calculate_F(P_left);
	std::vector<double> F_right = calculate_F(P_right);

	// Here we calculate the wave speeds and the intermediate states
	std::vector<double> lamdas = calculate_lamdas(P_left, P_right, p_g_left, p_g_right);
	double lamda_left = lamdas[0];
	double lamda_right = lamdas[1];

	std::vector<double> R_left = calculate_R(U_left, F_left, lamda_left);	  // temp lamda for now, we will replace this with the actual calculation later on
	std::vector<double> R_right = calculate_R(U_right, F_right, lamda_right); // temp lamda for now, we will replace this with the actual calculation later on

	double p_old = calculate_p_hll(P_left); // This is a placeholder for the actual calculation of p_hll, we will replace this with the actual calculation later on
	double p_new = calculate_p_hll(P_right); // This is a placeholder for the actual calculation of p_hll, we will replace this with the actual calculation later on

	double f_of_p_old = calculate_f_of_p(P_left, P_right, R_left, R_right, lamda_left, lamda_right, p_old);

	// This is a simple convergence criterion, we will replace this with the actual convergence criterion later on
	while (std::abs(p_old - p_new) > 1e-6) {
		
		// Here we use the secant method with f(p) and update p till our convergence criterion meets
		double f_of_p_new = calculate_f_of_p(P_left, P_right, R_left, R_right, lamda_left, lamda_right, p_new);

		double denominator = f_of_p_new - f_of_p_old;
		denominator = prevent_zero_division(denominator);
		double p_next = p_new -f_of_p_new * (p_new - p_old) / (denominator); // This is a placeholder for the actual update of p, we will replace this with the actual calculation later on
		
		p_old = p_new;
		f_of_p_old = f_of_p_new;

		p_new = p_next;
	}

	double p_found = p_new; // This is the value of p that we found after convergence

	// Now we calculate all the final values for the a left and a right and c state
	std::vector<double> v_a_left = calculate_v(R_left, P_left, lamda_left, p_found);	  
	std::vector<double> v_a_right = calculate_v(R_right, P_right, lamda_right, p_found); 

	std::vector<double> B_a_left = calculate_B(R_left, P_left, v_a_left, lamda_left);	  
	std::vector<double> B_a_right = calculate_B(R_right, P_right, v_a_right, lamda_right);

	double w_a_left = calculate_w(R_left, P_left, v_a_left, lamda_left, p_found);	  
	double w_a_right = calculate_w(R_right, P_right, v_a_right, lamda_right, p_found);

	double eta_a_left = calculate_eta(P_left, w_a_left, -1); // eta is -1 for the left state
	double eta_a_right = calculate_eta(P_right, w_a_right, 1); // eta is 1 for the right state

	std::vector<double> K_a_left = calculate_K(R_left, B_a_left, lamda_left, p_found, eta_a_left); 
	std::vector<double> K_a_right = calculate_K(R_right, B_a_right, lamda_right, p_found, eta_a_right);

	// B_c_left = B_c_right = B_c, since B_c is constant across the discontinuity
	std::vector<double> B_c = calculate_B_c(v_a_left, v_a_right, B_a_left, B_a_right, lamda_left, lamda_right, K_a_left[0], K_a_right[0]);

	double Y_a_left = calculate_Y_L(K_a_left, K_a_right, B_c, eta_a_left); // eta is -1 for the left state
	double Y_a_right = calculate_Y_R(K_a_left, K_a_right, B_c, eta_a_right); // eta is 1 for the right state

	double lamda_c = 0; // This is a placeholder for the actual calculation of lamda_c, we will replace this with the actual calculation later on
	double lamda_a_left = K_a_left[0];
	double lamda_a_right = K_a_right[0];

	std::vector<double> final_flux;
	if (lamda_left > 0) {
		final_flux = F_left; // flux_chooser = F_left
		
	} else if (lamda_left < 0 && lamda_a_left > 0) {
		double denominator = lamda_left - v_a_left[0];
		denominator = prevent_zero_division(denominator);
		double D_a_left = R_left[0] / (denominator); // D = R_D / (lamda - v_x), where R_D is the first component of R
		conserved_variables P_a_left = calculate_final_P(D_a_left, v_a_left, B_a_left, w_a_left, p_found); 
		final_flux = calculate_F(P_a_left); // F_left

	} else if (lamda_a_left < 0 && lamda_c > 0) {
		double denominator = lamda_left - v_a_left[0];
		denominator = prevent_zero_division(denominator);
		double D_a_left = R_left[0] / (denominator); // D = R_D / (lamda - v_x), where R_D is the first component of R
		
		conserved_variables P_a_left = calculate_final_P(D_a_left, v_a_left, B_a_left, w_a_left, p_found); 
		std::vector<double> U_a_left = calculate_U(P_a_left); // U_left
		std::vector<double> F_a_left = calculate_F(P_a_left); // F_left

		std::vector<double> v_c_left = calculate_v_c(K_a_left, B_a_left, eta_a_left); // v_c is a placeholder for the actual calculation of v_c, we will replace this with the actual calculation later on
		denominator = lamda_a_left - v_c_left[0];
		denominator = prevent_zero_division(denominator);
		
		double D_c_left = D_a_left * (lamda_a_left - v_a_left[0]) / (denominator); // D_c = D_a * (lamda_a - v_a) / (lamda_c - v_a)
		double E_a_left = U_a_left[4]; // E_a = U[4]
		double m_a_x_left = F_a_left[4]; // m_a_x = F[4]
		double E_c_left = (lamda_a_left * E_a_left - m_a_x_left + p_found * v_c_left[0] - (v_c_left[0] * B_c[0] + v_c_left[1] * B_c[1] + v_c_left[2] * B_c[2]) * B_a_left[0]) / (denominator); // E_c = lamda_a * E_a - m_a_x + p * v_c
		
		std::vector<double> m_c_left(3);
		double v_c_dot_B_c = v_c_left[0] * B_c[0] + v_c_left[1] * B_c[1] + v_c_left[2] * B_c[2];
		m_c_left[0] = (E_c_left + p_found) * v_c_left[0] - v_c_dot_B_c * B_c[0];
		m_c_left[1] = (E_c_left + p_found) * v_c_left[1] - v_c_dot_B_c * B_c[1];
		m_c_left[2] = (E_c_left + p_found) * v_c_left[2] - v_c_dot_B_c * B_c[2];

		std::vector<double> U_c_left = calculate_U_intermediate_region(D_c_left, E_c_left, B_c, m_c_left); 
		final_flux = calculate_flux_intermediate_region(U_a_left, U_c_left, F_a_left, lamda_c);

	} else if (lamda_c < 0 && lamda_a_right > 0) {
		double denominator = lamda_right - v_a_right[0];
		denominator = prevent_zero_division(denominator);
		double D_a_right = R_right[0] / (denominator); // D = R_D / (lamda - v_x), where R_D is the first component of R
		
		conserved_variables P_a_right = calculate_final_P(D_a_right, v_a_right, B_a_right, w_a_right, p_found); 
		std::vector<double> U_a_right = calculate_U(P_a_right); // U_right
		std::vector<double> F_a_right = calculate_F(P_a_right); // F_right

		std::vector<double> v_c_right = calculate_v_c(K_a_right, B_a_right, eta_a_right); 

		denominator = lamda_a_right - v_c_right[0];
		denominator = prevent_zero_division(denominator);
		
		double D_c_right = D_a_right * (lamda_a_right - v_a_right[0]) / (denominator); // D_c = D_a * (lamda_a - v_a) / (lamda_c - v_a)
		double E_a_right = U_a_right[4]; // E_a = U[4]
		double m_a_x_right = F_a_right[4]; // m_a_x = F[4]
		double E_c_right = (lamda_a_right * E_a_right - m_a_x_right + p_found * v_c_right[0] - (v_c_right[0] * B_c[0] + v_c_right[1] * B_c[1] + v_c_right[2] * B_c[2]) * B_a_right[0]) / (denominator); // E_c = lamda_a * E_a - m_a_x + p * v_c
		
		std::vector<double> m_c_right(3);
		double v_c_dot_B_c = v_c_right[0] * B_c[0] + v_c_right[1] * B_c[1] + v_c_right[2] * B_c[2];
		m_c_right[0] = (E_c_right + p_found) * v_c_right[0] - v_c_dot_B_c * B_c[0];
		m_c_right[1] = (E_c_right + p_found) * v_c_right[1] - v_c_dot_B_c * B_c[1];
		m_c_right[2] = (E_c_right + p_found) * v_c_right[2] - v_c_dot_B_c * B_c[2];

		std::vector<double> U_c_right = calculate_U_intermediate_region(D_c_right, E_c_right, B_c, m_c_right); 
		final_flux = calculate_flux_intermediate_region(U_a_right, U_c_right, F_a_right, lamda_c);

	} else if (lamda_a_right < 0 && lamda_right > 0) {
		double D_a_right = R_right[0] / (lamda_right - v_a_right[0]); // D = R_D / (lamda - v_x), where R_D is the first component of R
		conserved_variables P_a_right = calculate_final_P(D_a_right, v_a_right, B_a_right, w_a_right, p_found); 
		final_flux = calculate_F(P_a_right); // F_right
	
	} else if (lamda_right < 0) {
		final_flux = F_right; // flux_chooser = F_right
	}

	return final_flux;
}
