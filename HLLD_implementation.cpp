#include <iostream>
#include <vector>
#include <cmath>

const double GAMMA_EOS = 5.0f / 3.0f;
const double SPEED_OF_LIGHT = 299 792 458.0f;

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
	std::vector<double> result = HLLD_solver(P_left, P_right);

	return 0;
}

double calculate_gamma_factor(double v_x, double v_y, double v_z)
{
	double gamma_factor = 1.0 / sqrt(1 - (v_x * v_x + v_y * v_y + v_z * v_z) / (SPEED_OF_LIGHT * SPEED_OF_LIGHT)); // gamma = 1 / sqrt(1 - v^2 / c^2)
	return gamma_factor;
}

double calculate_b0(double B_x, double B_y, double B_z, double v_x, double v_y, double v_z, double gamma_factor)
{
	double b0 = (B_x * v_x + B_y * v_y + B_z * v_z) * gamma_factor;
	return b0;
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

	std::vector<double> F;
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
	P.w = w_g;
	P.p = p;

	return P;
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

// These are temporary functions to calculate the wave speeds, we will replace these with the actual calculations later on
double calculate_lamda_right()
{
	return 0.99 * 3e8f;
}

double calculate_lamda_left()
{
	return -0.99 * 3e8f;
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

	v[0] = (P.B_x * (A * P.B_x + C * lamda) - (A + G) * (p_guess + R[1])) / X;
	v[1] = (Q * R[2] + R[6] * (C + P.B_x* (lamda * R[1] - R[4]))) / X;
	v[2] = (Q * R[3] + R[7] * (C + P.B_x * (lamda * R[1] - R[4]))) / X;

	return v;
}

std::vector<double> calculate_B(std::vector<double> R, conserved_variables P, double lamda)
{
	std::vector<double> B(3);

	B[0] = P.B_x; // B_x is constant across the discontinuity
	B[1] = (R[6] - P.B_x * P.v_y) / (lamda - P.v_x); // B_y = R_b2 - B_x * v_y / (lamda - v_x)
	B[2] = (R[7] - P.B_x * P.v_z) / (lamda - P.v_x); // B_z = R_b3 - B_x * v_z / (lamda - v_x)

	return B;
}

double calculate_w(std::vector<double> R, conserved_variables P, double lamda, double p_guess)
{
	double w = p_guess + (R[4] - (P.v_x * R[1] + P.v_y * R[2] + P.v_z * R[3])) / (lamda - P.v_x); // w = w + R_4 - (v_x * R_1 + v_y * R_2 + v_z * R_3) / (lamda - v_x)
	return w;
}

std::vector<double> calculate_K(std::vector<double> R, conserved_variables P, double lamda)
{
	std::vector<double> K(3);

	for (int i = 0; i < 3; i++)
	{
		K[i] = (R[5 + i] - P.B_x * P.v_x) / (lamda - P.v_x); // K = R_bk - B_x * v_k / (lamda - v_x)
	}

	return K;
}

std::vector<double> calculate_B_c(std::vector<double> R, conserved_variables P, double lamda)
{
	std::vector<double> B_c(3);

	for (int i = 0; i < 3; i++)
	{
		B_c[i] = (R[5 + i] - P.B_x * P.v_x) / (lamda - P.v_x); // B_c = R_bk - B_x * v_k / (lamda - v_x)
	}

	return B_c;
}

double calculate_p_hll(conserved_variables P_left, conserved_variables P_right)
{
	double p1 = P_left.p;
	double p2 = P_right.p;
	return (p1 + p2) / 2; // This is a very simple approximation for the HLL pressure, we will replace this with the actual calculation later on
}

double calculate_Y_R(std::vector<double> K_left, std::vector<double> K_right, std::vector<double> B_c)
{
	double Y_R;
	// This is a placeholder for the actual calculation of Ys, we will replace this with the actual calculation later on
	Y_R = 0;
	return Y_R;
}

double calculate_Y_L(std::vector<double> K_left, std::vector<double> K_right, std::vector<double> B_c)
{
	double Y_L;
	// This is a placeholder for the actual calculation of Ys, we will replace this with the actual calculation later on
	Y_L = 0;
	return Y_L;
}


double calculate_f_of_p(conserved_variables P_left, conserved_variables P_right,
						std::vector<double> R_left, std::vector<double> R_right, double lamda_left, double lamda_right, double p_guess)
{

	std::vector<double> v_left = calculate_v(R_left, P_left, lamda_left, p_guess);	  
	std::vector<double> v_right = calculate_v(R_right, P_right, lamda_right, p_guess); 

	std::vector<double> B_left = calculate_B(R_left, P_left, lamda_left);	  
	std::vector<double> B_right = calculate_B(R_right, P_right, lamda_right);

	double w_left = calculate_w(R_left, P_left, lamda_left, p_guess);	  
	double w_right = calculate_w(R_right, P_right, lamda_right, p_guess);

	// step 3: Calculate K and the B_c field in the intermediate state
	std::vector<double> K_left = calculate_K(R_left, P_left, lamda_left);
	std::vector<double> K_right = calculate_K(R_right, P_right, lamda_right);

	std::vector<double> B_c = calculate_B_c(R_left, P_left, lamda_left);

	double Y_left = calculate_Y_L(K_left, K_right, B_c);
	double Y_right = calculate_Y_R(K_left, K_right, B_c);

	double delta_K_x = K_right[0] - K_left[0];
	double f_of_p = (delta_K_x * (1 - P_left.B_x * (Y_right - Y_left))); // This is a placeholder for the actual calculation of f(p), we will replace this with the actual calculation later on
	return f_of_p;
}

std::vector<double> HLLD_solver(conserved_variables P_left, conserved_variables P_right)
{

	// Here we calculate the conservative variables for the left and right states
	std::vector<double> U_left = calculate_U(P_left);
	std::vector<double> U_right = calculate_U(P_right);

	std::vector<double> F_left = calculate_F(P_left);
	std::vector<double> F_right = calculate_F(P_right);

	// Here we calculate the wave speeds and the intermediate states
	double lamda_left = calculate_lamda_left();
	double lamda_right = calculate_lamda_right();

	std::vector<double> R_left = calculate_R(U_left, F_left, lamda_left);	  // temp lamda for now, we will replace this with the actual calculation later on
	std::vector<double> R_right = calculate_R(U_right, F_right, lamda_right); // temp lamda for now, we will replace this with the actual calculation later on

	double p_old, p_new = calculate_p_hll(P_left, P_right);
	double f_of_p_old = calculate_f_of_p(P_left, P_right, R_left, R_right, lamda_left, lamda_right, p_old);

	// This is a simple convergence criterion, we will replace this with the actual convergence criterion later on
	while (p_old - p_new > 1e-6) {
		
		// Here we use the secant method with f(p) and update p till our convergence criterion meets
		double f_of_p_new = calculate_f_of_p(P_left, P_right, R_left, R_right, lamda_left, lamda_right, p_new);
		p_new += f_of_p_new * (p_new - p_old) / (f_of_p_new - f_of_p_old); // This is a placeholder for the actual update of p, we will replace this with the actual calculation later on
		p_old = p_new;
		f_of_p_old = f_of_p_new;
	}

	std::vector<double> result;
	return result;
}
