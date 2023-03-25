#pragma once
#include <vector>
/*CR3BP: Class definition used as container for system mass parameter, equations of motion, and other CR3BP-type stuff.
* You can use this as (one possible) template for your own code and expand upon it as necessary.
*/
typedef std::vector<double> state_type;

class CR3BP
{
public:
	// Member Variables
	double m_mu; // System Mass Parameter

	// Constructor Declaration
	CR3BP(double mu);

	// Declare Operator for Class (Direct Call to EOM's).  This allows us to 
	// easily pass the mass parameter when we integrate, since we can just pass
	// the CR3BP object.
	void operator()(const state_type& state, state_type& state_dot, const double t);

	// CR3BP EOM's, No Variational Equations
	void eom6(const state_type& state, state_type& state_dot, const double t);

	// CR3BP EOM's, w/ Variational Equations (1st Order)
	void eom42(const state_type& state, state_type& state_dot, const double t);

	// Evaluate Distance r to second primary
	double eval_r(const state_type& state);

	// Evaulate the Distance d from the first primary (larger)
	double eval_d(const state_type& state);

	// Evaluate uxx pseudopotential term
	double eval_uxx(const state_type& state);

	// Evaluate uxy pseudopotential term
	double eval_uxy(const state_type& state);

	// Evaluate uxz pseudopotential term
	double eval_uxz(const state_type& state);

	// Evaluate the uyy pseudopotential term
	double eval_uyy(const state_type& state);

	// Evaluate the uzz pseudopotential term
	double eval_uzz(const state_type& state);

	// Evaluate uyz pseudopotential term
	double eval_uyz(const state_type& state);

	// Evaluate the uyx Term
	double eval_uyx(const state_type& state);

	// Evaluate the uzx Term
	double eval_uzx(const state_type& state);

	// Evaluate the uzy Term
	double eval_uzy(const state_type& state);

	// Evaluate the jacobian A Matrix at the current state
	state_type A_mat(const state_type& state);
};

// Define CR3BP Constructor
CR3BP::CR3BP(double mu)
{
	m_mu = mu;
}

// Define CR3BP Operator
void CR3BP::operator()(const state_type& state, state_type& state_dot, const double t)
{
	if (state.size() == 6)
	{
		eom6(state, state_dot, t);
	}
	else if (state.size() == 42)
	{
		eom42(state, state_dot, t);
	}
}

// Define CR3BP Equations of Motion (State Only)
void CR3BP::eom6(const state_type& state, state_type& state_dot, const double t)
{
	// Unpack state for readability
	double x = state[0];
	double y = state[1];
	double z = state[2];
	double xd = state[3];
	double yd = state[4];
	double zd = state[5];

	// Define Distances from Primaries
	double d = sqrt(pow(x + m_mu, 2) + pow(y, 2) + pow(z, 2)); // Distance from larger primary
	double r = sqrt(pow(x - 1 + m_mu, 2) + pow(y, 2) + pow(z, 2)); // Distance from smaller primary

	// Equations of Motion
	state_dot[0] = xd;
	state_dot[1] = yd;
	state_dot[2] = zd;

	state_dot[3] = 2 * yd + x - ((1 - m_mu) * (x + m_mu)) / pow(d, 3) - (m_mu * (x - 1 + m_mu)) / pow(r, 3);
	state_dot[4] = -2 * xd + y - ((1 - m_mu) * y) / pow(d, 3) - (m_mu * y) / pow(r, 3);
	state_dot[5] = -((1 - m_mu) * z) / pow(d, 3) - (m_mu * z) / pow(r, 3);
}

// Define CR3BP Equations of Motion (w/ 1st Order Variations)
void CR3BP::eom42(const state_type& state, state_type& state_dot, const double t)
{
	// Initialize State Derivatives
	eom6(state, state_dot, t);

	// Get A Matrix
	state_type A = A_mat(state); // Row major layout


	// Linear Variational Derivatives
	int aug_state_index; // Index for current state element in state_dot
	int state_size6 = 6; // Size of pos+vel state

	// Perform multiplication Phidot = APhi, using vectors & linear indices
	// Elements 6:end of state_dot are in COLUMN MAJOR order
	for (int j = 0; j < state_size6; j++)
	{
		for (int i = 0; i < state_size6; i++)
		{
			aug_state_index = state_size6 + j * state_size6 + i;
			state_dot[aug_state_index] = 0.0;

			for (int k = 0; k < state_size6; k++)
			{
				state_dot[aug_state_index] += A[i * state_size6 + k] * state[state_size6 + j * state_size6 + k];
			}
		}
	}
}

/******************************** CR3BP Additional Member Function Definitions ***************************************/

// Compute distance from second primary
double CR3BP::eval_r(const state_type& state)
{
	return sqrt(pow(state[0] - 1 + m_mu, 2) + pow(state[1], 2) + pow(state[2], 2)); // Distance from smaller primary
}

// Compute distance from first primary
double CR3BP::eval_d(const state_type& state)
{
	return sqrt(pow(state[0] + m_mu, 2) + pow(state[1], 2) + pow(state[2], 2));
}

// Compute Pseudopotential Terms:
double CR3BP::eval_uxx(const state_type& state)
{
	double r = eval_r(state);
	double d = eval_d(state);
	double x = state[0];
	double y = state[1];
	double z = state[2];

	double uxx = 1 - (1 - m_mu) / pow(d, 3) - m_mu / pow(r, 3);
	uxx += (3 * (1 - m_mu) * pow(x + m_mu, 2)) / pow(d, 5);
	uxx += (3 * m_mu * pow(x - 1 + m_mu, 2)) / pow(r, 5);
	return uxx;
}

double CR3BP::eval_uxy(const state_type& state)
{
	double r = eval_r(state);
	double d = eval_d(state);
	double x = state[0];
	double y = state[1];
	double z = state[2];

	double uxy = (3 * (1 - m_mu) * (x + m_mu) * y) / pow(d, 5);
	uxy += (3 * m_mu * (x - 1 + m_mu) * y) / pow(r, 5);

	return uxy;
}

double CR3BP::eval_uxz(const state_type& state)
{
	double r = eval_r(state);
	double d = eval_d(state);
	double x = state[0];
	double y = state[1];
	double z = state[2];

	double uxz = (3 * (1 - m_mu) * (x + m_mu) * z) / pow(d, 5);
	uxz += (3 * m_mu * (x - 1 + m_mu) * z) / pow(r, 5);
	return uxz;
}

double CR3BP::eval_uyy(const state_type& state)
{
	double r = eval_r(state);
	double d = eval_d(state);
	double x = state[0];
	double y = state[1];
	double z = state[2];

	double uyy = 1 - (1 - m_mu) / pow(d, 3) - m_mu / pow(r, 3);
	uyy += (3 * (1 - m_mu) * pow(y, 2)) / pow(d, 5);
	uyy += 3 * m_mu * pow(y, 2) / pow(r, 5);

	return uyy;
}

double CR3BP::eval_uyz(const state_type& state)
{
	double r = eval_r(state);
	double d = eval_d(state);
	double x = state[0];
	double y = state[1];
	double z = state[2];

	double uyz = (3 * (1 - m_mu) * y * z) / pow(d, 5);
	uyz += 3 * m_mu * y * z / pow(r, 5);

	return uyz;
}

double CR3BP::eval_uzz(const state_type& state)
{
	double r = eval_r(state);
	double d = eval_d(state);
	double x = state[0];
	double y = state[1];
	double z = state[2];

	double uzz = -(1 - m_mu) / pow(d, 3) - m_mu / pow(r, 3);
	uzz += (3 * (1 - m_mu) * pow(z, 2)) / pow(d, 5);
	uzz += 3 * m_mu * pow(z, 2) / pow(r, 5);

	return uzz;
}

double CR3BP::eval_uyx(const state_type& state)
{
	return eval_uxy(state);
}

double CR3BP::eval_uzx(const state_type& state)
{
	return eval_uxz(state);
}

double CR3BP::eval_uzy(const state_type& state)
{
	return eval_uyz(state);
}

state_type CR3BP::A_mat(const state_type& state)
{
	state_type A(36, 0);

	// Position Derivatives
	A[3] = 1.0;
	A[10] = 1.0;
	A[17] = 1.0;

	// Pseudopotential Terms
	A[18] = eval_uxx(state);
	A[19] = eval_uxy(state);
	A[20] = eval_uxz(state);

	A[24] = eval_uyx(state);
	A[25] = eval_uyy(state);
	A[26] = eval_uyz(state);

	A[30] = eval_uzx(state);
	A[31] = eval_uzy(state);
	A[32] = eval_uzz(state);

	// Velocity Terms
	A[22] = 2.0;
	A[27] = -2.0;

	return A;
}