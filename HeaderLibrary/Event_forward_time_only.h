// Header file for definition of Event base Class
#pragma once
#include <vector>
#include <stdexcept>
#include <iostream>
#include "EventObserver.h"
#include <boost/numeric/odeint.hpp>
// #include "PrintUtils.h" // For debugging

/* General Helper Functions for Event */
// Return a Vector with sign of each of the elements in v
template <typename Numeric> std::vector<int> vsgn(const std::vector<Numeric>& v);

/*
* Event: Base class for defining an event
* 
* Event is an abstract class that contains the necessary functions for event finding. Subclasses
* of Event can be used to apply event finding to any of the Boost Runge-Kutta integrators.
* 
* Subclasses of Event must define the following three functions:
	* 1) event_fcn: Function that returns a vector corresponding to the value of the event function
	* 2) direction_fcn: Function that defines the direction of sign change for event finding
	* 3) terminate_fcn: Function that determines whether an event is terminal
* Note that event_fcn, direction_fcn, and terminate_fcn MUST return vectors of the same length.  
* An error will be thrown if this is not the case
* 
* This setup is somewhat similar to MATLAB's implementation of event functions where vectors event, isterminal, and direction
* are defined.  However, this can be limiting (e.g., it can be difficult to find the first n instances of an event).  As such,
* Event defines protected variables:
*		1) n_max_vals: Vector of maximum number of events of each type to find
*		2) n_curr_vals: Vector of current number of events found of each type
*		3) n_event_types: Number of distinct event types (size of n_max_vals and n_curr_vals)
* 
* Subclasses of Event can set n_max_vals and then use this value in defining their overwrite of event_fcn or terminate_fcn.
* n_max_vals must always be specified as a vector of integers in the constructor of a subclass of event.  However, this variable is not used
* UNLESS the terminate_fcn in the subclass uses it to define termination conditions.
* 
* Event uses a bisection algorithm to locate events between integration steps.  The bisection tolerance is defined by time.  By default,
* the tolerance is set to machine precision.  HOWEVER, this can be overwridden by specifying t_tol in the superclass constructor call of 
* a child class of Event.
* 
* For an example of a subclass of Event, see ApseEvent.h
*/

class Event
{
public:
	typedef std::vector<double> state_type; // Define type to reperesent states

	// Virtual Function to deine events such that event_fcn(t, state) = 0 at event.
	// Must return vector of length n_event_types
	virtual std::vector<double> event_fcn(const double t, const state_type& state) = 0;

	// Virtual termination specification function. Must return vector of length n_event_types
	virtual std::vector<int> terminate_fcn(const double t, const state_type& state) = 0;

	// Virtual direction specification function.  Must return vector of length n_event_types
	// Implementation is similar to MATLAB: -1 for decrease, 0 for either, +1 increase
	virtual std::vector<int> direction_fcn() = 0;

	// Print Protected Variables to Console
	void print_values();
	
	// Event finding function called within integration loop.
	// Also updates EventObserver eo by adding events eo.ie, eo.te, and eo.xe
	// Returns true if a terminal event is found.  Used in integration loop 
	template <class ControlledStepper, class System>
	int check_event_step(ControlledStepper stepper, const System& sys, const state_type& prev_state, const state_type& curr_state,
		double t_pre, double t_post, EventObserver& eo);

protected:

	// Variables
	std::vector<int> n_max_vals;  // Max number of events of each type
	std::vector<int> n_curr_vals; // Current number of events of each type found
	size_t n_event_types;			  // Number of event types to look for (size of n_max_vals, n_curr_vals)

	// Constructor
	Event(std::vector<int> n_max_vals)
	{
		this->n_max_vals = n_max_vals;
		this->n_event_types = n_max_vals.size();

		for (int i = 0; i < this->n_event_types; i++)
		{
			n_curr_vals.push_back(0);
		}

	}

private:
	typedef std::vector<int> IntVec;

	// Check if number of returned elements for each vector is correct. Used in check_event_step
	void check_sizing(const IntVec& events, const IntVec& directions, const IntVec& isterminal);

	// Check direction of event sign change.  Used in check_event_step
	int check_direction(int sgn_pre, int sgn_post, int dir);

	// Bisection algorithm to locate t_event in interval (t_pre, t_post).  Used in check_event_step
	template <class ControlledStepper, class System>
	double bisection_fcn(ControlledStepper stepper, System& sys, const state_type& prev_state, const state_type& curr_state,
		const double t_pre, const double t_post, int ie);
};


template <class ControlledStepper, class System>
int Event::check_event_step(ControlledStepper stepper, const System& sys, const state_type& prev_state, const state_type& curr_state,
	double t_pre, double t_post, EventObserver& eo)
{
	using namespace std;

	// Termination Flag
	bool terminate = false;

	// Define vectors for locating events
	vector<int> sgn_pre = vsgn(event_fcn(t_pre, prev_state));
	vector<int> sgn_post = vsgn(event_fcn(t_post, curr_state));
	vector<int> direction = direction_fcn();
	vector<int> isterminal = terminate_fcn(t_pre, prev_state); // For sizing.  Update again later after anything incremented

	// Ensure that Event Function is Valid
	check_sizing(sgn_pre, direction, isterminal);

	// Iterate Through To Check if Event Exists within Interval.  If yes, mark and save the
	// time - index pair in time_and_index
	vector<pair<double, int>> time_and_index;

	for (int ie = 0; ie < n_event_types; ie++)
	{
		// Event Has occurred in interval
		if ((sgn_pre[ie] != sgn_post[ie] &&
			check_direction(sgn_pre[ie], sgn_post[ie], direction[ie])) ||
			(sgn_pre[ie] == 0 && check_direction(sgn_pre[ie], sgn_post[ie], direction[ie])))
		{
			//std::cout << "Event Found!!!" << std::endl; // For debugging
			n_curr_vals.at(ie)++; // Increment number of found events for ie
			if (sgn_pre[ie] == 0)
			{
				//std::cout << "No Bisection Necessary!!" << std::endl // For debugging;
				time_and_index.push_back(make_pair(t_pre, ie));
			}
			else
			{
				//std::cout << "Bisection Required!!!" << std::endl; // For debugging
				double ie_time = bisection_fcn(stepper, sys, prev_state, curr_state, t_pre, t_post, ie);
				time_and_index.push_back(make_pair(ie_time, ie));
			}
		}
	}

	// Sort Time-Event Pair By Time.  This way, if an event is terminal we don't
	// save events after the termination time
	if (!time_and_index.empty())
	{
		sort(time_and_index.begin(), time_and_index.end());
		//print_pair_vec(time_and_index); // For debugging
	}

	// Iterate Through Time + Event Pair and Step to Event Time
	// Break from loop if a terminal event is found
	for (int i = 0; i < time_and_index.size(); i++)
	{
		// Specify Integration Params
		double t_event = time_and_index[i].first;
		double dt_event = time_and_index[i].first - t_pre;
		state_type xe_step = prev_state;

		// Step to event location
		boost::numeric::odeint::integrate_adaptive(stepper, sys, xe_step, t_pre, t_event, dt_event);

		// Save Event Results
		eo.ie.push_back(time_and_index[i].second + 1); // Add 1 for consistency with MATLAB event numbering
		eo.te.push_back(t_event);
		eo.xe.push_back(xe_step);

		// Check if Event is Terminal
		isterminal = terminate_fcn(t_event, xe_step); // Recheck terminate_fcn() with updated n_curr_events and time

		// If event is terminal, break and save results into the regular t and x vectors
		if (isterminal[time_and_index[i].second])
		{
			terminate = true;
			eo.t.push_back(t_event); // Becomes final time in regular integration results
			eo.x.push_back(xe_step); // Becomes final state in regular integration results
			break;
		}
	}

	return terminate;
}

// Check to make sure that all returned resultls from each _fcn are the same (and correct) size
void Event::check_sizing(const IntVec& events, const IntVec& directions, const IntVec& isterminal)
{
	if (events.size() != n_event_types || directions.size() != n_event_types || isterminal.size() != n_event_types)
	{
		throw std::invalid_argument("This event function has an improperely defined event_fcn, terminate_fcn, or direction_fcn");
	}
}

// Check to see if direction is correct
int Event::check_direction(int sgn_pre, int sgn_post, int dir)
{
	int iscorrect = 0;

	if (sgn_post > sgn_pre && dir >= 0)
	{
		iscorrect = 1;
	}

	if (sgn_post < sgn_pre && dir <= 0)
	{
		iscorrect = 1;
	}



	return iscorrect;
}

// Bisection Function in time to locate event
template <class ControlledStepper, class System>
double Event::bisection_fcn(ControlledStepper stepper, System& sys, const state_type& prev_state, const state_type& curr_state,
	const double t_pre, const double t_post, int ie)
{
	using namespace std;
	double t_event = t_pre; // Value to return (finding time of event)

	double t_a = t_pre; // Left Bracketing Time
	double t_b = t_post; // Right Bracketing Time
	double t_m = t_pre + (t_post - t_pre) / 2.0; // Midpoint time

	state_type a_state = prev_state; // Left Bracketing State
	state_type b_state = curr_state; // Right Bracketing State

	// Iterate Bisection.  Goes to machine precision for size of t_a and t_b
	// This means accuracy technically decreases as t_a and t_b become larger
	while ((t_a < (t_b+t_a)/2.0) && ((t_b+t_a)/2.0 < t_b))
	{
		state_type m_state = a_state; // Initial state is at a
		double dt_i = t_m - t_a; // Step for this bisection

		// Step from t_a -> t_m
		boost::numeric::odeint::integrate_adaptive(stepper, sys, m_state, t_a, t_m, dt_i);

		// Compute Signs
		vector<int> sgn_a = vsgn(event_fcn(t_a, a_state));
		vector<int> sgn_b = vsgn(event_fcn(t_b, b_state));
		vector<int> sgn_m = vsgn(event_fcn(t_m, m_state));

		// Changing Bounds (std. Bisection Method)
		if (sgn_a[ie] * sgn_m[ie] < 0)
		{
			t_b = t_m;
			b_state = m_state;
		}
		else if (sgn_a[ie] == 0)
		{
			t_event = t_a;
			break;
		}
		else
		{
			t_a = t_m;
			a_state = m_state; // So that we now integrate from new a_state
		}

		t_m = t_a + (t_b - t_a) / 2.0; // Update t_m
		t_event = t_m; // Save t_m in t_event (value to return)
	}
	return t_event;

}

// Print protected variables to terminal
void Event::print_values()
{
	using namespace std;
	cout << "n_event_types: ";
	cout << n_event_types << endl;
	cout << "" << endl;

	// Print max vals
	cout << "n_max_vals: " << endl;
	cout << "[";
	for (int i : n_max_vals)
	{
		cout << i;
		cout << " ";
	}
	cout << "]" << endl;

	// Print current vals
	cout << "n_curr_vals: " << endl;
	cout << "[";
	for (int i : n_curr_vals)
	{
		cout << i;
		cout << " ";
	}
	cout << "]" << endl;


}

// Get a Vector with sign of each of original elements
template <typename Numeric> std::vector<int> vsgn(const std::vector<Numeric>& v)
{
	std::vector<int> sgns;

	for (Numeric i : v)
	{
		if (i > 0)
		{
			sgns.push_back(1);
		}
		else if (i < 0)
		{
			sgns.push_back(-1);
		}
		else
		{
			sgns.push_back(0);
		}
	}
	return sgns;
}