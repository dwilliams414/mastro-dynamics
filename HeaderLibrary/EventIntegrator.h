// Class Definition for EventIntegrator

#pragma once
#pragma warning(push)
#pragma warning(disable : 4996)
#include <boost/numeric/odeint.hpp>
#include "integrate_adaptive_event.h"
#include <vector>
#include <iostream>
#include "Event.h"
#pragma warning(pop)

// EventIntegrator Class Description:
// 
// The EventIntegrator Class allows for easy management of adaptive-step integration with event functions. All integration details
// are handled internally after the constructor is called.  
// 
// An EventIntegrator can be created with any uncontrolled Runge-Kutta stepper from the Boost odeint library.  When the constructor
// is called, the appropriate controlled stepper is created based on the provided abs_tol and rel_tol values.
//

template <typename UncontrolledRKStepper>
class EventIntegrator
{
public:
	// Type Definitions
	typedef boost::numeric::odeint::controlled_runge_kutta<UncontrolledRKStepper> StepTypeControl;
	typedef std::vector<double> state_type;

	// Member Variables
	StepTypeControl stepper; // Controlled Stepper
	double abs_tol;    // Absolute Tolerance for Integration
	double rel_tol;    // Relative tolerance for integration
	double init_step;  // Initial Step size
	Event& event_obj;   // Instance of an Event-type object
	EventObserver integration_results;  // Event observer object

	// Member Functions
	EventIntegrator(Event& event_object, double abs_tolerance, double rel_tolerance, double dt_init) : event_obj(event_object)
	{
		abs_tol = abs_tolerance;
		rel_tol = rel_tolerance;
		init_step = dt_init;
		stepper = boost::numeric::odeint::make_controlled(abs_tol, rel_tol, UncontrolledRKStepper());
	}

	// Integration Routines
	template <typename System> void integrate_adaptive(System sys, const state_type& x0, const state_type& tspan);
};

// Integration Function Definition
template <typename UncontrolledRKStepper>
template <typename System>
void EventIntegrator<UncontrolledRKStepper>::integrate_adaptive(System sys, const state_type& x0,
	const state_type& tspan)
{
	state_type x0_cpy = x0;
	integration_results.clear_integration();

	//integrate_adaptive_event_nochange(stepper, sys, x0_cpy, tspan.front(), tspan.back(), init_step, std::ref(eo));
	integrate_adaptive_event(stepper, sys, x0_cpy, tspan.front(), tspan.back(), init_step, std::ref(integration_results), event_obj);
}