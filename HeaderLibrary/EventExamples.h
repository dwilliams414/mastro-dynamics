#pragma once
#include "Event.h"

/* Example Event: Locate the First n Apses along a trajectory */
class ApseEvent : public Event
{
public:
	// Constructor
	ApseEvent(int n_apses) : Event(std::vector<int>(1, n_apses)) {}

	// Event Function
	state_type event_fcn(const double t, const state_type& state) override
	{
		state_type event = { state[0] * state[3] + state[1] * state[4] + state[2] * state[5] };
		return event;
	}

	// Terminate Function
	std::vector<int> terminate_fcn(double t) override
	{
		std::vector<int> terminate = { 0 };
		if (n_curr_vals[0] >= n_max_vals[0])
		{
			terminate[0] = 1;
		}
		return terminate;
	}

	// Direction Function
	std::vector<int> direction_fcn() override
	{
		std::vector<int> direction = { 0 };
		return direction;
	}
};