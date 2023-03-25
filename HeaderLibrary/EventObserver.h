// Class definition for EventObserver.  This essentially works as a struct that saves
// integration results, including events, similar to what is returned by MATLAB from an integration
#pragma once
#include <vector>

// Class for Event Observation (Stores Integration Results)
class EventObserver
{
public:
	std::vector<double> t; // Integration times
	std::vector<std::vector<double>> x; // Integration States
	std::vector<double> te;  // Event Times
	std::vector<std::vector<double>> xe;  // Event states
	std::vector<int> ie;  // Event indices

	// Clear all saved vector data
	void clear_integration()
	{
		t.clear();
		x.clear();
		te.clear();
		xe.clear();
		ie.clear();
		ie.clear();
	}

	// Operator to Update Integration States (only)
	void operator()(const std::vector<double>& x_curr, double t_curr)
	{
		t.push_back(t_curr);
		x.push_back(x_curr);
	}
};