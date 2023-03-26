// This file defines the Adaptive-Step integration routines that allow for event finding.  
// These functions are slightly modified versions of the Boost integrate_adaptive functions, but
// implement event finding.
//
// Note: To use this header library, you MUST specify the locations of your BOOST_ROOT directory
#pragma once

#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/integrate/max_step_checker.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_const.hpp>
#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <iostream>
#include "Event.h"

/*
 * Adaptive-Step integration routine with event-finding added.  This is the actual loop
 * that executes the integration.
 */
template< class Stepper, class System, class State, class Time>
size_t integrate_adaptive_event_loop(
    Stepper stepper, System system, State& start_state,
    Time& start_time, Time end_time, Time& dt,
    EventObserver& observer, Event& event_obj, boost::numeric::odeint::controlled_stepper_tag
)
{
    using namespace boost::numeric;
    using namespace boost::numeric::odeint;
    using namespace boost::numeric::odeint::detail;

    EventObserver& obs = observer;
    typename odeint::unwrap_reference< Stepper >::type& st = stepper;

    // Copy Inputs for Readability in Loop w/ Error Checking.  
    State state_post_step = start_state;
    Time time_post_step = start_time;

    // Pre Step States and Times
    State state_pre_step;
    Time time_pre_step;

    // Termination Flag
    int event_terminate = 0;

    failed_step_checker fail_checker;  // to throw a runtime_error if step size adjustment fails
    size_t count = 0;
    while (less_with_sign(time_post_step, end_time, dt))
    {
        obs(state_post_step, time_post_step);
        if (less_with_sign(end_time, static_cast<Time>(time_post_step + dt), dt))
        {
            dt = end_time - time_post_step;
        }

        state_pre_step = state_post_step;
        time_pre_step = time_post_step;
        controlled_step_result res;
        do
        {
            res = st.try_step(system, state_post_step, time_post_step, dt);
            fail_checker();  // check number of failed steps
        } while (res == fail);
        fail_checker.reset();  // if we reach here, the step was successful -> reset fail checker

        // Implement Event Checking Code Here
        event_terminate = event_obj.check_event_step(st, system, state_pre_step, state_post_step, time_pre_step, time_post_step, obs);
        if (event_terminate) { break; }
        ++count;
    }
    if (!event_terminate) // Terminal time handled in check_event_step if termination due to event
    {
        obs(state_post_step, time_post_step); // Record final observation if final time reached
    }
    // For classic functionality with boost, copy back to original variables
    start_state = state_post_step;
    start_time = time_post_step;
    return count;
}

/*
 * Adaptive-Step integration routine with event-finding added.  This is the function that
 * should be called by any code that you implement yourself (so you don't have to deal with the 
 * controlled_stepper_tag stuff in the above).
 */
template< class Stepper, class System, class State, class Time>
size_t integrate_adaptive_event(
    Stepper stepper, System system, State& start_state,
    Time start_time, Time end_time, Time dt,
    EventObserver& observer, Event& event_obj)
{
    using namespace boost::numeric;
    using namespace boost::numeric::odeint;

    typedef typename odeint::unwrap_reference< Stepper >::type::stepper_category stepper_category;
    return integrate_adaptive_event_loop(
        stepper, system, start_state,
        start_time, end_time, dt,
        observer, event_obj, stepper_category());
}