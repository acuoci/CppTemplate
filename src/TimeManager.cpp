/************************************************************************************
 * @file        TimeManager.cpp
 * @brief       The TimeManager class handles time-stepping operations and manages
 *              output intervals for logging, mapping, and screen printing. It can
 *              also stop the simulation based on time or step limits.
 * 
 * @author      Alberto Cuoci
 * @email       alberto.cuoci@polimi.it
 * @affiliation CRECK Modeling Lab (http://creckmodeling.chem.polimi.it/)
 *              Department of Chemistry, Materials, and Chemical Engineering
 *              Politecnico di Milano
 *              P.zza Leonardo da Vinci 32, 20133 Milano, Italy
 * @date        YYYY-MM-DD (automatically updated by Git)
 * @copyright   This file is part of the software written by Alberto Cuoci
 *              at CRECK Modeling Lab. All rights reserved.
 * @license     MIT 
 ************************************************************************************/

#include "TimeManager.h"
#include <iostream>
#include <cassert>

/**
 * @brief Constructor for the TimeManager class.
 * Initializes all parameters for time-stepping and output intervals to default values.
 */
TimeManager::TimeManager() :
    tStart_(0.),
    tEnd_(100.),
    Co_(0.25),
    tMax_(100.),
    stepsMax_(1000),
    stop_(false),
    dtPrintMap_(10.),
    dtPrintLog_(10.),
    dtPrintScreen_(10.),
    printMap_(false),
    printLog_(false),
    printScreen_(false),
    step_(0),
    t_(0.),
    dt_(1.)
{
    tHistory_.reserve(10000);       // Pre-allocate memory for the time history vector
    tHistory_.push_back(tStart_);   // Add the start time to the history
}

/**
 * @brief Sets the start and end time of the simulation.
 * Updates the time window and calculates the maximum step time.
 * @param tStart Start time of the simulation.
 * @param tEnd End time of the simulation.
 */
void TimeManager::SetTimeWindow(const double tStart, const double tEnd)
{
    assert( (tEnd-tStart > 0) && "The final time must be larger than the initial time");

    tStart_ = tStart;
    tEnd_ = tEnd;
    tMax_ = (tEnd_ - tStart_)/100.; // Set tMax as 1% of the time window
    tHistory_[0] = tStart_;         // Update the first entry of the time history
}

/**
 * @brief Sets the Courant number for controlling time-stepping.
 * @param Co Courant number.
 */
void TimeManager::SetCourantNumber(const double Co)
{
    assert(Co > 0. && "Courant number must be positive");

    Co_ = Co;
}

/**
 * @brief Sets the maximum number of steps allowed for the simulation.
 * @param stepsMax Maximum number of steps.
 */
void TimeManager::SetStepsMax(const int stepsMax)
{
    assert(stepsMax > 0 && "The maximum number of steps must be positive");

    stepsMax_ = stepsMax;
}

/**
 * @brief Sets the intervals for map, log, and screen printing.
 * @param dtPrintMap Interval for map output.
 * @param dtPrintLog Interval for log output.
 * @param dtPrintScreen Interval for screen output.
 */
void TimeManager::SetPrintIntervals(const double dtPrintMap, const double dtPrintLog, const double dtPrintScreen)
{
    assert(dtPrintMap > 0. && "The Interval for map output must be positive");
    assert(dtPrintLog > 0. && "The Interval for log output must be positive");
    assert(dtPrintScreen > 0. && "The Interval for screen output must be positive");
    
    dtPrintMap_ = dtPrintMap;
    dtPrintLog_ = dtPrintLog;
    dtPrintScreen_ = dtPrintScreen;    
}

/**
 * @brief Advances the simulation by a given time step and checks output conditions.
 * Identifies the next output time and updates internal state accordingly.
 * @param dt The time step to advance.
 */
void TimeManager::Advance(const double dt)
{
    // Identification of target time for output operations
    constexpr double eps = 1.e-6;
    const double tPrintMap = (std::floor((t_+eps)/dtPrintMap_)+1)*dtPrintMap_;
    const double tPrintLog = (std::floor((t_+eps)/dtPrintLog_)+1)*dtPrintLog_;
    const double tPrintScreen = (std::floor((t_+eps)/dtPrintScreen_)+1)*dtPrintScreen_;
    const double tTarget = std::min(tPrintMap,std::min(tPrintLog,std::min(tPrintScreen, tEnd_)));

    // Advance time
    step_++;
    const double told = t_;
    t_ = std::min(t_+dt, tTarget);
    dt_ = t_-told;  
    tHistory_.push_back(t_);

    // Identification of output targets
    printMap_ = false;
    printLog_ = false;
    printScreen_ = false;
    if ( std::fabs(t_-tPrintMap)/tPrintMap < eps)
        printMap_ = true;
    if ( std::fabs(t_-tPrintLog)/tPrintLog < eps)
        printLog_ = true;
    if ( std::fabs(t_-tPrintScreen)/tPrintScreen < eps)
        printScreen_ = true;  

    // Identification of stop conditions
    if (step_ >= stepsMax_ || t_ >= tEnd_)
        stop_ = true;
}