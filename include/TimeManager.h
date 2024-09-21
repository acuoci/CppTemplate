/************************************************************************************
 * @file        TimeManager.h
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

#ifndef TimeManager_H
#define TimeManager_H

#include <vector>
#include "boost/filesystem.hpp"

/**
 * @class TimeManager
 * @brief Manages time-stepping and output intervals for simulations.
 *
 * The TimeManager class handles time-stepping operations and manages
 * output intervals for logging, mapping, and screen printing. It can
 * also stop the simulation based on time or step limits.
 */
class TimeManager
{
    public:

        /**
         * @brief Constructor for the TimeManager class.
         * Initializes time-stepping parameters and output intervals.
         */
        TimeManager();

        /**
         * @brief Sets the time window for the simulation.
         * @param tStart Start time of the simulation.
         * @param tEnd End time of the simulation.
         */
        void SetTimeWindow(const double tStart, const double tEnd);

        /**
         * @brief Sets the Courant number for the simulation.
         * @param Co Courant number to control time-stepping.
         */
        void SetCourantNumber(const double Co);

        /**
         * @brief Sets the maximum number of steps for the simulation.
         * @param stepsMax Maximum number of time steps allowed.
         */
        void SetStepsMax(const int stepsMax);

        /**
         * @brief Sets the print intervals for different output operations.
         *
         * @param dtPrintMap Interval for map printing.
         * @param dtPrintLog Interval for log printing.
         * @param dtPrintScreen Interval for screen printing.
         */
        void SetPrintIntervals(const double dtPrintMap, const double dtPrintLog, const double dtPrintScreen);

        /**
         * @brief Advances the simulation by a given time step.
         * @param dt The time step to advance.
         */
        void Advance(const double dt);

        /**
         * @brief Gets the current Courant number.
         * @return The Courant number.
         */
        double Co() const { return Co_; }

        /**
         * @brief Gets the current simulation time.
         * @return The current time.
         */
        double t() const { return t_; }

        /**
         * @brief Gets the current time step size.
         * @return The time step size.
         */
        double dt() const { return dt_; }

        /**
         * @brief Gets the current step number.
         * @return The current step.
         */
        int step() const { return step_; }

        /**
         * @brief Checks if the simulation should stop.
         * @return True if the simulation should stop, false otherwise.
         */
        bool stop() const { return stop_; }

        /**
         * @brief Checks if map printing is required.
         * @return True if map printing is needed, false otherwise.
         */
        bool printMap() const { return printMap_; }

        /**
         * @brief Checks if screen printing is required.
         * @return True if screen printing is needed, false otherwise.
         */
        bool printScreen() const { return printScreen_; }

        /**
         * @brief Checks if log printing is required.
         * @return True if log printing is needed, false otherwise.
         */
        bool printLog() const { return printLog_; }

        /**
         * @brief Gets the history of simulation times.
         * @return A vector containing the time history.
         */
        const std::vector<double>& tHistory() const { return tHistory_; }

    private:

        double tStart_;                 //!< Start time of the simulation
        double tEnd_;                   //!< End time of the simulation
        std::vector<double> tHistory_;  //!< History of time steps
        double Co_;                     //!< Courant number
        double tMax_;                   //!< Maximum time for a single step
        int stepsMax_;                  //!< Maximum number of steps
        bool stop_;                     //!< Flag to stop the simulation

        double dtPrintMap_;             //!< Interval for map output
        double dtPrintLog_;             //!< Interval for log output
        double dtPrintScreen_;          //!< Interval for screen output
    
        bool printMap_;                 //!< Flag to indicate map printing
        bool printLog_;                 //!< Flag to indicate log printing
        bool printScreen_;              //!< Flag to indicate screen printing

        int step_;                      //!< Current step number
        double t_;                      //!< Current simulation time
        double dt_;                     //!< Current time step size
};

#endif 