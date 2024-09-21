/************************************************************************************
 * @file        main.cpp
 * @brief       Wildfire Simulator based on the Cellular Automata technique
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
#include "RectangularMap2D.h"
#include "KarafyllidisModel.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

/**
 * @brief Main function for simulating fire spread using the Karafyllidis model.
 * 
 * This program initializes a 2D rectangular map, configures the time manager and simulation
 * parameters, and runs the fire spread simulation based on the input terrain features,
 * elevation map, and wind conditions. The output is generated in VTK format for visualization.
 * 
 * @return int Returns 0 upon successful completion.
 */
int main(int argc, char* argv[])
{
    // ---------------------------------------------------------------------------------------------------------- //
    // Handle command-line options
    // ---------------------------------------------------------------------------------------------------------- //

    // Default values for input variables
    std::string mapFile = "map.xml";          // Default map file
    std::string outFolder = "Output";         // Default output folder name
    double tEnd = 1000.0;                     // Default final time (s)
    double dtMap = 10.0;                      // Default interval for map output (s)
    double dtLog = 10.0;                      // Default interval for log output (s)
    double dtScreen = 1.0;                    // Default interval for screen output (s)
    double Co = 0.25;                         // Default Courant number

    try
    {
        // Define available options
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "Produce help message")
            ("map-file", po::value<std::string>(&mapFile)->default_value(mapFile), "Path to the XML map file")
            ("output", po::value<std::string>(&outFolder)->default_value(outFolder), "Name of the output folder")
            ("tEnd", po::value<double>(&tEnd)->default_value(tEnd), "Final time of the simulation (s)")
            ("dtMap", po::value<double>(&dtMap)->default_value(dtMap), "Time interval for map output (s)")
            ("dtLog", po::value<double>(&dtLog)->default_value(dtLog), "Time interval for log output (s)")
            ("dtScreen", po::value<double>(&dtScreen)->default_value(dtScreen), "Time interval for screen output (s)")
            ("Co", po::value<double>(&Co)->default_value(Co), "Courant number for the simulation");

        // Parse the command-line arguments
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        // If help option is provided, display the help message and exit
        if (vm.count("help")) 
        {
            std::cout << desc << "\n";
            return 0;
        }
    } 
    catch (const po::error &ex) 
    {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

// ---------------------------------------------------------------------------------------------------------- //
    // Ensure the output folder doesn't already exist, return an error if it does
    // ---------------------------------------------------------------------------------------------------------- //
    boost::filesystem::path outputPath(outFolder);
    if (boost::filesystem::exists(outputPath)) 
    {
        std::cerr << "Error: The output folder '" << outFolder << "' already exists. "
                  << "Please remove the folder or specify a different output folder." << std::endl;
        return 1;
    } 
    else 
    {
        // Create the output folder if it doesn't exist
        boost::filesystem::create_directories(outputPath);
    }

    // ---------------------------------------------------------------------------------------------------------- //
    // Time manager
    // ---------------------------------------------------------------------------------------------------------- //
    
    /// @brief Initializes the time manager to control the simulation.
    TimeManager runTime;
    
    // Set initial and final times (in seconds).
    runTime.SetTimeWindow(0., tEnd);           // Initial and final times (s)
    
    // Set the Courant number for determining the time step.
    runTime.SetCourantNumber(Co);             // Courant number
    
    // Set the maximum number of steps for the simulation.
    runTime.SetStepsMax(10000);               // Maximum number of steps
    
    // Set time intervals (in seconds) for printing output: maps, log file, screen outputs.
    runTime.SetPrintIntervals(dtMap, dtLog, dtScreen);    // Time steps for output operations: maps, log-file, screen


    // ---------------------------------------------------------------------------------------------------------- //
    // Mesh definition
    // ---------------------------------------------------------------------------------------------------------- //
    
    // Print the name of the XML map file being used
    std::cout << "Using map file: " << mapFile << std::endl;

    /// @brief Initializes a rectangular 2D map based on an input XML file.
    RectangularMap2D map(mapFile);


    // ---------------------------------------------------------------------------------------------------------- //
    // Initialize model (Karafyllidis)
    // ---------------------------------------------------------------------------------------------------------- //
    
    /// @brief Initializes the Karafyllidis fire spread model with the defined map and time manager.
    KarafyllidisModel model(runTime, map.nx(), map.ny(), map.deltax(), map.deltay());

    // Set the output log file for the model.
    model.SetLogFile((outputPath / "Solution.out").string());
    
    // Set the folder for VTK output files.
    model.SetVTKFolder(outputPath / "VTK");

    // Set terrain features maps (land use and vegetation density).
    model.SetTerrainFeatures(map.landUse(), map.density());

    // Set elevation map from the terrain data.
    model.SetElevation(map.h());

    // Set wind map (horizontal and vertical wind components).
    model.SetWind(map.wx(), map.wy());


    // ---------------------------------------------------------------------------------------------------------- //
    // Set ignition point(s)
    // ---------------------------------------------------------------------------------------------------------- //
    
    /// @brief Sets an initial ignition point at the center of the terrain.
    model.SetIgnitionCell(map.nx()/2, map.ny()/2);


    // ---------------------------------------------------------------------------------------------------------- //
    // Solve the model
    // ---------------------------------------------------------------------------------------------------------- //
    
    // Calculate the fire spread rate based on the input parameters.
    model.CalculateSpreadRate();

    /// @brief Main simulation loop.
    // Loop over time until the stop condition is met.
    for (;;)
    {
        // Advance the simulation by one time step.
        model.Advance();

        // Check if the stop condition is met
        if (runTime.stop() == true)
            break;
    }

    return 0;
}