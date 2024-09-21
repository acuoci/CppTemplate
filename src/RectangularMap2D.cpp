/************************************************************************************
 * @file        RectangularMap2D.cpp
 * @brief       A class representing a rectangular 2D map, including elevation, wind, vegetation density, and land use.
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

#include "RectangularMap2D.h"
#include "Utilities.h"
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

RectangularMap2D::RectangularMap2D(const boost::filesystem::path &xmlfile)
{
    // Create a property tree object to hold the XML data
    boost::property_tree::ptree ptree;

    // Load the XML file into the property tree
    boost::property_tree::read_xml(xmlfile.string(), ptree);

    // Read the values from the XML using the property tree
    try 
    {
        // Extract the values from the <map> node
        nx_ = ptree.get<int>("map.nx");
        ny_ = ptree.get<int>("map.ny");
        deltax_ = ptree.get<double>("map.deltax");
        deltay_ = ptree.get<double>("map.deltay");
    } 
    catch (const boost::property_tree::ptree_error& e) 
    {
        std::cerr << "Error parsing the XML file: " << e.what() << std::endl;
        throw; // Re-throw the exception after logging it
    }

    // Validate nx and ny
    if (nx_ <= 2 || ny_ <= 2)
        throw std::invalid_argument("Invalid values for nx or ny. Must be greater than 2.");

    // Processing data
    lengthx_ = (nx_-2)*deltax_;
    lengthy_ = (ny_-2)*deltay_;
    areatot_ = lengthx_*lengthy_;
    areacell_ = deltax_*deltay_;
    ncells_ = nx_*ny_;

    // Memory allocation
    h_.setZero(nx_, ny_);
    wx_.setZero(nx_, ny_);
    wy_.setZero(nx_, ny_);
    density_.setConstant(nx_, ny_, 1.);
    landUse_.setConstant(nx_, ny_, 0);

    // Reading elevation
    std::cout << "* Reading elevation from xml file..." << std::endl;
    ImportXMLMap(ptree, "map.elevation", h_);

    // Reading wind
    std::cout << "* Reading wind (x component) from xml file..." << std::endl;
    ImportXMLMap(ptree, "map.windx", wx_);
    ImportXMLMap(ptree, "map.windy", wy_);

    // Reading vegetation density
    std::cout << "* Reading vegetation density from xml file..." << std::endl;
    ImportXMLMap(ptree, "map.density", density_);

    // Reading land use
    std::cout << "* Reading land use from xml file..." << std::endl;
    ImportXMLMap(ptree, "map.landuse", landUse_);
}





