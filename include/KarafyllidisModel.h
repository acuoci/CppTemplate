/************************************************************************************
 * @file        KarafyllidisModel.h
 * @brief       Implementation of a class representing a 2D simulation model based 
 *              on the Karafyllidis method for fire spread modeling.
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

#ifndef KarafyllidisModel_H
#define KarafyllidisModel_H

#include <Eigen/Dense>
#include "boost/filesystem.hpp"
#include <boost/filesystem/fstream.hpp>
#include "TimeManager.h"
#include "Utilities.h"

/**
 * @class KarafyllidisModel
 * @brief A class representing a 2D simulation model based on the Karafyllidis method for fire spread modeling.
 */
class KarafyllidisModel
{
    public:

        /// Enumeration for various land use types in the model.
        enum LandUseType {  
            WATER = 0, 
            BUILDINGS = 1, 
            ROCKS = 2, 
            BROADLEAVES = 10,
            SHRUBS = 11, 
            GRASSLAND = 12, 
            TROPICALGRASSLAND = 13, 
            MEDITERRANEANWOODLAND = 14, 
            FIREPRONECONIFERS = 20, 
            AGROFORESTRY = 21, 
            NOTFIREPRONEFORESTRY = 22, 
            ALPINEFOREST = 23, 
            SAVANNA = 30, 
            TAIGA = 31, 
            TUNDRA = 32, 
            PEATLANDS = 40
        };

        /**
         * @brief Constructor: initialize a 2D map (the ghost layers are included in nx and ny)
         * @param runTime Reference to a TimeManager instance.
         * @param nx Number of cells along the horizontal direction (from West to East)
         * @param ny Number of cells along the vertical direction (from North to South)
         * @param deltax Size of cells along the horizontal direction (in meters)
         * @param deltay Size of cells along the vertical direction (in meters)
         */
        KarafyllidisModel(TimeManager& runTime, const int nx, const int ny, const double deltax, const double deltay);

        /**
         * @brief Set the log file name.
         * @param name Name of the log file.
         */
        void SetLogFile(const std::string name);

        /**
         * @brief Set the folder where Tecplot files will be written.
         * @param folder Path to the Tecplot folder.
         */
        void SetTecplotFolder(const boost::filesystem::path folder);

        /**
         * @brief Set the folder where VTK files will be written.
         * @param folder Path to the VTK folder.
         */
        void SetVTKFolder(const boost::filesystem::path folder);

        /**
         * @brief Set a single ignition cell based on its indices.
         * @param ix Index of the ignition cell along the horizontal direction.
         * @param iy Index of the ignition cell along the vertical direction.
         */
        void SetIgnitionCell(const int ix, const int iy);

        /**
         * @brief Set a single ignition cell based on physical coordinates.
         * @param x X-coordinate of the ignition point.
         * @param y Y-coordinate of the ignition point.
         */
        void SetIgnitionCell(const double x, const double y);

        /**
         * @brief Set multiple ignition cells using a matrix of indices.
         * @param indices Matrix of ignition cell indices.
         */
        void SetIgnitionCells(const Eigen::MatrixXi& indices);

        /**
         * @brief Set multiple ignition cells within a rectangular area.
         * @param xc Center x-coordinate of the rectangle.
         * @param yc Center y-coordinate of the rectangle.
         * @param dx Width of the rectangle.
         * @param dy Height of the rectangle.
         */
        void SetIgnitionCells(const double xc, const double yc, const double dx, const double dy);

        /**
         * @brief Set terrain features such as land use and vegetation density.
         * @param landUse Matrix representing the land use types.
         * @param density Matrix representing the vegetation density.
         */
        void SetTerrainFeatures(const Eigen::MatrixXi& landUse, const Eigen::MatrixXd& density);

        /**
         * @brief Set the elevation map.
         * @param h Matrix representing the terrain elevation (in meters).
         */
        void SetElevation(const Eigen::MatrixXd& h);

        /**
         * @brief Set the wind vector fields for the simulation.
         * @param wx Matrix representing the wind speed in the horizontal direction.
         * @param wy Matrix representing the wind speed in the vertical direction.
         */
        void SetWind(const Eigen::MatrixXd& wx, const Eigen::MatrixXd& wy);

        /**
         * @brief Calculate the fire spread rate based on terrain and wind conditions.
         */
        void CalculateSpreadRate();

        /**
         * @brief Advance the simulation by one timestep.
         */
        void Advance();

        /**
         * @brief Write the simulation log to a file.
         */
        void PrintLogFile();


    private:

        /**
         * @brief Extract a 3x3 submatrix of the elevation map centered at cell (i, j).
         * @param i Row index of the center cell.
         * @param j Column index of the center cell.
         * @return A 3x3 matrix of elevation data.
         */
        Eigen::Matrix<double,3,3> ElevationSubMatrix(const int i, const int j);

        /**
         * @brief Extract a 3x3 submatrix of the wind map centered at cell (i, j).
         * @param i Row index of the center cell.
         * @param j Column index of the center cell.
         * @return A 3x3 matrix of wind data.
         */
        Eigen::Matrix<double,3,3> WindSubMatrix(const int i, const int j);

    private:

        TimeManager& runTime_;           //!< Reference to time manager

        Eigen::MatrixXd a_;              //!< State map (values from 0 = unburnt to 1 = burnt)

        int nx_;                         //!< Number of cells (including ghost layers) along the horizontal direction  
        int ny_;                         //!< Number of cells (including ghost layers) along the vertical direction
        double deltax_;                  //!< Size of cells along the horizontal direction (in meters)
        double deltay_;                  //!< Size of cells along the vertical direction (in meters)

        double lengthx_;                 //!< Total length of the physical domain (without ghost layers) (in meters)
        double lengthy_;                 //!< Total length of the physical domain (without ghost layers) (in meters)
        double areatot_;                 //!< Total area of physical domain (without ghost layers) (in square meters)
        
        double areacell_;                //!< Area of a single cell (in square meters)

        int ncells_;                     //!< Number of physical cells (without ghost layers)
        int ncellstot_;                  //!< Total number of cells (including ghost layers)

        Eigen::MatrixXd h_;              //!< Elevation map (in meters)
        Eigen::MatrixXd wx_;             //!< Wind map (horizontal component) (in meters per second)
        Eigen::MatrixXd wy_;             //!< Wind map (vertical component) (in meters per second)
        Eigen::MatrixXd density_;        //!< Vegetation density map (dimensionless, from 0 to +Inf, 1 means regular density)
        Eigen::MatrixXi landUse_;        //!< Land use map (see LandUseType)

        Eigen::MatrixXd R_;                             //!< Fire spread rate map (in meters per second)

        std::vector<Eigen::Vector<double,4>> hadj_;     //!< Local (adjacent cells) elevation map correction coefficient
        std::vector<Eigen::Vector<double,4>> hdiag_;    //!< Local (diagonal cells) elevation map correction coefficient

        std::vector<Eigen::Vector<double,4>> wadj_;     //!< Local (adjacent cells) wind map correction coefficient
        std::vector<Eigen::Vector<double,4>> wdiag_;    //!< Local (diagonal cells) wind map correction coefficient

        boost::filesystem::ofstream fLog_;              //!< Log file stream
        boost::filesystem::path tecplotFolder_;         //!< Folder path for Tecplot output
        boost::filesystem::path vtkFolder_;             //!< Folder path for VTK output
        bool tecplotOutput_;                            //!< Flag to enable/disable Tecplot output
        bool vtkOutput_;                                //!< Flag to enable/disable VTK output
};

#endif