/************************************************************************************
 * @file        RectangularMap2D.h
 * @brief       A class representing a rectangular 2D map, including elevation, wind, 
 *              vegetation density, and land use.
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

#ifndef RectangularMap2D_H
#define RectangularMap2D_H

#include <Eigen/Dense>
#include "boost/filesystem.hpp"

/**
 * @class RectangularMap2D
 * @brief A class representing a rectangular 2D map, including elevation, wind, vegetation density, and land use.
 *
 * This class reads map data from an XML file, which includes attributes such as the number of cells (`nx` and `ny`),
 * cell sizes (`deltax` and `deltay`), and different maps such as elevation, wind (horizontal and vertical components),
 * vegetation density, and land use.
 */
class RectangularMap2D
{
    public:

        /**
         * @brief Constructor that initializes the RectangularMap2D from an XML file.
         * This constructor loads the map data, including elevation, wind components, vegetation density,
         * and land use, from an XML file using Boost's property tree parser.
         * @param xmlfile The path to the XML file containing the map data.
         */
        RectangularMap2D(const boost::filesystem::path& xmlfile);

        /**
         * @brief Gets the number of cells in the x-direction.
         * @return The number of cells in the x-direction (including ghost layers).
         */
        int nx() const { return nx_; }

        /**
         * @brief Gets the number of cells in the y-direction.
         * @return The number of cells in the y-direction (including ghost layers).
         */
        int ny() const { return ny_; }

        /**
         * @brief Gets the cell size in the x-direction.
         * @return The size of a cell in the x-direction (in meters).
         */
        double deltax() const { return deltax_; }

        /**
         * @brief Gets the cell size in the y-direction.
         * @return The size of a cell in the y-direction (in meters).
         */
        double deltay() const { return deltay_; }

        /**
         * @brief Gets the elevation map.
         * @return A reference to the elevation matrix (in meters).
         */
        const Eigen::MatrixXd& h() const { return h_; }

        /**
         * @brief Gets the wind map (horizontal component).
         * @return A reference to the matrix containing the wind's horizontal component (in meters per second).
         */
        const Eigen::MatrixXd& wx() const { return wx_; }

        /**
         * @brief Gets the wind map (vertical component).
         * @return A reference to the matrix containing the wind's vertical component (in meters per second).
         */
        const Eigen::MatrixXd& wy() const { return wy_; }

        /**
         * @brief Gets the vegetation density map.
         * @return A reference to the matrix containing the vegetation density.
         *         Values range from 0 to +Inf (1 means regular density).
         */
        const Eigen::MatrixXd& density() const { return density_; }

        /**
         * @brief Gets the land use map.
         * @return A reference to the matrix containing the land use type (as integers).
         */
        const Eigen::MatrixXi& landUse() const { return landUse_; }

    private:

        int nx_;                //!< Number of cells (including ghost layers) along the horizontal direction (x).
        int ny_;                //!< Number of cells (including ghost layers) along the vertical direction (y).
        double deltax_;         //!< Size of the cells in the horizontal direction (in meters).
        double deltay_;         //!< Size of the cells in the vertical direction (in meters).

        double lengthx_;        //!< Total length of the physical domain in the horizontal direction (without ghost layers, in meters).
        double lengthy_;        //!< Total length of the physical domain in the vertical direction (without ghost layers, in meters).
        double areatot_;        //!< Total area of the physical domain (without ghost layers, in square meters).
        double areacell_;       //!< Area of a single cell (in square meters).

        int ncells_;            //!< Total number of physical cells.

        Eigen::MatrixXd h_;             //!< Elevation map (in meters).
        Eigen::MatrixXd wx_;            //!< Wind map (horizontal component, in meters per second).
        Eigen::MatrixXd wy_;            //!< Wind map (vertical component, in meters per second).
        Eigen::MatrixXd density_;       //!< Vegetation density map (dimensionless, with 1 meaning regular density).
        Eigen::MatrixXi landUse_;       //!< Land use map (integer types representing land use categories, see `LandUseType`).
};

#endif
