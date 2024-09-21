/************************************************************************************
 * @file        Utilities.cpp
 * @brief       Utilities
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

#include "Utilities.h"
#include <stdexcept>
#include <map>

template<typename Matrix>
void WriteVTKFileASCII(const boost::filesystem::path fileName, const std::string fieldName, const Matrix& v, const double deltax, const double deltay)
{
    // Number of points
    const int nx = v.rows();
    const int ny = v.cols();

    // Open file
    boost::filesystem::ofstream fVTK(fileName, std::ios::out);
    if (!fVTK.is_open())
        throw std::runtime_error("Failed to open VTK file: " + fileName.string());
    fVTK.setf(std::ios::scientific);

    // Create an ostringstream buffer to accumulate output
    std::ostringstream buffer;

    // Accumulate header and metadata
    buffer << "# vtk DataFile Version 2.0" << std::endl;
    buffer << "vtk output" << std::endl;
    buffer << "ASCII" << std::endl;

    buffer << "DATASET STRUCTURED_POINTS" << std::endl;
    buffer << "DIMENSIONS " << nx-2 << " " << ny-2 << " 1" << std::endl;
    buffer << "SPACING " << deltax << " " << deltay << " " << 0. << std::endl;
    buffer << "ORIGIN " << 0. << " " << 0. << " " << 0. << std::endl;

    buffer << "POINT_DATA " << (nx - 2) * (ny - 2) << std::endl;
    buffer << "SCALARS " << fieldName << " double 1" << std::endl;
    buffer << "LOOKUP_TABLE default" << std::endl;

    // Accumulate the actual data
    for (int j = ny - 2; j > 0; j--) 
        for (int i = 1; i < nx - 1; i++) 
        {
            const double value = (v(i, j) < 1.e-16) ? 0. : v(i, j);
            buffer << value << std::endl;
        }

    // Write the buffered content to the file in one go
    fVTK << buffer.str();

    // Optionally, check for any write errors
    if (fVTK.fail())
        throw std::runtime_error("Failed to write data to VTK file: " + fileName.string());

    // No need to explicitly close the file, as fVTK will close automatically
    fVTK.close();
}

template<typename Matrix>
void WriteVTKFileASCII(const boost::filesystem::path fileName, const std::string fieldName, 
                       const Matrix& vx, const Matrix& vy, const double deltax, const double deltay)
{
    // Number of points
    const int nx = vx.rows();
    const int ny = vx.cols();

    // Open file
    boost::filesystem::ofstream fVTK(fileName, std::ios::out);
    if (!fVTK.is_open())
        throw std::runtime_error("Failed to open VTK file: " + fileName.string());
    fVTK.setf(std::ios::scientific);

    // Create an ostringstream buffer to accumulate output
    std::ostringstream buffer;

    // Accumulate header and metadata
    buffer << "# vtk DataFile Version 2.0" << std::endl;
    buffer << "vtk output" << std::endl;
    buffer << "ASCII" << std::endl;

    buffer << "DATASET STRUCTURED_POINTS" << std::endl;
    buffer << "DIMENSIONS " << nx - 2 << " " << ny - 2 << " 1" << std::endl;
    buffer << "SPACING " << deltax << " " << deltay << " " << 0. << std::endl;
    buffer << "ORIGIN " << 0. << " " << 0. << " " << 0. << std::endl;

    buffer << "POINT_DATA " << (nx - 2) * (ny - 2) << std::endl;

    // Specify that we are outputting vector data
    buffer << "VECTORS " << fieldName << " double" << std::endl;

    // Accumulate the vector field data (x and y components)
    for (int j = ny - 2; j > 0; j--) 
        for (int i = 1; i < nx - 1; i++) 
        {
            const double value_x = (vx(i, j) < 1.e-16) ? 0. : vx(i, j);
            const double value_y = (vy(i, j) < 1.e-16) ? 0. : vy(i, j);
            const double value_z = 0.0;  // z component is zero for 2D vectors

            // Output vector (value_x, value_y, value_z)
            buffer << value_x << " " << value_y << " " << value_z << std::endl;
        }

    // Write the buffered content to the file in one go
    fVTK << buffer.str();

    // Optionally, check for any write errors
    if (fVTK.fail())
        throw std::runtime_error("Failed to write data to VTK file: " + fileName.string());

    // No need to explicitly close the file, as fVTK will close automatically
    fVTK.close();
}


template<typename Matrix>
void WriteTecplotFileASCII(const boost::filesystem::path fileName, const std::string fieldName, const Matrix& v, const double deltax, const double deltay)
{
    // Number of points
    const int nx = v.rows();
    const int ny = v.cols();
    const double lengthy = (ny-2)*deltay;

    // Open file
    boost::filesystem::ofstream fTecplot(fileName, std::ios::out);
    if (!fTecplot.is_open())
        throw std::runtime_error("Failed to open Tecplot file: " + fileName.string());
    fTecplot.setf(std::ios::scientific);

    // Create an ostringstream buffer to accumulate the output
    std::ostringstream buffer;

    // Write the header and metadata into the buffer
    buffer << "Title=Solution" << std::endl;
    buffer << "Variables=" << "\"x[m]\", \"y[m]\", \""<< fieldName << "\"" << std::endl;
    buffer << "Zone I=" << nx - 2 << ", J=" << ny - 2 << ", F=POINT" << std::endl;

    // Accumulate data into the buffer
    for (int i = 1; i < nx - 1; i++)
    {
        const double x = (i - 0.5) * deltax;
        for (int j = 1; j < ny - 1; j++)
        {
            const double y = lengthy - (j - 0.5) * deltay;
            const double value = (v(i, j) < 1.e-16) ? 0. : v(i, j);
            buffer << x << " " << y << " " << value << std::endl;
        }
    }

    // Write the entire buffer to the file in one go
    fTecplot << buffer.str();

    // Optionally, check for any write errors
    if (fTecplot.fail())
        throw std::runtime_error("Failed to write data to Tecplot file: " + fileName.string());

    // No need to explicitly close the file; it will close automatically when out of scope
    fTecplot.close();
}

template<typename Matrix>
void WriteTecplotFileASCII(const boost::filesystem::path fileName, const std::string fieldName, 
                           const Matrix& vx, const Matrix& vy, const double deltax, const double deltay)
{
    // Number of points
    const int nx = vx.rows();
    const int ny = vx.cols();
    const double lengthy = (ny - 2) * deltay;

    // Open file
    boost::filesystem::ofstream fTecplot(fileName, std::ios::out);
    if (!fTecplot.is_open())
        throw std::runtime_error("Failed to open Tecplot file: " + fileName.string());
    fTecplot.setf(std::ios::scientific);

    // Create an ostringstream buffer to accumulate the output
    std::ostringstream buffer;

    // Write the header and metadata into the buffer
    buffer << "Title=Solution" << std::endl;
    // Modify variables to include vector components vx and vy
    buffer << "Variables=" << "\"x[m]\", \"y[m]\", \"" << fieldName << "x\", \"" << fieldName << "y\"" << std::endl;
    buffer << "Zone I=" << nx - 2 << ", J=" << ny - 2 << ", F=POINT" << std::endl;

    // Accumulate data into the buffer
    for (int i = 1; i < nx - 1; i++)
    {
        const double x = (i - 0.5) * deltax;
        for (int j = 1; j < ny - 1; j++)
        {
            const double y = lengthy - (j - 0.5) * deltay;
            const double value_x = (vx(i, j) < 1.e-16) ? 0. : vx(i, j);
            const double value_y = (vy(i, j) < 1.e-16) ? 0. : vy(i, j);
            buffer << x << " " << y << " " << value_x << " " << value_y << std::endl;
        }
    }

    // Write the entire buffer to the file in one go
    fTecplot << buffer.str();

    // Optionally, check for any write errors
    if (fTecplot.fail())
        throw std::runtime_error("Failed to write data to Tecplot file: " + fileName.string());

    // No need to explicitly close the file; it will close automatically when out of scope
    fTecplot.close();
}

template<typename T>
void ImportXMLMap(boost::property_tree::ptree& ptree, const std::string name, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& M)
{
    // Locate the "elevation" section
        boost::property_tree::ptree elevationNode = ptree.get_child(name);

        // Read the "nx" and "ny" attributes
        const int nx = elevationNode.get<int>("<xmlattr>.nx");
        const int ny = elevationNode.get<int>("<xmlattr>.ny");

        // Read the content (sequence of nx * ny numbers)
        std::string content = elevationNode.get_value<std::string>();
        
        // Create a string stream from the content to extract numbers
        std::istringstream iss(content);

        // In case of uniform values
        if (nx == 1 && ny == 1)
        {
            T value;
            iss >> value; 
            M.setConstant(value);
        }
        // In case of non uniform values
        else
        {
            T value;
            for (int i=0;i<nx;i++)
                for (int j=0;j<ny;j++)
                {
                    iss >> value; 
                    M(i,j) = value;
                }
        }
}

// Bilinear interpolation for double-sized matrix
Eigen::MatrixXd BilinearInterpolate(const Eigen::MatrixXd& matrix)
{
    const int originalRows = matrix.rows();
    const int originalCols = matrix.cols();
    const int newRows = originalRows * 2 - 1;
    const int newCols = originalCols * 2 - 1;

    Eigen::MatrixXd result(newRows, newCols);

    // Copy the original values to the even-indexed positions
    for (int i = 0; i < originalRows; ++i)
        for (int j = 0; j < originalCols; ++j)
            result(2 * i, 2 * j) = matrix(i, j);

    // Interpolate rows
    for (int i = 0; i < originalRows - 1; ++i)
        for (int j = 0; j < originalCols; ++j)
            result(2 * i + 1, 2 * j) = 0.5 * (matrix(i, j) + matrix(i + 1, j));

    // Interpolate columns
    for (int i = 0; i < originalRows; ++i)
        for (int j = 0; j < originalCols - 1; ++j)
            result(2 * i, 2 * j + 1) = 0.5 * (matrix(i, j) + matrix(i, j + 1));

    // Interpolate centers
    for (int i = 0; i < originalRows - 1; ++i)
        for (int j = 0; j < originalCols - 1; ++j)
            result(2 * i + 1, 2 * j + 1) = 0.25 * (matrix(i, j) + matrix(i + 1, j) +
                                                   matrix(i, j + 1) + matrix(i + 1, j + 1));

    return result;
}

// Nearest-neighbor interpolation for integer matrices (e.g., land use map)
Eigen::MatrixXi NearestNeighborInterpolate(const Eigen::MatrixXi& matrix)
{
    const int originalRows = matrix.rows();
    const int originalCols = matrix.cols();
    const int newRows = originalRows * 2 - 1;
    const int newCols = originalCols * 2 - 1;

    Eigen::MatrixXi result(newRows, newCols);

    // Copy the original values to the even-indexed positions
    for (int i = 0; i < originalRows; ++i)
        for (int j = 0; j < originalCols; ++j)
            result(2 * i, 2 * j) = matrix(i, j);

    // Nearest-neighbor interpolation for rows
    for (int i = 0; i < originalRows - 1; ++i)
        for (int j = 0; j < originalCols; ++j)
            result(2 * i + 1, 2 * j) = matrix(i, j); // Nearest-neighbor from above

    // Nearest-neighbor interpolation for columns
    for (int i = 0; i < originalRows; ++i)
        for (int j = 0; j < originalCols - 1; ++j)
            result(2 * i, 2 * j + 1) = matrix(i, j); // Nearest-neighbor from left

    // Nearest-neighbor interpolation for centers
    for (int i = 0; i < originalRows - 1; ++i)
        for (int j = 0; j < originalCols - 1; ++j)
            result(2 * i + 1, 2 * j + 1) = matrix(i, j); // Nearest-neighbor from top-left

    return result;
}

// Coarsen the matrix by averaging 2x2 blocks
Eigen::MatrixXd CoarsenAverage(const Eigen::MatrixXd& matrix)
{
    const int originalRows = matrix.rows();
    const int originalCols = matrix.cols();

    // Check if the number of rows and columns is even
    if (originalRows % 2 != 0 || originalCols % 2 != 0)
        throw std::invalid_argument("The number of rows (nx) and columns (ny) must be even for coarsening.");

    const int newRows = originalRows / 2;
    const int newCols = originalCols / 2;

    Eigen::MatrixXd result(newRows, newCols);

    for (int i = 0; i < newRows; ++i)
        for (int j = 0; j < newCols; ++j)
        {
            // Calculate the average of the 2x2 block
            result(i, j) = 0.25 * (matrix(2 * i, 2 * j) + matrix(2 * i + 1, 2 * j) +
                                   matrix(2 * i, 2 * j + 1) + matrix(2 * i + 1, 2 * j + 1));
        }

    return result;
}

// Coarsen the integer matrix using the most frequent value in each 2x2 block
Eigen::MatrixXi CoarsenMostFrequent(const Eigen::MatrixXi& matrix)
{
    const int originalRows = matrix.rows();
    const int originalCols = matrix.cols();

    // Check if the number of rows and columns is even
    if (originalRows % 2 != 0 || originalCols % 2 != 0)
        throw std::invalid_argument("The number of rows (nx) and columns (ny) must be even for coarsening.");
            
    const int newRows = originalRows / 2;
    const int newCols = originalCols / 2;

    Eigen::MatrixXi result(newRows, newCols);

    for (int i = 0; i < newRows; ++i)
        for (int j = 0; j < newCols; ++j)
        {
            // Get the 2x2 block
            int values[4] = {
                matrix(2 * i, 2 * j),
                matrix(2 * i + 1, 2 * j),
                matrix(2 * i, 2 * j + 1),
                matrix(2 * i + 1, 2 * j + 1)
            };

            // Count the occurrences of each value
            std::map<int, int> freq;
            for (int k = 0; k < 4; ++k) {
                freq[values[k]]++;
            }

            // Find the most frequent value in the 2x2 block
            int mostFrequent = values[0];
            int maxCount = freq[values[0]];
            for (const auto& [value, count] : freq) {
                if (count > maxCount) {
                    mostFrequent = value;
                    maxCount = count;
                }
            }

            result(i, j) = mostFrequent;
        }

    return result;
}


// Explicit instantiations

template void WriteVTKFileASCII<Eigen::MatrixXi>(const boost::filesystem::path fileName, const std::string fieldName, const Eigen::MatrixXi& v, const double deltax, const double deltay);
template void WriteVTKFileASCII<Eigen::MatrixXd>(const boost::filesystem::path fileName, const std::string fieldName, const Eigen::MatrixXd& v, const double deltax, const double deltay);

template void WriteTecplotFileASCII<Eigen::MatrixXi>(const boost::filesystem::path fileName, const std::string fieldName, const Eigen::MatrixXi& v, const double deltax, const double deltay);
template void WriteTecplotFileASCII<Eigen::MatrixXd>(const boost::filesystem::path fileName, const std::string fieldName, const Eigen::MatrixXd& v, const double deltax, const double deltay);

template void WriteVTKFileASCII<Eigen::MatrixXd>(const boost::filesystem::path fileName, const std::string fieldName, const Eigen::MatrixXd& vx, const Eigen::MatrixXd& vy, const double deltax, const double deltay);
template void WriteTecplotFileASCII<Eigen::MatrixXd>(const boost::filesystem::path fileName, const std::string fieldName, const Eigen::MatrixXd& vx, const Eigen::MatrixXd& vy, const double deltax, const double deltay);

template void ImportXMLMap<int>(boost::property_tree::ptree& ptree, const std::string name, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& M);
template void ImportXMLMap<double>(boost::property_tree::ptree& ptree, const std::string name, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& M);