/************************************************************************************
 * @file        KarafyllidisModel.cpp
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

#include "KarafyllidisModel.h"
#include <iostream>
#include <boost/math/constants/constants.hpp>

/**
 * @brief Constructor: initialize the KarafyllidisModel instance and allocate memory for required matrices.
 */
KarafyllidisModel::KarafyllidisModel(TimeManager& runTime, const int nx, const int ny, const double deltax, const double deltay) :
    runTime_(runTime),
    nx_(nx),
    ny_(ny),
    deltax_(deltax),
    deltay_(deltay)
{
    lengthx_ = (nx_-2)*deltax_;
    lengthy_ = (ny_-2)*deltay_;
    
    areatot_ = lengthx_ * lengthy_;
    areacell_ = deltax_ * deltay_;

    ncells_ = (nx_-2) * (ny_-2);
    ncellstot_ = nx_ * ny_;

    a_.setZero(nx_, ny_);
    h_.setZero(nx_, ny_);
    wx_.setZero(nx_, ny_);
    wy_.setZero(nx_, ny_);
    density_.setConstant(nx_,ny_,1.);
    landUse_.setConstant(nx_, ny_, KarafyllidisModel::GRASSLAND);
    R_.setZero(nx_, ny_);

    tecplotFolder_ = "Tecplot";
    vtkFolder_ = "VTK";
    tecplotOutput_ = false;
    vtkOutput_ = false;
}

/**
 * @brief Set the log file name and open the file for output.
 */
void KarafyllidisModel::SetLogFile(const std::string name)
{
    fLog_.open(name, std::ios::out);
    if (!fLog_.is_open())
        throw std::runtime_error("Failed to open log file: " + name);
    fLog_.setf(std::ios::scientific);

    // Write header line
    fLog_ << std::setw(12) << std::left << "Step";
    fLog_ << std::setw(16) << std::left << "Time(s)";
    fLog_ << std::setw(16) << std::left << "Interval(s)";
    fLog_ << std::setw(16) << std::left << "BurntArea(m2)";
    fLog_ << std::setw(16) << std::left << "BurntArea(%)";
    fLog_ << std::setw(16) << std::left << "Diameter(m)";
    fLog_ << std::setw(16) << std::left << "BurntCells";
    fLog_ << std::setw(16) << std::left << "BurntCells(%)";
    fLog_ << std::endl;
}

/**
 * @brief Set the folder for Tecplot output and ensure it exists.
 */
void KarafyllidisModel::SetTecplotFolder(const boost::filesystem::path folder)
{
    tecplotFolder_ = folder;
    tecplotOutput_ = true;

    // Create directory and handle errors
    try 
    {
        if (!boost::filesystem::create_directory(tecplotFolder_))
            throw std::runtime_error("Failed to create Tecplot folder: " + tecplotFolder_.string());
    } 
    catch (const boost::filesystem::filesystem_error& ex) 
    {
        throw std::runtime_error("Filesystem error creating Tecplot folder: " + tecplotFolder_.string() + " - " + ex.what());
    }
}

/**
 * @brief Set the folder for VTK output and ensure it exists.
 */
void KarafyllidisModel::SetVTKFolder(const boost::filesystem::path folder)
{
    vtkFolder_ = folder;
    vtkOutput_ = true;

    // Create directory and handle errors
    try 
    {
        if (!boost::filesystem::create_directory(vtkFolder_)) 
        {
            throw std::runtime_error("Failed to create VTK folder: " + vtkFolder_.string());
        }
    } 
    catch (const boost::filesystem::filesystem_error& ex) 
    {
        throw std::runtime_error("Filesystem error creating VTK folder: " + vtkFolder_.string() + " - " + ex.what());
    }
}

/**
 * @brief Set a single ignition cell by its matrix indices.
 */
void KarafyllidisModel::SetIgnitionCell(const int ix, const int iy)
{
    a_(ix, iy) = 1.;
}

/**
 * @brief Set a single ignition cell using physical coordinates.
 */
void KarafyllidisModel::SetIgnitionCell(const double x, const double y)
{
    const int ix = std::ceil(x / deltax_);
    const int iy = std::ceil((lengthy_ - y) / deltay_);
    a_(ix, iy) = 1.;
}

/**
 * @brief Set multiple ignition cells using a matrix of indices.
 */
void KarafyllidisModel::SetIgnitionCells(const Eigen::MatrixXi& indices)
{
    for (int i = 0; i < indices.rows(); i++)
        a_(indices(i, 0), indices(i, 1)) = 1.;
}

/**
 * @brief Set multiple ignition cells within a rectangular region defined by a center and dimensions.
 */
void KarafyllidisModel::SetIgnitionCells(const double xc, const double yc, const double dx, const double dy)
{
    const int xm = std::ceil((xc - dx / 2.) / deltax_);
    const int ym = std::ceil((yc - dy / 2.) / deltay_);
    const int xp = std::ceil((xc + dx / 2.) / deltax_);
    const int yp = std::ceil((yc + dy / 2.) / deltay_);

    for (int i = xm; i < xp; i++)
        for (int j = ym; j < yp; j++)
            a_(i, j) = 1.;
}

/**
 * @brief Set terrain features, including land use and vegetation density.
 */
void KarafyllidisModel::SetTerrainFeatures(const Eigen::MatrixXi& landUse, const Eigen::MatrixXd& density)
{
    // Assign maps
    landUse_ = landUse;
    density_ = density;

    // Write maps
    if (tecplotOutput_ == true)
    {
        const boost::filesystem::path fileNameLandUse = tecplotFolder_ / ("LandUse.tec");
        WriteTecplotFileASCII(fileNameLandUse, "landUse", landUse_, deltax_, deltay_);
        const boost::filesystem::path fileNameDensity = tecplotFolder_ / ("Density.tec");
        WriteTecplotFileASCII(fileNameDensity, "density", density_, deltax_, deltay_);
    }
    if (vtkOutput_ == true)
    {
        const boost::filesystem::path fileName = vtkFolder_ / ("LandUse.vtk");
        WriteVTKFileASCII(fileName, "landUse", landUse_, deltax_, deltay_);
        const boost::filesystem::path fileNameDensity = vtkFolder_ / ("Density.vtk");
        WriteVTKFileASCII(fileNameDensity, "density", density_, deltax_, deltay_);
    }
}

/**
 * @brief Set the elevation map and compute the necessary correction coefficients for adjacent and diagonal cells.
 */
void KarafyllidisModel::SetElevation(const Eigen::MatrixXd& h)
{
    // Assign map
    h_ = h;

    // Calculate submatrices
    hadj_.resize(nx_ * ny_);
    hdiag_.resize(nx_ * ny_);
    for (int i = 1; i < nx_ - 1; i++)
        for (int j = 1; j < ny_ - 1; j++)
        {
            const Eigen::MatrixXd htot = ElevationSubMatrix(i, j);
            const int index = j + i * nx_;

            hadj_[index]  << htot(1, 0), htot(0, 1), htot(1, 2), htot(2, 1);
            hdiag_[index] << htot(2, 0), htot(0, 0), htot(0, 2), htot(2, 2);
        }

    // Write maps
    if (tecplotOutput_ == true)
    {
        const boost::filesystem::path fileName = tecplotFolder_ / ("Elevation.tec");
        WriteTecplotFileASCII(fileName, "h", h_, deltax_, deltay_);
    }
    if (vtkOutput_ == true)
    {
        const boost::filesystem::path fileName = vtkFolder_ / ("Elevation.vtk");
        WriteVTKFileASCII(fileName, "h", h_, deltax_, deltay_);
    }
}

/**
 * @brief Extract a 3x3 submatrix centered at cell (i, j) from the elevation matrix.
 * @param i Row index of the center cell.
 * @param j Column index of the center cell.
 * @return A 3x3 matrix representing the slopes between adjacent and diagonal cells.
 */
Eigen::Matrix<double,3,3> KarafyllidisModel::ElevationSubMatrix(const int i, const int j)
{
    const double L = std::min(deltax_, deltay_);
    const double Ld = L * std::sqrt(2.);
    
    Eigen::Matrix<double,3,3> slope;
    slope.setZero();

    // Adjacent cells
    slope(1, 0) = std::atan((h_(i, j) - h_(i, j - 1)) / L);
    slope(0, 1) = std::atan((h_(i, j) - h_(i - 1, j)) / L);
    slope(1, 2) = std::atan((h_(i, j) - h_(i, j + 1)) / L);
    slope(2, 1) = std::atan((h_(i, j) - h_(i + 1, j)) / L);

    // Diagonal cells
    slope(2, 0) = std::atan((h_(i, j) - h_(i + 1, j - 1)) / Ld);
    slope(0, 0) = std::atan((h_(i, j) - h_(i - 1, j - 1)) / Ld);
    slope(0, 2) = std::atan((h_(i, j) - h_(i - 1, j + 1)) / Ld);
    slope(2, 2) = std::atan((h_(i, j) - h_(i + 1, j + 1)) / Ld);

    // Apply exponential transformation
    const double pi = boost::math::constants::pi<double>();
    const double a = 0.0693 * (360. / (2. * pi));

    Eigen::Matrix<double,3,3> h;
    for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
            h(ii, jj) = std::exp(a * slope(ii, jj));

    return h;
}

/**
 * @brief Set the wind velocity components and compute the wind correction coefficients for adjacent and diagonal cells.
 * @param wx Matrix representing the wind speed in the horizontal direction (x-axis).
 * @param wy Matrix representing the wind speed in the vertical direction (y-axis).
 */
void KarafyllidisModel::SetWind(const Eigen::MatrixXd& wx, const Eigen::MatrixXd& wy)
{
    // Assign maps
    wx_ = wx;
    wy_ = wy;

    // Calculate submatrices
    wadj_.resize(nx_*ny_);
    wdiag_.resize(nx_*ny_);
    for (int i = 1; i < nx_ - 1; i++)
        for (int j = 1; j < ny_ - 1; j++)
        {
            const Eigen::MatrixXd wtot = WindSubMatrix(i, j);
            const int index = j + i * nx_;

            wadj_[index]  << wtot(1, 0), wtot(0, 1), wtot(1, 2), wtot(2, 1);
            wdiag_[index] << wtot(2, 0), wtot(0, 0), wtot(0, 2), wtot(2, 2);
        }

    // Write maps
    if (tecplotOutput_ == true)
    {
        const boost::filesystem::path fileName = tecplotFolder_ / ("Wind.tec");
        WriteTecplotFileASCII(fileName, "wind", wx_, wy_, deltax_, deltay_);
    }
    if (vtkOutput_ == true)
    {
        const boost::filesystem::path fileName = vtkFolder_ / ("Wind.vtk");
        WriteVTKFileASCII(fileName, "wind", wx_, wy_, deltax_, deltay_);
    }
}

/**
 * @brief Compute the wind correction coefficient based on wind magnitude and angle.
 * @param wmag Wind magnitude (m/s).
 * @param cosTeta Cosine of the angle between the wind vector and the desired direction.
 * @return The wind correction coefficient.
 */
double WindCorrectionCoefficient(const double wmag, const double cosTeta)
{   
    const double c1 = 0.045;
    const double c2 = 0.131;
    return std::exp(wmag * (c1 + c2 * (cosTeta - 1.)));
}

/**
 * @brief Extract a 3x3 submatrix for the wind correction coefficients around cell (i, j).
 * @param i Row index of the center cell.
 * @param j Column index of the center cell.
 * @return A 3x3 matrix containing wind correction factors for adjacent and diagonal cells.
 */
Eigen::Matrix<double,3,3> KarafyllidisModel::WindSubMatrix(const int i, const int j)
{
    // Wind velocity components and magnitude (m/s)
    const double wx = wx_(i, j);
    const double wy = wy_(i, j);
    const double wmag = std::sqrt(wx * wx + wy * wy);

    Eigen::Matrix<double, 3, 3> w;
    w.setConstant(1.);

    // Skip if there is no wind
    if (wmag > 0.)
    {
        // Adjacent cells
        w(2, 1) = WindCorrectionCoefficient(wmag, -wx / wmag);
        w(1, 2) = WindCorrectionCoefficient(wmag,  wy / wmag);
        w(0, 1) = WindCorrectionCoefficient(wmag,  wx / wmag);
        w(1, 0) = WindCorrectionCoefficient(wmag, -wy / wmag);

        // Diagonal cells
        w(2, 0) = WindCorrectionCoefficient(wmag, (-wx - wy) / (wmag * std::sqrt(2.)));
        w(2, 2) = WindCorrectionCoefficient(wmag, (-wx + wy) / (wmag * std::sqrt(2.)));
        w(0, 2) = WindCorrectionCoefficient(wmag, ( wx + wy) / (wmag * std::sqrt(2.)));
        w(0, 0) = WindCorrectionCoefficient(wmag, ( wx - wy) / (wmag * std::sqrt(2.)));
    }

    return w;
}

/**
 * @brief Calculate the spread rate of fire across the terrain for each cell based on land use type and vegetation density.
 */
void KarafyllidisModel::CalculateSpreadRate()
{
    const std::unordered_map<LandUseType, double> landUseSpreadRates = 
    {
        {WATER, 0.0}, {BUILDINGS, 0.0}, {ROCKS, 0.0}, {BROADLEAVES, 0.05},
        {SHRUBS, 0.15}, {GRASSLAND, 1.0}, {FIREPRONECONIFERS, 0.25}, 
        {AGROFORESTRY, 0.15}, {NOTFIREPRONEFORESTRY, 0.01}, {SAVANNA, 0.50},
        {TAIGA, 0.05}, {TROPICALGRASSLAND, 0.75}, {MEDITERRANEANWOODLAND, 0.20},
        {ALPINEFOREST, 0.05}, {TUNDRA, 0.01}, {PEATLANDS, 0.005}
    };

    #pragma omp parallel for collapse(2)  // Parallelize over both i and j
    for (int i = 0; i < nx_; i++) 
        for (int j = 0; j < ny_; j++) 
        {
            LandUseType landUseType = static_cast<LandUseType>(landUse_(i, j));
            R_(i, j) = landUseSpreadRates.at(landUseType) * density_(i, j);
        }
}

/**
 * @brief Advance the simulation by updating the fire spread over time and writing logs and output files if required.
 */
void KarafyllidisModel::Advance()
{
    // Determine maximum time step based on the physics
    const double RmaxAbs = R_.maxCoeff();
    const double L = std::min(deltax_,deltay_);
    const double deltatMax = L/RmaxAbs * runTime_.Co();

    // Update run time
    runTime_.Advance(deltatMax);

    // Correct the time step accounting for the output requests
    const double deltat = runTime_.dt();
    const double Rmax = L/deltat;

    // Print on the screen
    if (runTime_.printScreen())
        std::cout << std::scientific << runTime_.step() << " " << runTime_.t() << " " << runTime_.dt() << " " << RmaxAbs << " " << RmaxAbs/Rmax << std::endl;

    // Update cells
    Eigen::MatrixXd aold = a_;

    #pragma omp parallel for collapse(2)  // Parallelize over both i and j
    for (int i=1;i<nx_-1;i++)
        for (int j=1;j<ny_-1;j++)
        {
            // Global index
            const int index = j + i*nx_;

            // Look at cell partially burned only
            if (aold(i,j)<1.)
            {
                // Adjacent cells
                double adj = 0.;
                {
                    const std::vector<int> ii { i, i-1, i, i+1 };
                    const std::vector<int> jj { j-1, j, j+1, j };

                    for (int kk=0;kk<4;kk++)
                    {
                        const int ik = ii[kk];
                        const int jk = jj[kk];

                        const double l = wadj_[index](kk)*hadj_[index](kk) * R_(ik,jk) * deltat;
                        const double coefficient = l/L;
                        adj += coefficient * aold(ik,jk);    

                    }
                }

                // Diagonal cells
                double diag = 0.;
                {
                    const std::vector<int> ii { i+1, i-1, i-1, i+1 };
                    const std::vector<int> jj { j-1, j-1, j+1, j+1 };

                    for (int kk=0;kk<4;kk++)
                    {
                        const int ik = ii[kk];
                        const int jk = jj[kk];

                        const double l = wdiag_[index](kk)*hdiag_[index](kk) * R_(ik,jk) * deltat;
                        double coefficient = l*l/(L*L);
                        if (l>=std::sqrt(2.)*L/2.)
                            coefficient = (2.*std::sqrt(2.)*L*l-l*l-L*L)/(L*L);

                        diag += coefficient * aold(ik,jk);    
                    }
                }

                // Overall combination
                a_(i,j) = aold(i,j) + adj + diag;
                a_(i,j) = std::min(1.,a_(i,j));
            }
        }

    // Print on log file
    if (runTime_.printLog())
        PrintLogFile();   

    // Print maps
    if (runTime_.printMap() == true)
    {
        if (tecplotOutput_ == true)
        {
            const boost::filesystem::path fileName = tecplotFolder_ / ("Solution.tec." + std::to_string(runTime_.step()));
            WriteTecplotFileASCII(fileName, "a", a_, deltax_, deltay_);
        }

        if (vtkOutput_ == true)
        {
            const boost::filesystem::path fileName = vtkFolder_ / ("Solution.vtk." + std::to_string(runTime_.step()));
            WriteVTKFileASCII(fileName, "a", a_, deltax_, deltay_);
        }
    }
}

void KarafyllidisModel::PrintLogFile()
{
    // Check if the log file is open
    if (!fLog_.is_open())
        throw std::runtime_error("Log file is not open for writing.");

    // Total burned area (m2)
    const double burnedArea = a_.sum() * areacell_;
    
    // Total number of fully burned cells
    int fullBurnedCells = 0;
    for (int i=1;i<nx_-1;i++)
        for (int j=1;j<ny_-1;j++)
            if (a_(i,j) == 1.)
                fullBurnedCells++;

    // Equivalent diameter
    auto pi = boost::math::constants::pi<double>();
    const double De = std::sqrt(4.*burnedArea/pi);

    // Create an ostringstream to buffer the log output
    std::ostringstream logBuffer;
    logBuffer << std::setw(12) << std::left << runTime_.step();
    logBuffer << std::setw(16) << std::left << std::scientific << runTime_.t();
    logBuffer << std::setw(16) << std::left << std::scientific << runTime_.dt();
    logBuffer << std::setw(16) << std::left << std::scientific << burnedArea;
    logBuffer << std::setw(16) << std::left << std::scientific << 100.*burnedArea / areatot_;
    logBuffer << std::setw(16) << std::left << std::scientific << De;
    logBuffer << std::setw(16) << std::left << std::scientific << fullBurnedCells;
    logBuffer << std::setw(16) << std::left << std::scientific << 100.*static_cast<double>(fullBurnedCells) / static_cast<double>(ncells_);
    logBuffer << std::endl;

    // Now flush the buffered content to the log file
    fLog_ << logBuffer.str();

    // Check for I/O errors after flushing
    if (fLog_.fail())
        throw std::runtime_error("Failed to write to log file.");
}
