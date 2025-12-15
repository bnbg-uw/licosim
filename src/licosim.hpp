#pragma once

#include "Eigen/Dense"
#include "eigen/src/Eigenvalues/EigenSolver.h"
#include "utilities.hpp"
#include "projectarea.hpp"
#include "allometry.hpp"
#include "projectsettings.hpp"
#include "output.hpp"
#include "pca.hpp"

namespace licosim {
    class Licosim {
    public:
        rxtools::ProjectArea projectArea;
        rxtools::Treatment treater;
        rxtools::allometry::UnivariateLinearModel dbhModel;
        rxtools::allometry::FastFuels ffa;
        rxtools::Output output;
        std::default_random_engine dre;

        bool shp = true;
        lapis::VectorDataset<lapis::MultiPolygon> unitPoly;
        lapis::Raster<int> unitRaster;
        dbhFunction dbhFunc;

        std::vector<std::string> names;
        std::vector<int> ids;
        std::vector<rxtools::LmuType> types;
        std::vector<int> climateClass;
        std::vector<double> aet;
        std::vector<double> cwd;
        std::vector<double> tmn;
        std::vector<rxtools::StructureSummary> refStructures;
        double maxDist = 0;
        rxtools::utilities::Pca pca;
        rxtools::utilities::PCA_KDTree<3> allLmu;
        rxtools::utilities::PCA_KDTree<3> ridgeLmu;
        rxtools::utilities::PCA_KDTree<3> valleyLmu;
        rxtools::utilities::PCA_KDTree<3> swFacingLmu;
        rxtools::utilities::PCA_KDTree<3> neFacingLmu;

        lapis::Raster<int> paired;

        Licosim();
        std::vector<int> getKnn(rxtools::Lmu& lmu, const rxtools::LmuType& type, const double& maxDist, const int k);
        //void assignTargets(int nThread);
        void doTreatmentThreaded(int nThread, double dbhMin, double dbhMax);
        //void writeOutputs(std::string path);

    private:
        Licosim();
        void treatmentThread(int& sofar, std::mutex& mut, double dbhMin, double dbhMax, lapis::Raster<int>& unitZonal, lapis::TaoListMP& treatedTaos, const int thisThread);
        void assignTargetThread(int& sofar, int& nLmu, rxtools::Lmu& lmu, const int thisThread);
        //void createRxUnitsPolyThread(int& sofar, std::mutex& mut, const int thisThread, spatial::SpVectorDataset<spatial::SpMultiPolygon>& unitPoly, dbhFunction& dbhFunc, std::string& column);
        //void createRxUnitsRasterThread(int& sofar, std::mutex& mut, const int thisThread, spatial::Raster<int>& unitR, dbhFunction& dbhFunc);

    };
}