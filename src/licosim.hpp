#pragma once

#include "Eigen/Dense"
#include "eigen/src/Eigenvalues/EigenSolver.h"
#include "licosim/utilities.hpp"
#include "licosim/projectarea.hpp"
#include "licosim/allometry.hpp"
#include "licosim/projectsettings.hpp"
#include "licosim/output.h"

namespace licosim {
    class Licosim {
    public:
        ProjectArea projectArea;
        Treatment treater;
        Allometry allometry;
        fastFuelsAllometry ffa;
        ProjectSettings projectSettings;
        Output output;
        std::default_random_engine dre;

        bool shp = true;
        spatial::SpVectorDataset<spatial::SpMultiPolygon> unitPoly;
        spatial::Raster<int> unitRaster;
        dbhFunction dbhFunc;

        std::vector<std::string> names;
        std::vector<int> ids;
        std::vector<LmuType> types;
        std::vector<int> climateClass;
        std::vector<double> aet;
        std::vector<double> cwd;
        std::vector<double> tmn;
        std::vector<StructureSummary> refStructures;
        double maxDist = 0;
        Pca pca;
        PCA_KDTree<3> allLmu;
        PCA_KDTree<3> ridgeLmu;
        PCA_KDTree<3> valleyLmu;
        PCA_KDTree<3> swFacingLmu;
        PCA_KDTree<3> neFacingLmu;

        spatial::Raster<int> paired;

        Licosim(ProjectSettings& ps);
        std::vector<int> getKnn(Lmu& lmu, const LmuType& type, const double& maxDist, const int k);
        //void assignTargets(int nThread);
        void doTreatmentThreaded(int nThread, double dbhMin, double dbhMax);
        //void writeOutputs(std::string path);

    private:
        Licosim();
        void treatmentThread(int& sofar, std::mutex& mut, double dbhMin, double dbhMax, spatial::Raster<int>& unitZonal, lico::TaoList& treatedTaos, const int thisThread);
        void assignTargetThread(int& sofar, int& nLmu, Lmu& lmu, const int thisThread);
        //void createRxUnitsPolyThread(int& sofar, std::mutex& mut, const int thisThread, spatial::SpVectorDataset<spatial::SpMultiPolygon>& unitPoly, dbhFunction& dbhFunc, std::string& column);
        //void createRxUnitsRasterThread(int& sofar, std::mutex& mut, const int thisThread, spatial::Raster<int>& unitR, dbhFunction& dbhFunc);

    };
}