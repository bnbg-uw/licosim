//Singleton struct for holding options parsed from cmd line. 

//TODO: handle landscape indices.

#pragma once

#ifndef licosim_projectsettings_h
#define licosim_projectsettings_h

#include "boost/program_options.hpp"

namespace licosim {
    struct ProjectSettings {
    public:
        std::string exeParentPath;
        std::string commandLine;
        std::string outputPath;
        bool writeUnits = false;
        bool fastFuels = false;
        std::string lidarDatasetPath;
        std::string projectPolygonPath;
        std::string unitPolygonPath;
        std::string lmuRasterPath;
        std::string priorityMethod;
        std::string priorityColumn;
        std::string terrain = "moderate";
        std::string referenceDatasetPath;
        std::string fireRefPath = "/resources/reference/fire.csv";
        std::string hydroRefPath = "/resources/reference/hydro.csv";
        std::string habitatRefPath = "/resources/reference/habitat.csv";
        std::string aetPath = "/resources/biophysicalanalogs/aet.img";
        std::string cwdPath = "/resources/biophysicalanalogs/cwd.img";
        std::string tmnPath = "/resources/biophysicalanalogs/tmn.img";
        std::string climateClassPath = "/resources/biophysicalanalogs/ClimateClasses.img";
        std::string fiaPath = "/resources/fia/";
        //std::string mcsPath = "/resources/treatment/mcs_prop.csv";
        std::string defaultRefPath;
        int nThread = 1;
        int seed = -1;
        double dbhMin = 15.24;
        double dbhMax = 53.34;
        double slope = 0;
        double intercept = 0;
        double rsq = 0;
        rxtools::allometry::UnivariateLinearModel::Transform transform;
        lapis::LinearUnit inUnit;
        lapis::LinearUnit outUnit;

        bool useFireRef = true;
        bool useHydroRef = false;
        bool useHabitatRef = false;
        bool doFireModel = false;

        bool subdivideLmus = false;
        bool overrideTargets = false;

        //Treatment behavior indices:
        //double percentLandscape = std::nan("");
        //double totalLandscape = std::nan("");
        //double clump1 = 0.5;
        //double clump2 = 0.5;  
        std::vector<double> forestType = { 1, 1, 1 };

        static ProjectSettings& get() {
            static ProjectSettings instance;
            return instance;
        }

        void loadFromParentPath(const std::string exeParentPath) {
            this->exeParentPath = exeParentPath;

            fiaPath = exeParentPath + fiaPath;
            fireRefPath = exeParentPath + fireRefPath;
            hydroRefPath = exeParentPath + hydroRefPath;
            habitatRefPath = exeParentPath + habitatRefPath;
            defaultRefPath = fireRefPath;
            aetPath = exeParentPath + aetPath;
            cwdPath = exeParentPath + cwdPath;
            tmnPath = exeParentPath + tmnPath;
            climateClassPath = exeParentPath + climateClassPath;
            //mcsPath = exeParentPath + mcsPath;
        }

        //TODO: Implement landscape indices.
        void loadFromOptionsAndParentPath(const std::string exeParentPath, const boost::program_options::variables_map& vm) {
            this->exeParentPath = exeParentPath;

            fiaPath = exeParentPath + fiaPath;
            fireRefPath = exeParentPath + fireRefPath;
            hydroRefPath = exeParentPath + hydroRefPath;
            habitatRefPath = exeParentPath + habitatRefPath;
            //mcsPath = exeParentPath + mcsPath;
            defaultRefPath = fireRefPath;

            lidarDatasetPath = vm["lidar"].as<std::string>();
            if (vm.count("output"))
                outputPath = vm["output"].as<std::string>();
            else
                outputPath = lidarDatasetPath;
            if (std::filesystem::exists(outputPath)) {
                auto fs = std::filesystem::path(outputPath);
                fs = fs / "licosim";
                std::filesystem::create_directory(fs);
                outputPath = fs.string();
            }
            else {
                throw std::runtime_error("Output path does not exist.");
            }

            if (vm.count("writeunits")) writeUnits = true;
            if (vm.count("fastfuels")) fastFuels = true;

            projectPolygonPath = vm["projectpoly"].as<std::string>();
            unitPolygonPath = vm["unitpoly"].as<std::string>();

            lmuRasterPath = vm.count("lmu") ? vm["lmu"].as<std::string>() : "";

            priorityMethod = vm["priority"].as<std::string>();
            priorityColumn = vm.count("column") ? vm["column"].as<std::string>() : "column";
            try { if (std::stoi(priorityColumn) < 0) throw std::invalid_argument("Invalid column index"); } catch (...) {}
            if (vm.count("terrain")) terrain = vm["terrain"].as<std::string>();
            referenceDatasetPath = vm.count("reference") ? vm["reference"].as<std::string>() : "";
            if (vm.count("aet"))
                aetPath = vm["aet"].as<std::string>();
            else
                aetPath = exeParentPath + aetPath;

            if (vm.count("cwd"))
                cwdPath = vm["cwd"].as<std::string>();
            else
                cwdPath = exeParentPath + cwdPath;
            if (vm.count("janmin"))
                tmnPath = vm["janmin"].as<std::string>();
            else
                tmnPath = exeParentPath + tmnPath;
            climateClassPath = exeParentPath + climateClassPath;


            if (vm.count("allom_slope") &&
                vm.count("allom_intercept") &&
                vm.count("allom_rsq") &&
                vm.count("allom_transform") &&
                vm.count("allom_inunits") &&
                vm.count("allom_outunits")) {
                slope = vm["allom_slope"].as<double>();
                intercept = vm["allom_intercept"].as<double>();
                rsq = vm["allom_rsq"].as<double>();

                auto tStr = vm["allom_transform"].as<std::string>();
                for (char& c : tStr) {
                    c = std::tolower(static_cast<unsigned char>(c));
                }
                if(tStr == "none") {
                    transform = rxtools::allometry::UnivariateLinearModel::Transform::None;
                }
                else if (tStr == "square") {
                    transform = rxtools::allometry::UnivariateLinearModel::Transform::Square;
                }
                else if (tStr == "cube") {
                    transform = rxtools::allometry::UnivariateLinearModel::Transform::Cube;
                }
                else if (tStr == "power") {
                    transform = rxtools::allometry::UnivariateLinearModel::Transform::Power;
                }
                else {
                    std::cout << "\"" + tStr + "\" is not a valid tranform. Choose one of (none, square, cube, power)";
                    throw std::invalid_argument("invalid allom transform");
                }

                tStr = vm["allom_inunits"].as<std::string>();
                for (char& c : tStr) {
                    c = std::tolower(static_cast<unsigned char>(c));
                }
                if (tStr == "feet") {
                    inUnit = lapis::linearUnitPresets::internationalFoot;
                }
                else if (tStr == "meters") {
                    inUnit = lapis::linearUnitPresets::meter;
                }
                else {
                    std::cout << "\"" + tStr + "\" is not a valid input unit. Choose one of (feet, meters)";
                    throw std::invalid_argument("invalid allom inunit");
                }


                tStr = vm["allom_outunits"].as<std::string>();
                for (char& c : tStr) {
                    c = std::tolower(static_cast<unsigned char>(c));
                }
                if (tStr == "inches") {
                    inUnit = lapis::linearUnitPresets::internationalFoot;
                }
                else if (tStr == "centimeters") {
                    inUnit = lapis::linearUnitPresets::meter;
                }
                else if (tStr == "meters") {
                    inUnit = lapis::linearUnitPresets::meter;
                }
                else {
                    std::cout << "\"" + tStr + "\" is not a valid input unit. Choose one of (inches, centimeters, meters)";
                    throw std::invalid_argument("invalid allom inunit");
                }
            }
            nThread = vm["thread"].as<int>();

            if (vm.count("seed")) seed = vm["seed"].as<int>();

            useFireRef = vm.count("usefire") ? vm["usefire"].as<bool>() : false;
            useHydroRef = vm.count("usehydro") ? vm["usehydro"].as<bool>() : false;
            useHabitatRef = vm.count("usehabitat") ? vm["usehabitat"].as<bool>() : false;
            doFireModel = vm.count("usefire") ? vm["units"].as<bool>() : false;

            //if(options.count("pland"))
            //    percentLandscape = options["pland"].as<double>();
            //if (options.count("tland") && !options.count("pland"))
            //    totalLandscape = options["tland"].as<double>();
            
            dbhMin = vm.count("dbhmin") ? vm["dbhmin"].as<double>() : dbhMin;
            dbhMax = vm.count("dbhmax") ? vm["dbhmax"].as<double>() : dbhMax;

            if (vm.count("cover"))
                forestType = vm["cover"].as<std::vector<double>>();
            if (forestType.size() != 3)
                throw std::invalid_argument("forest type vector should be of length 3.");
            
            double forestSum = 0;
            for (auto v : forestType)
                forestSum += v;
            for (int i = 0; i < forestType.size(); i++)
                forestType[i] /= forestSum;

            overrideTargets = vm.count("overridetargets") ? true : false;
            subdivideLmus = vm.count("subdivideLmus") ? true : false;

        }

    private:
        ProjectSettings() {}
        ProjectSettings(const ProjectSettings&) = delete;
        ProjectSettings& operator=(const ProjectSettings&) = delete;
    };
}  // namespace licosim

#endif  // !licosim_projectsettings_h