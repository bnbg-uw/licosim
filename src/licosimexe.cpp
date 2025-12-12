// licosim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <boost/program_options.hpp>
#include "licosim.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    cxxopts::Options opt{"LicoSim", "Given lidar data and a lot of area parameters, calculate proposed LMU based treatment alternatives based on given reference conditions."};
    opt.add_options()
        ("p,projectpoly", "REQUIRED: Path to project area polygon (esri shapefile).", cxxopts::value<std::string>())
        ("u,unitpoly", "REQUIRED: Path to unit polygons (esri shapefile, or raster).", cxxopts::value<std::string>())
        ("o,output", "Path to location outputs should be written (licosim folder will be created). Defaults to lidar data directory.", cxxopts::value<std::string>())
        ("t,thread", "How many threads to process on.", cxxopts::value<int>()->default_value("1"))
        ("l,lidar", "Path to lidar dataset, in standard uw processing format. Defaults to current directory.", cxxopts::value<std::string>()->default_value("."))
        ("lmu", "Path to lmu raster (\".img\").  See documentation for expected format.", cxxopts::value<std::string>())
        ("subdivideLmus", "If present, subdivide lmus first by climate class and the by kmeans to a certain size")
        ("priority", "What method to prioritize units? Options are: column.", cxxopts::value<std::string>()->default_value("column"))
        ("column", "column number or name in shapefile that is the priority attribute (0 indexed. Will look for \"priority\" if blank", cxxopts::value<std::string>())
        ("r, reference", "Path to reference database, if not using one of the internal databases.", cxxopts::value<std::string>())
        ("allom", "Comma separted coefficients to a linear allometric equation (b0,b1,b2...), will derive them from FIA data if missing.", cxxopts::value<std::vector<double>>())
        ("a,aet", "Path to AET layer, to override default", cxxopts::value<std::string>())
        ("c,cwd", "Path to CWD layer, to override default", cxxopts::value<std::string>())
        ("j,janmin", "Path to Janmin layer, to override default", cxxopts::value<std::string>())
        //("usefire", "Use built in fire reference database??", cxxopts::value<bool>()->default_value("true"))
        //("usehydro", "Use built in hydro reference database??", cxxopts::value<bool>()->default_value("false"))
        //("usehabitat", "Use built in habitat reference database?", cxxopts::value<bool>()->default_value("false"))
        //("dofire", "Run fire modeling after treatment?", cxxopts::value<bool>()->default_value("false"))
        ("terrain", "What terrain type to generate LMUs based on? (steep or moderate)", cxxopts::value<std::string>()->default_value("moderate"))
        //("pland", "What proportion of the landscape should be treated? (0-1) Overrides tland.", cxxopts::value<double>())
        //("tland", "What is the total landscape area that should be treated (ha if units are meters, ac if units are feet)", cxxopts::value<double>())
        ("dbhmin", "All trees smaller than this size will be cut in cm. Default 15.24cm (6in)", cxxopts::value<double>())
        ("dbhmax", "All trees larger than this will be retained in cm. Default 53.34cm (21in)", cxxopts::value<double>())
        //("cover", "Three comma separated values between 0-1, representing proportion of cover left in each cover class: 0-40 (open), 40-60 (moderate), 60-100 (dense).", cxxopts::value<std::vector<double>>())
        ("s,seed", "Positive int value to seed the randomness.", cxxopts::value<int>())
        ("writeunits", "Write out lmus and rxunits to the output folder (treelists). WARNING can take up A LOT of space.")
        ("fastfuels", "write csvs in fastfuels format too if writeunits is set.")
        ("overridetargets", "Use Ba and DBH cutoffs defined in unit file to override internal targets")
        ("h,help", "Display this help message and exit.");

    std::string commandLine = "";
    for (int i = 0; i < argc; ++i) {
        commandLine += std::string(argv[i]) + " ";
    }

    if (argc == 1) {
        std::cout << opt.help() << '\n';
        exit(EXIT_SUCCESS);
    }

    auto options = opt.parse(argc, argv);

    if (options.count("help")) {
        std::cout << opt.help() << '\n';
        exit(EXIT_SUCCESS);
    }

    std::cout << "\n";
    std::cout << "---------------------------------------------------------\n";
    std::cout << "-------------- Licosim v0.5.1   03/01/2021 --------------\n";
    std::cout << "---------------------------------------------------------\n";
    std::cout << "\n";

    //auto searchpath = "";
    //proj_context_set_search_paths(nullptr, 1, &searchpath);

    auto x = proj_context_get_database_path(nullptr);
    std::cout << x << "\n";

    std::cout << "Beginning processing:\n";
    std::cout << "\tReading projectsettings...";
    auto ps = licosim::ProjectSettings(licosim::getExeLocation(argc, argv).string(), options);
    ps.commandLine = commandLine;
    std::cout << " Done!\n";
    
    std::cout << "\tConstructing primary data objects\n";
    auto ls = licosim::Licosim(ps);
    std::cout << "\tPrimary data objects done!\n";

    std::cout << "Beginning treatments on " + std::to_string(ps.nThread) + " threads:\n ";
    ls.doTreatmentThreaded(ps.nThread, ps.dbhMin, ps.dbhMax);
    std::cout << "Treatment done.\n";

    std::cout << "Writing outputs...";
    ls.output.write(ps.outputPath);
    std::cout << " Done!\n";
    std::cout << "Processing stopped!\n";

    return(0);
}
