// licosimexe.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "boost/program_options.hpp"
#include "licosim.hpp"
#include "ProjWrappers.hpp"

int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description desc("Licosim: Given lidar data and a lot of area parameters, calculate proposed LMU based treatment alternatives based on given reference conditions.");
    desc.add_options()
        ("projectpoly,p", po::value<std::string>()->required(), "REQUIRED: Path to project area polygon (esri shapefile).")
        ("unitpoly,u", po::value<std::string>()->required(), "REQUIRED: Path to unit polygons (esri shapefile, or raster).")
        ("output,o", po::value<std::string>(), "Path to location outputs should be written (licosim folder will be created). Defaults to lidar data directory.")
        ("thread,t", po::value<int>()->default_value(1), "How many threads to process on.")
        ("lidar,l", po::value<std::string>()->default_value("."), "Path to lidar dataset, in standard uw processing format. Defaults to current directory.")
        ("lmu", po::value<std::string>(), "Path to lmu raster (\".img\").  See documentation for expected format.")
        ("subdivideLmus", "If present, subdivide lmus first by climate class and the by kmeans to a certain size")
        ("priority", po::value<std::string>()->default_value("column"), "What method to prioritize units? Options are: column.")
        ("column", po::value<std::string>(), "column number or name in shapefile that is the priority attribute (0 indexed. Will look for \"priority\" if blank")
        ("reference,r", po::value<std::string>(), "Path to reference database, if not using one of the internal databases.")
        ("allom_slope", po::value<double>(), "Slope for a univariate linear model. Derived from FIA data if any parameter missing.")
        ("allom_intercept", po::value<double>(), "Intercept for a univariate linear model. Derived from FIA data if any parameter missing.")
        ("allom_rsq", po::value<double>(), "R-squared for a univariate linear model. Derived from FIA data if any parameter missing.")
        ("allom_transform", po::value<std::string>(), "One of (none, square, cube, power). Derived from FIA data if any parameter missing.")
        ("allom_inunits", po::value<double>(), "Input units (feet, meters). Derived from FIA data if any parameter missing.")
        ("allom_outunits", po::value<double>(), "Output units (inches, centimeters, meters). Derived from FIA data if any parameter missing.")
        ("aet,a", po::value<std::string>(), "Path to AET layer, to override default")
        ("cwd,c", po::value<std::string>(), "Path to CWD layer, to override default")
        ("janmin,j", po::value<std::string>(), "Path to Janmin layer, to override default")
        //("usefire", "Use built in fire reference database??", cxxopts::value<bool>()->default_value("true"))
        //("usehydro", "Use built in hydro reference database??", cxxopts::value<bool>()->default_value("false"))
        //("usehabitat", "Use built in habitat reference database?", cxxopts::value<bool>()->default_value("false"))
        //("dofire", "Run fire modeling after treatment?", cxxopts::value<bool>()->default_value("false"))
        ("terrain", po::value<std::string>()->default_value("moderate"), "What terrain type to generate LMUs based on? (steep or moderate)")
        //("pland", "What proportion of the landscape should be treated? (0-1) Overrides tland.", cxxopts::value<double>())
        //("tland", "What is the total landscape area that should be treated (ha if units are meters, ac if units are feet)", cxxopts::value<double>())
        ("dbhmin", po::value<double>(), "All trees smaller than this size will be cut in cm. Default 15.24cm (6in)")
        ("dbhmax", po::value<double>(), "All trees larger than this will be retained in cm. Default 53.34cm (21in)")
        //("cover", "Three comma separated values between 0-1, representing proportion of cover left in each cover class: 0-40 (open), 40-60 (moderate), 60-100 (dense).", cxxopts::value<std::vector<double>>())
        ("seed,s", po::value<int>(), "Positive int value to seed the randomness.")
        ("writeunits", "Write out lmus and rxunits to the output folder (treelists). WARNING can take up A LOT of space.")
        ("fastfuels", "write csvs in fastfuels format too if writeunits is set.")
        ("overridetargets", "Use Ba and DBH cutoffs defined in unit file to override internal targets")
        ("help,h", "Display this help message and exit.");

    std::string commandLine = "";
    for (int i = 0; i < argc; ++i) {
        commandLine += std::string(argv[i]) + " ";
    }

    if (argc == 1) {
        std::cout << desc << '\n';
        return 0;
    }

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    }
    catch (const po::error& e) {
        std::cerr << "Error: " << e.what() << "\n";
        std::cerr << desc << "\n";
        return 1;
    }

    if (vm.count("help")) {
        std::cout << desc << '\n';
        return 0;
    }

    std::cout << "\n";
    std::cout << "---------------------------------------------------------\n";
    std::cout << "-------------- Licosim v0.8.1   12/12/2025 --------------\n";
    std::cout << "---------------------------------------------------------\n";
    std::cout << "\n";

    std::cout << "Beginning processing:\n";
    std::cout << "\tReading projectsettings...";
    licosim::ProjectSettings::get().loadFromOptionsAndParentPath(lapis::executableFilePath(), vm);
    licosim::ProjectSettings::get().commandLine = commandLine;
    std::cout << " Done!\n";

    std::cout << "\tConstructing primary data objects\n";
    auto ls = licosim::Licosim();
    std::cout << "\tPrimary data objects done!\n";

    std::cout << "Beginning treatments on " + std::to_string(licosim::ProjectSettings::get().nThread) + " threads:\n ";
    ls.doTreatmentThreaded(licosim::ProjectSettings::get().nThread, licosim::ProjectSettings::get().dbhMin, licosim::ProjectSettings::get().dbhMax);
    std::cout << "Treatment done.\n";

    std::cout << "Writing outputs...";
    ls.output.write(licosim::ProjectSettings::get().outputPath);
    std::cout << " Done!\n";
    std::cout << "Processing stopped!\n";

    return(0);
}
