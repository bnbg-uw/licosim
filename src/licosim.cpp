#include "licosim.hpp"
#include <chrono>

namespace licosim {

    Licosim::Licosim() {
        if (ps.seed > -1) {
            srand(ps.seed);
            dre = std::default_random_engine(ps.seed);
        }
        else {
            dre = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());
        }

        std::cout << "\tConstructing project area:\n";
        projectArea = rxtools::ProjectArea(ps.lidarDatasetPath, ps.projectPolygonPath, ps.aetPath, ps.cwdPath, ps.tmnPath, ps.nThread, ps.lmuRasterPath, ps.terrain);
        treater = rxtools::Treatment(dre);
        std::cout << "\tProject Area done!\n";

        if (ps.subdivideLmus) {
            projectArea.subdivideLmus(ps.climateClassPath, ps.nThread);
        }

        output.lmus = projectArea.lmuRaster;
        output.lmuIds = projectArea.lmuIds;
        output.commandLine = projectSettings.commandLine;

        std::cout << "\tCreating allometry from ";
        if(ps.allom_coefficients.size()) {
            std::cout << "user provided coefficients...";
            allometry = Allometry(ps.allom_coefficients);
            if(ps.allomPower)
                dbhFunc = [this](lico::adapt_type<spatial::unit_t> ht) {return allometry.getDbhFromHeightAuto(ht); };
            else
                dbhFunc = [this](lico::adapt_type<spatial::unit_t> ht) {return allometry.getDbhFromHeightLinear(ht); };

            std::cout << "  Done!\n";
        }
        else {
            std::cout << "FIA plots...";
            allometry = Allometry(ps.fiaPath);
            auto backup = allometry;
            if(!allometry.model.limitByExtent(projectArea.projectPoly)) {
                std::cout << "No fia plots in study area, falling back to 5km buffer around study area...";
                auto conv = projectArea.lidarDataset->getConvFactor();
                spatial::Extent e(projectArea.projectPoly.xmin() - 5000 * conv, projectArea.projectPoly.xmax() + 5000 * conv,
                                  projectArea.projectPoly.ymin() - 5000 * conv, projectArea.projectPoly.ymax() + 5000 * conv,
                                  projectArea.projectPoly.projection());
                if (!backup.model.limitByExtent(e))
                    throw std::runtime_error("no fia in this place");
                allometry = backup;
            }
            allometry.model.initializeModel();

            dbhFunc = [this](lico::adapt_type<spatial::unit_t> ht) { return allometry.getDbhFromHeightAuto(ht); };
            std::cout << " Done!\n";
            std::cout << allometry.model.intercept << " " << allometry.model.slope << " " << int(allometry.model.transform) << "\n";
        }

        if (projectSettings.writeUnits && projectSettings.fastFuels) {
            ffa = rxtools::allometry::FastFuels(); wrong input to this function;
        }

        std::cout << "\tReading and reprojecting unit layer...";
        if (std::filesystem::path(projectSettings.unitPolygonPath).extension() != ".shp")
            shp = false;

        if (shp) {
            unitPoly = lapis::VectorDataset<lapis::MultiPolygon>(projectSettings.unitPolygonPath);
            if (!lapis::consistentProjection(unitPoly.crs(), projectArea.projectPoly.crs()))
                unitPoly.project(projectArea.projectPoly.projection());
            if (!unitPoly.overlaps(projectArea.projectPoly)) {
                std::cerr << "Unit polygon does not overlap project polygon\n";
                throw lapis::OutsideExtentException();
            }

        }
        else {
            unitRaster = lapis::Raster<int>(projectSettings.unitPolygonPath);
            lapis::Alignment projectAlign = lapis::Alignment((lapis::Extent)projectSettings.unitPolygonPath, unitRaster.nrow(), unitRaster.ncol());
            unitRaster = lapis::resampleNGB(unitRaster, projectAlign);
            unitRaster.trim();
            if (!unitRaster.overlaps(projectArea.projectPoly)) {
                std::cerr << "Unit raster does not overlap project polygon\n";
                throw lapis::OutsideExtentException();
            }
        }
        std::cout << " Done!\n";


        projectArea.createCoreGapAndReadTaos(projectSettings.nThread, projectSettings.dbhMax, dbhFunc);

        std::cout << "\t Performing PCA on climate data...";
        std::filebuf fb;
        if (ps.referenceDatasetPath != "") {
            if (!fb.open(ps.referenceDatasetPath, std::ios::in)) throw std::runtime_error("Cannot open reference table.");
        }
        else {
            if (!fb.open(ps.defaultRefPath, std::ios::in)) throw std::runtime_error("Cannot open internal reference table.");
        }
        std::istream is{ &fb };

        rxtools::utilities::readCSVLine(is); //skip colnames.
        while (!is.eof()) {
            auto row = rxtools::utilities::readCSVLine(is);
            if (row.size() <= 1)
                continue;
            names.push_back(row[0]);
            ids.push_back(std::stoi(row[1]));
            types.push_back(static_cast<rxtools::LmuType>(std::stoi(row[2])));
            refStructures.push_back(rxtools::StructureSummary(std::stod(row[3]), std::stod(row[4]), std::stod(row[5]), std::stod(row[6]), std::stod(row[7]),
                std::vector<double>{std::stod(row[8]), std::stod(row[9]), std::stod(row[10]), std::stod(row[11]), std::stod(row[12]), std::stod(row[13])}));
            aet.push_back(std::stod(row[14]));
            cwd.push_back(std::stod(row[15]));
            tmn.push_back(std::stod(row[16]));
        }
        fb.close();

        //Performing PCA.
        Eigen::MatrixXd data(aet.size(), 3);
        for (int i = 0; i < aet.size(); ++i) {
            data(i, 0) = aet.at(i);
            data(i, 1) = cwd.at(i);
            data(i, 2) = tmn.at(i);
        }
        pca = rxtools::utilities::Pca(data);

        //Calculate distance weighted by axis importance.
        for (int i = 0; i < pca.data.rows(); i++) {
            for (int j = i + 1; j < pca.data.rows(); j++) {
                double thisDist = 0;
                for (int k = 0; k < pca.data.cols(); k++) {
                    thisDist += (pca.data(i, k) - pca.data(j, k)) * (pca.data(i, k) - pca.data(j, k)) * pca.importance(k);
                }
                if (thisDist > maxDist) {
                    maxDist = thisDist;
                }
            }
        }

        maxDist = std::sqrt(maxDist);
        maxDist /= 5;

        allLmu = rxtools::utilities::PCA_KDTree<3>(pca);

        std::vector<bool> use(pca.data.rows(), true);
        for (int i = 0; i < types.size(); ++i)
            if (types[i] != rxtools::LmuType::ridgeTop)
                use[i] = false;
        ridgeLmu = rxtools::utilities::PCA_KDTree<3>(pca, use);

        for (int i = 0; i < types.size(); ++i) {
            if (types[i] == rxtools::LmuType::valleyBottom)
                use[i] = true;
            else
                use[i] = false;
        }
        valleyLmu = rxtools::utilities::PCA_KDTree<3>(pca, use);

        for (int i = 0; i < types.size(); ++i) {
            if (types[i] == rxtools::LmuType::swFacing)
                use[i] = true;
            else
                use[i] = false;
        }
        swFacingLmu = rxtools::utilities::PCA_KDTree<3>(pca, use);

        for (int i = 0; i < types.size(); ++i) {
            if (types[i] == rxtools::LmuType::neFacing)
                use[i] = true;
            else
                use[i] = false;
        }
        neFacingLmu = rxtools::utilities::PCA_KDTree<3>(pca, use);

        std::cout << "\tPCA done and reference area distance threshold done!\n";
    }

    std::vector<int> Licosim::getKnn(rxtools::Lmu& lmu, const rxtools::LmuType& type, const double& maxDist, const int k) {
        auto rTmpClim = projectArea.aet;
        rTmpClim.crop(lmu.mask);
        rTmpClim.extend(lmu.mask);
        rTmpClim.mask(lmu.mask);
        auto thisaet = xt::mean(xt::filter(rTmpClim.values().value(), rTmpClim.values().has_value()))();

        rTmpClim = lapis::Raster(projectArea.cwd);
        rTmpClim.crop(lmu.mask);
        rTmpClim.extend(lmu.mask);
        rTmpClim.mask(lmu.mask);
        auto thiscwd = xt::mean(xt::filter(rTmpClim.values().value(), rTmpClim.values().has_value()))();

        rTmpClim = lapis::Raster(projectArea.tmn);
        rTmpClim.crop(lmu.mask);
        rTmpClim.extend(lmu.mask);
        rTmpClim.mask(lmu.mask);
        auto thistmn = xt::mean(xt::filter(rTmpClim.values().value(), rTmpClim.values().has_value()))();

        Eigen::RowVector3d pcClimate; 
        pcClimate << thisaet, thiscwd, thistmn;

        if (type == rxtools::LmuType::all)
            return allLmu.detailedQuery(pcClimate, k, false, maxDist);
        else if (type == rxtools::LmuType::ridgeTop)
            return ridgeLmu.detailedQuery(pcClimate, k, false, maxDist);
        else if (type == rxtools::LmuType::valleyBottom)
            return valleyLmu.detailedQuery(pcClimate, k, false, maxDist);
        else if (type == rxtools::LmuType::swFacing)
            return swFacingLmu.detailedQuery(pcClimate, k, false, maxDist);
        else
            return neFacingLmu.detailedQuery(pcClimate, k, false, maxDist);
    }

    void Licosim::assignTargetThread(int& sofar, int& nLmu, rxtools::Lmu& lmu, const int thisThread) {
        std::cout << "\t Assigning targets to Lmu " + std::to_string(sofar) + "/" + std::to_string(nLmu) + " on thread " + std::to_string(thisThread) + "\n";
        auto k = getKnn(lmu, lmu.type, maxDist, 20);
        if (k.size() == 0) {
            k = getKnn(lmu, rxtools::LmuType::all, maxDist, 20);
        }
        if (k.size() == 0) {
            k = getKnn(lmu, rxtools::LmuType::all, 0, 20); //defaults to max dist = max double, 10 samples.
        }

        auto thisOsiNum = lapis::crop(projectArea.bbOsiNum, lmu.mask);
        auto thisOsiDen = lapis::crop(projectArea.bbOsiDen, lmu.mask);
        for (int i : k) {
            lmu.targetLmuNames.push_back(names[i] + "_" + std::to_string(ids[i]));
            lmu.structures.push_back(refStructures[i]);
        }
        lmu.assignUnitTargets(dre, projectSettings.dbhMax, projectSettings.overrideTargets);
    }

    void Licosim::doTreatmentThreaded(int nThread, double dbhMin, double dbhMax) {
        std::filesystem::path p(projectSettings.outputPath);
        if (projectSettings.writeUnits) {
            p /= "lmus";
            if (!std::filesystem::exists(p)) {
                std::filesystem::create_directories(p);
            }
        }

        lapis::Raster<int> unitZonal{ (lapis::Alignment)projectArea.lmuIds };
        paired = lapis::Raster<int>{ (lapis::Alignment)projectArea.lmuIds };

        rxtools::TaoListMP treatedTaos{};
        std::mutex mut{};
        int sofar = 0;
        std::vector<std::thread> threads{};
        auto threadFunc = [&](int i) { treatmentThread(sofar, mut, dbhMin, dbhMax, unitZonal, treatedTaos, i); };
        for (int i = 0; i < nThread; ++i) {
            threads.push_back(std::thread(threadFunc, i));
        }
        for (int i = 0; i < nThread; ++i) {
            threads[i].join();
        }

        std::cout << "Calculating post treatment OSI\n";
        //Post osi
        auto postNum = lapis::Raster<int>{ (lapis::Alignment)projectArea.lmuIds };
        auto postDen = lapis::Raster<int>{ (lapis::Alignment)projectArea.lmuIds };
        std::pair<lapis::coord_t, lapis::coord_t> expectedRes{};
        if (projectArea.lidarDataset->getConvFactor() == 1) {
            expectedRes.first = 0.75;
        }
        else {
            expectedRes.first = 2.4606;
        }
        expectedRes.second = 0;
        threads.clear();
        sofar = 0;
        auto postGapFunc = [&](int i) { projectArea.postGapThread(postNum, postDen, treatedTaos, projectSettings.nThread, i, mut, sofar, projectArea.lmuRaster, 2, 6, expectedRes); };
        for (int i = 0; i < nThread; ++i) {
            threads.push_back(std::thread(postGapFunc, i));
        }
        for (int i = 0; i < nThread; ++i) {
            threads[i].join();
        }

        std::cout << "Calculating post OSI zones...\n";
        lsmetrics::zonal_function<xtl::xoptional<int>, int> zSum = lsmetrics::zonalNoDataSum<int>;
        auto numZones = lsmetrics::zonalStatisticsByRaster(postNum, unitZonal, zSum);
        auto denZones = lsmetrics::zonalStatisticsByRaster(postDen, unitZonal, zSum);

        lapis::Raster<double> postOsi{ (spatial::Alignment)unitZonal };
        for (auto z : numZones) {
            bool hv = z.second.has_value();
            if (hv) {
                double osi = static_cast<double>(z.second.value()) / static_cast<double>(denZones.find(z.first)->second.value()) * 100;
                xt::filtration(postOsi.values().value(), xt::equal(unitZonal.values().value(), z.first)) = osi;
                xt::filtration(postOsi.values().has_value(), xt::equal(unitZonal.values().value(), z.first)) = true;
            }
        }
        postOsi.crop(output.pre[3]);
        postOsi.extend(output.pre[3]);
        postOsi.mask(output.pre[3]);
        output.post[3] = postOsi;

        auto ft = output.shp.getFeaturesPtr();
        for (size_t i = 0; i < ft->size(); ++i) {
            auto v = std::stoi(ft->at(i).getAttribute("id"));
            if (numZones.find(v) == numZones.end()) {
                std::cout << "failed to find: " << v << "\n";
                continue;
            }
            double osi = static_cast<double>(numZones.at(v).value()) / static_cast<double>(denZones.at(v).value()) * 100;
            ft->at(i).setAttribute(18, osi);
        }
        std::cout << "Treatment done\n";
    }

    void Licosim::treatmentThread(int& sofar, std::mutex& mut, double dbhMin, double dbhMax, lapis::Raster<int>& unitZonal, lico::TaoList& treatedTaos, const int thisThread) {
        int nLmu = projectArea.regionType.size();
        std::filesystem::path p(projectSettings.outputPath);
        p /= "lmus";

        int failedToAdd = 0;
        std::set<int> failedIds;

        while (true) {

            mut.lock();
            int i = sofar;
            ++sofar;    
            mut.unlock();
            //if (i != 4507) continue;

            if (i >= nLmu) {
                break;
            }

            try {
                auto before = std::chrono::high_resolution_clock::now();
                //if (i != 6461) continue;
                auto lmu = projectArea.createLmuThread(i, thisThread);
                auto e = lapis::Extent(projectArea.aet.xmin() + projectArea.aet.xres() / 2,
                    projectArea.aet.xmax() - projectArea.aet.xres() / 2,
                    projectArea.aet.ymin() + projectArea.aet.yres() / 2,
                    projectArea.aet.ymax() - projectArea.aet.yres() / 2);
                if (!lmu.mask.overlaps(e)) continue;

                rxtools::TaoListMP taos{};
                for (size_t j = 0; j < projectArea.allTaos.size(); ++j) {
                    if (lmu.mask.extract(projectArea.allTaos.x(j), projectArea.allTaos.y(j)).has_value())
                        taos.addTAO(projectArea.allTaos[j]);
                }

                auto after = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
                std::cout << "Thread " + std::to_string(thisThread) + " completed in " + std::to_string(duration.count()) + " seconds.\n";

                std::cout << "\tCreating Units for Lmu " + std::to_string(i) + "/" + std::to_string(nLmu) + " on thread " + std::to_string(thisThread) + "\n";
                before = std::chrono::high_resolution_clock::now();
                auto thisOsiNum = lapis::crop(projectArea.osiNum, lmu.mask);
                auto thisOsiDen = lapis::crop(projectArea.osiDen, lmu.mask);

                if (shp)
                    lmu.makeUnits(unitPoly, taos, thisOsiNum, thisOsiDen, projectArea.lidarDataset->getConvFactor(), dbhFunc, projectSettings.overrideTargets);
                else
                    lmu.makeUnits(unitRaster, taos, thisOsiNum, thisOsiDen, projectArea.lidarDataset->getConvFactor(), dbhFunc);

                after = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
                std::cout << "Thread " + std::to_string(thisThread) + " completed in " + std::to_string(duration.count()) + " seconds.\n";

                //assign targets here.
                before = std::chrono::high_resolution_clock::now();
                assignTargetThread(i, nLmu, lmu, thisThread);
                after = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
                std::cout << "Thread " + std::to_string(thisThread) + " completed in " + std::to_string(duration.count()) + " seconds.\n";

                //treat like below
                std::cout << "\tTreating Lmu " + std::to_string(i) + "/" + std::to_string(nLmu) + " on thread " + std::to_string(thisThread) + "\n";
                before = std::chrono::high_resolution_clock::now();
                for (int j = 0; j < lmu.units.size(); ++j) {
                    if (projectSettings.overrideTargets) {
                        if (lmu.units[j].dbhMin < 0) {
                            lmu.units[j].dbhMin = dbhMin;
                            lmu.units[j].dbhMax = dbhMax;
                        }
                    }
                    else {
                        lmu.units[j].dbhMin = dbhMin;
                        lmu.units[j].dbhMax = dbhMax;
                    }

                    if (lmu.units[j].currentStructure.ba > lmu.units[j].targetStructure.ba) {
                        int a = lmu.units[j].taos.size();
                        try {
                            auto trt = treater.doTreatment(lmu.units[j], lmu.units[j].dbhMin, lmu.units[j].dbhMax, 3);
                            lmu.units[j].treatedTaos = std::get<0>(trt);
                        }
                        catch (std::exception e) {
                            std::cout << e.what();
                            std::cout << " duplicates in treat " + std::to_string(i) + "  " + std::to_string(j) + "\n";
                            std::filesystem::create_directory(p / std::to_string(i));
                            lmu.write((p / std::to_string(i)).string(), ffa);
                            //throw std::runtime_error("you've been dooped");
                        }
                        int b = lmu.units[j].treatedTaos.size();

                        lmu.units[j].treatedStructure = lmu.units[j].summarizeStructure(lmu.units[j].treatedTaos, 0, lmu.units[j].dbhFunc);
                        int c = lmu.units[j].treatedTaos.size();

                        if (lmu.units[j].targetStructure.ba - lmu.units[j].treatedStructure.ba > 1) {
                            std::cout << "Post treat <<< target ba " + std::to_string(i) + "  " + std::to_string(j) + " " + std::to_string(a) + " " + std::to_string(b) + " " + std::to_string(c) + "\n";
                            std::filesystem::create_directory(p / std::to_string(i));
                            lmu.write((p / std::to_string(i)).string(), ffa);
                            std::cout << (p / std::to_string(i)).string() << "\n";
                            mut.unlock();
                            throw std::runtime_error("break this shit");
                        }
                    }
                    else {
                        lmu.units[j].treatedTaos = lmu.units[j].taos;
                        lmu.units[j].treatedStructure = lmu.units[j].currentStructure;
                    }
                    lmu.units[j].treated = true;
                    lmu.units[j].unitMask.values().value().fill(i * 10000 + j);

                    std::vector<lapis::Raster<int>*> v = { &lmu.units[j].unitMask, &unitZonal };

                    mut.lock();
                    std::cout << "Thread " + std::to_string(thisThread) + " adding rxunit. i = " + std::to_string(i) + "\n";
                    //try {
                        output.addRxUnit(lmu.units[j], static_cast<int>(lmu.type));
                    /* }
                    catch (std::exception e) {
                        failedToAdd++;
                        failedIds.emplace(i * 10000 + j);
                        mut.unlock();
                        continue;
                    }*/
                    for (size_t k = 0; k < lmu.units[j].treatedTaos.size(); ++k)
                        treatedTaos.addTAO(lmu.units[j].treatedTaos[k]);
                    unitZonal = lapis::rasterMergeInterior(v);

                    auto tmp = lmu.units[j].unitMask;
                    tmp.values().value().fill(lmu.units[j].paired);
                    v = { &tmp, &paired };
                    paired = lapis::rasterMergeInterior(v);

                    mut.unlock();
                }
                after = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
                std::cout << "Thread " + std::to_string(thisThread) + " completed in " + std::to_string(duration.count()) + " seconds.\n";

                if (projectSettings.writeUnits) {
                    std::filesystem::create_directory(p / std::to_string(i));
                    lmu.write((p / std::to_string(i)).string(), ffa);
                }
            }
            catch (std::exception e) {
                std::cerr << e.what() << " on thread " << std::to_string(thisThread);
                throw e;
            }
        }
        /*std::cout << "Failed to add " << failedToAdd << "Treatment units\n";
        for (int i : failedIds) {
            std::cout << i << "\n";
        }*/
    }
}