#include "licosim.hpp"
#include <chrono>

namespace licosim {

    Licosim::Licosim() {
        ProjectSettings* ps = &ProjectSettings::get();
        if (ps->seed > -1) {
            srand(ps->seed);
            dre = std::default_random_engine(ps->seed);
        }
        else {
            dre = std::default_random_engine(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
        }

        std::cout << "\tConstructing project area:\n";
        projectArea = rxtools::ProjectArea(ps->lidarDatasetPath, ps->projectPolygonPath, ps->aetPath, ps->cwdPath, ps->tmnPath, ps->nThread, ps->lmuRasterPath, ps->terrain);
        treater = rxtools::Treatment(dre);
        std::cout << "\tProject Area done!\n";

        if (ps->subdivideLmus) {
            projectArea.subdivideLmus(ps->climateClassPath, ps->nThread, ps->noClimateClasses);
        }

        output = rxtools::Output(projectArea.lmuRaster);
        output.lmus = projectArea.lmuRaster;
        output.lmuIds = projectArea.lmuIds;
        output.commandLine = ps->commandLine;

        std::cout << "\tReading and reprojecting unit layer...";
        if (std::filesystem::path(ps->unitPolygonPath).extension() != ".shp")
            shp = false;

        if (shp) {
            try {
                unitPoly = lapis::VectorDataset<lapis::MultiPolygon>(ps->unitPolygonPath);
            }
            catch (std::exception& e) {
                std::cerr << "Error reading unit polygon shapefile: " << e.what() << "\n";
                throw;
            }
            if (!unitPoly.crs().isConsistent(projectArea.projectPoly.crs()))
                unitPoly.projectInPlace(projectArea.projectPoly.crs());
            if (!unitPoly.extent().overlaps(projectArea.projectPoly.extent())) {
                std::cerr << "Unit polygon does not overlap project polygon\n";
                throw lapis::OutsideExtentException();
            }

        }
        else {
            try {
                unitRaster = lapis::Raster<int>(ps->unitPolygonPath);
            }
            catch (std::exception& e) {
                std::cerr << "Error reading unit polygon raster: " << e.what() << "\n";
                throw;
            }
            lapis::Alignment projectAlign = lapis::Alignment((lapis::Extent)ps->unitPolygonPath, unitRaster.nrow(), unitRaster.ncol());
            unitRaster = lapis::resampleRaster(unitRaster, projectAlign, lapis::ExtractMethod::near);
            unitRaster = lapis::trimRaster(unitRaster);
            if (!unitRaster.overlaps(projectArea.projectPoly.extent())) {
                std::cerr << "Unit raster does not overlap project polygon\n";
                throw lapis::OutsideExtentException();
            }
        }
        std::cout << " Done!\n";

        std::cout << "\tCreating allometry from ";
        auto before = std::chrono::high_resolution_clock::now();

        if (ps->slope > 0) {
            std::cout << "user provided coefficients...";
            //dbhModel = rxtools::allometry::DbhModel(ps->slope, ps->intercept, ps->transform, ps->rsq, ps->inUnit, ps->outUnit);
            std::cout << "  Deprecated function for now as as allometry rasters have been added\n";
        }
        else {
            std::cout << "FIA plots...";
            std::unique_ptr<rxtools::allometry::FIAReader> reader;
            try {
                reader = std::make_unique<rxtools::allometry::FIAReader>(ps->fiaPath);
            }
            catch (std::exception& e) {
                std::cerr << "Error reading FIA data: " << e.what() << "\n";
                throw;
            }
            auto dist = projectArea.lidarDataset->units().value().convertOneToThis(10000, lapis::linearUnitPresets::meter);
            lapis::Extent e(projectArea.projectPoly.extent().xmin() - dist, projectArea.projectPoly.extent().xmax() + dist,
                projectArea.projectPoly.extent().ymin() - dist, projectArea.projectPoly.extent().ymax() + dist,
                projectArea.projectPoly.crs());
            if (!reader->limitByExtent(e))
                throw std::runtime_error("no fia in this place");
            reader->makePlotTreeMap(std::vector<std::string>{ "DIA" });
            auto allTrees = reader->collapsePlotTreeMap();

            try {
                allTrees.writeCsv(ProjectSettings::get().outputPath + "/fia_tree_data.csv");
            }
            catch (std::exception& e) {
                std::cerr << "Error writing FIA tree data to CSV: " << e.what() << "\n";
                throw;
            }
            
            dbhAllomRaster = std::make_shared<rxtools::allometry::AllometryRaster>(rxtools::allometry::calculateAllometryOverAOI<rxtools::allometry::DbhModel>(
                projectArea.lmuRaster,
                *reader,
                "DIA",
                lapis::linearUnitPresets::internationalInch,
                std::nullopt,
                std::nullopt,
                ps->get().nThread
            ));

            if (ps->writeUnits && ps->fastFuels) {
                ffa = rxtools::allometry::FastFuels(allTrees);
            }
        }
        std::cout << " Done!\n";
        auto after = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = after - before;
        std::cout << " Done! Time taken: " << elapsed.count() << " seconds\n";

        lapis::Raster<int> tmp((lapis::Alignment)*dbhAllomRaster);
        for (lapis::cell_t c = 0; c < tmp.ncell(); ++c) {
            if ((*dbhAllomRaster)[c].has_value()) {
                tmp[c].has_value() = true;
                tmp[c].value() = 0;
            }
        }
        lapis::Raster<int> lmu((lapis::Alignment)projectArea.lmuRaster);
        for (lapis::cell_t c = 0; c < lmu.ncell(); ++c) {
            if (projectArea.lmuRaster[c].has_value()) {
                lmu[c].has_value() = true;
                lmu[c].value() = 0;
            }

        }
        lmu.writeRaster("F:/lmu_raster.tif");

        std::cout << "Setting taos...";
        before = std::chrono::high_resolution_clock::now();
        auto dbhG = [
            dr = dbhAllomRaster,
            hg = projectArea.lidarDataset->heightGetterPolygon(),
            cg = projectArea.lidarDataset->coordGetterPolygon()
        ](const lapis::ConstFeature<lapis::MultiPolygon>& ft)->double {
            auto xy = cg(ft);
            lapis::cell_t c;
            try { //since we load whole tiles, some taos will be outside of allometry
                c = dr->cellFromXY(xy.x, xy.y);
            }
            catch (lapis::OutsideExtentException e) {
                return(0);
            }
            if ((*dr)[c].has_value()) {
                return (*dr)[c].value()->predict(hg(ft), lapis::linearUnitPresets::meter, lapis::linearUnitPresets::centimeter);
            }
            else {
                return(0);
            }
        };
        projectArea.setTaos(dbhG, licosim::ProjectSettings::get().nThread, licosim::ProjectSettings::get().fixedRadiusMeters);
        std::cout << " Done!\n";
        after = std::chrono::high_resolution_clock::now();
        elapsed = after - before;
        std::cout << " Done! Time taken: " << elapsed.count() << " seconds\n";

        std::cout << "\t Performing PCA on climate data...";
        std::unique_ptr<csv::CSVReader> csv;
        try {
            if (ps->referenceDatasetPath != "") {
                csv = std::make_unique<csv::CSVReader>(ps->referenceDatasetPath);
            }
            else {
                csv = std::make_unique<csv::CSVReader>(ps->defaultRefPath);
            }
        }
        catch (std::exception& e) {
            std::cerr << "Error reading reference dataset CSV: " << e.what() << "\n";
            throw;
        }
        for (auto& row : *csv) {
            names.push_back(row["name"].get<>());
            ids.push_back(row["id"].get<int>());
            types.push_back(static_cast<rxtools::LmuType>(row["type"].get<int>()));
            refStructures.push_back(
                rxtools::StructureSummary(
                    row["ba"].get<double>(), row["tph"].get<double>(), row["mcs"].get<double>(), row["cc"].get<double>(),
                    std::vector<double>{ 
                        row["c1"].get<double>(), row["c2to4"].get<double>(), row["c5to9"].get<double>(),
                        row["c10to14"].get<double>(), row["c15to29"].get<double>(), row["c30to35"].get<double>()
                    }
                )
            );
            aet.push_back(row["aet"].get<double>());
            cwd.push_back(row["cwd"].get<double>());
            tmn.push_back(row["tmn"].get<double>());
        }

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
        rTmpClim = lapis::cropRaster(rTmpClim, lmu.mask, lapis::SnapType::out);
        rTmpClim = lapis::extendRaster(rTmpClim, lmu.mask, lapis::SnapType::out);
        rTmpClim.mask(lmu.mask);
        double thisaet = 0;
        int aetcount = 0;
        for (lapis::cell_t c = 0; c < rTmpClim.ncell(); ++c) {
            if (rTmpClim[c].has_value()) {
                thisaet += rTmpClim[c].value();
                aetcount++;
            }
        }
        thisaet /= aetcount;

        rTmpClim = lapis::Raster(projectArea.cwd);
        rTmpClim = lapis::cropRaster(rTmpClim, lmu.mask, lapis::SnapType::out);
        rTmpClim = lapis::extendRaster(rTmpClim, lmu.mask, lapis::SnapType::out);
        rTmpClim.mask(lmu.mask);
        double thiscwd = 0;
        int cwdcount = 0;
        for (lapis::cell_t c = 0; c < rTmpClim.ncell(); ++c) {
            if (rTmpClim[c].has_value()) {
                thiscwd += rTmpClim[c].value();
                cwdcount++;
            }
        }
        thiscwd /= cwdcount;

        rTmpClim = lapis::Raster(projectArea.tmn);
        rTmpClim = lapis::cropRaster(rTmpClim, lmu.mask, lapis::SnapType::out);
        rTmpClim = lapis::extendRaster(rTmpClim, lmu.mask, lapis::SnapType::out);
        rTmpClim.mask(lmu.mask);
        double thistmn = 0;
        int tmncount = 0;
        for (lapis::cell_t c = 0; c < rTmpClim.ncell(); ++c) {
            if (rTmpClim[c].has_value()) {
                thistmn += rTmpClim[c].value();
                tmncount++;
            }
        }
        thistmn /= tmncount;

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
        else if (type == rxtools::LmuType::neFacing)
            return neFacingLmu.detailedQuery(pcClimate, k, false, maxDist);
        else 
            return allLmu.detailedQuery(pcClimate, k, false, maxDist);
    }

    void Licosim::assignTargetThread(size_t& sofar, size_t& nLmu, rxtools::Lmu& lmu, const int thisThread) {
        std::cout << "\t Assigning targets to Lmu " + std::to_string(sofar+1) + "/" + std::to_string(nLmu) + " on thread " + std::to_string(thisThread) + "\n";
        auto k = getKnn(lmu, lmu.type, maxDist, 20);
        if (k.size() == 0) {
            k = getKnn(lmu, rxtools::LmuType::all, maxDist, 20);
        }
        if (k.size() == 0) {
            k = getKnn(lmu, rxtools::LmuType::all, 0, 20); //defaults to max dist = max double, 10 samples.
        }

        for (int i : k) {
            lmu.targetLmuNames.push_back(names[i] + "_" + std::to_string(ids[i]));
            lmu.structures.push_back(refStructures[i]);
        }
        lmu.assignUnitTargets(dre, ProjectSettings::get().dbhMax, ProjectSettings::get().overrideTargets);
    }

    void Licosim::doTreatmentThreaded(int nThread, double dbhMin, double dbhMax) {
        std::filesystem::path p(ProjectSettings::get().outputPath);
        if (ProjectSettings::get().writeUnits) {
            p /= "lmus";
            if (!std::filesystem::exists(p)) {
                std::filesystem::create_directories(p);
            }
        }

        lapis::Raster<lapis::cell_t> unitZonal{ (lapis::Alignment)projectArea.lmuIds };
        paired = lapis::Raster<lapis::cell_t>{ (lapis::Alignment)projectArea.lmuIds };

        rxtools::TaoList treatedTaos{ projectArea.getTaos().at(0).second.crs()};
        std::mutex mut{};
        size_t sofar = 0;
        std::vector<std::thread> threads{};
        auto threadFunc = [&](int i) { treatmentThread(sofar, mut, dbhMin, dbhMax, unitZonal, treatedTaos, i); };
        for (int i = 0; i < nThread; ++i) {
            threads.push_back(std::thread(threadFunc, i));
        }
        for (int i = 0; i < nThread; ++i) {
            threads[i].join();
        }
        std::cout << "Treatment done\n";
    }

    void Licosim::treatmentThread(size_t& sofar, std::mutex& mut, double dbhMin, double dbhMax, lapis::Raster<lapis::cell_t>& unitZonal, rxtools::TaoList& treatedTaos, const int thisThread) {
        size_t nLmu = projectArea.regionType.size();
        std::filesystem::path p(ProjectSettings::get().outputPath);
        p /= "lmus";

        int failedToAdd = 0;
        std::set<int> failedIds;

        while (true) {

            mut.lock();
            size_t i = sofar;
            ++sofar;    
            mut.unlock();
            //if (i != 4507) continue;

            if (i >= nLmu) {
                break;
            }

            //try {
                auto before = std::chrono::high_resolution_clock::now();
                //if (i != 6461) continue;
                auto lmu = projectArea.createLmuThread(i, thisThread);
                auto e = lapis::Extent(projectArea.aet.xmin() + projectArea.aet.xres() / 2,
                    projectArea.aet.xmax() - projectArea.aet.xres() / 2,
                    projectArea.aet.ymin() + projectArea.aet.yres() / 2,
                    projectArea.aet.ymax() - projectArea.aet.yres() / 2);
                if (!lmu.mask.overlaps(e)) continue;

                rxtools::TaoList taos{ projectArea.getTaos().at(0).second.crs()};
                for (size_t j = 0; j < projectArea.getTaos().size(); ++j) {
                    auto& thisTile = projectArea.getTaos().at(j);
                    if (thisTile.first.overlaps(lmu.mask)) {
                        for (size_t k = 0; k < thisTile.second.size(); ++k) {
                            if (lmu.mask.extract(thisTile.second.x(k), thisTile.second.y(k), lapis::ExtractMethod::near).has_value())
                                taos.addTao(thisTile.second.xy(k),
                                    thisTile.second.height(k),
                                    thisTile.second.radius(k),
                                    thisTile.second.area(k),
                                    thisTile.second.dbh(k));
                        }
                    }
                }
                std::cout << "LMU taos size: " << taos.size() << "\n";

                auto after = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
                std::cout << "Thread " + std::to_string(thisThread) + " completed in " + std::to_string(duration.count()) + " seconds.\n";

                std::cout << "\tCreating Units for Lmu " + std::to_string(i) + "/" + std::to_string(nLmu) + " on thread " + std::to_string(thisThread) + "\n";
                before = std::chrono::high_resolution_clock::now();

                if (shp)
                    lmu.makeUnits(unitPoly, taos, ProjectSettings::get().overrideTargets);
                else
                    lmu.makeUnits(unitRaster, taos);

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
                std::cout << "\tTreating Lmu " + std::to_string(i+1) + "/" + std::to_string(nLmu) + " on thread " + std::to_string(thisThread) + "\n";
                before = std::chrono::high_resolution_clock::now();
                for (int j = 0; j < lmu.units.size(); ++j) {
                    if (ProjectSettings::get().overrideTargets) {
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
                        size_t a = lmu.units[j].taos.size();
                        //try {
                            auto trt = treater.doTreatment(lmu.units[j], lmu.units[j].dbhMin, lmu.units[j].dbhMax, 3);
                            lmu.units[j].treatedTaos = std::get<0>(trt);
                        /* }
                        catch (std::exception e) {
                            std::cout << e.what();
                            std::cout << " duplicates in treat " + std::to_string(i) + "  " + std::to_string(j) + "\n";
                            std::filesystem::create_directory(p / std::to_string(i));
                            lmu.write((p / std::to_string(i)).string(), ffa);
                            throw std::runtime_error("you've been dooped");
                        }*/
                        size_t b = lmu.units[j].treatedTaos.size();

                        lmu.units[j].treatedStructure = rxtools::StructureSummary(lmu.units[j].treatedTaos, lmu.units[j].unitMask, lmu.units[j].areaHa);
                        size_t c = lmu.units[j].treatedTaos.size();

                        if (lmu.units[j].targetStructure.ba - lmu.units[j].treatedStructure.ba > 1) {
                            mut.lock();
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
                    for (lapis::cell_t c = 0; c < lmu.units[j].unitMask.ncell(); ++c) {
                        lmu.units[j].unitMask[c].value() = i * 10000 + j;
                    }

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
                    for (size_t k = 0; k < lmu.units[j].treatedTaos.size(); ++k) {
                        treatedTaos.addTao(lmu.units[j].treatedTaos.xy(k),
                            lmu.units[j].treatedTaos.height(k),
                            lmu.units[j].treatedTaos.radius(k),
                            lmu.units[j].treatedTaos.area(k),
                            lmu.units[j].treatedTaos.dbh(k));
                    }

                    unitZonal.overlayInside(lmu.units[j].unitMask);

                    auto tmp = lmu.units[j].unitMask;
                    for (lapis::cell_t c = 0; c < tmp.ncell(); ++c) {
                        tmp[c].value() = lmu.units[j].paired;
                    }
                    paired.overlayInside(tmp);

                    mut.unlock();
                }
                after = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
                std::cout << "Thread " + std::to_string(thisThread) + " completed in " + std::to_string(duration.count()) + " seconds.\n";

                if (ProjectSettings::get().writeUnits) {
                    std::filesystem::create_directory(p / std::to_string(i));
                    lmu.write((p / std::to_string(i)).string(), ffa);
                }
            /*}
            catch (std::exception e) {
                std::cerr << e.what() << " on thread " << std::to_string(thisThread);
                throw e;
            }*/
        }
        /*std::cout << "Failed to add " << failedToAdd << "Treatment units\n";
        for (int i : failedIds) {
            std::cout << i << "\n";
        }*/
    }
}