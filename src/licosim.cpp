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
            projectArea.subdivideLmus(ps->climateClassPath, ps->nThread);
        }

        output = rxtools::Output(projectArea.lmuRaster);
        output.lmus = projectArea.lmuRaster;
        output.lmuIds = projectArea.lmuIds;
        output.commandLine = ps->commandLine;

        std::cout << "\tCreating allometry from ";
        if(ps->slope > 0) {
            std::cout << "user provided coefficients...";
            dbhModel = rxtools::allometry::UnivariateLinearModel(ps->slope, ps->intercept, ps->transform, ps->rsq, ps->inUnit, ps->outUnit);
            std::cout << "  Done!\n";
        }
        else {
            std::cout << "FIA plots...";
            auto reader = rxtools::allometry::FIAReader(ps->fiaPath);
            auto dist = projectArea.lidarDataset->units().value().convertOneToThis(10000, lapis::linearUnitPresets::meter);
            lapis::Extent e(projectArea.projectPoly.extent().xmin() - dist, projectArea.projectPoly.extent().xmax() + dist,
                projectArea.projectPoly.extent().ymin() - dist, projectArea.projectPoly.extent().ymax() + dist,
                projectArea.projectPoly.crs());
            if (!reader.limitByExtent(e))
                throw std::runtime_error("no fia in this place");
            reader.makePlotTreeMap(std::vector<std::string>{ "DIA" });
            auto allTrees = reader.collapsePlotTreeMap();
            allTrees.writeCsv(ProjectSettings::get().outputPath + "/fia_tree_data.csv");
            dbhModel = rxtools::allometry::UnivariateLinearModel(allTrees, "DIA", rxtools::linearUnitPresets::inch);
            std::cout << " Done!\n";
            std::cout << "intercept: " << dbhModel.parameters.intercept << "\n";
            std::cout << "slope: " << dbhModel.parameters.slope << "\n";
            std::cout << "transform: " << int(dbhModel.parameters.transform) << "\n";
            std::cout << "rsq: " << dbhModel.parameters.rsq << "\n";
            std::cout << "inUnit: " << dbhModel.inputUnit.name() << "\n";
            std::cout << "outUnit: " << dbhModel.outputUnit.name() << "\n";

            if (ps->writeUnits && ps->fastFuels) {
                ffa = rxtools::allometry::FastFuels(allTrees);
            }
        }

        std::cout << "\tReading and reprojecting unit layer...";
        if (std::filesystem::path(ps->unitPolygonPath).extension() != ".shp")
            shp = false;

        if (shp) {
            unitPoly = lapis::VectorDataset<lapis::MultiPolygon>(ps->unitPolygonPath);
            if (!unitPoly.crs().isConsistent(projectArea.projectPoly.crs()))
                unitPoly.projectInPlacePreciseExtent(projectArea.projectPoly.crs());
            if (!unitPoly.extent().overlaps(projectArea.projectPoly.extent())) {
                std::cerr << "Unit polygon does not overlap project polygon\n";
                throw lapis::OutsideExtentException();
            }

        }
        else {
            unitRaster = lapis::Raster<int>(ps->unitPolygonPath);
            lapis::Alignment projectAlign = lapis::Alignment((lapis::Extent)ps->unitPolygonPath, unitRaster.nrow(), unitRaster.ncol());
            unitRaster = lapis::resampleRaster(unitRaster, projectAlign, lapis::ExtractMethod::near);
            unitRaster = lapis::trimRaster(unitRaster);
            if (!unitRaster.overlaps(projectArea.projectPoly.extent())) {
                std::cerr << "Unit raster does not overlap project polygon\n";
                throw lapis::OutsideExtentException();
            }
        }
        std::cout << " Done!\n";

        rxtools::TaoGettersMP getters = rxtools::TaoGettersMP(
            lapis::lico::alwaysAdd<lapis::VectorDataset<lapis::MultiPolygon>>,
            projectArea.lidarDataset->coordGetter(),
            projectArea.lidarDataset->heightGetter(),
            projectArea.lidarDataset->radiusGetter(),
            projectArea.lidarDataset->areaGetter(),
            //TODO: at some point it would be good to not implicitly assume height units are not in meters.
            [dm = dbhModel, hg = projectArea.lidarDataset->heightGetter()](const lapis::ConstFeature<lapis::MultiPolygon>& ft)->double {
                return dm.predict(hg(ft), lapis::linearUnitPresets::meter, rxtools::linearUnitPresets::centimeter);
            }
        );

        projectArea.createCoreGapAndReadTaos(ps->nThread, ps->dbhMax, getters);
        std::cout << "Alltaos size: " << projectArea.allTaos.size() << "\n";

        std::cout << "\t Performing PCA on climate data...";
        std::filebuf fb;
        if (ps->referenceDatasetPath != "") {
            if (!fb.open(ps->referenceDatasetPath, std::ios::in)) throw std::runtime_error("Cannot open reference table.");
        }
        else {
            if (!fb.open(ps->defaultRefPath, std::ios::in)) throw std::runtime_error("Cannot open internal reference table.");
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
        else
            return neFacingLmu.detailedQuery(pcClimate, k, false, maxDist);
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

        auto thisOsiNum = lapis::cropRaster(projectArea.bbOsiNum, lmu.mask, lapis::SnapType::out);
        auto thisOsiDen = lapis::cropRaster(projectArea.bbOsiDen, lmu.mask, lapis::SnapType::out);
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

        rxtools::TaoListMP treatedTaos{ projectArea.allTaos, true };
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

        std::cout << "Calculating post treatment OSI\n";
        //Post osi
        auto postNum = lapis::Raster<int>{ (lapis::Alignment)projectArea.lmuIds };
        auto postDen = lapis::Raster<int>{ (lapis::Alignment)projectArea.lmuIds };
        std::pair<lapis::coord_t, lapis::coord_t> expectedRes{};
        if (projectArea.lidarDataset->type() == processedfolder::RunType::fusion) {
            if (projectArea.lidarDataset->units() == lapis::linearUnitPresets::meter) {
                expectedRes.first = 0.75;
            }
            else {
                expectedRes.first = 2.4606;
            }
        }
        else {
            expectedRes.first = -1;
            for (size_t i = 0; i < projectArea.lidarDataset->nTiles(); ++i) {
                auto csmFile = projectArea.lidarDataset->csmRaster(i);
                if (csmFile) {
                    expectedRes.first = lapis::Raster<lapis::coord_t>(csmFile.value().string()).xres();
                    break;
                }
            }
            if (expectedRes.first < 0) {
                throw processedfolder::FileNotFoundException("no tiles with data found.");
            }
        }
        expectedRes.second = 0;

        threads.clear();
        std::queue<rxtools::ProjectArea::CoreGapWorkItem> workQ;
        std::condition_variable cvWork;   // I/O → Workers
        std::condition_variable cvQ;  // Workers → I/O
        bool ioComplete = false;
        double coreGapDist = 6;

        std::thread ioThread([&]() {
            projectArea.coreGapIOThread(workQ, mut, cvWork, cvQ, ioComplete, 
                projectArea.lmuRaster, coreGapDist, expectedRes, 
                treatedTaos.getters, true, nThread);
            });

        auto postGapFunc = [&](int i) { 
            projectArea.postGapWorkThread(
                postNum, postDen, treatedTaos,
                i, mut, workQ, cvWork, cvQ, ioComplete,
                projectArea.lmuRaster, 2, coreGapDist
            ); };

        for (int i = 0; i < nThread; ++i) {
            threads.push_back(std::thread(postGapFunc, i));
        }

        ioThread.join();
        for (int i = 0; i < nThread; ++i) {
            threads[i].join();
        }

        std::cout << "Calculating post OSI zones...\n";
        auto numZones = lapis::zonalSum(postNum, unitZonal);
        auto denZones = lapis::zonalSum(postDen, unitZonal);

        lapis::Raster<double> postOsi{ (lapis::Alignment)unitZonal };
        for (auto z : numZones) {
            double osi = static_cast<double>(z.second) / static_cast<double>(denZones.find(z.first)->second) * 100;
            for (lapis::cell_t c = 0; c < postOsi.ncell(); ++c) {
                if (unitZonal[c].value() == z.first) {
                    postOsi[c].value() = osi;
                    postOsi[c].has_value() = true;
                }
            }
        }
        postOsi = lapis::cropRaster(postOsi, output.pre[3], lapis::SnapType::out);
        postOsi = lapis::extendRaster(postOsi, output.pre[3], lapis::SnapType::out);
        postOsi.mask(output.pre[3]);
        output.post[3] = postOsi;

        for (size_t i = 0; i < output.atts.nFeature(); ++i) {
            auto v = output.atts.getIntegerField(i, "ID");
            if (numZones.find(v) == numZones.end()) {
                std::cout << "failed to find: " << v << "\n";
                continue;
            }
            double osi = static_cast<double>(numZones.at(v)) / static_cast<double>(denZones.at(v)) * 100;
            output.atts.setRealField(i, "trtOSI", osi);
        }
        std::cout << "Treatment done\n";
    }

    void Licosim::treatmentThread(size_t& sofar, std::mutex& mut, double dbhMin, double dbhMax, lapis::Raster<lapis::cell_t>& unitZonal, rxtools::TaoListMP& treatedTaos, const int thisThread) {
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

                rxtools::TaoListMP taos{ projectArea.allTaos, true };
                for (size_t j = 0; j < projectArea.allTaos.size(); ++j) {
                    if (lmu.mask.extract(projectArea.allTaos.x(j), projectArea.allTaos.y(j), lapis::ExtractMethod::near).has_value())
                        taos.taoVector.addFeature(projectArea.allTaos.taoVector.getFeature(j));
                }
                std::cout << "LMU taos size: " << taos.size() << "\n";

                auto after = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
                std::cout << "Thread " + std::to_string(thisThread) + " completed in " + std::to_string(duration.count()) + " seconds.\n";

                std::cout << "\tCreating Units for Lmu " + std::to_string(i) + "/" + std::to_string(nLmu) + " on thread " + std::to_string(thisThread) + "\n";
                before = std::chrono::high_resolution_clock::now();
                auto thisOsiNum = lapis::cropRaster(projectArea.osiNum, lmu.mask, lapis::SnapType::out);
                auto thisOsiDen = lapis::cropRaster(projectArea.osiDen, lmu.mask, lapis::SnapType::out);

                if (shp)
                    lmu.makeUnits(unitPoly, taos, thisOsiNum, thisOsiDen, ProjectSettings::get().overrideTargets);
                else
                    lmu.makeUnits(unitRaster, taos, thisOsiNum, thisOsiDen);

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

                        lmu.units[j].treatedStructure = rxtools::StructureSummary(lmu.units[j].treatedTaos, lmu.units[j].unitMask, lmu.units[j].areaHa, 0);
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
                    for (size_t k = 0; k < lmu.units[j].treatedTaos.size(); ++k)
                        treatedTaos.taoVector.addFeature(lmu.units[j].treatedTaos.taoVector.getFeature(k));

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