#include "roofitCombined.h"

// Simultaneous fit of the PbPb years of one tree (ntmix_X3872, ntmix_PSI2S, ntKp, ntphi, ntKstar), bin by bin in
// VAR. Every parameter is common to all years except the yields: the MC signal shape (from the MC of all years, then
// locked), the peak position and the data/MC resolution scale, the background shapes (combinatorial and, for ntKp,
// the partially reconstructed B). The single-year fit (../roofitB.C) stays the fit of ppRef and the check of each year.
//
// Inputs, per year: the scored DATA and MC samples of selectionER (scored_samples/<tree>/<year>/) and the score
// threshold. The selection of a year is the pre-cut and score branch of its scored DATA metadata plus the threshold.
//
// Products, in the folders of roofitB.C (doCombinedRoofit.sh, README.md). Every year gets the files of a single-year
// fit, under the same names, replacing the ones there:
//   fitER/ROOTfiles/<YEAR>/          nominalFitModel_<TREE>_<YEAR>.root (FULL=1, validation), fitResults_ and
//                                    mcFitResults_<TREE>_<VAR>_<YEAR>.root (FULL=0, Acc x Eff), systematicFILES/
//                                    (SYST=1, resultER)
//   fitER/results/<YEAR>/            fit plots (<TREE>/<VAR>/), tables/, Graphs/ of the year
//   fitER/ROOTfiles/<COMBO>/         the simultaneous workspace, the yields of every year and their sum
//   fitER/results/<COMBO>/           the add-up plots (<TREE>/<VAR>/: data and fits of all years summed), Graphs/
// COMBO is the years joined by "_" (PbPb23_PbPb24).
void roofitCombined(TString TREE = "ntKp", int FULL = 1,
                    std::vector<TString> YEARS = {}, std::vector<TString> INPUTDATA = {}, std::vector<TString> INPUTMC = {},
                    std::vector<double> THRESHOLDS = {}, TString VAR = "Bpt", std::vector<double> VAR_BINS = {},
                    int SYST = 0)
{
	const int nYears = YEARS.size();
	const int nBins = VAR_BINS.size() - 1;
	TString combo = YEARS[0], comboLabel = YEARS[0] + "'", yearList = YEARS[0];
	for (int y = 1; y < nYears; ++y) {
		combo += "_" + YEARS[y];
		comboLabel += " + " + YEARS[y] + "'";
		yearList += "," + YEARS[y];
	}

	cout << endl << "Combined fit of " << TREE << " in " << VAR << " (" << nBins << " bins), years " << yearList
	     << ", every parameter but the yields common to all years" << endl;

	// The folders of roofitB.C, from fitER/combinedANA
	const TString comboRoot = "../ROOTfiles/" + combo;
	const TString comboPlots = Form("../results/%s/%s/%s", combo.Data(), TREE.Data(), VAR.Data());
	const TString comboGraphs = "../results/" + combo + "/Graphs";
	gSystem->mkdir(comboRoot, true);
	gSystem->mkdir(comboPlots, true);
	gSystem->mkdir(comboGraphs, true);
	for (const TString& year : YEARS) {
		gSystem->mkdir(Form("../ROOTfiles/%s/systematicFILES", year.Data()), true);
		gSystem->mkdir(Form("../results/%s/%s/%s", year.Data(), TREE.Data(), VAR.Data()), true);
		gSystem->mkdir(Form("../results/%s/tables", year.Data()), true);
		gSystem->mkdir(Form("../results/%s/Graphs", year.Data()), true);
	}

	// VARIABLES (ranges of roofitB.C)
	const TreeSetup& setup = GetTreeSetup(TREE);
	const double massMin = setup.massMin, massMax = setup.massMax;
	const int nMassBins = setup.nMassBins;
	RooRealVar* mass = new RooRealVar("Bmass", "Bmass", massMin, massMax);
	mass->setRange("all", massMin, massMax);
	mass->setRange("signal", setup.peak - setup.signalHalfWidth, setup.peak + setup.signalHalfWidth);
	RooRealVar* pt = new RooRealVar("Bpt", "Bpt", 0, 300);
	RooRealVar* rapidity = new RooRealVar("By", "By", -2.4, 2.4);
	RooRealVar* nChargedTracks = new RooRealVar("nChargedTracks", "nChargedTracks", 0, 2000000000);
	RooRealVar* centBin = new RooRealVar("CentBin", "CentBin", 0, 100);
	RooRealVar* pThatWeight = new RooRealVar("pThatreweight", "pThatreweight", 0, 1.0e12);
	RooArgSet dataVars(*mass, *pt, *rapidity, *nChargedTracks, *centBin);
	RooArgSet mcVars(dataVars);
	mcVars.add(*pThatWeight);
	RooCategory* cat = new RooCategory("year", "year");
	for (const TString& year : YEARS) cat->defineType(year.Data());

	// DATA and MC SAMPLES of every year, with the selection of the year's scored samples
	std::vector<TString> cuts(nYears);
	std::vector<RooDataSet*> data(nYears), mc(nYears);
	for (int y = 0; y < nYears; ++y) {
		TFile* scored = TFile::Open(INPUTDATA[y]);
		cuts[y] = Form("(%s) && (%s > %g)", scored->Get<TObjString>("metadata/pre_cut")->GetString().Data(),
		               scored->Get<TObjString>("metadata/score_branch")->GetString().Data(), THRESHOLDS[y]);
		scored->Close();
		data[y] = ReadSample(INPUTDATA[y], setup.dataTree, cuts[y], dataVars, "data_" + YEARS[y], massMin, massMax);
		mc[y] = ReadSample(INPUTMC[y], TREE, cuts[y], mcVars, "mc_" + YEARS[y], massMin, massMax);
	}

	// Per-year export workspaces (the contract of roofitB.C) and the combined workspace
	std::vector<RooWorkspace*> wsYear(nYears), wsYearMC(nYears);
	for (int y = 0; y < nYears; ++y) {
		wsYear[y] = new RooWorkspace("ws");
		wsYearMC[y] = new RooWorkspace("ws_mc");
		wsYear[y]->import(*data[y], Rename("data"));
	}
	RooWorkspace* wsComb = new RooWorkspace("ws_combined");

	// FIT INFO, [year][bin]
	std::vector<std::vector<double>> yield(nYears, std::vector<double>(nBins)), yieldErr = yield, chi2 = yield;
	std::vector<std::vector<double>> scaleVal = yield, scaleErr = yield, resol = yield, resolErr = yield, signif = yield;
	std::vector<double> sumYield(nBins), sumYieldErr(nBins), sumChi2(nBins), combSignif(nBins, -1.0);
	std::vector<double> xMean(nBins), xLow(nBins), xHigh(nBins);
	std::vector<TH1D*> hPt(nYears), hPtMC(nYears);
	for (int y = 0; y < nYears; ++y) {
		hPt[y] = new TH1D(Form("hPt_%s", YEARS[y].Data()), "", nBins, VAR_BINS.data());
		hPtMC[y] = new TH1D(Form("hPtMC_%s", YEARS[y].Data()), "", nBins, VAR_BINS.data());
	}
	TH1D* hPtSum = new TH1D("hPtSum", "", nBins, VAR_BINS.data());

	// MODELS for syst studies; relative yield deviations [year][bin][variation], summary [year][bin] = {bkg, sig, total}
	const std::vector<SystVariationConfig> background = GetBackgroundSystematicModels(TREE);
	const std::vector<SystVariationConfig> signal = GetSignalSystematicModels(TREE);
	std::vector<std::vector<std::vector<double>>> backRel(nYears), sigRel(nYears), general(nYears), stat(nYears);

	//BIN ANALYSIS START
	for (int i = 0; i < nBins; ++i) {
		const int bin = i + 1;
		const double low = VAR_BINS[i], high = VAR_BINS[i + 1], width = high - low;
		const TString binTag = VAR == "By" ? Form("abs(By)_%0.1f_%0.1f", low, high) : Form("%s_%i_%i", VAR.Data(), (int)low, (int)high);
		const TString binCut = Form("(abs(%s)>=%f && abs(%s)<%f)", VAR.Data(), low, VAR.Data(), high);

		std::vector<RooDataSet*> dataBin(nYears), mcBin(nYears);
		std::map<std::string, RooDataSet*> dataMap, mcMap;
		for (int y = 0; y < nYears; ++y) {
			dataBin[y] = static_cast<RooDataSet*>(data[y]->reduce(binCut));
			dataBin[y]->SetName(Form("data%d_%s", bin, YEARS[y].Data()));
			mcBin[y] = static_cast<RooDataSet*>(mc[y]->reduce(binCut));
			mcBin[y]->SetName(Form("mc%d", bin));
			dataMap[YEARS[y].Data()] = dataBin[y];
			mcMap[YEARS[y].Data()] = mcBin[y];
			cout << YEARS[y] << " bin " << bin << ": " << dataBin[y]->numEntries() << " DATA, " << mcBin[y]->numEntries() << " MC entries" << endl;
		}
		RooDataSet* combData = new RooDataSet(Form("combData%d", bin), "", dataVars, Index(*cat), Import(dataMap));
		RooDataSet* combMC = new RooDataSet(Form("combMC%d", bin), "", mcVars, Index(*cat), Import(mcMap), WeightVar(*pThatWeight));

		double sumVar = 0.0;
		for (int e = 0; e < combData->numEntries(); ++e) sumVar += std::abs(combData->get(e)->getRealValue(VAR));
		xMean[i] = sumVar / combData->numEntries();
		xLow[i] = xMean[i] - low;
		xHigh[i] = high - xMean[i];

		////////// FITFITFITFITFITFITFITFITFITFITFITFIT //////////
		cout << "Starting the combined fit for " << TREE << " " << VAR << " [" << low << ", " << high << "]" << endl;
		CombinedFit nominal = FitCombined(TREE, bin, "", "", YEARS, mass, cat, dataBin, mcBin, combData, combMC);
		const RooArgList& floating = nominal.result->floatParsFinal();
		////////// FITFITFITFITFITFITFITFITFITFITFITFIT //////////

		// Per year: export (the products of roofitB.C), plots, fit info
		for (int y = 0; y < nYears; ++y) {
			const YearModel& m = nominal.years[y];
			YearModel exportMC = ExportYearModel(TREE, bin, "", "", m, *nominal.mcValues, mass, dataBin[y], mcBin[y]);
			TCanvas* cMC = new TCanvas(Form("cMC%d_%s", bin, YEARS[y].Data()), "", 700, 700);
			DrawMCFit(cMC, TREE, mass, mcBin[y], exportMC, nMassBins);
			cMC->SaveAs(Form("../results/%s/%s/%s/mc_%s_%s_%s.pdf", YEARS[y].Data(), TREE.Data(), VAR.Data(), YEARS[y].Data(),
			                 binTag.Data(), TREE.Data()));
			wsYearMC[y]->import(*mcBin[y]);
			wsYearMC[y]->import(*exportMC.modelMC, RecycleConflictNodes());
			std::unique_ptr<RooArgSet> mcPars(exportMC.modelMC->getParameters(*mcBin[y]));
			wsYearMC[y]->saveSnapshot(Form("mcFitPars_bin%d", bin), *mcPars, true);

			YearModel exported = ExportYearModel(TREE, bin, "", "", m, *nominal.dataValues, mass, dataBin[y], mcBin[y]);
			RooArgSet floatingYear;
			for (size_t k = 0; k < m.pars.size(); ++k) {
				if (floating.find(m.pars[k]->GetName())) floatingYear.add(*exported.pars[k]);
			}
			wsYear[y]->import(*exported.model, RecycleConflictNodes());
			wsYear[y]->import(*exported.nsigMC);
			wsYear[y]->saveSnapshot(Form("nominalPars_bin%d", bin), floatingYear, true);
			wsYear[y]->loadSnapshot(Form("nominalPars_bin%d", bin));

			std::unique_ptr<RooArgSet> yearPars(exported.model->getParameters(*dataBin[y]));
			int nFloat = 0;
			for (RooAbsArg* p : *yearPars) nFloat += !static_cast<RooRealVar*>(p)->isConstant();
			TCanvas* c = new TCanvas(Form("c%d_%s", bin, YEARS[y].Data()), "", 700, 700);
			RooArgSet erfc;
			if (exported.erfc) erfc.add(*exported.erfc);
			chi2[y][i] = DrawDataFit(c, TREE, mass, dataBin[y], exported.model, RooArgSet(*exported.sig), RooArgSet(*exported.bkg), erfc,
			                         nFloat, nMassBins);
			RooRealVar chi2Var(Form("chi2_data_norm%d_", bin), "", chi2[y][i]);
			wsYear[y]->import(chi2Var);

			yield[y][i] = m.nsig->getVal();
			yieldErr[y][i] = m.nsig->getError();
			scaleVal[y][i] = m.scale->getVal();
			scaleErr[y][i] = m.scale->getError();
			resol[y][i] = std::sqrt(m.sig1frac->getVal() * std::pow(m.sigma1->getVal(), 2) +
			                        (1 - m.sig1frac->getVal()) * std::pow(m.sigma2->getVal(), 2)) * scaleVal[y][i];
			resolErr[y][i] = scaleErr[y][i] / scaleVal[y][i] * resol[y][i];
			hPt[y]->SetBinContent(bin, yield[y][i] / width);
			hPt[y]->SetBinError(bin, yieldErr[y][i] / width);
			hPtMC[y]->SetBinContent(bin, m.nsigMC->getVal() / width);
			hPtMC[y]->SetBinError(bin, m.nsigMC->getError() / width);

			TString significance = "";
			if (setup.significance) {
				signif[y][i] = GetSignificance(wsYear[y], bin, mass, 2.0);
				significance = Form("S = %.2f", signif[y][i]);
			}
			DrawFitLabels(TREE, VAR, low, high, chi2[y][i], yield[y][i], yieldErr[y][i], "", significance);
			DrawCmsHeader(c, YEARS[y]);
			c->SaveAs(Form("../results/%s/%s/%s/data_%s_%s_%s.pdf", YEARS[y].Data(), TREE.Data(), VAR.Data(), YEARS[y].Data(),
			               binTag.Data(), TREE.Data()));
		}

		// Combined workspace
		wsComb->import(*nominal.sim, RecycleConflictNodes());
		wsComb->import(*combData);
		wsComb->saveSnapshot(Form("nominalPars_bin%d", bin), floating, true);

		// Sum of the yields, its error from the covariance of the combined fit
		const TMatrixDSym cov = nominal.result->covarianceMatrix();
		for (int a = 0; a < nYears; ++a) {
			sumYield[i] += yield[a][i];
			for (int b = 0; b < nYears; ++b) {
				sumYieldErr[i] += cov(floating.index(nominal.years[a].nsig->GetName()), floating.index(nominal.years[b].nsig->GetName()));
			}
		}
		sumYieldErr[i] = std::sqrt(sumYieldErr[i]);
		hPtSum->SetBinContent(bin, sumYield[i] / width);
		hPtSum->SetBinError(bin, sumYieldErr[i] / width);

		// X(3872), B0: significance of the signal of all years together, likelihood ratio against no signal in any year
		// (k = number of years bounded signal yields: chi-bar-squared, sum_j C(k,j)/2^k Prob(q0, j)). Past q0 ~ 1400 the
		// p-value is below the smallest double: Z is then the asymptotic sqrt(q0) of one bounded yield.
		TString combSignificance = "";
		if (setup.significance) {
			for (YearModel& m : nominal.years) {
				m.nsig->setVal(0.0);
				m.nsig->setConstant(true);
			}
			RooFitResult* noSignal = nominal.sim->fitTo(*combData, Save(), Extended(true), Range("all"), PrintLevel(-1));
			const double q0 = std::max(0.0, 2.0 * (noSignal->minNll() - nominal.result->minNll()));
			double pValue = 0.0;
			for (int j = 1; j <= nYears; ++j) pValue += TMath::Binomial(nYears, j) / std::pow(2.0, nYears) * TMath::Prob(q0, j);
			combSignif[i] = pValue > 0.0 ? ROOT::Math::normal_quantile_c(pValue, 1.0) : std::sqrt(q0);
			combSignificance = Form("Z_{LR} = %.2f", combSignif[i]);
			cout << "Combined significance: q0 = " << q0 << ", p = " << pValue << ", Z = " << combSignif[i] << endl;
			std::unique_ptr<RooArgSet> simPars(nominal.sim->getParameters(*combData));
			simPars->assignValueOnly(*nominal.dataValues);
			for (YearModel& m : nominal.years) m.nsig->setConstant(false);
		}

		// Add-up plot: the data of all years summed, the combined fit, the signal of every year, the background of all years
		std::vector<RooAbsPdf*> yearModels, yearSignals;
		RooArgSet backgrounds;
		for (const YearModel& m : nominal.years) {
			yearModels.push_back(m.model);
			yearSignals.push_back(m.sig);
			backgrounds.add(*m.bkg);
			if (m.erfc) backgrounds.add(*m.erfc);
		}
		TCanvas* cSum = new TCanvas(Form("cSum%d", bin), "", 700, 700);
		sumChi2[i] = DrawAddUpFit(cSum, TREE, mass, combData, yearModels, yearSignals, backgrounds, YEARS, floating.getSize(), nMassBins);
		DrawFitLabels(TREE, VAR, low, high, sumChi2[i], sumYield[i], sumYieldErr[i], "", combSignificance);
		DrawCmsHeader(cSum, "");
		TLatex header;
		header.SetNDC();
		header.SetTextFont(42);
		header.SetTextSize(0.030);
		header.SetTextAlign(31);
		header.DrawLatex(1.06 - cSum->GetRightMargin(), 1.0 - 0.7 * cSum->GetTopMargin(), comboLabel + " #sqrt{s_{NN}}=5.36 TeV");
		cSum->SaveAs(Form("%s/data_%s_%s_%s.pdf", comboPlots.Data(), combo.Data(), binTag.Data(), TREE.Data()));

		// Fit-model variations: the whole combined fit again with the variation in every year
		if (SYST == 1 && FULL == 0) {
			for (int y = 0; y < nYears; ++y) {
				backRel[y].push_back({});
				sigRel[y].push_back({});
				stat[y].push_back({yieldErr[y][i] / yield[y][i] * 100.0});
			}
			for (int isSignal = 0; isSignal < 2; ++isSignal) {
				for (const SystVariationConfig& v : isSignal ? signal : background) {
					const TString sigCode = isSignal ? v.code.c_str() : "";
					const TString bkgCode = isSignal ? "" : v.code.c_str();
					CombinedFit f = FitCombined(TREE, bin, sigCode, bkgCode, YEARS, mass, cat, dataBin, mcBin, combData, combMC);
					for (int y = 0; y < nYears; ++y) {
						const YearModel& m = f.years[y];
						if (isSignal && v.code != "fixed") {
							YearModel exportMC = ExportYearModel(TREE, bin, sigCode, bkgCode, m, *f.mcValues, mass, dataBin[y], mcBin[y]);
							TCanvas* cMC = new TCanvas(Form("cMC%d_%s_%s", bin, v.code.c_str(), YEARS[y].Data()), "", 700, 700);
							DrawMCFit(cMC, TREE, mass, mcBin[y], exportMC, nMassBins);
							cMC->SaveAs(Form("../results/%s/%s/%s/mc_%s_%s_%s_%s.pdf", YEARS[y].Data(), TREE.Data(), VAR.Data(),
							                 YEARS[y].Data(), binTag.Data(), v.code.c_str(), TREE.Data()));
						}
						YearModel exported = ExportYearModel(TREE, bin, sigCode, bkgCode, m, *f.dataValues, mass, dataBin[y], mcBin[y]);
						std::unique_ptr<RooArgSet> yearPars(exported.model->getParameters(*dataBin[y]));
						int nFloat = 0;
						for (RooAbsArg* p : *yearPars) nFloat += !static_cast<RooRealVar*>(p)->isConstant();
						TCanvas* c = new TCanvas(Form("c%d_%s_%s", bin, v.code.c_str(), YEARS[y].Data()), "", 700, 700);
						RooArgSet erfc;
						if (exported.erfc) erfc.add(*exported.erfc);
						const double chi2Variation = DrawDataFit(c, TREE, mass, dataBin[y], exported.model, RooArgSet(*exported.sig),
						                                         RooArgSet(*exported.bkg), erfc, nFloat, nMassBins);
						DrawFitLabels(TREE, VAR, low, high, chi2Variation, m.nsig->getVal(), m.nsig->getError(), v.label, "");
						DrawCmsHeader(c, YEARS[y]);
						c->SaveAs(Form("../results/%s/%s/%s/data_%s_%s_%s_%s.pdf", YEARS[y].Data(), TREE.Data(), VAR.Data(), YEARS[y].Data(),
						               binTag.Data(), v.code.c_str(), TREE.Data()));
						const double deviation = std::abs((yield[y][i] - m.nsig->getVal()) / yield[y][i]) * 100.0;
						(isSignal ? sigRel : backRel)[y][i].push_back(deviation);
					}
				}
			}
			for (int y = 0; y < nYears; ++y) {
				const double maxBack = *std::max_element(backRel[y][i].begin(), backRel[y][i].end());
				const double maxSig = *std::max_element(sigRel[y][i].begin(), sigRel[y][i].end());
				general[y].push_back({maxBack, maxSig, std::hypot(maxBack, maxSig)});
			}
		}
	}
	//BIN ANALYSIS END

	// PRODUCTS
	// FULL=1 writes the inclusive model used by validation, FULL=0 the analysis-bin models used by Acc x Eff.
	const TString combinedFile = FULL == 1 ? Form("%s/combinedFitModel_%s_%s.root", comboRoot.Data(), TREE.Data(), combo.Data())
	                                       : Form("%s/combinedFitResults_%s_%s_%s.root", comboRoot.Data(), TREE.Data(), VAR.Data(), combo.Data());
	auto writeYearMetadata = [&](TFile* f, int y, const TString& fitMode, const TString& fitRole) {
		f->cd();
		WriteFitMetadata(f, TREE, VAR, YEARS[y], cuts[y], VAR_BINS, massMin, massMax, nMassBins);
		TObjString(INPUTMC[y]).Write("inputMC");
		TObjString(INPUTDATA[y]).Write("inputData");
		TObjString("pThatreweight").Write("mcWeight");
		TObjString(fitMode).Write("fitMode");
		TObjString(fitRole).Write("fitRole");
		TObjString("combined").Write("fitStrategy");
		TObjString(yearList).Write("combinedYears");
		TObjString("all but the yields").Write("sharedParameters");
		TObjString(combinedFile).Write("combinedFile");
		TParameter<double>("scoreThreshold", THRESHOLDS[y]).Write("scoreThreshold");
	};
	for (int y = 0; y < nYears; ++y) {
		const TString yearRoot = "../ROOTfiles/" + YEARS[y];
		TFile* out = new TFile(FULL == 1 ? Form("%s/nominalFitModel_%s_%s.root", yearRoot.Data(), TREE.Data(), YEARS[y].Data())
		                                 : Form("%s/fitResults_%s_%s_%s.root", yearRoot.Data(), TREE.Data(), VAR.Data(), YEARS[y].Data()),
		                       "RECREATE");
		wsYear[y]->Write("ws_nominal");
		hPt[y]->Write("hPt");
		writeYearMetadata(out, y, FULL == 1 ? "inclusive" : "binned", FULL == 1 ? "validation sPlot" : "Acc x Eff data correction");
		out->Close();
		if (FULL == 0) {
			TFile* outMC = new TFile(Form("%s/mcFitResults_%s_%s_%s.root", yearRoot.Data(), TREE.Data(), VAR.Data(), YEARS[y].Data()), "RECREATE");
			wsYearMC[y]->Write("ws_mc");
			hPtMC[y]->Write("hPtMC");
			writeYearMetadata(outMC, y, "binnedMC", "Acc x Eff closure");
			outMC->Close();
		}
	}

	TH1D* hSignif = new TH1D("hCombinedSignificance", ";" + FitVarLabel(VAR) + ";Z_{LR}", nBins, VAR_BINS.data());
	for (int i = 0; i < nBins; ++i) hSignif->SetBinContent(i + 1, combSignif[i]);
	TFile* outComb = new TFile(combinedFile, "RECREATE");
	wsComb->Write("ws_combined");
	hSignif->Write();
	for (int y = 0; y < nYears; ++y) hPt[y]->Write();
	hPtSum->Write();
	WriteFitMetadata(outComb, TREE, VAR, combo, "per year: selectionCut_<year>", VAR_BINS, massMin, massMax, nMassBins);
	TObjString(yearList).Write("combinedYears");
	TObjString("all but the yields").Write("sharedParameters");
	TObjString(FULL == 1 ? "inclusive" : "binned").Write("fitMode");
	for (int y = 0; y < nYears; ++y) {
		TObjString(cuts[y]).Write("selectionCut_" + YEARS[y]);
		TObjString(INPUTDATA[y]).Write("inputData_" + YEARS[y]);
		TObjString(INPUTMC[y]).Write("inputMC_" + YEARS[y]);
	}
	outComb->Close();

	if (FULL == 1) return;

	// Summary graphs: every year and their sum
	std::vector<TString> labels = YEARS;
	std::vector<std::vector<double>> zero(nYears + 1, std::vector<double>(nBins, 0.0));
	std::vector<std::vector<double>> yields = yield, yieldErrs = yieldErr, chi2s = chi2;
	labels.push_back("Sum");
	yields.push_back(sumYield);
	yieldErrs.push_back(sumYieldErr);
	chi2s.push_back(sumChi2);
	DrawBinGraphs(Form("%s/raw_yield_%s_%s_%s.pdf", comboGraphs.Data(), TREE.Data(), combo.Data(), VAR.Data()), FitVarLabel(VAR), "Y_{S}",
	              xMean, xLow, xHigh, labels, yields, yieldErrs);
	DrawBinGraphs(Form("%s/chi2_%s_%s_%s.pdf", comboGraphs.Data(), TREE.Data(), combo.Data(), VAR.Data()), FitVarLabel(VAR), "#chi^{2}/NDF",
	              xMean, xLow, xHigh, labels, chi2s, zero);
	DrawBinGraphs(Form("%s/resolution_%s_%s_%s.pdf", comboGraphs.Data(), TREE.Data(), combo.Data(), VAR.Data()), FitVarLabel(VAR), "Resolution [GeV/c^{2}]",
	              xMean, xLow, xHigh, YEARS, resol, resolErr);
	DrawBinGraphs(Form("%s/scale_variation_%s_%s_%s.pdf", comboGraphs.Data(), TREE.Data(), combo.Data(), VAR.Data()), FitVarLabel(VAR), "Scale factor",
	              xMean, xLow, xHigh, YEARS, scaleVal, scaleErr);

	if (SYST == 0) return;

	// Systematic tables, leading uncertainty (resultER) and graphs of every year
	std::vector<std::string> labelsBack, labelsSignal;
	std::vector<std::string> colBack = {"Background Model"}, colSignal = {"Signal Model"}, colGeneral = {"Systematic Source"};
	for (const auto& v : background) labelsBack.push_back(v.label);
	for (const auto& v : signal) labelsSignal.push_back(v.label);
	for (int i = 0; i < nBins; ++i) {
		const std::string column = GetSystematicColumnLabel(VAR, VAR_BINS[i], VAR_BINS[i + 1]);
		colBack.push_back(column);
		colSignal.push_back(column);
		colGeneral.push_back(column);
	}
	for (int y = 0; y < nYears; ++y) {
		WriteSystematicsTablesDocument(Form("../results/%s/tables/systematics_tables_%s_%s_%s", YEARS[y].Data(), TREE.Data(), YEARS[y].Data(), VAR.Data()),
		                               colSignal, colBack, colGeneral, labelsSignal, labelsBack, GetGeneralSystematicLabels(),
		                               sigRel[y], backRel[y], BuildGeneralSystematicNumbers(general[y], stat[y]));

		TH1D hLeading("hLeadingTotalUncPercent", ";" + FitVarLabel(VAR) + ";Leading total uncertainty (%)", nBins, VAR_BINS.data());
		hLeading.SetStats(0);
		for (int i = 0; i < nBins; ++i) hLeading.SetBinContent(i + 1, general[y][i][2]);
		TFile fLeading(Form("../ROOTfiles/%s/systematicFILES/fit_totalSystematic_leadingUnc_%s_%s_%s.root", YEARS[y].Data(),
		                    TREE.Data(), YEARS[y].Data(), VAR.Data()), "RECREATE");
		hLeading.Write("hLeadingTotalUncPercent");
		fLeading.Close();

		std::vector<TString> systLabels;
		std::vector<std::vector<double>> systValues;
		for (size_t j = 0; j < background.size(); ++j) {
			systLabels.push_back(background[j].label);
			systValues.push_back({});
			for (int i = 0; i < nBins; ++i) systValues.back().push_back(backRel[y][i][j]);
		}
		for (size_t j = 0; j < signal.size(); ++j) {
			systLabels.push_back(signal[j].label);
			systValues.push_back({});
			for (int i = 0; i < nBins; ++i) systValues.back().push_back(sigRel[y][i][j]);
		}
		systLabels.insert(systLabels.end(), {"Total Systematic", "Statistical"});
		systValues.push_back({});
		systValues.push_back({});
		for (int i = 0; i < nBins; ++i) {
			systValues[systValues.size() - 2].push_back(general[y][i][2]);
			systValues.back().push_back(stat[y][i][0]);
		}
		DrawBinGraphs(Form("../results/%s/Graphs/general_systematics_plot_%s_%s_%s.pdf", YEARS[y].Data(), TREE.Data(), YEARS[y].Data(), VAR.Data()),
		              FitVarLabel(VAR), "Uncertainty (%)", xMean, xLow, xHigh, systLabels, systValues,
		              std::vector<std::vector<double>>(systValues.size(), std::vector<double>(nBins, 0.0)));
	}
}
