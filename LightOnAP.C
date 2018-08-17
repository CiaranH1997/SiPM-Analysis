# SiPM-Analysis
C++ and ROOT code for analysing SiPM data
#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraphPainter.h>
#include <TFile.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <cmath>

// SiPM Analysis Code for reading SiPM waveforms 
// This code finds the trigger peaks and creates a histogram of SER peaks
// Code also contains an alsogrithm to count potential afterpulses

// Author: Ciaran Hasnip
// Last updated: 17/8/2018

using namespace std;

//**********************************************************************

// write smoothing as seperate functions to be applied anywhere as many times as I like
// unweighted moving-average smoothing function; 
// it simply replaces each point in the signal with the average of M adjacent points

vector<double> SmoothingRaw(vector<double> &vecV) {
	// make M 80 as this is approximatly the width of that weird noise fequency
	int M = 80;
	for (int v = 0; v < (vecV.size()-M); v++) {
		double SmV = 0;
		int t = 0;
		do {
			t++;
			SmV += vecV[v+t];			
		} while(t < M);
	vecV[v] = SmV/M;
	}
}

// Function to get the RMS of the baseline noise with an uncertainty on the calculation

double GetRMS(vector<double> &AverageV) {
	
	double N = AverageV.size();
	double X = 0;
	double E = 0;

	for (int k = 0; k < AverageV.size(); k++) {
		X += pow(AverageV[k],2);
	}
	
	double Vrms = sqrt(X/N);
	
	// now find error on Vrms
	for (int q = 0; q < AverageV.size(); q++) {
		E += pow((Vrms - AverageV[q]),2);
	}
	
	double Erms = sqrt(E/N); 
	
	cout << "\nThe RMS noise value is: " << Vrms << " +/- " << Erms << "\n\n";
	return Vrms;
}
	
	


//**********************************************************************************************


// function to plot the raw waveforms
void RawWaveform() {
	
	// make a root file for the waveforms
	TFile* fWave = new TFile("Waveforms.root","RECREATE");	
	// make a canvas
	TCanvas* c = new TCanvas("c1","total waveform",0,0,500,500);
	c->SetGrid();
	
	// multigraph to plot graphs over one another
	TMultiGraph *mg = new TMultiGraph();
	
	// loop over all the data files wanted
	// WARNING too many and this takes a long time and mg gets messy
	for (int i=1; i<41; i++) {
		// read in data
		stringstream filename;
		filename << "oldV54.5a3.5_" << i << ".csv"; 
		
		// new graph for each file
		TGraph *g = new TGraph(filename.str().c_str(),"%lg %lg", ",");
		
		// set style for each graph
		gROOT->SetStyle("Plain");
		g->SetMarkerStyle(kFullDotSmall);
		g->SetMarkerColor(i+2);
		g->SetLineColor(i+2);

		// name each graph and write to file
		stringstream ss;
		ss << "graph" << i;
		g->SetName(ss.str().c_str());
		g->Write();

		g->Draw("AL");
		// add each graph to multigraph
		mg->Add(g);
		// draw multigraph
		mg->Draw("AL");
		mg->GetXaxis()->SetTitle("Time [s]");
		mg->GetYaxis()->SetTitle("Amplitude [V]");
		
		gPad->Update();
		gPad->Modified();
   		mg->SetMinimum(-0.005);
   		mg->SetMaximum(0.02);

		c->Print("SiPM total waveform.pdf");
		}
	// write multigraph to file
	mg->SetName("mg");
	mg->Write();
	fWave->Write();
	}

//********************************************************************************


// Function to extract data from text file, find peak value, plot a histogram and plot the average waveform.
void HistMacro() {

	// make the canvas to draw on
	TCanvas* CHist = new TCanvas("CHist","Peak Amplitude Hist",0,0,500,500);
	CHist->SetGrid();

	TCanvas* CAveW = new TCanvas("CAveW","Average Waveform",0,0,500,500);
	CAveW->SetGrid();

	// create root files for histogram and average waveform
	TFile* fHist = new TFile("Histos.root","RECREATE");
	TFile* fAverage = new TFile("AverageW.root","RECREATE");

	// create histogram
	TH1F* ampHist = new TH1F("Amplitude Hist","amplitude",80,-0.002,0.01);

	// create fitting functions
	// HARDCODE to fit specific histogram
	fit1 = new TF1("m1","gaus",-0.0008,0.0006);
	fit2 = new TF1("m2","gaus",0.0006,0.0019);
	fit3 = new TF1("m3","gaus",0.0019,0.0033);

	total = new TF1("totalFit","gaus(0) + gaus(3) + gaus(6)",-0.0009,0.004);
	
	// declare vectors for average waveform data and afterpulses
	vector<double> AverageV, SumVec, AfterP;
	double Vrms = 0;

	// perform this section of programm twice
	// first loop calculates RMS noise value, the second run subtracts the basiline and plots the histogram
	for (int c=0; c<=1; c++) {
		// loop over all the data files from the oscilloscope
		for (int i=1; i<=900; i++) {
			// define vectors for voltage and time readings from oscilloscope
			vector<double> vecS, vecV;
			double S, V;
	
			// read in the data
			ifstream data;
			
			string filename;
			stringstream ss;
			ss << i;
			// always name files e.g. "newfile_1.csv, newfile_2.csv, etc."
			filename = "oldV54.5a3.5_" + ss.str() + ".csv"; // hard coded for specific file format
			data.open(filename.c_str());
		
			// check file is open, else call error
			if (data.is_open()) {
				cout << filename << endl;
			}
			else { 
				cout << "Error loading file!!!" << "\n";
				break;
			}

			// skip first 2 lines of file as these are headers
			string line; 
			getline(data, line);
			getline(data, line);
		
			// now split the lines of data and put in relevant vectors
			while (getline(data,line)) {
				string sublineS, sublineV;
				int len = line.length();
				size_t comma = line.find(",");
			
				sublineS = line.substr(0,comma);
				S = strtod(sublineS.c_str(),NULL);

				sublineV = line.substr(comma + 1,len);
				V = strtod(sublineV.c_str(),NULL);
			
				vecS.push_back(S);
				vecV.push_back(V);
			}
			

			//***************************************************************		

			// Find AMPLITUDE HISTOGRAM section
		
			// apply smoothing algorithm to raw amplitude data
			SmoothingRaw(vecV);

			// find the max V value for each file in order to get the amplitude
			double temp = -1;
			double maxV2 = 0;

			// ONLY do this if Vrms has been calculated
			if (Vrms > 0 ) {
		
				for (int j=0; j < vecV.size();j++) {
					// take away RMS noise baseline
					// Vrms calculated at end of loop
					double Vnew = vecV[j] - Vrms;
			
					// define a region of interest around the trigger and only include data within that region for the peak
					if (-1.2e-7 < vecS[j] && vecS[j] < -0.8e-7) {
						// use temporary variable to find peak
						if (Vnew > temp) {
		                   	temp = Vnew;
							// variable to fill histogram with
							double maxV = temp;				
						}		
					}	
				}
				// fill histogram with peak values
				ampHist->Fill(maxV);

				/********************************
				* Algorithm to find AFTERPULSES *
				********************************/

				int start = 2000; // make this 2000 so only record start once
				int end = 0;
				// only an afterpulse if it follows a light pulse trigger
				if (maxV > 0.0007) {
					for (int p=0; p < vecV.size();p++) {
						// do this as signal gets weird at the end of data & make sure trigger has decayed a bit
						if (vecS[p] < 2.0e-7 && vecS[p] > -0.5e-7) {
							// take away RMS noise baseline
							// Vrms calculated at end of loop
							double Vnew2 = vecV[p] - Vrms;

							// make a variable threshold
							// I know that p=800 is were time of -0.05us always is and I want the voltage at this time to be my threshold
							double Thresh = vecV[880] - Vrms;
							// look at the regions outside the trigger for peaks
							if (Vnew2 > Thresh) {
								// record start only once
								if (start > p) {
									start = p;
								}
								// use temporary variable to find peak		
								if (Vnew2 > temp) {
									temp = Vnew2;
									// variable to fill histogram with
									maxV2 = temp; 													
								}
							}
							// if maxV2 has been given a value and the waveform has gone below the threshold						
							else if (Vnew2 < Thresh && maxV2 > 0) {
								end = p;
								break;
							}
						}		
					}
				}
				// small oscillations approx. 100 iterations wide 
				// so width must be greater than 50 to be counted as afterpulse not electrical noise
				if (start < 2000) {
					int width = end - start;
					if (maxV2 > 0 && width > 50) {
						AfterP.push_back(maxV2);	
					}
				}
			}
		

			//******************************************************
		
			// Find AVERAGE WAVEFORM section
		
			// only do this if Vrms has not been calculated yet
			if (Vrms == 0) {
				double sumV = 0;
				// declare temporary Sum vector
				vector<double> SumVec2;

				// loop over all values in vecV
				for (int t = 0; t < vecV.size(); t++) {
					if (vecS[t] < -1.4e-7) {  //use this if statement to know the rms noise baseline before trigger
						// if sum vector is not empty
						if (!SumVec.empty()) {
							sumV = SumVec[t] + vecV[t];
							// push back on to temporary vector
							SumVec2.push_back(sumV);
						}
						// if the vector hasnt been filled yet (first file run)
						else {
							SumVec2.push_back(vecV[t]);
						}
					}			
				}
				// swap temporary sum vector into total sum vector
				swap(SumVec,SumVec2);	
				// Write the graphs of the smoothed waveforms to a graph if you want to check the smoothing
				TGraph *RawG = new TGraph(vecV.size(),&vecS[0],&vecV[0]);
				RawG->Write();
			}			
		}
		// END DATA FOR LOOP

		//**************************************************************
	
		// only do this if SumVec has been filled for averaging (i.e. first loop only)
		if (!SumVec.empty() && AverageV.empty()) {
			// number of files
			double n = 900;
			// get average value vector by dividing sum vector by number of files
			for (int z = 0; z < SumVec.size(); z++) {
				AverageV.push_back(SumVec[z]/n);				
			}
			// use GetRMS function to find the RMS noise value and its error
			Vrms = GetRMS(AverageV);
		}
	}
	
	// state how many afterpulses have been found
	cout << "\nThe number of After Pulses identified is:\n\n" << AfterP.size() << " out of 900 waveforms" << "\n\n";

	//**************************************************************

	// DRAW GRAPHS section
	gStyle->SetOptStat("nemruo");
	gROOT->SetStyle("Plain");

	// draw average waveform graph 
	TGraph *gA = new TGraph(AverageV.size(), &(vecS[0]), &(AverageV[0]));
	CAveW->cd();
	gA->SetTitle("Average Waveform");
	gA->GetXaxis()->SetTitle("Time [s]");
	gA->GetYaxis()->SetTitle("Amplitude [V]");
	gA->Draw("AL");
	fAverage->Write();
	gA->Write();
	CAveW->Print("Average Waveform");

	CHist->cd();
	// FIT gaussian to the histogram
	Double_t par[9];

	ampHist->Fit(fit1,"NR");
	ampHist->Fit(fit2,"NR+");		
	ampHist->Fit(fit3,"NR+");
	//ampHist->Fit(fit4,"R+");

	fit1->GetParameters(&par[0]);
	fit2->GetParameters(&par[3]);
	fit3->GetParameters(&par[6]);
	//fit3->GetParameters(&par[9]);
	
	total->SetParameters(par);
	gStyle->SetOptFit(1011); 
	ampHist->Fit(total,"R+");

	// histogram drawing options and write file
	
	ampHist->SetTitle("Amplitude Pk-Pk");
	ampHist->SetName("AmpHistPk-Pk");
	ampHist->GetXaxis()->SetTitle("Amplitude [V]");
	ampHist->GetYaxis()->SetTitle("Frequency [cnts]");
    ampHist->Draw();
	fHist->Write();
	ampHist->Write();
	CHist->Print("SiPM Hist");
	gPad->Update();
	gPad->Modified();
}	


//*******************************************************************************




// Main function for everything

int main() {
	macro2loop();
	HistMacro();
	return 0;
}

// FINISH
