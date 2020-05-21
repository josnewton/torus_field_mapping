#include "Riostream.h"
#include "TROOT.h"
#include "TClass.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TVirtualFitter.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TBrowser.h"
#include "TObjString.h"
#include "TError.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TString.h"
#include "TNtuple.h"
#include "TFile.h"
#include <fstream>
#include "TGraph.h"
#include <string> 
#include "TLegend.h"

using namespace std;

struct fieldcollection {
  double bvalue1;
  double bvalue2;
  double bvalue3;
  double bvalue4;
  double bvalue5;
  double bvalue6;

};


fieldcollection getDistortionValues(int iteration, double arr[], double arr2[], int a1, int a2, int a3, int a4, int a5, int a6);

void filldata(TString file1, TString file2, TString file3, TString file4, TString ljj, TString file5, TString file6);



void data_format_script()
{
  
  TString root_file = "fittingdata.root";
  TFile *file = new TFile(root_file,"RECREATE");

  TString measured = "measured.txt";
  //TString measured = "newfield.txt";
  //TString nominal = "geometry.txt";
  TString nominal = "geometry_aug.txt";
  TString dist1 = "distortion1.txt";   
  TString dist2 = "distortion2.txt";
  TString dist3 = "distortion3.txt";
  TString dist4 = "lastdistortion.txt";
  TString outputfile = "alldata.dat";
  
 
  filldata(nominal,dist1,dist2,dist3,dist4,outputfile,measured);



}


void filldata(TString file1, TString file2, TString file3, TString file4, TString ljj, TString file5, TString file6){
  ifstream inFile;
  ifstream inFile2;
  ifstream inFile3;
  ifstream inFile4;
  ofstream outFile;
  ifstream inFile6;
  ifstream jFile;


  inFile.open(file1);
  inFile2.open(file2);
  inFile3.open(file3);
  inFile4.open(file4);
  outFile.open(file5);
  inFile6.open(file6);
  jFile.open(ljj);




  Int_t sectornom[1302],sectordist1[1302],sectordist2[1302], sectordist3[1302], sectordist4[1302], sectormeas[1302], sector;
  Double_t znom[1302], zdist1[1302], zdist2[1302], zdist3[1302],zdist4[1302],zmeas[1302], bxnom[1302], bynom[1302], bznom[1302], z, bx, by, bz;
  Double_t bxdist1[1302], bxdist2[1302], bxdist3[1302], bxdist4[1302];
  Double_t bydist1[1302], bydist2[1302], bydist3[1302], bydist4[1302];
  Double_t bzdist1[1302], bzdist2[1302], bzdist3[1302], bzdist4[1302];
  Double_t bxmeas[1302], bymeas[1302], bzmeas[1302];
  Double_t transform = 0.0;
  Double_t current = 0.0;
  current = 3000.0/3770.0;
  //  cout << "*******************************" << endl;

  if(!inFile)
    {
      cout << "Unable to open file" << endl;
      exit(0);
    }
  else
    {
 for(Int_t i1 = 0 ; i1 <= 1007; i1++){
	      inFile >> sector >> z >> bx >> by >> bz;
	      transform = 60.0*(sector-1);
        sectornom[i1] = sector;
        znom[i1] = z;
        bxnom[i1] = (bx*cos(-transform/57.2958) - by*sin(-transform/57.2958))*10000*current;
        bynom[i1] = (bx*sin(-transform/57.2958) + by*cos(-transform/57.2958))*10000*current;
	    bznom[i1] = bz*10000*current;
	
      }
    }

  inFile.close();
 
  if(!inFile2)
    {
      cout << "Unable to open file" << endl;
      exit(0);
    }
  else
    {
      for(Int_t i1 = 0 ; i1 <= 1007; i1++){
	      inFile2 >> sector >> z >> bx >> by >> bz;
        transform = 60.0*(sector-1);
        sectordist1[i1] = sector;
        zdist1[i1] = z;
        bxdist1[i1] = (bx*cos(-transform/57.2958) - by*sin(-transform/57.2958))*10000*current;
        bydist1[i1] = (bx*sin(-transform/57.2958) + by*cos(-transform/57.2958))*10000*current;
	      bzdist1[i1] = bz*10000*current;

      }
    }

  inFile2.close();

  if(!inFile3)
    {
      cout << "Unable to open file" << endl;
      exit(0);
    }
  else
    {
      for(Int_t i1 = 0 ; i1 <= 1007; i1++){
	    inFile3 >> sector >> z >> bx >> by >> bz;
	      transform = 60.0*(sector-1);
        sectordist2[i1] = sector;
        zdist2[i1] = z;
        bxdist2[i1] = (bx*cos(-transform/57.2958) - by*sin(-transform/57.2958))*10000*current;
        bydist2[i1] = (bx*sin(-transform/57.2958) + by*cos(-transform/57.2958))*10000*current;
	      bzdist2[i1] = bz*10000*current;

      }
    }

  inFile3.close();

  if(!inFile4)
    {
      cout << "Unable to open file" << endl;
      exit(0);
    }
  else
    {
      for(Int_t i1 = 0 ; i1 <= 1007; i1++){
	inFile4 >> sector >> z >> bx >> by >> bz;
	transform = 60*(sector-1);
        sectordist3[i1] = sector;
        zdist3[i1] = z;
        bxdist3[i1] = (bx*cos(-transform/57.2958) - by*sin(-transform/57.2958))*10000*current;
        bydist3[i1] = (bx*sin(-transform/57.2958) + by*cos(-transform/57.2958))*10000*current;
	bzdist3[i1] = bz*10000*current;

      }
    }

  inFile4.close();

  if(!jFile)
    {
      cout << "Unable to open file" << endl;
      exit(0);
    }
  else
    {
      for(Int_t i1 = 0 ; i1 <= 1007; i1++){
        jFile >> sector >> z >> bx >> by >> bz;
        transform = 60*(sector-1);
        sectordist4[i1] = sector;
        zdist4[i1] = z;
        bxdist4[i1] = (bx*cos(-transform/57.2958) - by*sin(-transform/57.2958))*10000*current;
        bydist4[i1] = (bx*sin(-transform/57.2958) + by*cos(-transform/57.2958))*10000*current;
        bzdist4[i1] = bz*10000*current;

      }
    }

  jFile.close();


  if(!inFile6)
    {
      cout << "Unable to open file" << endl;
      exit(0);
    }
  else
    {
      for(Int_t i1 = 0 ; i1 <= 1007; i1++){
        inFile6 >> sector >> z >> bx >> by >> bz;
        sectormeas[i1] = sector;
        zmeas[i1] = z;
        transform = 60*(sector-1);
        bxmeas[i1] = bx*10000;
        bymeas[i1] = by*10000;
        bzmeas[i1] = bz*10000;

      }
    }

  inFile6.close();


  //dBx Distortions due to effect 1, 2, 3, ..., 18

  
  TH1D *h1 = new TH1D("h1", "dBmag", 200, -500.0, 500.0);
  h1->SetXTitle("Difference Between Total Magnitude of Measured Field and Newly Re-Fitted Field (Gauss)");

  TH1D *hx = new TH1D("hx", "dBx", 200, -500.0, 500.0);
  hx->SetXTitle("Difference Between X Value of Measured Field and Newly Re-Fitted Field (Gauss)");

  TH1D *hy = new TH1D("hy", "dBy", 200, -500.0, 500.0);
  hy->SetXTitle("Difference Between Y Value of Measured Field and Newly Re-Fitted Field (Gauss)");

  TH1D *hz = new TH1D("hz", "dBz", 200, -500.0, 500.0);
  hz->SetXTitle("Difference Between Z Value of Measured Field and Newly Re-Fitted Field (Gauss)");

  Double_t z1[42], z2[42], z3[42], z4[42];
  Double_t dbmod1[42], dbmod2[42], dbmod3[42], dbmod4[42];

  for(Int_t i1 = 0 ; i1 <= 1007; i1++) {

        Int_t run_size = 168;
        Double_t thesector = sectornom[i1];
        Double_t thez = znom[i1];

	Double_t dBx[18], dBy[18], dBz[18];
      
	
	fieldcollection dbx1to6, dbx7to12, dbx13to18, dbx19to24; 
	fieldcollection dby1to6, dby7to12, dby13to18, dby19to24; 
	fieldcollection dbz1to6, dbz7to12, dbz13to18, dbz19to24;



	

	if(thesector==1) {

	  
	  dbx1to6   = getDistortionValues(i1, bxnom, bxdist1, 0, 1, 2, 3, 4, 5);
	  dbx7to12  = getDistortionValues(i1, bxnom, bxdist2, 0, 1, 2, 3, 4, 5);
	  dbx13to18 = getDistortionValues(i1, bxnom, bxdist3, 0, 1, 2, 3, 4, 5);
	  dbx19to24 = getDistortionValues(i1, bxnom, bxdist4, 0, 1, 2, 3, 4, 5);
	  dby1to6   = getDistortionValues(i1, bynom, bydist1, 0, 1, 2, 3, 4, 5);
	  dby7to12  = getDistortionValues(i1, bynom, bydist2, 0, 1, 2, 3, 4, 5);
	  dby13to18  = getDistortionValues(i1, bynom, bydist3, 0, 1, 2, 3, 4, 5);
	  dby19to24 = getDistortionValues(i1, bynom, bydist4, 0, 1, 2, 3, 4, 5);
	  dbz1to6   = getDistortionValues(i1, bznom, bzdist1, 0, 1, 2, 3, 4, 5);
	  dbz7to12  = getDistortionValues(i1, bznom, bzdist2, 0, 1, 2, 3, 4, 5);
	  dbz13to18  = getDistortionValues(i1, bznom, bzdist3, 0, 1, 2, 3, 4, 5);
	  dbz19to24 = getDistortionValues(i1, bznom, bzdist4, 0, 1, 2, 3, 4, 5);
	  

	  }

  
          if(thesector==2) {

	  dbx1to6   = getDistortionValues(i1, bxnom, bxdist1, 0, 1, 2, 3, 4, -1);
	  dbx7to12  = getDistortionValues(i1, bxnom, bxdist2, 0, 1, 2, 3, 4, -1);
	  dbx13to18 = getDistortionValues(i1, bxnom, bxdist3, 0, 1, 2, 3, 4, -1);
	  dbx19to24 = getDistortionValues(i1, bxnom, bxdist4, 0, 1, 2, 3, 4, -1);
	  dby1to6   = getDistortionValues(i1, bynom, bydist1, 0, 1, 2, 3, 4, -1);
	  dby7to12  = getDistortionValues(i1, bynom, bydist2, 0, 1, 2, 3, 4, -1);
	  dby13to18  = getDistortionValues(i1, bynom, bydist3, 0, 1, 2, 3, 4, -1);
	  dby19to24 = getDistortionValues(i1, bynom, bydist4, 0, 1, 2, 3, 4, -1);
	  dbz1to6   = getDistortionValues(i1, bznom, bzdist1, 0, 1, 2, 3, 4, -1);
	  dbz7to12  = getDistortionValues(i1, bznom, bzdist2, 0, 1, 2, 3, 4, -1);
	  dbz13to18  = getDistortionValues(i1, bznom, bzdist3, 0, 1, 2, 3, 4, -1);
	  dbz19to24 = getDistortionValues(i1, bznom, bzdist4, 0, 1, 2, 3, 4, -1);

	  }

          if(thesector==3) {
	  
	  dbx1to6   = getDistortionValues(i1, bxnom, bxdist1, 0, 1, 2, 3, -2, -1);
	  dbx7to12  = getDistortionValues(i1, bxnom, bxdist2, 0, 1, 2, 3, -2, -1);
	  dbx13to18 = getDistortionValues(i1, bxnom, bxdist3, 0, 1, 2, 3, -2, -1);
	  dbx19to24 = getDistortionValues(i1, bxnom, bxdist4, 0, 1, 2, 3, -2, -1);
	  dby1to6   = getDistortionValues(i1, bynom, bydist1, 0, 1, 2, 3, -2, -1);
	  dby7to12  = getDistortionValues(i1, bynom, bydist2, 0, 1, 2, 3, -2, -1);
	  dby13to18  = getDistortionValues(i1, bynom, bydist3, 0, 1, 2, 3, -2, -1);
	  dby19to24 = getDistortionValues(i1, bynom, bydist4, 0, 1, 2, 3, -2, -1);
	  dbz1to6   = getDistortionValues(i1, bznom, bzdist1, 0, 1, 2, 3, -2, -1);
	  dbz7to12  = getDistortionValues(i1, bznom, bzdist2, 0, 1, 2, 3, -2, -1);
	  dbz13to18  = getDistortionValues(i1, bznom, bzdist3, 0, 1, 2, 3, -2, -1);
	  dbz19to24 = getDistortionValues(i1, bznom, bzdist4, 0, 1, 2, 3, -2, -1);

	  }

          if(thesector==4) {

	  dbx1to6   = getDistortionValues(i1, bxnom, bxdist1, 0, 1, 2, -3, -2, -1);
	  dbx7to12  = getDistortionValues(i1, bxnom, bxdist2, 0, 1, 2, -3, -2, -1);
	  dbx13to18 = getDistortionValues(i1, bxnom, bxdist3, 0, 1, 2, -3, -2, -1);
	  dbx19to24 = getDistortionValues(i1, bxnom, bxdist4, 0, 1, 2, -3, -2, -1);
	  dby1to6   = getDistortionValues(i1, bynom, bydist1, 0, 1, 2, -3, -2, -1);
	  dby7to12  = getDistortionValues(i1, bynom, bydist2, 0, 1, 2, -3, -2, -1);
	  dby13to18  = getDistortionValues(i1, bynom, bydist3, 0, 1, 2, -3, -2, -1);
	  dby19to24 = getDistortionValues(i1, bynom, bydist4, 0, 1, 2, -3, -2, -1);
	  dbz1to6   = getDistortionValues(i1, bznom, bzdist1, 0, 1, 2, -3, -2, -1);
	  dbz7to12  = getDistortionValues(i1, bznom, bzdist2, 0, 1, 2, -3, -2, -1);
	  dbz13to18  = getDistortionValues(i1, bznom, bzdist3, 0, 1, 2, -3, -2, -1);
	  dbz19to24 = getDistortionValues(i1, bznom, bzdist4, 0, 1, 2, -3, -2, -1);

	  }

          if(thesector==5) {
	  
	  dbx1to6   = getDistortionValues(i1, bxnom, bxdist1, 0, 1, -4, -3, -2, -1);
	  dbx7to12  = getDistortionValues(i1, bxnom, bxdist2, 0, 1, -4, -3, -2, -1);
	  dbx13to18 = getDistortionValues(i1, bxnom, bxdist3, 0, 1, -4, -3, -2, -1);
	  dbx19to24 = getDistortionValues(i1, bxnom, bxdist4, 0, 1, -4, -3, -2, -1);
	  dby1to6   = getDistortionValues(i1, bynom, bydist1, 0, 1, -4, -3, -2, -1);
	  dby7to12  = getDistortionValues(i1, bynom, bydist2, 0, 1, -4, -3, -2, -1);
	  dby13to18  = getDistortionValues(i1, bynom, bydist3, 0, 1, -4, -3, -2, -1);
	  dby19to24 = getDistortionValues(i1, bynom, bydist4, 0, 1, -4, -3, -2, -1);
	  dbz1to6   = getDistortionValues(i1, bznom, bzdist1, 0, 1,-4, -3, -2, -1);
	  dbz7to12  = getDistortionValues(i1, bznom, bzdist2, 0, 1, -4, -3, -2, -1);
	  dbz13to18  = getDistortionValues(i1, bznom, bzdist3, 0, 1, -4, -3, -2, -1);
	  dbz19to24 = getDistortionValues(i1, bznom, bzdist4, 0, 1, -4, -3, -2, -1);


	  }

          if(thesector==6) {

	  dbx1to6   = getDistortionValues(i1, bxnom, bxdist1, 0, -5, -4, -3, -2, -1);
	  dbx7to12  = getDistortionValues(i1, bxnom, bxdist2, 0, -5, -4, -3, -2, -1);
	  dbx13to18 = getDistortionValues(i1, bxnom, bxdist3, 0, -5, -4, -3, -2, -1);
	  dbx19to24 = getDistortionValues(i1, bxnom, bxdist4, 0, -5, -4, -3, -2, -1);
	  dby1to6   = getDistortionValues(i1, bynom, bydist1, 0, -5, -4, -3, -2, -1);
	  dby7to12  = getDistortionValues(i1, bynom, bydist2, 0, -5, -4, -3, -2, -1);
	  dby13to18  = getDistortionValues(i1, bynom, bydist3, 0, -5, -4, -3, -2, -1);
	  dby19to24 = getDistortionValues(i1, bynom, bydist4, 0, -5, -4, -3, -2, -1);
	  dbz1to6   = getDistortionValues(i1, bznom, bzdist1, 0, -5,-4, -3, -2, -1);
	  dbz7to12  = getDistortionValues(i1, bznom, bzdist2, 0, -5, -4, -3, -2, -1);
	  dbz13to18  = getDistortionValues(i1, bznom, bzdist3, 0, -5, -4, -3, -2, -1);
	  dbz19to24 = getDistortionValues(i1, bznom, bzdist4, 0, -5, -4, -3, -2, -1);

	      }
	
	  if(thez >= 2929.24 && thez <= 4779.24) {
	    




	    double bx1, by1, bz1, bx2, by2, bz2;
	    bx1 = bxnom[i1];
	    by1 = bynom[i1];
	    bz1 = bznom[i1];

	    bx2 = bxmeas[i1];
	    by2 = bymeas[i1];
	    bz2 = bzmeas[i1];

	    
	    double dbmag = sqrt(bx1*bx1 + by1*by1 + bz1*bz1) - sqrt(bx2*bx2 + by2*by2 + bz2*bz2);
	    double dbxmag = bx1 - bx2;
	    double dbymag = by1 - by2;
	    double dbzmag = bz1 - bz2;

	  

	    h1->Fill(dbmag);
	    hx->Fill(dbxmag);
	    hy->Fill(dbymag);
	    hz->Fill(dbzmag);

	
	    outFile << (bxnom[i1]-bxmeas[i1]) << "  " << (bynom[i1]-bymeas[i1]) << "  " << (bznom[i1]-bzmeas[i1]) << " ";
	    outFile << dbx1to6.bvalue1 << " " << dbx1to6.bvalue2 << " ";
	    outFile <<  dbx1to6.bvalue3 << " " << dbx1to6.bvalue4 << " " << dbx1to6.bvalue5 << " " << dbx1to6.bvalue6 << " "; 
	    outFile << dbx7to12.bvalue1 <<  " " << dbx7to12.bvalue2 << " " << dbx7to12.bvalue3 << " " << dbx7to12.bvalue4 << " "; 
	    outFile << dbx7to12.bvalue5 << " " << dbx7to12.bvalue6 << " " << dbx13to18.bvalue1 << " " << dbx13to18.bvalue2 << " "; 
	    outFile << dbx13to18.bvalue3 << " " << dbx13to18.bvalue4 << " " << dbx13to18.bvalue5 << " " << dbx13to18.bvalue6 << " "; 
	    outFile << dbx19to24.bvalue1 << " " << dbx19to24.bvalue2 << " " << dbx19to24.bvalue3 << " " << dbx19to24.bvalue4 << " " ;
	    outFile << dbx19to24.bvalue5 << " " << dbx19to24.bvalue6 << " ";
	    outFile <<  dby1to6.bvalue1 << " " << dby1to6.bvalue2 << " " << dby1to6.bvalue3 << " " << dby1to6.bvalue4 << " "; 
	    outFile << dby1to6.bvalue5 <<  " " << dby1to6.bvalue6 << " " << dby7to12.bvalue1 << " " << dby7to12.bvalue2 << " "; 
	    outFile << dby7to12.bvalue3 << " " << dby7to12.bvalue4 << " " << dby7to12.bvalue5 << "  " << dby7to12.bvalue6 << " "; 
	    outFile << dby13to18.bvalue1 << " " << dby13to18.bvalue2 << " " << dby13to18.bvalue3 << " " << dby13to18.bvalue4 << " ";
	    outFile << dby13to18.bvalue5 << " " << dby13to18.bvalue6 << " ";
            outFile << dby19to24.bvalue1 << " "<< dby19to24.bvalue2 <<" " << dby19to24.bvalue3 << " "<< dby19to24.bvalue4 <<" " ;
            outFile << dby19to24.bvalue5 << " "<< dby19to24.bvalue6 <<" ";
	    outFile <<  dbz1to6.bvalue1 << " " << dbz1to6.bvalue2 << " " << dbz1to6.bvalue3 << " " << dbz1to6.bvalue4 << " "; 
	    outFile << dbz1to6.bvalue5 <<  " " << dbz1to6.bvalue6 << " " << dbz7to12.bvalue1 << " " << dbz7to12.bvalue2 << " "; 
	    outFile << dbz7to12.bvalue3 << " " << dbz7to12.bvalue4 << " " << dbz7to12.bvalue5 << "  " << dbz7to12.bvalue6 << " "; 
	    outFile << dbz13to18.bvalue1 << " " << dbz13to18.bvalue2 << " " << dbz13to18.bvalue3 << " " << dbz13to18.bvalue4 << " "; 	            
	    outFile << dbz13to18.bvalue5 << " " << dbz13to18.bvalue6 << " ";
            outFile << dbz19to24.bvalue1 << " "<< dbz19to24.bvalue2 <<" " << dbz19to24.bvalue3 << " "<< dbz19to24.bvalue4 <<" " ;
            outFile << dbz19to24.bvalue5 << " "<< dbz19to24.bvalue6 <<" ";

	    }

	  
  }


  TCanvas *l = new TCanvas("l","Sampling Fraction",200,10,700,500);
  l->Divide(2,2);
  l->cd(1);
  h1->Draw();
  l->cd(2);
  hx->Draw();
  l->cd(3);
  hy->Draw();
  l->cd(4);
  hz->Draw();

}




fieldcollection getDistortionValues(int iteration, double arr[], double arr2[], int a1, int a2, int a3, int a4, int a5, int a6) {
  double b1, b2, b3, b4, b5, b6;
  fieldcollection bgroup;
  int sector_size = 0;

  if(iteration<=1007) {
    sector_size = 168;
  }
   
  if(iteration>=1008) {
    sector_size = 42;
  }

  bgroup.bvalue1 = arr[iteration] - arr2[iteration + sector_size*a1];
  bgroup.bvalue2 = arr[iteration] - arr2[iteration + sector_size*a2];
  bgroup.bvalue3 = arr[iteration] - arr2[iteration + sector_size*a3];
  bgroup.bvalue4 = arr[iteration] - arr2[iteration + sector_size*a4];
  bgroup.bvalue5 = arr[iteration] - arr2[iteration + sector_size*a5];
  bgroup.bvalue6 = arr[iteration] - arr2[iteration + sector_size*a6];

  return bgroup;
} 


