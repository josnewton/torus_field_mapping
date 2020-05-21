


double myFunction(double x1, double x2, double x3,double x4,double x5, double x6, double x7, double x8, double x9, double x10, double x11, double x12,
		  double x13, double x14, double x15, double x16, double x17, double x18, double x19, double x20, double x21, double x22, double x23,
		  double x24) {
 


  ifstream inFile;
  inFile.open("alldata.dat");
    
  Int_t fitting_size = 912;

   Int_t ncols;
   Int_t nlines = 0;

   
   Float_t field[fitting_size][75];
   Float_t final;
   nlines = 0;



   if(!inFile)
     {
       cout << "Unable to open file" << endl;
       exit(0);
     }
   else
     {
       for(Int_t i1 = 0 ; i1 <= (fitting_size - 1); i1++){
	 double bpoints[75];
	 for(Int_t b = 0; b <= 74 ; b++) {
	   inFile >> bpoints[b];
	   field[i1][b] = bpoints[b]/100;
	   }
       }
     }

   inFile.close();


  
   final = 0;
   int i;
   for(i=0;i<=(fitting_size - 1);i++){
    
     
     Float_t xcomponent, ycomponent, zcomponent;

     double dxmeas = field[i][0];

     double dxcalc = x1*field[i][3] + x2*field[i][4] + x3*field[i][5] + x4*field[i][6] + x5*field[i][7] + x6*field[i][8] + x7*field[i][9]
       + x8*field[i][10] + x9*field[i][11] + x10*field[i][12] + x11*field[i][13] + x12*field[i][14]
       + x13*field[i][15] + x14*field[i][16] + x15*field[i][17] + x16*field[i][18] + x17*field[i][19] + x18*field[i][20]
       + x19*field[i][21] + x20*field[i][22] + x21*field[i][23] + x22*field[i][24] + x23*field[i][25] + x24*field[i][26];//Add 4th DOF


     double dymeas = field[i][1];

     double dycalc = x1*field[i][27] + x2*field[i][28] + x3*field[i][29] + x4*field[i][30] + x5*field[i][31] + x6*field[i][32] + x7*field[i][33]
       + x8*field[i][34] + x9*field[i][35] + x10*field[i][36] + x11*field[i][37] + x12*field[i][38]
       + x13*field[i][39] + x14*field[i][40] + x15*field[i][41] + x16*field[i][42] + x17*field[i][43] + x18*field[i][44]
       + x19*field[i][45] + x20*field[i][46] + x21*field[i][47] + x22*field[i][48] + x23*field[i][49] + x24*field[i][50];//Add 4th DOF


     double dzmeas = field[i][2];

     double dzcalc = x1*field[i][51] + x2*field[i][52] + x3*field[i][53] + x4*field[i][54] + x5*field[i][55] + x6*field[i][56] + x7*field[i][57]
       + x8*field[i][58] + x9*field[i][59] + x10*field[i][60] + x11*field[i][61] + x12*field[i][62] 
       + x13*field[i][63] + x14*field[i][64] + x15*field[i][65] + x16*field[i][66] + x17*field[i][67] + x18*field[i][68]
       + x19*field[i][69] + x20*field[i][70] + x21*field[i][71] + x22*field[i][72] + x23*field[i][73] + x24*field[i][74];//Add 4th DOF


     final = final + pow((dxmeas-dxcalc)/50,2) + pow((dymeas-dycalc)/5,2) + pow((dzmeas-dzcalc)/50,2);	
     //final = final + pow((dymeas-dycalc),2);
     //cout << (dxmeas-dxcalc) << " " << (dymeas-dycalc) << " " << (dzmeas-dzcalc) << endl;
	}

   return final;
}

void minuitFunction(int& nDim, double* gout, double&result, double par[], int flg) {
  result = myFunction(par[0],par[1],par[2],par[3],par[4],par[5],par[6], par[7], par[8] , par[9], par[10], par[11],
		      par[12],par[13],par[14],par[15],par[16],par[17],par[18],par[19],par[20],par[21], par[22], par[23]);

  // result = myFunction(par[0], par[1], par[2], par[3], par[4], par[5]);

}

void twelve2() {
  //TFitter* minimizer = new TFitter(18);
  TFitter* minimizer = new TFitter(24);


  double bestx1 = 0.0, bestx2 = 0.0, bestx3 = 0.0, bestx4 = 0.0, bestx5 = 0.0, bestx6 = 0.0, bestx7 = 0.0;
  double bestx8 = 0.0, bestx9 = 0.0, bestx10 = 0.0, bestx11 = 0.0, bestx12 = 0.0, bestx13 = 0.0, bestx14 = 0.0;
  double bestx15 = 0.0, bestx16 = 0.0, bestx17 = 0.0, bestx18 = 0.0, bestx19 = 0.0, bestx20 = 0.0, bestx21 = 0.0;
  double bestx22 = 0.0, bestx23 = 0.0, bestx24 = 0.0;

  for(int i = 0; i < 20; i++) {

  {
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
  }

  minimizer->SetFCN(minuitFunction);



  minimizer->SetParameter(0,"x1",bestx1,1,0,0);
  minimizer->SetParameter(1,"x2",bestx2,1,0,0);
  minimizer->SetParameter(2,"x3",bestx3,1,0,0);
  minimizer->SetParameter(3,"x4",bestx4,1,0,0);
  minimizer->SetParameter(4,"x5",bestx5,1,0,0);
  minimizer->SetParameter(5,"x6",bestx6,1,0,0);
    
  minimizer->SetParameter(6,"x7",bestx7,1,0,0); 
  minimizer->SetParameter(7,"x8",bestx8,1,0,0);
  minimizer->SetParameter(8,"x9",bestx9,1,0,0);
  minimizer->SetParameter(9,"x10",bestx10,1,0,0);
  minimizer->SetParameter(10,"x11",bestx11,1,0,0);
  minimizer->SetParameter(11,"x12",bestx12,1,0,0);
  
  minimizer->SetParameter(12,"x13",bestx13,1,0,0);
  minimizer->SetParameter(13,"x14",bestx14,1,0,0);
  minimizer->SetParameter(14,"x15",bestx15,1,0,0);
  minimizer->SetParameter(15,"x16",bestx16,1,0,0);
  minimizer->SetParameter(16,"x17",bestx17,1,0,0);
  minimizer->SetParameter(17,"x18",bestx18,1,0,0);

  minimizer->SetParameter(18,"x19",bestx19,1,0,0);
  minimizer->SetParameter(19,"x20",bestx20,1,0,0);
  minimizer->SetParameter(20,"x21",bestx21,1,0,0);
  minimizer->SetParameter(21,"x22",bestx22,1,0,0);
  minimizer->SetParameter(22,"x23",bestx23,1,0,0);
  minimizer->SetParameter(23,"x24",bestx24,1,0,0);


 

  //  minimizer->ExecuteCommand("MIGRAD",0,0);
    minimizer->ExecuteCommand("SIMPLEX",0,0); 
    minimizer->ExecuteCommand("MIGRAD",0,0); 
  


   bestx1 = minimizer->GetParameter(0);
   bestx2 = minimizer->GetParameter(1);
   bestx3 = minimizer->GetParameter(2);
   bestx4 = minimizer->GetParameter(3);
   bestx5 = minimizer->GetParameter(4);
   bestx6 = minimizer->GetParameter(5);
    
   bestx7 = minimizer->GetParameter(6);
   bestx8 = minimizer->GetParameter(7);
   bestx9 = minimizer->GetParameter(8);
   bestx10 = minimizer->GetParameter(9);
   bestx11 = minimizer->GetParameter(10);
   bestx12 = minimizer->GetParameter(11);
  
   bestx13 = minimizer->GetParameter(12);
   bestx14 = minimizer->GetParameter(13);
   bestx15 = minimizer->GetParameter(14);
   bestx16 = minimizer->GetParameter(15);
   bestx17 = minimizer->GetParameter(16);
   bestx18 = minimizer->GetParameter(17);

   bestx19 = minimizer->GetParameter(18);
   bestx20 = minimizer->GetParameter(19);
   bestx21 = minimizer->GetParameter(20);
   bestx22 = minimizer->GetParameter(21);
   bestx23 = minimizer->GetParameter(22);
   bestx24 = minimizer->GetParameter(23);

  
  cout << bestx1 << endl;
  cout << bestx2 << endl;
  cout << bestx3 << endl;
  cout << bestx4 << endl;
  cout << bestx5 << endl;
  cout << bestx6 << endl;
  cout << bestx7 << endl;
  cout << bestx8 << endl;
  cout << bestx9 << endl;
  cout << bestx10 << endl;
  cout << bestx11 << endl;
  cout << bestx12 << endl;
  cout << bestx13 << endl;
  cout << bestx14 << endl;
  cout << bestx15 << endl;
  cout << bestx16 << endl;
  cout << bestx17 << endl;
  cout << bestx18 << endl;
  cout << bestx19 << endl;
  cout << bestx20 << endl;
  cout << bestx21 << endl;
  cout << bestx22 << endl;
  cout << bestx23 << endl;
  cout << bestx24 << endl;


  cout << "Weighted Results " << endl;
  
  cout << "Coil 1 Radial Shift in mm     " << bestx1*2 << endl;
  cout << "Coil 2 Radial Shift in mm     " << bestx6*2 << endl;
  cout << "Coil 3 Radial Shift in mm     " << bestx5*2 << endl;
  cout << "Coil 4 Radial Shift in mm     "  << bestx4*2 << endl;
  cout << "Coil 5 Radial Shift in mm     "  << bestx3*2 << endl;
  cout << "Coil 6 Radial Shift in mm     "  << bestx2*2 << endl;
  
  cout << "Coil 1 Downstream Shift in mm     " << bestx7*2 << endl;
  cout << "Coil 2 Downstream Shift in mm     " << bestx12*2 << endl;
  cout << "Coil 3 Downstream Shift in mm     " << bestx11*2 << endl;
  cout << "Coil 4 Downstream Shift in mm     "  << bestx10*2 << endl;
  cout << "Coil 5 Downstream Shift in mm     "  << bestx9*2 << endl;
  cout << "Coil 6 Downstream Shift in mm     "  << bestx8*2 << endl;

  cout << "Coil 1 Azimuthal Shift in mm     " << bestx13*2 << endl;
  cout << "Coil 2 Azimuthal Shift in mm     " << bestx18*2 << endl;
  cout << "Coil 3 Azimuthal Shift in mm     " << bestx17*2 << endl;
  cout << "Coil 4 Azimuthal Shift in mm     "  << bestx16*2 << endl;
  cout << "Coil 5 Azimuthal Shift in mm     "  << bestx15*2 << endl;
  cout << "Coil 6 Azimuthal Shift in mm     "  << bestx14*2 << endl;

  cout << "Coil 1 Shape Change in mm     " << bestx19*4 << endl;
  cout << "Coil 2 Shape Change in mm     " << bestx24*4 << endl;
  cout << "Coil 3 Shape Change in mm     " << bestx23*4 << endl;
  cout << "Coil 4 Shape Change in mm     "  << bestx22*4 << endl;
  cout << "Coil 5 Shape Change in mm     "  << bestx21*4 << endl;
  cout << "Coil 6 Shape Change in mm     "  << bestx20*4 << endl;

  }
 
}
