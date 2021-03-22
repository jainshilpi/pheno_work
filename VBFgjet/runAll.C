#include <string>
//#include "TChain.h
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "Analyse.C+"
#else
#include "Analyse.C"
#endif

#include <fstream>

void runAll()
{

  cout<<"inside runAll.C"<<endl;


  
  //double lumi = 41*1000;  //pb-1
  double lumi = 1;  //pb-1
  cout<<"now making class object "<<endl;

  Analyse t;
  
  cout<<"done making class object "<<endl;
  
  
  cout<<"====Now running over samples===="<<endl;

  ifstream is("files.list");
  //ifstream is("files_data.list");
  string str;
  while(getline(is, str))
    {

      string fname, fnameout;
      int isBkg, itype;
      string xsec_st;
      is >> fname >> fnameout >> isBkg >> xsec_st >> itype;
      int fstr = fname.find("#",0);
      cout<<"fstr "<<fstr<<endl;
      cout<<"fname : fnameout: isBkg : xsec_st : itype : "<<fname<<" "<<fnameout<<" "<<isBkg<<" "<<xsec_st<<" "<<itype<<endl;
      if(fstr!=string::npos)
	{
	  continue;
	}
            
      
      double xsec = atof(xsec_st.c_str());
      cout<<" ==== xsec === "<<xsec<<endl;
      xsec = xsec*lumi;
      
      //t.Loop(fname, fnameout, xsec, isBkg);      
      t.Loop(fname, fnameout);      
    }
  


  
}//void runAll()
