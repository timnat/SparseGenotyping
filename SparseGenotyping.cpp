#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>

/*
  Changes in v3. In v1 and v2 assigned genotype 1 only if number of AT reads is 0
  but it looks like a big portion of those "~" (out of scope) should be actually 1
                 In v3 assign genotype 1 if ratio (M/(M+T) > e
*/

#define e 0.8  //low limit on ratio to assign homozigous genotype e must be > maxr

using namespace std;
int extractAM_AT(int i, string g, int* AM, int *AT)
{string s_genotype,s_AM, s_AT;
 //cout << "g_"<<i<<": " << g << endl;
 stringstream stream_g(g);

 getline(stream_g,s_genotype, ':'); 
 //cout << "g_"<<i<<": " << s_genotype << endl;
 getline(stream_g,s_AM, ','); 
 //cout << "s_AM " << s_AM << endl;
 getline(stream_g,s_AT, ':'); 
 //cout << "s_AT " << s_AT << endl;

 AM[i]+=atoi(s_AM.c_str()); 
 AT[i]+=atoi(s_AT.c_str());
 
 if(s_genotype=="0/0") return 1;
 else return 0;
}
//--------------------------------------------
float gtype_float(int AM_rd, int AT_rd, int min_read_number,float min_r,float max_r)
{
 int sum=0;
 sum=AM_rd+AT_rd;
 if (sum<min_read_number) return -1;
 
 float ratio_AM_to_sum=(float)AM_rd/sum;
 if(ratio_AM_to_sum==1) return 3; 
 if(ratio_AM_to_sum>=min_r && ratio_AM_to_sum<=max_r) return 1;

 cout<<ratio_AM_to_sum<<"\t";
 return ratio_AM_to_sum;

}

//--------------------------------------------
char gtype_char(int AM_rd, int AT_rd, int min_read_number,float min_r,float max_r, float &ratio_AM_to_sum)
{
 int sum=0;
 ratio_AM_to_sum=-1;
 sum=AM_rd+AT_rd;
 if (sum<min_read_number) return '-';
 
 ratio_AM_to_sum=(float)AM_rd/sum;
 /* if(ratio_AM_to_sum==1) return '0'; // was in V1 and V2 */
 if(ratio_AM_to_sum>e) return '0'; //  NEW in V3
 if(ratio_AM_to_sum>=min_r && ratio_AM_to_sum<=max_r) return '1';

// cout<<ratio_AM_to_sum<<"\t";
 return '~';

}

//-------------------------------------------- 
int main( int argc , char** argv ) {

   if(argc<5) 
    {
     cout<<"wrong number of arguments"; 
     cout<<"Usage: program input_thin.vcf min_read_number min_ratio max_ratio" <<endl;
     cout<<"Read depth Output will be in input_thin.vcf.rd";
     cout<<"Read depth ratio Output will be in input_thin.vcf.ratios"; //NEW in v2!!!!!!!!!!
     cout<<"Genotype Output will be in input_thin.vcf.gtv3";
     return 0;
    }
   

   int min_read_number, Znumber=0, Unumber=0, dash_number=0, tilda_number=0;
   float min_r,max_r, ratio_AM_to_sum;
   std::string line1;

   char sss[1000], sum_out[1004]="", s_out[1003]="", gs_out[1003]="", rs_out[1003]="",  G;
   
   string Name,Name_pred="", ID1, genotype[50];
   int i, AM[50],AT[50], linecount=0, BEGIN, sum_AM, sum_AMT, sum_low, sum_lim;
    sprintf(sss,"%s",argv[1]);
//cout<<"len "<<strlen(sss) <<endl;
    cout<<"Input file: "<< sss << endl;
    std::ifstream infile1( sss );
    if (!infile1) {cout<<"Can't open file "<<sss<<endl; return 0;}

    min_read_number=atoi(argv[2]);	
    min_r=atof(argv[3]);
    max_r=atof(argv[4]);
    cout<<"parameters: min_read_number="<<min_read_number<<" min_r="<<min_r<<" max_r="<<max_r<<"\n";


//-- outputs----------
    strncat(sss,".r",2);  strncat(sss,argv[2],strlen(argv[2]));
    strncat(sss,".",1);    strncat(sss,argv[3],strlen(argv[3]));
    strncat(sss,"_",1);    strncat(sss,argv[4],strlen(argv[4]));
    //cout<<"sss "<< sss <<endl;	
    //cout<<"Genotype output will be in "<< sum_out <<endl;	
    strncpy(sum_out,sss,strlen(sss));	
    //cout<<"len "<<strlen(sss) <<endl;
    //cout<<"Genotype output will be in "<< sum_out <<endl;
    strncat(sum_out,".sumv3", 6);
    std::ofstream sumf_out(sum_out) ;
    cout<<"Genotype output will be in "<< sum_out <<endl;

    strncpy(s_out,sss,strlen(sss));	
    strncat(s_out,".rdv3", 5);
    std::ofstream f_out(s_out) ;
    cout<<"Read depth output will be in "<< s_out <<endl;

    strncpy(rs_out,sss,strlen(sss));	
    strncat(rs_out,".ratiosv3", 9);
    std::ofstream r_out(rs_out) ;
    cout<<"Read depth ratios output will be in "<< rs_out <<endl;
  
    strncpy(gs_out,sss,strlen(sss));
    strncat(gs_out,".gtv3", 5);
    std::ofstream g_out(gs_out) ;
    cout<<"Genotype output will be in "<< gs_out <<endl;
    cout<<"0 - AM/AM; 1 - AM/AT; - - low count of reads; ~ ratio is out of scope \n";
 


   BEGIN=0;
   sum_AM=sum_AMT=sum_low=sum_lim=0;
   for(i=0;i<50;i++) AM[i]=AT[i]=0;

   if ( infile1 ) {
      while ( getline( infile1 , line1 ) ) {
 	 linecount++;
	// f_out << linecount << ": " << line1 << '\n' ;//supposing '\n' to be line end

	 stringstream ss1(line1);
         getline(ss1,Name, '\t');	
	 //cout << "ss1: " << Name << endl; 

	 /* skip header lines */
	if(BEGIN==0)
	{
	 if(Name=="#CHROM")
	   {for (i=0;i<8;i++) getline(ss1,ID1, '\t'); 
            f_out<<"Contig_name\t"<<"number_of_SNPsites\t";
 	    g_out<<"Contig_name\t"<<"number_of_SNPsites\t";
 	    r_out<<"Contig_name\t"<<"number_of_SNPsites\t";
	    sumf_out<<"Contig_name\t"<<"number_of_SNPsites\t"<<"number_of_0(AM/AM)\t"<<"number_of_1(AM/AT)\t"<<"number_of_-(low_count_of_reads)\t"<<"number_of_~(out of ["<<min_r<<","<<max_r<<"])\n";
            for (i=0;i<50;i++)
		{ getline(ss1,ID1, '\t'); 
		  f_out<<ID1<<"\t";
		  g_out<<ID1<<"\t";
		  r_out<<ID1<<"\t";
		  //sumf_out<<ID1<<"\t";
		}
		f_out << endl;
		g_out << endl;
		r_out << endl;
		sumf_out << endl;

	    BEGIN=1;
  	    linecount=0; 
	    continue;
	   }
         else continue;	
	}

	 if((Name_pred!="")&&(Name!=Name_pred))
	  {
	   g_out<<Name_pred<<"\t"<<linecount-1<<"\t";
	   r_out<<Name_pred<<"\t"<<linecount-1<<"\t";
	   sumf_out<<Name_pred<<"\t"<<linecount-1<<"\t";
	  
	   for(i=0;i<50;i++) 
	    {G = gtype_char(AM[i],AT[i],min_read_number,min_r,max_r,ratio_AM_to_sum);
	     g_out<< G <<"\t";  
	     r_out<< (int(100*(ratio_AM_to_sum)))/100.0 <<"\t";  

	  //  if(i>0 && i<49)
	  // if(i>=0 && i<48)
	     if(i>=2) //AM in pos 0, ATT in pos 1
	     {
	      if (G=='0') sum_AM++, Znumber++;
	      if (G=='1') sum_AMT++, Unumber++;
	      if (G=='-') sum_low++, dash_number++;
	      if (G=='~') sum_lim++, tilda_number++;
	     }

	    }

	   sumf_out<<sum_AM<<"\t"<<sum_AMT<<"\t"<<sum_low<<"\t"<<sum_lim<<"\n";

	   f_out<<Name_pred<<"\t"<<linecount-1<<"\t";
           for(i=0;i<50;i++) 
	    {
	     f_out<<AM[i]<<"\t"<<AT[i]<<"\t";
 	     AM[i]=AT[i]=0;
	    }	
	   linecount=1;
	   sum_AM=sum_AMT=sum_low=sum_lim=0;
           f_out<<"\n";
           g_out<<"\n";
           r_out<<"\n";
          }
         
         getline(ss1,ID1, 'L'); 
         	//cout << "ID1: " << ID1 << endl;
         getline(ss1,ID1, '\t'); 
         	//cout << "ID1: " << ID1 << endl;
	 for(int i=0;i<50;i++)
	   {
            getline(ss1,genotype[i], '\t'); 
         	//cout << "genotype"<<i<<": " << genotype[i] << endl;
	    extractAM_AT(i,genotype[i],AM,AT);	    
	   }

	Name_pred=Name;
     }
	//last contig
	   linecount++;
	   g_out<<Name_pred<<"\t"<<linecount-1<<"\t";
	   r_out<<Name_pred<<"\t"<<linecount-1<<"\t";
	   sumf_out<<Name_pred<<"\t"<<linecount-1<<"\t";
	   //for(i=0;i<50;i++) g_out<< gtype_char(AM[i],AT[i],min_read_number,min_r,max_r)<<" "; 
	   for(i=0;i<50;i++) 
             {
	       G = gtype_char(AM[i],AT[i],min_read_number,min_r,max_r,ratio_AM_to_sum);
	       g_out<< G <<"\t";  
	      r_out<< (int(100*(ratio_AM_to_sum)))/100.0 <<"\t"; 

              if(i>=2) 
	      // if(i>=0 && i<48)
	       {
	      if (G=='0') sum_AM++, Znumber++;
	      if (G=='1') sum_AMT++, Unumber++;
	      if (G=='-') sum_low++, dash_number++;
	      if (G=='~') sum_lim++, tilda_number++;
	       }
             }
	   f_out<<Name_pred<<"\t"<<linecount-1<<"\t";

           for(i=0;i<50;i++) 
	    {
	     f_out<<AM[i]<<"\t"<<AT[i]<<"\t";
 	     AM[i]=AT[i]=0;
	    }	
           f_out<<"\n";
	   sumf_out<<sum_AM<<"\t"<<sum_AMT<<"\t"<<sum_low<<"\t"<<sum_lim<<"\n";

            cout<<"Number of 0 (i.e. Mex, M(M+T)>e)="<<Znumber<<"\n";
            cout<<"Number of 1 (i.e. MexTig)="<<Unumber<<"\n";
            cout<<"Number of - (i.e. too few reads)="<<dash_number<<"\n";
            cout<<"Number of 0 (i.e. out of scope[] and M(M+T)<e)="<<tilda_number<<"\n";

   infile1.close();
   f_out.close();   g_out.close();   r_out.close(); sumf_out.close();
  }
 else {
  /* could not open directory */
  perror ("can't open files");
 }

    return 0 ;
}


