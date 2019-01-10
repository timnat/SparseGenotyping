/*
SparseGenotyping - is a program that allows to define a consensus genotype for scaffolds based on low-coverage SNP calls by accumulating number of reads that support one or another genotype across the length of a scaffold

Compilation
g++ SparseGenotyping.cpp -o SparseGenotyping

Usage details, test data and updates can be found at 
https://github.com/timnat/SparseGenotyping

Notice: 
in the study, this program was designed for, genotypes of first two (F=2) individuals were predefined (first was A.Mexicanum, second was A.Tigranum) and used to filter vcf file to keep only SNPs that truly differentiate between A.Mexicanum and A.Tigranum. So the genotypes of interest were genotypes of only 48 individuals (excluding first two given in vcf file). TO DO task is to generalize the program for any number of samples in vcf file. In this version one have to change NInd and F according to the given vcf file and purposes of their study and recompile the program.

*/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <sstream>
using namespace std;

#define TT 1000
#define NInd 50  //NInd is a number of individuals/samples in vcf file. 
#define F 2

//-------------------------------------------
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
//-------------------------------------------
char gtype_char(int AM_rd, int AT_rd, int min_read_number,float min_r,float max_r, float e, float &ratio_AM_to_sum)
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

int main(int argc, char *argv[])
{
    int opt,min_read_number;
    float e,min_r, max_r, ratio_AM_to_sum;
    char fname[TT]="", sss[TT], sum_out[TT]="", s_out[TT]="", gs_out[TT]="", rs_out[TT]="";

   min_read_number=4;
   min_r=0.25;
   max_r=0.75;
   e=0.8;
   /* Parsing command line*/
    while ((opt = getopt(argc, argv, "r:m:M:e:")) != -1) {
        switch (opt) {
        case 'r':
	    min_read_number = atoi(optarg);
            break;
        case 'm':
	    min_r = atof(optarg);
            break;
        case 'M':
	    max_r = atof(optarg);
            break;
        case 'e':
            e = atof(optarg);
            break;
        default: 
            fprintf(stderr, "Usage: %s [-r min_read_number [4]] [-m min_ratio [0.25]] [-M max_ratio [0.75]] [-e homoz_min_ratio [0.8]] vcf_file\n",
                    argv[0]);
            exit(EXIT_FAILURE);
        }
    }


   if (optind >= argc) {
        fprintf(stderr, "Expected vcf_file\n");
	fprintf(stderr, "Usage: %s [-r min_read_number [4]] [-m min_ratio [0.25]] [-M max_ratio [0.75]] [-e homoz_min_ratio [0.8]] vcf_file\n",
                    argv[0]);
        exit(EXIT_FAILURE);
    }

   strncat(fname,argv[optind],strlen(argv[optind]));
   printf("vcf_file: %s\nmin_read_number=%d\nmin_ratio=%.2f\nmax_ratio=%.2f\ne=%.2f\n", fname,min_read_number,min_r,max_r,e);
   
   std::ifstream infile( fname );
    if (!infile) {
	fprintf(stderr, "Can't open vcf_file %s\n",fname);	
        exit(EXIT_FAILURE);
    }

   //puts("prepear output files");
  //-- outputs----------
    sprintf(sss,"%d",min_read_number); strncat(fname,".r",2); strncat(fname,sss,strlen(sss));
    sprintf(sss,"%.2f",min_r);         strncat(fname,".",1);  strncat(fname,sss,strlen(sss));
    sprintf(sss,"%.2f",max_r);         strncat(fname,"_",1);  strncat(fname,sss,strlen(sss));
    strncpy(sum_out,fname,strlen(fname));
    cout<<"Output will be in files with prefix: "<< sum_out <<endl;	

    strncat(sum_out,".sum", 4);
    std::ofstream sumf_out(sum_out);
    //cout<<"sumf_out will be in "<< sum_out <<endl;

    strncpy(s_out,fname,strlen(fname));	strncat(s_out,".rd", 3);
    std::ofstream f_out(s_out);
    //cout<<"Read depth output will be in "<< s_out <<endl;

    strncpy(rs_out,fname,strlen(fname));	
    strncat(rs_out,".ratios", 7);
    std::ofstream r_out(rs_out) ;
    //cout<<"Read depth ratios output will be in "<< rs_out <<endl;

    strncpy(gs_out,fname,strlen(fname));
    strncat(gs_out,".gt", 3);
    std::ofstream g_out(gs_out) ;
    cout<<"Genotypes output will be in "<< gs_out <<endl;

    //cout<<"0 - AM/AM; 1 - AM/AT; - low count of reads; ~ ratio is out of scope \n";
    cout<<"Genotypes labeled with\n 0 - homoz;\n 1 - heteroz;\n - number of reads is too low;\n ~ read numbers ratio is out of scope\n";


   /* Main code*/
    int begin=0, sum_AM=0, sum_AMT=0, sum_low=0, sum_lim=0, linecount=0, i, AM[NInd],AT[NInd];
    int Znumber=0, Unumber=0, dash_number=0, tilda_number=0;
    string line, Name,Name_pred="", ID1, genotype[NInd];
    char G;
   

    for(i=0;i<NInd;i++) AM[i]=AT[i]=0;

    while ( getline( infile , line ) ) {

 	 linecount++;
	 std::stringstream ss1(line);
         getline(ss1,Name,'\t');	

	 /* skip header lines */
	 if(begin==0)
	 {
	  if(Name=="#CHROM")
	  {
	    for (i=0;i<8;i++) getline(ss1,ID1, '\t'); 
            f_out<<"Contig_name\t"<<"number_of_SNPsites\t";
 	    g_out<<"Contig_name\t"<<"number_of_SNPsites\t";
 	    r_out<<"Contig_name\t"<<"number_of_SNPsites\t";
	    sumf_out<<"Contig_name\t"<<"number_of_SNPsites\t"<<"number_of_0(AM/AM)\t"<<"number_of_1(AM/AT)\t"<<"number_of_-(low_number_of_reads)\t"<<"number_of_~(out of ["<<min_r<<","<<max_r<<"])\n";
            for (i=0;i<NInd;i++)
		{ getline(ss1,ID1, '\t'); 
		  f_out<<ID1<<"_M\t"<<ID1<<"_T\t";
		  g_out<<ID1<<"\t";
		  r_out<<ID1<<"\t";
		}
		f_out << endl;
		g_out << endl;
		r_out << endl;
		sumf_out << endl;

	    begin=1;
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
	  
	   for(i=0;i<NInd;i++) 
	    {
	     G = gtype_char(AM[i],AT[i],min_read_number,min_r,max_r,e,ratio_AM_to_sum);
	     g_out<< G <<"\t";  
	     r_out<< (int(100*(ratio_AM_to_sum)))/100.0 <<"\t";  

	     if(i>=F) //skip AM in column 0, and ATT in column 1
	      {
	       if (G=='0') sum_AM++, Znumber++;
	       if (G=='1') sum_AMT++, Unumber++;
	       if (G=='-') sum_low++, dash_number++;
	       if (G=='~') sum_lim++, tilda_number++;
	      }
            }

	   sumf_out<<sum_AM<<"\t"<<sum_AMT<<"\t"<<sum_low<<"\t"<<sum_lim<<"\n";

	   f_out<<Name_pred<<"\t"<<linecount-1<<"\t";
           for(i=0;i<NInd;i++) 
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
	 for(int i=0;i<NInd;i++)
	   {
            getline(ss1,genotype[i], '\t'); 
         	//cout << "genotype"<<i<<": " << genotype[i] << endl;
	    extractAM_AT(i,genotype[i],AM,AT);	    
	   }

	Name_pred=Name;
     }

   	//process last line
	   linecount++;
	   g_out<<Name_pred<<"\t"<<linecount-1<<"\t";
	   r_out<<Name_pred<<"\t"<<linecount-1<<"\t";
	   sumf_out<<Name_pred<<"\t"<<linecount-1<<"\t";

	   for(i=0;i<NInd;i++) 
             {
	       G = gtype_char(AM[i],AT[i],min_read_number,min_r,max_r,e,ratio_AM_to_sum);
	       g_out<< G <<"\t";  
	       r_out<< (int(100*(ratio_AM_to_sum)))/100.0 <<"\t"; 

              if(i>=F) 
	       {
	        if (G=='0') sum_AM++, Znumber++;
	        if (G=='1') sum_AMT++, Unumber++;
	        if (G=='-') sum_low++, dash_number++;
	        if (G=='~') sum_lim++, tilda_number++;
	       }
             }
	   f_out<<Name_pred<<"\t"<<linecount-1<<"\t";

           for(i=0;i<NInd;i++) 
	    {
	     f_out<<AM[i]<<"\t"<<AT[i]<<"\t";
 	     AM[i]=AT[i]=0;
	    }	
           f_out<<"\n";
	   sumf_out<<sum_AM<<"\t"<<sum_AMT<<"\t"<<sum_low<<"\t"<<sum_lim<<"\n";


	    cout<<"Stats summary\n";
            cout<<"Number of 0 (i.e. Mex, M/(M+T)>e)="<<Znumber<<"\n";
            cout<<"Number of 1 (i.e. MexTig)="<<Unumber<<"\n";
            cout<<"Number of - (i.e. too few reads, min_read_number= "<< min_read_number <<" )="<<dash_number<<"\n";
            cout<<"Number of ~ (i.e. out of scope[] and M/(M+T)<e)="<<tilda_number<<"\n";

   infile.close();
   f_out.close();   g_out.close();   r_out.close(); sumf_out.close();

   exit(EXIT_SUCCESS);
}
