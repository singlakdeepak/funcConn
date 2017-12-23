#ifdef BUILD1
    extern "C" {

		#include <cblas.h>

	}
#elif BUILD2
	extern "C" {

		#include <nvblas.h>
	} 
#else
    #include "mkl.h"
    #include "mkl_cblas.h"
#endif

bool measure = false;
// #ifdef MAT
// 	measure = false;
// #else
// 	measure = true;
// #endif

#include <omp.h>
#include <iostream>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h> 
#include <sys/stat.h>
#include "image/image.hpp"
#include <vector>
#include <fstream>
#include <time.h>
#include <map>
#include <cstdlib>
#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y))

struct Valid{
int x,y,z;
};

std::string type(".csv"),ofname,ipfilename;
bool input = false,roi1 = false,roi2 = false,output=false,gzip = false;
bool mask = false;
std::string roifname,maskfilename;
int thresh=0;
bool seedcmp = false,all=false,flops = false;
int seedx ,seedy,seedz;
int ROI_MAX;
long long no_of_oper=0;
double time_taken =0.0;




void showhelpinfo();
void correl();
void seed_correl();
void all_pair_corr();
void avg_corr_roi();
void avg_roi_time_corr();
void getattributes(int argc,char *argv[]);

int main (int argc,char *argv[])
{
  // omp_set_num_threads(8);
  /*if the program is ran witout options ,it will show the usage and exit*/
  if(argc < 2)
  {
	showhelpinfo();
	exit(1);
  }
  
  if (flops==true)
  {
  		if (measure==true)
  		{
  			std::cout<<"Measuring flops for all ops"<<std::endl;
  		}
  		else
  		{
  			std::cout<<"Measuring flops for only Dgemm"<<std::endl;
  		}
  }
  
 // std::cout<<" 1"<<std::endl;

  getattributes(argc,argv);
  
  if(input==false||output==false ||(roi2 == false && roi1 == false && all == false && seedcmp == false))
  {
	showhelpinfo();
	exit(1);
  }
//std::cout<<" 2"<<std::endl;
  if( ( roi1 == true&&roi2 == true&& all == true ) || ( roi1 == true&&roi2 == true&& seedcmp == true ) || ( all == true && seedcmp == true ))
  {
	showhelpinfo();
	exit(1);
  }
  //std::cout<<" 3"<<std::endl;
	///######################### Creating the project folder #######################################################
   
	mode_t nMode = 0733; // UNIX style permissions
	int nError = 0;
	#if defined(_WIN32)
	  nError = _mkdir(ofname.c_str()); // can be used on Windows
	#else 
	  nError = mkdir(ofname.c_str(),nMode); // can be used on non-Windows
	#endif
	if (nError != 0) {
	  std::cout<<"Problem in creating directory "<<std::endl;
	  exit(0);
	}
  if(roi1==false&&roi2==false){
	std::string sPath = "/voxel";
	sPath = ofname + sPath;
	mode_t nMode = 0733; // UNIX style permissions
	int nError = 0;
	#if defined(_WIN32)
	  nError = _mkdir(sPath.c_str()); // can be used on Windows
	#else 
	  nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
	#endif
	if (nError != 0) {
	  std::cout<<"Problem in creating directory "<<std::endl;
	  exit(0);
	}
  }
  else{
	std::string sPath = "/roi";
  sPath = ofname + sPath;
  mode_t nMode = 0733; // UNIX style permissions
  int nError = 0;
	#if defined(_WIN32)
	  nError = _mkdir(sPath.c_str()); // can be used on Windows
	#else 
	  nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
	#endif
	if (nError != 0) {
	  std::cout<<"Problem in creating directory "<<std::endl;
	  exit(0);
	}
  }

	//############################## END ###########################################################################
  correl();
  return 0;
}


void seed_correl(){

	//DEFINING VARS
  image::io::nifti nifti_parser,nifti_parser2;
  image::basic_image<double,4> image_data;
  image::basic_image<double,4> image_data_cpy;
  image::basic_image<int,3> mask_image;
  std::vector<Valid> valid;

  image_data_cpy = image_data;

  std::string t_ipfilename,t_maskfilename;
  t_ipfilename = ipfilename;
  t_maskfilename = maskfilename;

  if(ipfilename.find(".gz")!=std::string::npos){

  		t_ipfilename = t_ipfilename.substr(0,(t_ipfilename.length()-3));

  }
  if(mask && maskfilename.find(".gz")!=std::string::npos){

  		t_maskfilename = t_maskfilename.substr(0,(t_maskfilename.length()-3));

  }

	//LOADING THE IMAGE
  if(nifti_parser.load_from_file(t_ipfilename));
	  nifti_parser >> image_data;

	// LOADING THE MASK
  if(mask&&nifti_parser2.load_from_file(t_maskfilename))
  	   nifti_parser2 >> mask_image;


if(ipfilename.find(".gz")!=std::string::npos){

		  std::string cmd = "rm ";
  		  cmd+=t_ipfilename;
  		  system(cmd.c_str());

  		

  }
  if(mask && maskfilename.find(".gz")!=std::string::npos){

  	  	  std::string cmd = "rm ";
  		  cmd+=t_maskfilename;
  		 // system(cmd.c_str());

  }


  
  
  image::geometry<4> g = image_data.geometry();
  image::geometry<4> g_copy;
  image::geometry<3> g_mask;

	//CHECK IF THE GEOMETRY OF THE MASK AND IMAGE IS SAME
  if(mask){
  	g_mask = mask_image.geometry();
  	if(g_mask[0]!=g[0]||g_mask[1]!=g[1]||g_mask[2]!=g[2]){
  		std::cout<<":::: ERROR INVALID MASK ::::"<<std::endl;
  		return;
  	}
  }


	//############################# DEFINING MORE VARIABLES ####################################

  double *EX = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  double *EX_sq = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  std::string filename2("/voxelmap");
  filename2 = ofname + filename2;
  const char *s1 = filename2.c_str();
  std::ofstream f(s1,std::ios::binary);

	// #########################CALCULATING THE MEAN AND STD DEVIATION OF DATA AND CREATING A VALID VOXEL MAPPINGS ################################
  
  clock_t tStart;

  for (int z = 0; z < g[2]; ++z){
	for (int y = 0; y < g[1]; ++y)
	  {
		for (int x = 0; x < g[0]; ++x)
		{

		  tStart = clock();
		  double temp  = 0,temp2 = 0;
		  for (int t = 0; t < g[3]; ++t)
		  {
			temp += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3] ;
			temp2 += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]*(double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3];
		  }

		  time_taken+=(double)(clock()-tStart)/CLOCKS_PER_SEC;

		  if(!mask){
			  if(temp > max(20,thresh)){
				Valid tv ;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				int temp = valid.size();
				f.write(reinterpret_cast<char *>(&temp),sizeof(int));
				valid.push_back(tv);
			  }else{
				int temp = -1;
				f.write(reinterpret_cast<char *>(&temp),sizeof(int));
				// direc.insert( ( ( z*g[1] + y )*g[0] + x ) ,-1);
			  }
			}else{
				if(mask_image[(z*g[1]+y)*g[0]+x]==1){
				Valid tv ;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				int temp = valid.size();
				f.write(reinterpret_cast<char *>(&temp),sizeof(int));
				valid.push_back(tv);
			  }else{
				int temp = -1;
				f.write(reinterpret_cast<char *>(&temp),sizeof(int));
				
			  }


			}
		  
		  EX[ (x*g[1] +y)*g[2] + z] = temp;
		  EX_sq[ (x*g[1] +y)*g[2] + z] = sqrt(temp2- temp*temp);

		}
	  } 

  }

  f.close();


	// ########################################### MAKING A MATRIX OF NORMALIZED DATA###########################################
  double * Valid_matrix = (double * )malloc( sizeof(double)*valid.size()*g[3]);

  tStart = clock();
  # pragma omp parallel for
  for (int i = 0; i < valid.size(); ++i)
  {
	for (int t = 0; t < g[3]; ++t){
		double currentdev = EX_sq[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z];
		double tempNormdata;
		if (currentdev < 0.005){
			tempNormdata = 0;
			image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
			Valid_matrix[i*g[3]+t] = tempNormdata;
		}
		else{
			tempNormdata = (double)(image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] - EX[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z])/currentdev;
	  		image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
	  		Valid_matrix[i*g[3]+t] = tempNormdata;
		}

	}
  }
  time_taken+=(double)(clock()-tStart)/CLOCKS_PER_SEC;
	// #############INTIALIZE#############################################################################
   // clock_t tStart = clock();
  // free(EXY);
    free(EX_sq);
    
  	std::vector<Valid> seeds;
	double * seed = (double * )malloc( sizeof(double)*g[3]);
	double * Corr = (double * )malloc( sizeof(double)*valid.size());
	
	
	Valid tv;
	tv.x = seedx;
	tv.y = seedy;
	tv.z = seedz;
	seeds.push_back(tv);

	int tempi = 0;		

	for (int t = 0; t < g[3]; ++t){

		if(t == 0)
			seed[t]= (image_data[((t*g[2] + seeds[tempi].z)*g[1]+seeds[tempi].y)*g[0]+seeds[tempi].x]);
		else
			seed[t]+= (image_data[((t*g[2] + seeds[tempi].z)*g[1]+seeds[tempi].y)*g[0]+seeds[tempi].x]);

	}


	g_copy.dim[0] = g[0];
	g_copy.dim[1] = g[1];
	g_copy.dim[2] = g[2];
	g_copy.dim[3] = 1;
	image_data_cpy.resize(g_copy);

	#pragma omp parallel for collapse(3)
	for (int x = 0; x < g[0]; ++x)
	  for (int y = 0; y < g[1]; ++y)
		 for (int z = 0; z < g[2]; ++z)
			image_data_cpy[((z)*g[1]+y)*g[0]+x] = 0;

	// #############################FINDING MATRIX FOR CORRELATION OF A SEED#####################################

	// clock_t tStart = clock();

	if(EX[(seedx*g[1] +seedy)*g[2] + seedz] == 0){
	  
	   for (int i = 0; i < valid.size(); ++i){
			for (int t = 0; t < g[3]; ++t){
		
		image_data_cpy[((valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = 0;

			}
	  
		}
	}
	else{
	// vector matrix multiplication
	// free(EXY);
	free(EX);
	int m = valid.size();
	int dime = g[3];
	float timeseries = 1/(float)dime;
	// std::cout<<m<<std::endl;
	// cblas_dgemv(CblasRowMajor,CblasNoTrans,m,dime,timeseries,Valid_matrix,dime,seed,1,0.0,Corr,1);
	#ifdef BUILD2
	  	std::cout<<"DEPRECATED RUN ON CEREBRUM....."<<std::endl;
	  	exit(0);
	  #else
	    cblas_dgemv(CblasRowMajor,CblasNoTrans,m,dime,timeseries,Valid_matrix,dime,seed,1,0.0,Corr,1);
	  #endif


	for (int i = 0; i < valid.size(); ++i)
		image_data_cpy[((valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = Corr[i];

  

	}

	//image::io::nifti nifti_parser2;
	nifti_parser2.load_from_image(image_data_cpy);
	std::string filename("output.nii");
	filename = ofname + filename;
	nifti_parser2.save_to_file(filename.c_str());
	// std::cout<< ((double)(clock() - tStart)/CLOCKS_PER_SEC)<<",";
}

void all_pair_corr(){
 
	//DEFINING VARIABLES  
  image::io::nifti nifti_parser,nifti_parser2;
  image::basic_image<double,4> image_data;
  image::basic_image<double,4> image_data_cpy;
  image::basic_image<int,3> mask_image;
  std::vector<Valid> valid;

  image_data_cpy = image_data;
  std::string t_ipfilename,t_maskfilename;
  t_ipfilename = ipfilename;
  t_maskfilename = maskfilename;

  if(ipfilename.find(".gz")!=std::string::npos){

  		t_ipfilename = t_ipfilename.substr(0,(t_ipfilename.length()-3));

  }
  if(mask && maskfilename.find(".gz")!=std::string::npos){

  		t_maskfilename = t_maskfilename.substr(0,(t_maskfilename.length()-3));

  }

	//LOADING THE IMAGE
  if(nifti_parser.load_from_file(t_ipfilename));
	  nifti_parser >> image_data;

	// LOADING THE MASK
  if(mask&&nifti_parser2.load_from_file(t_maskfilename))
  	   nifti_parser2 >> mask_image;
 
 if(ipfilename.find(".gz")!=std::string::npos){

		  std::string cmd = "rm ";
  		  cmd+=t_ipfilename;
  		  system(cmd.c_str());

  		

  }
  if(mask && maskfilename.find(".gz")!=std::string::npos){

  	  	  std::string cmd = "rm ";
  		  cmd+=t_maskfilename;
  		  //system(cmd.c_str());

  }



  image::geometry<4> g = image_data.geometry();
  image::geometry<4> g_copy;
  image::geometry<3> g_mask;

	//CHECK IF THE GEOMETRY OF THE MASK AND IMAGE IS SAME
  if(mask){
  	g_mask = mask_image.geometry();
  	if(g_mask[0]!=g[0]||g_mask[1]!=g[1]||g_mask[2]!=g[2]){
  		std::cout<<":::: ERROR INVALID MASK ::::"<<std::endl;
  		return;
  	}
  }

  // std::string cmd = "gzip ";
  // cmd+=ipfilename;
  // system(cmd.c_str());

	//############### DEFINING VARIABLES ####################################

  double *EX = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  double *EX_sq = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  std::string filename2("/voxelmap");
  filename2 = ofname + filename2;
  const char *s1 = filename2.c_str();
  std::ofstream f(s1,std::ios::binary);

	//#CALCULATING THE MEAN AND STD DEVIATION OF DATA AND CREATING A VALID VOXEL MAPPINGS ################################
  clock_t tStart ;
  for (int z = 0; z < g[2]; ++z){
	for (int y = 0; y < g[1]; ++y)
	  {
		for (int x = 0; x < g[0]; ++x)
		{
		  tStart = clock();
		  double temp  = 0,temp2 = 0;
		  for (int t = 0; t < g[3]; ++t)
		  {
			temp += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3] ;
			temp2 += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]*(double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3];
		  }
		  if (measure)
		  {
		  	time_taken+=(double)(clock()-tStart)/CLOCKS_PER_SEC;	
		  }
		  

		  if(!mask){
			  if(temp > max(20,thresh)){
				Valid tv ;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				int temp = valid.size();
				f.write(reinterpret_cast<char *>(&temp),sizeof(int));
				valid.push_back(tv);
			  }else{
				int temp = -1;
				f.write(reinterpret_cast<char *>(&temp),sizeof(int));
				// direc.insert( ( ( z*g[1] + y )*g[0] + x ) ,-1);
			  }
			}else{
				if(mask_image[(z*g[1]+y)*g[0]+x]==1){
				Valid tv ;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				int temp = valid.size();
				f.write(reinterpret_cast<char *>(&temp),sizeof(int));
				valid.push_back(tv);
			  }else{
				int temp = -1;
				f.write(reinterpret_cast<char *>(&temp),sizeof(int));
				
			  }


			}
		  
		  EX[ (x*g[1] +y)*g[2] + z] = temp;
		  EX_sq[ (x*g[1] +y)*g[2] + z] = sqrt(temp2- temp*temp);

		}
	  } 

  }

  f.close();

	//################## MAKING A MATRIX OF NORMALIZED DATA###########################################
  double * Valid_matrix = (double * )malloc( sizeof(double)*valid.size()*g[3]);
  tStart = clock();
  # pragma omp parallel for
  for (int i = 0; i < valid.size(); ++i)
  {
	for (int t = 0; t < g[3]; ++t){
		double currentdev = EX_sq[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z];
		double tempNormdata;
		if (currentdev == 0){
			tempNormdata = 0;
			image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
			Valid_matrix[i*g[3]+t] = tempNormdata;
		}
		else{
			tempNormdata = (double)(image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] - EX[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z])/currentdev;
	  		image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
	  		Valid_matrix[i*g[3]+t] = tempNormdata;
		}

	}
  }
  if (measure)
  {
  	time_taken+=(double)(clock()-tStart)/CLOCKS_PER_SEC;	
  }
  
	free(EX_sq);
	// clock_t tStart = clock();
	// free(EXY);


	// #######FINDING MATRIX FOR CORRELATION OF COMBINATION########################################
	free(EX);
	int dime = g[3];
	int n =  8192;
	double * res = (double * )malloc( sizeof(double)*n*n);
	double * Valid_matrix_chunk= (double * )malloc( sizeof(double)*g[3]*n);
	double * Valid_matrix_chunk_trans= (double * )malloc( sizeof(double)*n*g[3]);
  if(res==NULL){
	  free(Valid_matrix);
	  free(Valid_matrix_chunk);
	  free(res);
	  std::cout<<"res Allocating failed"<<std::endl;
	  exit(0);
	}

  if(Valid_matrix_chunk==NULL){
	  free(Valid_matrix);
	  free(Valid_matrix_chunk);
	  free(res);
	  std::cout<<"valid Allocating failed"<<std::endl;
	  exit(0);
	}
   
	double timeTk =0;
	//clock_t tStart = clock();
  for(int starti = 0;starti<valid.size();starti+=n)
	for(int startj = starti;startj<valid.size();startj+=n)
	{


	  // tStart = clock();
	 
	  // std::cout<<starti<<","<<startj<<std::endl;
	  int n2 = min(n,valid.size()-starti);
	  int n3 = min(n,valid.size()-startj);
	  int n_toget=max(n2,n3);

	  #pragma omp parallel for collapse(2)
	  for (int i = 0; i < n_toget; ++i)
	  {
		for (int j = 0; j < g[3]; ++j)
		{
		  //######### AS BOTH THE MATRICES ARE OF DIFFERENT SIZES SO SELECT ACCORDINGLY######################################### 
		  if(i<n2)
		  Valid_matrix_chunk[i*g[3] + j]=Valid_matrix[(starti+i)*g[3] + (j)];
		  if(i<n3)
		  Valid_matrix_chunk_trans[j*n + i]=Valid_matrix[(startj+i)*g[3] + (j)];
		}
	  }


	  double timeseries = 1/float(g[3]);

	  tStart = clock();
	  // cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ROI_MAX,valid_size,dime,timeseries,roi_avg,dime,Valid_matrix,valid_size,0.0,res,valid_size);
	  #ifdef BUILD2
	  	 char trans = 'T';
	  	 double gamma = 0.0;
	  	dgemm_(&trans,&trans,&n,&n,&dime,&timeseries,Valid_matrix_chunk,&dime,Valid_matrix_chunk_trans,&n,&gamma,res,&n);
	  #else
	    cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,n,n,dime,timeseries,Valid_matrix_chunk,dime,Valid_matrix_chunk_trans,n,0.0,res,n);
	  #endif

	  // cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n2,n3,dime,timeseries,Valid_matrix_chunk,dime,Valid_matrix_chunk_trans,n3,0.0,res,n3);
	  // cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,n2,n3,dime,timeseries,Valid_matrix_chunk,n2,Valid_matrix_chunk_trans,n3,0.0,res,n2);
	  // cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,n,n,dime,timeseries,Valid_matrix_chunk,dime,Valid_matrix_chunk_trans,n,0.0,res,n); FINALLY
	  time_taken+=(double)(clock()-tStart)/CLOCKS_PER_SEC;
	  //std::cout<<"Time taken:  for CORRELATION "<< ((double)(clock() - tStart)/CLOCKS_PER_SEC)<<std::endl;
	  //timeTk += (double)(clock() - tStart)/CLOCKS_PER_SEC ;
	  # pragma omp parallel for
	  for (int resi = 0; resi < n2; ++resi)
	  {
		// double *arr = (double *)malloc(sizeof(double)*n3);
		std::string filename("/voxel/voxel_");
		filename = ofname + filename;
		filename+= std::to_string(starti+resi);
		filename+=type;
		//direc.insert( ( ( valid[starti + resi].z*g[1] + valid[starti + resi].y )*g[0] + valid[starti + resi]x ) ,starti+resi);
		const char *s1 = filename.c_str();
		//std::ofstream f(s1,std::ios::binary|std::ios::app);
		std::ofstream f(s1,std::ios::binary|std::ios::app);
		for (int resj = 0; resj < n3; ++resj){
		   // f<< res[resi*n+resj]<<'\n';
		   // double tmp = res[resi*n+resj];
			double tmp = res[resj*n+resi];
		  f.write(reinterpret_cast<char *>(&tmp),sizeof(double));

		}
		  //f.write(reinterpret_cast<char *>(&arr),sizeof(double)*n3);
		  //fprintf(f, "%f ",res[resi*n+resj]);

		f.close();
		// free(arr);

	  }
	if(starti!=startj){
		#pragma omp parallel for
		for (int resj = 0; resj < n3; ++resj)
		{
		  //double *arr = (double *)malloc(sizeof(double)*n3);
		  std::string filename("/voxel/voxel_");
		  filename = ofname + filename;
		  filename+= std::to_string(startj+resj);
		  const char *s1 = filename.c_str();
		  filename+=type;
		  std::ofstream f(s1,std::ios::binary|std::ios::app);
		  for (int resi = 0; resi < n2; ++resi){
		  	// double tmp = res[resi*n+resj];
		  	double tmp = res[resj*n+resi];
		  	f.write(reinterpret_cast<char *>(&tmp),sizeof(double));
		  }
			// f << res[resi*n+resj]<<'\n';
			//f.write(reinterpret_cast<char *>(&arr),sizeof(double)*n3);
			// f.write((char*)&res[resi*n+resj],sizeof(double));
			// fprintf(f, "%f ",res[resi*n+resj]);
		  f.close();
		  //free(arr);
		}  
	  }

	}


	long long valid_size = valid.size();
	no_of_oper = 5*(g[0]*g[1]*g[2]*g[3])+(g[0]*g[1]*g[2]) + valid_size*valid_size*2*g[3];
	//std::cout<<"Time taken:  for CORRELATION "<< timeTk<<std::endl;
	std::cout<<"NO OF OPER :"<<no_of_oper<<std::endl;
	std::cout<<"time taken :"<<time_taken<<std::endl;

}


void avg_roi_time_corr(){
	//DEFINING VARIABLES  
  image::io::nifti nifti_parser,nifti_parser2;
  image::basic_image<double,4> image_data;
  image::basic_image<double,4> image_data_cpy;
  image::basic_image<int,3> roi_image;
  std::vector<Valid> valid;

  image_data_cpy = image_data;
  std::string t_ipfilename,t_maskfilename,t_roifname;
  t_ipfilename = ipfilename;
  t_maskfilename = maskfilename;
  t_roifname = roifname;

  if(ipfilename.find(".gz")!=std::string::npos){

  		t_ipfilename = t_ipfilename.substr(0,(t_ipfilename.length()-3));

  }
  if(roifname.find(".gz")!=std::string::npos){

  		t_roifname = t_roifname.substr(0,(t_roifname.length()-3));

  }
  if(mask && maskfilename.find(".gz")!=std::string::npos){

  		t_maskfilename = t_maskfilename.substr(0,(t_maskfilename.length()-3));

  }

	//LOADING THE IMAGE
  if(nifti_parser.load_from_file(t_ipfilename));
	  nifti_parser >> image_data;
 
   if(ipfilename.find(".gz")!=std::string::npos){

			  std::string cmd = "rm ";
	  		  cmd+=t_ipfilename;
	  		  system(cmd.c_str());

	  		

	  }

	// LOADING THE ROI
  if(nifti_parser2.load_from_file(t_roifname))
  	   nifti_parser2 >> roi_image;

  	if(roifname.find(".gz")!=std::string::npos){

			  std::string cmd = "rm ";
	  		  cmd+=t_roifname;
	  		  //system(cmd.c_str());

	  		

	  }
	// std::cout<<roifname<<std::endl;  
  image::geometry<4> g = image_data.geometry();
  //image::geometry<4> g_copy;
  image::geometry<3> g_roi;

	//CHECK IF THE GEOMETRY OF THE MASK AND IMAGE IS SAME
   	g_roi = roi_image.geometry();
   	

  	if(g_roi[0]!=g[0]||g_roi[1]!=g[1]||g_roi[2]!=g[2]){
  		std::cout<<":::: ERROR INVALID MASK ::::"<<std::endl;
  		return;
  	}
  image::basic_image<int,3> mask_image;
	if(mask&&nifti_parser.load_from_file(t_maskfilename))
	  	   nifti_parser >> mask_image;
	image::geometry<3> g_mask;

	 
	  if(mask && maskfilename.find(".gz")!=std::string::npos){

	  	  	  std::string cmd = "rm ";
	  		  cmd+=t_maskfilename;
	  		  //system(cmd.c_str());

	  }



		//CHECK IF THE GEOMETRY OF THE MASK AND IMAGE IS SAME
	  if(mask){
	  	g_mask = mask_image.geometry();
  	
	  	if(g_mask[0]!=g[0]||g_mask[1]!=g[1]||g_mask[2]!=g[2]){
	  		std::cout<<":::: ERROR INVALID MASK ::::"<<std::endl;
	  		return;
	  	}
	  }


  // std::string cmd = "gzip ";
  // cmd+=ipfilename;
  // system(cmd.c_str());

	//############### DEFINING VARIABLES ####################################

  double *EX = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  double *EX_sq = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  
	//#CALCULATING THE MEAN AND STD DEVIATION OF DATA AND CREATING A VALID VOXEL MAPPINGS ################################
  clock_t tStart;
  for (int z = 0; z < g[2]; ++z){
	for (int y = 0; y < g[1]; ++y)
	  {
		for (int x = 0; x < g[0]; ++x)
		{

		  double temp  = 0,temp2 = 0;
		  tStart = clock();
		  for (int t = 0; t < g[3]; ++t)
		  {
			temp += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3] ;
			temp2 += image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]*image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3];
		  }
		  if (measure)
		  {
		  	time_taken += (double)(clock()-tStart)/CLOCKS_PER_SEC;	
		  }
		  
			
		  // else{
				// int temp = -1;
				
				
		  // }
				  
		  EX[ (x*g[1] +y)*g[2] + z] = temp;
		  EX_sq[ (x*g[1] +y)*g[2] + z] = sqrt(temp2- temp*temp);

		  if((mask&&mask_image[(z*g[1]+y)*g[0]+x]==1)||(!mask&&EX_sq[ (x*g[1] +y)*g[2] + z] > 1e-2)){
				Valid tv ;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				// int temp = valid.size();
				valid.push_back(tv);
		  }

		}
	  } 

  }
std::cout<<"time taken 1:"<<std::endl;


	//################## MAKING A MATRIX OF NORMALIZED DATA###########################################
  long long valid_size = valid.size();
  double * Valid_matrix = (double * )malloc( sizeof(double)*valid_size*g[3]);
  // double * roi_result = (double * )malloc( sizeof(double)*ROI_MAX*valid_size);
  double * roi_avg = (double * ) calloc(ROI_MAX*g[3], sizeof(double));
  double * roi_var = (double * ) calloc(ROI_MAX, sizeof(double));
  double * roi_tot = (double * ) calloc(ROI_MAX, sizeof(double));
  double * roi_mean = (double * ) calloc(ROI_MAX, sizeof(double));	
  
  // std::cout<<"val "<<roi_avg[0]<<std::endl;
  tStart = clock();
  #pragma omp parallel for shared(roi_tot,roi_avg)
  for (int i = 0; i < valid.size(); ++i)
  {	
  	
  	int roi_no = roi_image[(valid[i].z*g[1] +valid[i].y)*g[0] + valid[i].x];
  	if(roi_no!=0)
  		roi_tot[roi_no-1]++;
	double currentdev = EX_sq[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z];
	for (int t = 0; t < g[3]; ++t){
		
		double tempNormdata;
		if (currentdev == 0){
			tempNormdata = 0;
			image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
			Valid_matrix[t*valid_size + i] = tempNormdata;

		}
		else{
			tempNormdata = (double)(image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] - EX[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z])/currentdev;
	  		if(roi_no!=0)
  				roi_avg[(roi_no-1)*g[3]+t ] += tempNormdata;
	  		image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
	  		Valid_matrix[t*valid_size + i] = tempNormdata;
	  		
		}

	}

	

  }
  if (measure)
  {
  	time_taken += (double)(clock()-tStart)/CLOCKS_PER_SEC;
  	std::cout<<"time taken 2:"<<time_taken<<std::endl;	
  }
  
  
  tStart = clock();

  // #pragma omp parallel for
  for (int i = 0; i < ROI_MAX; ++i)
  {	
  		double temp  = 0,temp2 = 0;
  		// std::cout<<i<<":"<<roi_tot[i]<<std::endl;
  		for (int t = 0; t < g[3]; ++t){
  			roi_avg[i*g[3]+t]/=roi_tot[i];
  	// 		temp += (double)roi_avg[i*g[3]+t]/g[3] ;
			// temp2 += roi_avg[i*g[3]+t]*roi_avg[i*g[3]+t]/g[3];
  		}

  		// roi_mean[i]= temp;
  		// roi_var[i] = sqrt(temp2- temp*temp);

  }
  if (measure)
  {

  time_taken += (double)(clock()-tStart)/CLOCKS_PER_SEC;
  std::cout<<"time taken 3:"<<time_taken<<std::endl;	
  }
  // #pragma omp parallel for
  // for (int i = 0; i < ROI_MAX; ++i){	
  // 		if(roi_var[i]!=0){
  // 			for (int t = 0; t < g[3]; ++t){
  // 				roi_avg[i*g[3]+t]=(roi_avg[i*g[3]+t]-roi_mean[i])/roi_var[i];
  // 				}
  // 		}else{
  // 			for (int t = 0; t < g[3]; ++t){
  // 				roi_avg[i*g[3]+t] = 0;
  // 			}

  // 		}	
  // }


  free(EX_sq);
  // free(roi_var);
  // free(roi_tot);
	// clock_t tStart = clock();
	// free(EXY);


	// #######FINDING MATRIX FOR CORRELATION OF COMBINATION########################################
	free(EX);
	int dime = g[3];
	int n =  8192;
	// double * roi_chunk = (double * )malloc( sizeof(double)*ROI_MAX*n);
	double * res = (double * )malloc( sizeof(double)*ROI_MAX*valid_size);

		image::geometry<4> g_copy;
		g_copy.dim[0] = g[0];
		g_copy.dim[1] = g[1];
		g_copy.dim[2] = g[2];
		g_copy.dim[3] = ROI_MAX;
		image_data_cpy.resize(g_copy);

	#pragma omp parallel for collapse(3)
		for (int x = 0; x < g[0]; ++x)
		  for (int y = 0; y < g[1]; ++y)
			 for (int z = 0; z < g[2]; ++z)
			 	for(int r = 0;r<ROI_MAX;r++)
					image_data_cpy[((r*g[2]+z)*g[1]+y)*g[0]+x] = 0;


	  double timeseries = 1/float(g[3]);

	  tStart = clock();

	   //cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ROI_MAX,valid_size,dime,timeseries,roi_avg,dime,Valid_matrix,valid_size,0.0,res,valid_size);

	  // cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ROI_MAX,valid_size,dime,timeseries,roi_avg,dime,Valid_matrix,valid_size,0.0,res,valid_size);
	  #ifdef BUILD2
	  	char trans = 'T';
	  	double gamma = 0.0;
	  	int valid_size_int = (int) valid_size;
	  	dgemm(&trans,&trans,&ROI_MAX,&valid_size_int,&dime,&timeseries,roi_avg,&ROI_MAX,Valid_matrix,&valid_size_int,&gamma,res,&ROI_MAX);
	  #else
	  cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,ROI_MAX,valid_size,dime,timeseries,roi_avg,ROI_MAX,Valid_matrix,valid_size,0.0,res,ROI_MAX);
	  // cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ROI_MAX,valid_size,dime,timeseries,roi_avg,dime,Valid_matrix,valid_size,0.0,res,valid_size);
	  #endif
	  time_taken += (double)(clock()-tStart)/CLOCKS_PER_SEC;
	  std::cout<<"time taken :"<<time_taken<<std::endl;
	#pragma omp parallel for
	for(int r = 0;r<ROI_MAX;r++)
		for (int i = 0; i < valid_size; ++i)
		 	image_data_cpy[((r*g[2]+valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = res[i*ROI_MAX+r];


	nifti_parser2.load_from_image(image_data_cpy);
	std::string filename("avg_roi_time_series.nii");
	filename = ofname+"/"+ filename;
	nifti_parser2.save_to_file(filename.c_str());

	if (measure)
	{
		no_of_oper = 2*ROI_MAX*g[3]*valid_size+5*(g[0]*g[1]*g[2]*g[3])+(g[0]*g[1]*g[2]);	
	}
	else
	{
		no_of_oper = 2*ROI_MAX*g[3]*valid_size;
	}
	
	std::cout<<"ROI_MAX, g[3], valid_size "<<ROI_MAX<<" "<<g[3]<<" "<<valid_size<<std::endl;
	std::cout<<"no_of_oper :"<<no_of_oper<<std::endl;
	std::cout<<"time taken :"<<time_taken<<std::endl;

}

void avg_corr_roi(){


 //DEFINING VARIABLES  
  image::io::nifti nifti_parser,nifti_parser2;
  image::basic_image<double,4> image_data;
  image::basic_image<double,4> image_data_cpy;
  image::basic_image<int,3> roi_image;
  std::vector<Valid> valid;

  image_data_cpy = image_data;
  std::string t_ipfilename,t_maskfilename,t_roifname;
  t_ipfilename = ipfilename;
  t_maskfilename = maskfilename;
  t_roifname = roifname;

  if(ipfilename.find(".gz")!=std::string::npos){

  		t_ipfilename = t_ipfilename.substr(0,(t_ipfilename.length()-3));

  }
  if(roifname.find(".gz")!=std::string::npos){

  		t_roifname = t_roifname.substr(0,(t_roifname.length()-3));

  }
  if(mask && maskfilename.find(".gz")!=std::string::npos){

  		t_maskfilename = t_maskfilename.substr(0,(t_maskfilename.length()-3));

  }

	//LOADING THE IMAGE
  if(nifti_parser.load_from_file(t_ipfilename));
	  nifti_parser >> image_data;
 
   if(ipfilename.find(".gz")!=std::string::npos){

			  std::string cmd = "rm ";
	  		  cmd+=t_ipfilename;
	  		  system(cmd.c_str());

	  		

	  }

	// LOADING THE ROI
  if(nifti_parser2.load_from_file(t_roifname))
  	   nifti_parser2 >> roi_image;

  	if(roifname.find(".gz")!=std::string::npos){

			  std::string cmd = "rm ";
	  		  cmd+=t_roifname;
	  		  //system(cmd.c_str());

	  		

	  }
	// std::cout<<roifname<<std::endl;  
  image::geometry<4> g = image_data.geometry();
  //image::geometry<4> g_copy;
  image::geometry<3> g_roi;

	//CHECK IF THE GEOMETRY OF THE MASK AND IMAGE IS SAME
   	g_roi = roi_image.geometry();
   	

  	if(g_roi[0]!=g[0]||g_roi[1]!=g[1]||g_roi[2]!=g[2]){
  		std::cout<<":::: ERROR INVALID MASK ::::"<<std::endl;
  		return;
  	}
  image::basic_image<int,3> mask_image;
	if(mask&&nifti_parser.load_from_file(t_maskfilename))
	  	   nifti_parser >> mask_image;
	image::geometry<3> g_mask;

	 
	  if(mask && maskfilename.find(".gz")!=std::string::npos){

	  	  	  std::string cmd = "rm ";
	  		  cmd+=t_maskfilename;
	  		  //system(cmd.c_str());

	  }



		//CHECK IF THE GEOMETRY OF THE MASK AND IMAGE IS SAME
	  if(mask){
	  	g_mask = mask_image.geometry();
  	
	  	if(g_mask[0]!=g[0]||g_mask[1]!=g[1]||g_mask[2]!=g[2]){
	  		std::cout<<":::: ERROR INVALID MASK ::::"<<std::endl;
	  		return;
	  	}
	  }


  // std::string cmd = "gzip ";
  // cmd+=ipfilename;
  // system(cmd.c_str());

	//############### DEFINING VARIABLES ####################################

  double *EX = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  double *EX_sq = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  
	//#CALCULATING THE MEAN AND STD DEVIATION OF DATA AND CREATING A VALID VOXEL MAPPINGS ################################
  clock_t tStart;
  for (int z = 0; z < g[2]; ++z){
	for (int y = 0; y < g[1]; ++y)
	  {
		for (int x = 0; x < g[0]; ++x)
		{

		  double temp  = 0,temp2 = 0;
		  tStart = clock();
		  for (int t = 0; t < g[3]; ++t)
		  {
			temp += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3] ;
			temp2 += image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]*image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3];
		  }
		  if (measure)
		  {
		  	time_taken += (double)(clock()-tStart)/CLOCKS_PER_SEC;	
		  }
		  
			
		  // else{
				// int temp = -1;
				
				
		  // }
				  
		  EX[ (x*g[1] +y)*g[2] + z] = temp;
		  EX_sq[ (x*g[1] +y)*g[2] + z] = sqrt(temp2- temp*temp);

		  if((mask&&mask_image[(z*g[1]+y)*g[0]+x]==1)||(!mask&&EX_sq[ (x*g[1] +y)*g[2] + z] > 1e-2)){
				Valid tv ;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				// int temp = valid.size();
				valid.push_back(tv);
		  }

		}
	  } 

  }



	//################## MAKING A MATRIX OF NORMALIZED DATA###########################################
  int valid_size = valid.size();
  double * Valid_matrix = (double * )malloc( sizeof(double)*valid_size*g[3]);
  // double * roi_result = (double * )malloc( sizeof(double)*ROI_MAX*valid_size);
  double * roi_avg = (double * ) calloc(ROI_MAX*g[3], sizeof(double));
  double * roi_var = (double * ) calloc(ROI_MAX, sizeof(double));
  double * roi_tot = (double * ) calloc(ROI_MAX, sizeof(double));
  double * roi_mean = (double * ) calloc(ROI_MAX, sizeof(double));	
  
  // std::cout<<"val "<<roi_avg[0]<<std::endl;
  tStart = clock();
  #pragma omp parallel for shared(roi_tot,roi_avg)
  for (int i = 0; i < valid.size(); ++i)
  {	
  	
  	int roi_no = roi_image[(valid[i].z*g[1] +valid[i].y)*g[0] + valid[i].x];
  	if(roi_no!=0)
  		roi_tot[roi_no-1]++;
	double currentdev = EX_sq[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z];
	for (int t = 0; t < g[3]; ++t){
		
		double tempNormdata;
		if (currentdev == 0){
			tempNormdata = 0;
			image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
			Valid_matrix[t*valid_size + i] = tempNormdata;

		}
		else{
			tempNormdata = (double)(image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] - EX[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z])/currentdev;
	  		if(roi_no!=0)
  				roi_avg[(roi_no-1)*g[3]+t ] += image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x];
	  		image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
	  		Valid_matrix[t*valid_size + i] = tempNormdata;
	  		
		}

	}

	

  }
  if (measure)
  {

  time_taken += (double)(clock()-tStart)/CLOCKS_PER_SEC;
  std::cout<<"time taken :"<<time_taken<<std::endl;	
  }
  tStart = clock();

  // #pragma omp parallel for
  for (int i = 0; i < ROI_MAX; ++i)
  {	
  		double temp  = 0,temp2 = 0;
  		// std::cout<<i<<":"<<roi_tot[i]<<std::endl;
  		for (int t = 0; t < g[3]; ++t){
  			roi_avg[i*g[3]+t]/=roi_tot[i];
  			temp += (double)roi_avg[i*g[3]+t]/g[3] ;
			temp2 += roi_avg[i*g[3]+t]*roi_avg[i*g[3]+t]/g[3];
  		}

  		roi_mean[i]= temp;
  		roi_var[i] = sqrt(temp2- temp*temp);

  }
  if (measure)
  {

  time_taken += (double)(clock()-tStart)/CLOCKS_PER_SEC;
  std::cout<<"time taken :"<<time_taken<<std::endl;	
  }
  #pragma omp parallel for
  for (int i = 0; i < ROI_MAX; ++i){	
  		if(roi_var[i]!=0){
  			for (int t = 0; t < g[3]; ++t){
  				roi_avg[i*g[3]+t]=(roi_avg[i*g[3]+t]-roi_mean[i])/roi_var[i];
  				}
  		}else{
  			for (int t = 0; t < g[3]; ++t){
  				roi_avg[i*g[3]+t] = 0;
  			}

  		}	
  }


  // free(EX_sq);
  // free(roi_var);
  // free(roi_tot);
	// clock_t tStart = clock();
	// free(EXY);


	// #######FINDING MATRIX FOR CORRELATION OF COMBINATION########################################
	// free(EX);
	int dime = g[3];
	int n =  8192;
	// double * roi_chunk = (double * )malloc( sizeof(double)*ROI_MAX*n);
	double * res = (double * )malloc( sizeof(double)*ROI_MAX*valid_size);

		image::geometry<4> g_copy;
		g_copy.dim[0] = g[0];
		g_copy.dim[1] = g[1];
		g_copy.dim[2] = g[2];
		g_copy.dim[3] = ROI_MAX;
		image_data_cpy.resize(g_copy);

	#pragma omp parallel for collapse(3)
		for (int x = 0; x < g[0]; ++x)
		  for (int y = 0; y < g[1]; ++y)
			 for (int z = 0; z < g[2]; ++z)
			 	for(int r = 0;r<ROI_MAX;r++)
					image_data_cpy[((r*g[2]+z)*g[1]+y)*g[0]+x] = 0;


	  double timeseries = 1/float(g[3]);

	  tStart = clock();
	  //cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ROI_MAX,valid_size,dime,timeseries,roi_avg,dime,Valid_matrix,valid_size,0.0,res,valid_size);
	  #ifdef BUILD2
	  	char trans = 'T';
	  	double gamma = 0.0;
	  	// std::cout<<"Start of dgemm"<<std::endl;
	  	dgemm(&trans,&trans,&ROI_MAX,&valid_size,&dime,&timeseries,roi_avg,&ROI_MAX,Valid_matrix,&valid_size,&gamma,res,&ROI_MAX);
	 	// std::cout<<"End of dgemm"<<std::endl; 
	  #else
	  cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,ROI_MAX,valid_size,dime,timeseries,roi_avg,ROI_MAX,Valid_matrix,valid_size,0.0,res,ROI_MAX);
	  #endif
	  time_taken += (double)(clock()-tStart)/CLOCKS_PER_SEC;
	  std::cout<<"time taken :"<<time_taken<<std::endl;
	#pragma omp parallel for
	for(int r = 0;r<ROI_MAX;r++)
		for (int i = 0; i < valid_size; ++i)
		 	image_data_cpy[((r*g[2]+valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = res[i*ROI_MAX + r];


	nifti_parser2.load_from_image(image_data_cpy);
	std::string filename("avg_roi_time_series.nii");
	filename = ofname+"/"+ filename;
	nifti_parser2.save_to_file(filename.c_str());

	if (measure)
	{
		no_of_oper = 2*ROI_MAX*g[3]*valid_size+5*(g[0]*g[1]*g[2]*g[3])+(g[0]*g[1]*g[2]);	
	}
	else
	{
		no_of_oper = 2*ROI_MAX*g[3]*valid_size;
	}
	
	std::cout<<"NO OF OPER :"<<no_of_oper<<std::endl;
	std::cout<<"time taken :"<<time_taken<<std::endl;




}

void correl()
{
  if(seedcmp==true){
  		seed_correl();
  }else if(roi1==true){
  		avg_roi_time_corr();
  		if(flops)
  			std::cout<<"FLOPS: "<<(double)(no_of_oper)/time_taken<<std::endl;
  }else if(roi2){
  		avg_corr_roi();
  		if(flops)
  			std::cout<<"FLOPS: "<<(double)(no_of_oper)/time_taken<<std::endl;
  }else if(all==true){
  		all_pair_corr();
  		if(flops)
  			std::cout<<"FLOPS: "<<(double)(no_of_oper)/time_taken<<std::endl;
  }
}


// #########################################get options #############################################################################
void getattributes(int argc,char *argv[])
{
  char tmp;
  int i =1;

  while(i<argc)
  {
  	char tmp = argv[i][1];
   // std::cout<<tmp<<std::endl;
	switch(tmp)
	{
	  /*option h show the help infomation*/
	  case 'h':
		showhelpinfo();
		exit(0);
		break;
	  /*option u pres ent the username*/
	  case 'i':
		input= true;
		// std::cout<<optarg<<std::endl;
		// strcpy(ipfilename,optarg);
		i++;
		ipfilename = argv[i];
		if(ipfilename.find(".gz")!=std::string::npos){
			try
			{
				gzip = true;
				std::string command = "gunzip -k ";
				command += ipfilename;
				system(command.c_str());
				// ipfilename = ipfilename.substr(0,(ipfilename.length()-3));
			}catch(...){

				std::cout<<"Gunzip problem caught"<<std::endl;
			}
		}

		// if(ipfilename.find(".dcm")!=std::string::npos)
		// 	dcm = true;

		break;
	  /*option p present the password*/ 
	  case 'r':
		roi1= true;
		i++;
		roifname = argv[i];
		if(roifname.find(".gz")!=std::string::npos){
			try{
				gzip = true;
				std::string command = "gunzip -k ";
				command += roifname;
				system(command.c_str());
				// roifname = roifname.substr(0,(roifname.length()-3));
			}catch(...){

				std::cout<<"Gunzip problem caught"<<std::endl;
			}
		}
		i++;
		ROI_MAX = std::stoi(argv[i]);
		
		break;
	  case 'R':
		roi2= true;
		i++;
		roifname = argv[i];
		if(roifname.find(".gz")!=std::string::npos){
			try{
				gzip = true;
				std::string command = "gunzip -k ";
				command += roifname;
				system(command.c_str());
				// roifname = roifname.substr(0,(roifname.length()-3));
			}catch(...){

				std::cout<<"Gunzip problem caught"<<std::endl;
			}
		}
		i++;
		ROI_MAX = std::stoi(argv[i]);
		break;
	  case 'o':
		output = true;
		//strcpy(ofname,optarg);
		i++;
		ofname = argv[i];
		break;
	  case 't':
		// std::cout<<"thresh"<<std::endl;
		// std::cout<<optarg<<std::endl;
	  	i++;
		thresh = std::stoi(argv[i]);
	  break;
	  case 's':
		 seedcmp = true;
		 i++;
		 seedx = std::stoi(argv[i]);
		 i++;
		 seedy = std::stoi(argv[i]);
		 i++;
		 seedz = std::stoi(argv[i]);
		break;
	  case 'a': all = true;
				// std::cout<<all<<std::endl;
				break;
	  /*invail input will get the heil infomation*/
	  case 'm': mask = true;
	  			i++;
	  			maskfilename = argv[i];
				if(maskfilename.find(".gz")!=std::string::npos){
					try{
						gzip = true;
						std::string command = "gunzip -k ";
						command += maskfilename;
						system(command.c_str());
						// maskfilename = maskfilename.substr(0,(maskfilename.length()-3));
					}catch(...){

						std::cout<<"Gunzip problem caught"<<std::endl;
					}
				}
	  			break;
	  case 'f': flops = true;
	  			
	  			break;

	  default:
				showhelpinfo();
				exit(0);
	  break;
	}

	i++;
  }
  // return 0;
}

/*################################funcion that show the help information################################################################*/
void showhelpinfo()
{
 std::cout<<"Usage :\n fconn -i <fmri.nii/nii.gz> -o <project name>  -[r <roi_filename> <N>/R <roi_filename> <N> /s <x> <y> <z>/a] -[f] -[t] -[h] \n";
  
  std::cout<<"Compulsory arguments(You have to specify the following)\n";
  std::cout<<"\t -i\t\t path of the input volume \n";
  std::cout<<"\t -o \t\t project name  \n ";
  std::cout<<"Only one of the arguments must be present \n";
  std::cout<<"\t -r <roi_filename> <N>\t\t path of the volume containing the desired ROI \n";
  std::cout<<"\t -R <roi_filename> <N> \t\t path of the volume containing the desired ROI \n";
  std::cout<<"\t -s x y z  \t\t for seed to all voxel mode(no argument) \n";
  std::cout<<"\t\t x \t\t x-coordinate for seed (compulsory in -s mode if -r option's not there)\n";
  std::cout<<"\t\t y \t\t y-coordinate for seed (compulsory in -s mode if -r option's not there)\n";
  std::cout<<"\t\t z \t\t z-coordinate for seed (compulspry in -s mode if -r option's not there)\n";
  std::cout<<"\t -a \t\t for all voxels to all voxels mode(no argument) \n";  
  std::cout<<"Optional arguments(You may optionally specify the following)\n";
  std::cout<<"\t -t \t an upper threshold  \n ";
  std::cout<<"\t -h \t display this message  \n ";
  std::cout<<"\t -f \t\t report FLOPS \n";
  // 
}
/* ###################################################### NON-OPTIMIZED CODE############################################ 

	//DEFINING VARIABLES  
  image::io::nifti nifti_parser,nifti_parser2;
  image::basic_image<double,4> image_data;
  image::basic_image<double,4> image_data_cpy;
  image::basic_image<int,3> roi_image;
  std::vector<Valid> valid;

  image_data_cpy = image_data;

	//LOADING THE IMAGE
  if(nifti_parser.load_from_file(ipfilename));
	  nifti_parser >> image_data;
  std::string cmd = "gzip ";
  cmd+=ipfilename;
  system(cmd.c_str());
	// LOADING THE ROI
  if(nifti_parser2.load_from_file(roifname))
  	   nifti_parser2 >> roi_image;
  
  image::geometry<4> g = image_data.geometry();
  image::geometry<4> g_copy;
  image::geometry<3> g_roi;

	//CHECK IF THE GEOMETRY OF THE MASK AND IMAGE IS SAME
   	g_roi = roi_image.geometry();
  	if(g_roi[0]!=g[0]||g_roi[1]!=g[1]||g_roi[2]!=g[2]){
  		std::cout<<":::: ERROR INVALID MASK ::::"<<std::endl;
  		return;
  	}
  cmd = "gzip ";
  cmd+=roifname;
  system(cmd.c_str());

   image::basic_image<int,3> mask_image;
	if(mask&&nifti_parser.load_from_file(maskfilename))
	  	   nifti_parser >> mask_image;
	image::geometry<3> g_mask;


		//CHECK IF THE GEOMETRY OF THE MASK AND IMAGE IS SAME
	  if(mask){
	  	g_mask = mask_image.geometry();
	  	cmd = "gzip ";
  		cmd+=maskfilename;
  		system(cmd.c_str());
	  	if(g_mask[0]!=g[0]||g_mask[1]!=g[1]||g_mask[2]!=g[2]){
	  		std::cout<<":::: ERROR INVALID MASK ::::"<<std::endl;
	  		return;
	  	}
	  }


	//############### DEFINING VARIABLES ####################################

  double *EX = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  double *EX_sq = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  
	//#CALCULATING THE MEAN AND STD DEVIATION OF DATA AND CREATING A VALID VOXEL MAPPINGS ################################
  for (int z = 0; z < g[2]; ++z){
	for (int y = 0; y < g[1]; ++y)
	  {
		for (int x = 0; x < g[0]; ++x)
		{

		  double temp  = 0,temp2 = 0;
		  for (int t = 0; t < g[3]; ++t)
		  {
			temp += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3] ;
			temp2 += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]*(double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3];
		  }
		  
			
				  
		  EX[ (x*g[1] +y)*g[2] + z] = temp;
		  EX_sq[ (x*g[1] +y)*g[2] + z] = sqrt(temp2- temp*temp);

		  if((mask&&mask_image[(z*g[1]+y)*g[0]+x]==1)||(!mask&&EX_sq[ (x*g[1] +y)*g[2] + z] > 1e-2)){
				Valid tv ;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				// int temp = valid.size();
				valid.push_back(tv);
		  }

		}
	  } 

  }

 // std::cout<<":::: mean std end ::::"<<std::endl;

	//################## MAKING A MATRIX OF NORMALIZED DATA###########################################
  int valid_size = valid.size();
  double * Valid_matrix = (double * )malloc( sizeof(double)*valid_size*g[3]);
  double * roi_result = (double * )malloc( sizeof(double)*ROI_MAX*valid_size);
  double * roi_result_temp = (double * )malloc( sizeof(double)*ROI_MAX*g[3]);
  double * roi_coeff = (double * )calloc( ROI_MAX*valid_size, sizeof(double));
  double * roi_avg = (double * )calloc(ROI_MAX*g[3], sizeof(double));
  double * roi_var = (double * ) calloc(ROI_MAX, sizeof(double));
  double * roi_tot = (double * ) calloc(ROI_MAX, sizeof(double));
	
 // std::cout<<":::: init start ::::"<<std::endl;
 # pragma omp parallel for shared(roi_avg,roi_tot)
  for (int i = 0; i < valid.size(); ++i)
  {	
  	int roi_no = roi_image[(valid[i].z*g[1] +valid[i].y)*g[2] + valid[i].x];
  	if(roi_no!=0)
  	  	roi_tot[roi_no-1]++;
  	double currentdev = EX_sq[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z];
	for (int t = 0; t < g[3]; ++t){
		
		double tempNormdata;
		if (currentdev == 0){
			tempNormdata = 0;
			image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
			Valid_matrix[i*g[3]+t] = tempNormdata;

		}
		else{
			tempNormdata = (double)(image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] - EX[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z])/currentdev;
	  		if(roi_no!=0)
	  			roi_avg[(roi_no-1)*g[3]+t ] += (double)(image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] - EX[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z]);
	  		image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = tempNormdata;
	  		Valid_matrix[i*g[3]+t] = tempNormdata;
	  		
		}

	}
	if(roi_no!=0)
		roi_coeff[(roi_no-1)*valid_size+i]=EX_sq[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z];

  }

  # pragma omp parallel for
  for (int i = 0; i < ROI_MAX; ++i)
  {	
  		double temp  = 0,temp2 = 0;
  		for (int t = 0; t < g[3]; ++t){
  			roi_avg[i*g[3]+t]/=roi_tot[i];
  			temp += (double)roi_avg[i*g[3]+t]/g[3] ;
			temp2 += (double)roi_avg[i*g[3]+t]*(double)roi_avg[i*g[3]+t]/g[3];
  		}

  		roi_var[i] = sqrt(temp2- temp*temp);

  }


  // # pragma omp parallel for
  for (int i = 0; i < valid.size(); ++i)
  {	
  	int roi_no = roi_image[(valid[i].z*g[1] +valid[i].y)*g[2] + valid[i].x];
  	if(roi_no==0)
  		continue;
  	roi_coeff[(roi_no-1)*valid_size+i]/=roi_var[roi_no];
  	roi_coeff[(roi_no-1)*valid_size+i]/=roi_tot[roi_no];
  }

   // std::cout<<":::: init end ::::"<<std::endl;

  free(EX_sq);
  // free(roi_var);
  // free(roi_tot);
  // free(roi_avg);
	// clock_t tStart = clock();
	// free(EXY);


	// #######FINDING MATRIX FOR CORRELATION OF COMBINATION########################################
	free(EX);
	int dime = g[3];
	float timeseries = 1/(float)dime;
	 //std::cout<<":::: blas start ::::"<<std::endl;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ROI_MAX,dime,valid_size,1.0,roi_coeff,valid_size,Valid_matrix,dime,0.0,roi_result_temp,dime);
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ROI_MAX,valid_size,dime,timeseries,roi_result_temp,dime,Valid_matrix,dime,0.0,roi_result,valid_size);

	// make result 
	g_copy.dim[0] = g[0];
	g_copy.dim[1] = g[1];
	g_copy.dim[2] = g[2];
	g_copy.dim[3] = ROI_MAX;
	image_data_cpy.resize(g_copy);

 // std::cout<<":::: output start ::::"<<std::endl;
	#pragma omp parallel for collapse(3)
	for (int x = 0; x < g[0]; ++x)
	  for (int y = 0; y < g[1]; ++y)
		 for (int z = 0; z < g[2]; ++z)
		 	for(int r = 0;r<ROI_MAX;r++)
				image_data_cpy[((r*g[2]+z)*g[1]+y)*g[0]+x] = 0;

	#pragma omp parallel for
	for(int r = 0;r<ROI_MAX;r++)
		for (int i = 0; i < valid_size; ++i)
		 	image_data_cpy[((r*g[2]+valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = roi_result[r*valid_size +i];

 	// save the result

	//image::io::nifti nifti_parser2;
	// cmd = "gzip ";
	// cmd+=ipfilename;
	// system(cmd.c_str());
	nifti_parser2.load_from_image(image_data_cpy);
	std::string filename("avg_roi_avg_time_series.nii");
	filename = ofname +"/"+ filename;
	nifti_parser2.save_to_file(filename.c_str());
	 //std::cout<<":::: output end  ::::"<<std::endl;
*/
