extern "C" {

#include <cblas.h>

}
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

std::string type(".bin"),ofname,ipfilename;
bool input = false,roi = false,output=false,gzip = false;
// dcm = false;
// char ipfilename[400]
char roifname[400],mask[400];
int thresh=0;
bool seedcmp = false,all=false;
int seedx ,seedy,seedz;


void showhelpinfo();
void correl();
void getattributes(int argc,char *argv[]);

int main (int argc,char *argv[])
{
  omp_set_num_threads(8);
  /*if the program is ran witout options ,it will show the usgage and exit*/
  if(argc < 2)
  {
	showhelpinfo();
	exit(1);
  }
  
  
  getattributes(argc,argv);
  
  if(input==false||output==false ||(roi == false&& all == false && seedcmp == false))
  {
	showhelpinfo();
	exit(1);
  }

  if( ( roi == true&& all == true ) || ( roi == true&& seedcmp == true ) || ( all == true && seedcmp == true ))
  {
	showhelpinfo();
	exit(1);
  }
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
  if(roi==false){
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




void correl()
{
  image::io::nifti nifti_parser;
  image::basic_image<double,4> image_data;
  image::basic_image<double,4> image_data_cpy;
  std::vector<Valid> valid;

  image_data_cpy = image_data;

  // if(dcm == true){
  // 	image::io::dicom dicom_parser;
  // 	if(dicom_parser.load_from_file(ipfilename))
  //   	dicom_parser >> image_data;
  // }
  // else{
  
  if(nifti_parser.load_from_file(ipfilename));
	  nifti_parser >> image_data;
  // }
  image::geometry<4> g = image_data.geometry();
  image::geometry<4> g_copy;

  std::string cmd = "gzip ";
  cmd+=ipfilename;
  system(cmd.c_str());

  //####################################### DEFINING VARIABLES ####################################

  double *EX = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  double *EX_sq = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  // double *EXY = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;



  // std::cout<<g[0]<<":"<<g[1]<<":"<<g[2]<<":"<<g[3]<<":"<<std::endl;

  std::string filename2("/voxelmap");
  filename2 = ofname + filename2;
  const char *s1 = filename2.c_str();
  std::ofstream f(s1,std::ios::binary);

  // #########################CALCULATING THE MEAN AND STD DEVIATION OF DATA AND CREATING A VALID VOXEL MAPPINGS #######################################################################################
  // for (int x = 0; x < g[0]; ++x){
  //   for (int y = 0; y < g[1]; ++y)
  //     {
  //       for (int z = 0; z < g[2]; ++z)
  //       {
  for (int z = 0; z < g[2]; ++z){
	for (int y = 0; y < g[1]; ++y)
	  {
		for (int x = 0; x < g[0]; ++x)
		{

		  double temp  = 0,temp2 = 0;
		  // std::cout<<temp<<std::endl;
		  for (int t = 0; t < g[3]; ++t)
		  {
			temp += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3] ;
			temp2 += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]*image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3];
		  }
		  // std::cout<<temp<<std::endl;
		  if(temp > max(20,thresh)){
			Valid tv ;
			tv.x = x;
			tv.y = y;
			tv.z = z;
			int temp = valid.size();
			f.write(reinterpret_cast<char *>(&temp),sizeof(int));
			valid.push_back(tv);
			// std::cout<<temp<<std::endl;
			// std::cout<<temp2<<std::endl;

			// std::cout<<x<<" "<<y<<" "<<z<<std::endl;
		  }else{
			int temp = -1;
			f.write(reinterpret_cast<char *>(&temp),sizeof(int));
			// direc.insert( ( ( z*g[1] + y )*g[0] + x ) ,-1);
		  }
		  
		  EX[ (x*g[1] +y)*g[2] + z] = temp;
		  EX_sq[ (x*g[1] +y)*g[2] + z] = sqrt(temp2);
		  // image_data_cpy[((z)*g[1]+y)*g[0]+x] = 0;
		}
	  } 

  }

  f.close();

  // ####################################################END##################################################################




  // ########################################### MAKING A MATRIX OF NORMALIZED DATA##########################################################################################################################
  double * Valid_matrix = (double * )malloc( sizeof(double)*valid.size()*g[3]);
  # pragma omp parallel for
  for (int i = 0; i < valid.size(); ++i)
  {
	for (int t = 0; t < g[3]; ++t){
	  image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] =(double)(image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] - EX[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z])/(EX_sq[(valid[i].x*g[1] +valid[i].y)*g[2] + valid[i].z]);
	  Valid_matrix[i*g[3]+t] = image_data[((t*g[2] + valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x];

	}
  }

  // ##########################################END######################################################################################################################################################################

   // clock_t tStart = clock();
  // free(EXY);
  free(EX_sq);
  if(seedcmp == true){
	// ########################################### MAKING A SEED VECTOR ##########################################################################################################################
  	int no_roi = 0;
  	std::vector<Valid> seeds;
	double * seed = (double * )malloc( sizeof(double)*g[3]);
	double * Corr = (double * )malloc( sizeof(double)*valid.size());
	if(roi == true){
		image::basic_image<double,3> image_data_cpy2;
		if(nifti_parser.load_from_file(roifname));
		nifti_parser >> image_data_cpy2;
		image::geometry<3> geo;
		geo = image_data_cpy2.geometry();



	 #pragma omp parallel for collapse(3) shared(no_roi)
	  for (int x = 0; x < geo[0]; ++x)
		for (int y = 0; y < geo[1]; ++y)
		   for (int z = 0; z < geo[2]; ++z)
			  if(image_data_cpy2[((z)*g[1]+y)*g[0]+x]>0){
			  	no_roi++;
			  	Valid tv;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				seeds.push_back(tv);

			  }
				

	}else{
		no_roi++;
		Valid tv;
		tv.x = seedx;
		tv.y = seedy;
		tv.z = seedz;
		seeds.push_back(tv);		

	}

	for (int tempi = 0;tempi<no_roi;tempi++)
		for (int t = 0; t < g[3]; ++t){

			if(tempi == 0)
				seed[t]= (image_data[((t*g[2] + seeds[tempi].z)*g[1]+seeds[tempi].y)*g[0]+seeds[tempi].x]);
			else
				seed[t]+= (image_data[((t*g[2] + seeds[tempi].z)*g[1]+seeds[tempi].y)*g[0]+seeds[tempi].x]);

		}

	for(int t =0;t<g[3] && no_roi!=1;t++){
		seed[t]/=no_roi;
	}

	// ##########################################END######################################################################################################################################################################

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

	// #############################FINDING MATRIX FOR CORRELATION OF A SEED######################################################################################################################

	clock_t tStart = clock();

	if(EX[(seedx*g[1] +seedy)*g[2] + seedz] == 0){
	  // return null values
	  // std::cout<<"NULL VALUES"<<std::endl;
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
	// std::cout<<m<<std::endl;
	cblas_dgemv(CblasRowMajor,CblasNoTrans,m,dime,1.0,Valid_matrix,dime,seed,1,0.0,Corr,1);
	 for (int i = 0; i < valid.size(); ++i){
		
		image_data_cpy[((valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = Corr[i];

		
	  }

	}

	image::io::nifti nifti_parser2;
	nifti_parser2.load_from_image(image_data_cpy);
	std::string filename("output.nii");
	filename = ofname + filename;
	nifti_parser2.save_to_file(ofname.c_str());
	//std::cout<<valid[0].x<<std::endl;
	//std::cout<<"nC2 begin"<<std::endl;
	 std::cout<< ((double)(clock() - tStart)/CLOCKS_PER_SEC)<<",";

	// #####################################################END#########################################################################################################################

  }else if(all == true){

	
	// #############################FINDING MATRIX FOR CORRELATION OF COMBINATION######################################################################################################################
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
	// std::cout<<valid.size()<<std::endl;
	 clock_t tStart = clock();
  // #pragma omp parallel for
  for(int starti = 0;starti<valid.size();starti+=n)
	for(int startj = starti;startj<valid.size();startj+=n)
	{


	  tStart = clock();
	 
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


	  // std::cout<<<<std::endl;
	  // std::cout<<"done something..."<<std::endl;
	  // std::cout<<"cblas_dgemm("<<CblasRowMajor<<","<<CblasNoTrans<<","<<CblasNoTrans<<","<<n<<","<<n<<","<<dime<<","<<1.0<<","<<"Valid_matrix"<<","<<dime<<","<<"Valid_matrix_trans,"<<n<<","<<0.0<<",res,"<<n<<")"<<std::endl;
	  
	  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n2,n3,dime,1.0,Valid_matrix_chunk,dime,Valid_matrix_chunk_trans,n3,0.0,res,n3);
	   //std::cout<<"Time taken:  for CORRELATION "<< ((double)(clock() - tStart)/CLOCKS_PER_SEC)<<std::endl;
	  timeTk += (double)(clock() - tStart)/CLOCKS_PER_SEC ;
	  # pragma omp parallel for
	  for (int resi = 0; resi < n2; ++resi)
	  {
		double *arr = (double *)malloc(sizeof(double)*n3);
		std::string filename("/voxel/voxel_");
		filename = ofname + filename;
		filename+= std::to_string(starti+resi);
		filename+=type;
		//direc.insert( ( ( valid[starti + resi].z*g[1] + valid[starti + resi].y )*g[0] + valid[starti + resi]x ) ,starti+resi);
		const char *s1 = filename.c_str();
		std::ofstream f(s1,std::ios::binary|std::ios::app);
		for (int resj = 0; resj < n3; ++resj)
		   arr[resj] = res[resi*n+resj];
		  f.write(reinterpret_cast<char *>(&arr),sizeof(double)*n3);
		  // fprintf(f, "%f ",res[resi*n+resj]);

		f.close();
		free(arr);

	  }
	if(starti!=startj){
		#pragma omp parallel for
		for (int resj = 0; resj < n3; ++resj)
		{
		  double *arr = (double *)malloc(sizeof(double)*n3);
		  std::string filename("/voxel/voxel_");
		  filename = ofname + filename;
		  filename+= std::to_string(startj+resj);
		  const char *s1 = filename.c_str();
		  filename+=type;
		  std::ofstream f(s1,std::ios::binary|std::ios::app);
		  for (int resi = 0; resi < n2; ++resi)
			arr[resj] = res[resi*n+resj];
			f.write(reinterpret_cast<char *>(&arr),sizeof(double)*n3);
			// f.write((char*)&res[resi*n+resj],sizeof(double));
			// fprintf(f, "%f ",res[resi*n+resj]);
		  f.close();
		  free(arr);
		}  
	  }

	}
	std::cout<<"Time taken:  for CORRELATION "<< timeTk<<std::endl;
  }
  else{
	if(nifti_parser.load_from_file(roifname));
		nifti_parser >> image_data_cpy;
	int no_roi = 0;
	g_copy = image_data_cpy.geometry();



	 #pragma omp parallel for collapse(3) shared(no_roi)
	  for (int x = 0; x < g_copy[0]; ++x)
		for (int y = 0; y < g_copy[1]; ++y)
		   for (int z = 0; z < g_copy[2]; ++z)
			  if(image_data_cpy[((z)*g[1]+y)*g[0]+x]>0)
				no_roi++;
	std::vector<Valid> roi_val;


	double* roi = (double * )malloc( sizeof(double)*no_roi*g[3]);
	int i =0 ;
	for (int x = 0; x < g_copy[0]; ++x)
		for (int y = 0; y < g_copy[1]; ++y)
		   for (int z = 0; z < g_copy[2]; ++z){
			  

			  if(image_data_cpy[((z)*g[1]+y)*g[0]+x]>0){
				Valid tv ;
				tv.x = x;
				tv.y = y;
				tv.z = z;
				roi_val.push_back(tv);
				for (int t = 0; t < g[3]; ++t)
				  roi[i*g[3]+t] = image_data[((t*g[2] + z)*g[1]+y)*g[0]+x];

				i++;

			  }
			}

	  free(EX);
	  int dime = g[3];
	  int n =  2048;
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
	 

	  // std::cout<<valid.size()<<std::endl;
	   clock_t tStart = clock();
	// #pragma omp parallel for
	for(int starti = 0;starti<roi_val.size();starti+=n)
	  for(int startj = 0;startj<valid.size();startj+=n){ 
		// std::cout<<starti<<","<<startj<<std::endl;
		int n2 = min(n,roi_val.size()-starti);
		int n3 = min(n,valid.size()-startj);
		int n_toget=max(n2,n3);

		#pragma omp parallel for collapse(2)
		for (int i = 0; i < n_toget; ++i)
		{
		  for (int j = 0; j < g[3]; ++j)
		  {
			if(i<n2)
			Valid_matrix_chunk[i*g[3] + j]=roi[(starti+i)*g[3] + (j)];
			if(i<n3)
			Valid_matrix_chunk_trans[j*n + i]=Valid_matrix[(startj+i)*g[3] + (j)];
		  }
		}


		// std::cout<<<<std::endl;
		// std::cout<<"done something..."<<std::endl;
		// std::cout<<"cblas_dgemm("<<CblasRowMajor<<","<<CblasNoTrans<<","<<CblasNoTrans<<","<<n<<","<<n<<","<<dime<<","<<1.0<<","<<"Valid_matrix"<<","<<dime<<","<<"Valid_matrix_trans,"<<n<<","<<0.0<<",res,"<<n<<")"<<std::endl;
		
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n2,n3,dime,1.0,Valid_matrix_chunk,dime,Valid_matrix_chunk_trans,n3,0.0,res,n3);
		 std::cout<<"Time taken:  for CORRELATION "<< ((double)(clock() - tStart)/CLOCKS_PER_SEC)<<std::endl;
		 tStart = clock();
		# pragma omp parallel for
		for (int resi = 0; resi < n2; ++resi)
		{
		  double *arr = (double *)malloc(sizeof(double)*n3);
		  std::string filename("/roi/roi_");
		  filename = ofname + filename;
		  filename+= std::to_string(starti+resi);
		  filename+=type;
		  const char *s1 = filename.c_str();
		  std::ofstream f(s1,std::ios::binary|std::ios::app);
		  for (int resj = 0; resj < n3; ++resj)
			 arr[resj] = res[resi*n+resj];
			f.write(reinterpret_cast<char *>(&arr),sizeof(double)*n3);
			// fprintf(f, "%f ",res[resi*n+resj]);

		  f.close();
		  free(arr);

		}


	  }

  }
}


// #########################################get options #############################################################################
void getattributes(int argc,char *argv[])
{
  char tmp;

  while((tmp=getopt(argc,argv,"h:i:r:o:t:s:x:y:z:a:"))!=-1)
  {
	switch(tmp)
	{
	  /*option h show the help infomation*/
	  case 'h':
		showhelpinfo();
		break;
	  /*option u present the username*/
	  case 'i':
		input= true;
		// std::cout<<optarg<<std::endl;
		// strcpy(ipfilename,optarg);
		ipfilename = optarg;
		if(ipfilename.find(".gz")!=std::string::npos){
			gzip = true;
			std::string command = "gunzip ";
			command += ipfilename;
			system(command.c_str());
			ipfilename = ipfilename.substr(0,(ipfilename.length()-3));
		}

		// if(ipfilename.find(".dcm")!=std::string::npos)
		// 	dcm = true;

		break;
	  /*option p present the password*/ 
	  case 'r':
		roi= true;
		strcpy(roifname,optarg);
		break;
	  case 'o':
		output = true;
		//strcpy(ofname,optarg);
		ofname =optarg;
		break;
	  case 't':
		// std::cout<<"thresh"<<std::endl;
		// std::cout<<optarg<<std::endl;
		thresh = std::stoi(optarg);
	  break;
	  case 's':
		 seedcmp = true;
		 tmp =getopt(argc,argv,"r");
		 
		 if(tmp!=':'){
		 roi= true;
		 strcpy(roifname,optarg);
		 break;

		}
		 tmp =getopt(argc,argv,"x");
		 switch(tmp)
		 {
		  case 'x':seedx = std::stoi(optarg); break;
		  default : showhelpinfo(); exit(1);


		 }
		 tmp =getopt(argc,argv,"y");
		 switch(tmp)
		 {
		  case 'y':seedy = std::stoi(optarg); break;
		  default : showhelpinfo(); exit(1);


		 }
		 tmp =getopt(argc,argv,"z");
		 switch(tmp)
		 {
		  case 'z':seedz = std::stoi(optarg); break;
		  default : showhelpinfo(); exit(1);


		 }
	  break;
	  case 'a': all = true;
				// std::cout<<all<<std::endl;
				break;
	  /*invail input will get the heil infomation*/
	  default:
		showhelpinfo();
	  break;
	}
  }
  // return 0;

}

/*################################funcion that show the help information################################################################*/
void showhelpinfo()
{
 std::cout<<"Usage :\n fconn -i <fmri.nii/nii.gz>  \n";
  
  std::cout<<"Compulsory arguments(You have to specify the following)\n";
  std::cout<<"\t -i\t\t filename of the input volume \n";
  std::cout<<"\t -o \t\t project name  \n ";
  std::cout<<" one of the arguments must be present \n";
  std::cout<<"\t -r \t\t filename of the volume containg the desired ROI \n";
  std::cout<<"\t -s 1 \t\t for seed to all voxel mode(no argument) \n";
  std::cout<<"\t\t -x \t\t x-coordinate for seed (compulosry in -s mode if -r options not there)\n";
  std::cout<<"\t\t -y \t\t y-coordinate for seed (compulosry in -s mode if -r options not there)\n";
  std::cout<<"\t\t -z \t\t z-coordinate for seed (compulosry in -s mode if -r options not there)\n";
  std::cout<<"\t\t -r \t\t roi whose average will be seed (compulosry in -s mode if -x -y -z options not there)\n";
  std::cout<<"\t -a 1 \t\t for all voxels to all voxels mode(no argument) \n";
  std::cout<<"Optional arguments(You may optionally specify the following)\n";
  std::cout<<"\t -t \t an upper threshold  \n ";
  std::cout<<"\t -h \t display this message  \n ";
  // 
}
