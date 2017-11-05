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

std::string type(".csv"),ofname,ipfilename;
bool input = false,roi = false,output=false,gzip = false;
bool mask = false;
char roifname[400],maskfilename[400];
int thresh=0;
bool seedcmp = false,all=false;
int seedx ,seedy,seedz;


void showhelpinfo();
void correl();
void seed_correl();
void all_pair_corr();
void avg_corr_roi();
void avg_roi_time_corr();
void getattributes(int argc,char *argv[]);

int main (int argc,char *argv[])
{
  omp_set_num_threads(8);
  /*if the program is ran witout options ,it will show the usage and exit*/
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


void seed_correl(){

//DEFINING VARS
  image::io::nifti nifti_parser,nifti_parser2;
  image::basic_image<double,4> image_data;
  image::basic_image<double,4> image_data_cpy;
  image::basic_image<int,4> mask_image;
  std::vector<Valid> valid;

  image_data_cpy = image_data;

//LOADING THE IMAGE
  if(nifti_parser.load_from_file(ipfilename));
	  nifti_parser >> image_data;

// LOADING THE MASK
  if(mask&&nifti_parser2.load_from_file(maskfilename))
  	   nifti_parser2 >> mask_image;
  
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

  std::string cmd = "gzip ";
  cmd+=ipfilename;
  system(cmd.c_str());

//############################# DEFINING MORE VARIABLES ####################################

  double *EX = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  double *EX_sq = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  std::string filename2("/voxelmap");
  filename2 = ofname + filename2;
  const char *s1 = filename2.c_str();
  std::ofstream f(s1,std::ios::binary);

// #########################CALCULATING THE MEAN AND STD DEVIATION OF DATA AND CREATING A VALID VOXEL MAPPINGS ################################
  
  for (int z = 0; z < g[2]; ++z){
	for (int y = 0; y < g[1]; ++y)
	  {
		for (int x = 0; x < g[0]; ++x)
		{

		  double temp  = 0,temp2 = 0;
		  for (int t = 0; t < g[3]; ++t)
		  {
			temp += (double)image_data[((t*g[2] + z)*g[1]+y)*g[0]+x]/g[3] ;
			temp2 += temp*temp*g[3];
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


// ########################################### MAKING A MATRIX OF NORMALIZED DATA###########################################
  double * Valid_matrix = (double * )malloc( sizeof(double)*valid.size()*g[3]);
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
	cblas_dgemv(CblasRowMajor,CblasNoTrans,m,dime,timeseries,Valid_matrix,dime,seed,1,0.0,Corr,1);
	 
	for (int i = 0; i < valid.size(); ++i)
		image_data_cpy[((valid[i].z)*g[1]+valid[i].y)*g[0]+valid[i].x] = Corr[i];

  

	}

	image::io::nifti nifti_parser2;
	nifti_parser2.load_from_image(image_data_cpy);
	std::string filename("output.nii");
	filename = ofname + filename;
	nifti_parser2.save_to_file(ofname.c_str());
	// std::cout<< ((double)(clock() - tStart)/CLOCKS_PER_SEC)<<",";




}

void all_pair_corr(){
 
//DEFINING VARIABLES  
  image::io::nifti nifti_parser,nifti_parser2;
  image::basic_image<double,4> image_data;
  image::basic_image<double,4> image_data_cpy;
  image::basic_image<int,4> mask_image;
  std::vector<Valid> valid;

  image_data_cpy = image_data;

//LOADING THE IMAGE
  if(nifti_parser.load_from_file(ipfilename));
	  nifti_parser >> image_data;

// LOADING THE MASK
  if(mask&&nifti_parser2.load_from_file(maskfilename))
  	   nifti_parser2 >> mask_image;
  
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

  std::string cmd = "gzip ";
  cmd+=ipfilename;
  system(cmd.c_str());

//############### DEFINING VARIABLES ####################################

  double *EX = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  double *EX_sq = (double *) malloc( sizeof(double) * g[0] * g[1]* g[2] ) ;
  std::string filename2("/voxelmap");
  filename2 = ofname + filename2;
  const char *s1 = filename2.c_str();
  std::ofstream f(s1,std::ios::binary);

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
			temp2 += temp*temp*g[3];
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

	  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n2,n3,dime,timeseries,Valid_matrix_chunk,dime,Valid_matrix_chunk_trans,n3,0.0,res,n3);
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
		   double tmp = res[resi*n+resj];
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
		  	double tmp = res[resi*n+resj];
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
	//std::cout<<"Time taken:  for CORRELATION "<< timeTk<<std::endl;

}


void avg_roi_time_corr(){

}

void avg_corr_roi(){

}

void correl()
{
  if(seedcmp==true){
  	seed_correl()
  }else if(roi==true){


  }else if(all==true){
  		all_pair_corr();
  }
}


// #########################################get options #############################################################################
void getattributes(int argc,char *argv[])
{
  char tmp;

  while((tmp=getopt(argc,argv,"h:i:r:o:t:s:x:y:z:a:m:"))!=-1)
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
	  case 'm': mask = true;
	  			strcpy(maskfilename,optarg);
	  			break;
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
