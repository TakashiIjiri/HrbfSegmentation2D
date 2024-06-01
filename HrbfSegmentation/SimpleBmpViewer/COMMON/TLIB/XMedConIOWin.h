#include<sstream>
#include<iostream>


#define MEDCON "medcon"



/*------------------------------------------------------------------------------
this class use "xmedcon.exe" for loading dcm file
the user of the final executable has to install xmedcon and setup path
------------------------------------------------------------------------------*/




// for UNIX
//#define MEDCON "plugin/medcon.exe"
/* the String "ext" should be one of
	  XMedCon Formats
	  "acr"	Acr/Nema 2.0
	  "anlz"	Analyze (SPM)
	  "conc"	Concorde/microPET
	  "dicom"	DICOM 3.0 (NM)
	  "ecat6" or "ecat"	CTI ECAT 6.4
	  "ecat7"	CTI ECAT 7.2
	  // not supported "gif"	Gif89a
	   "intf"	InterFile 3.3
	   "inw"	INW (RUG)
	   "nifti"	NIfTI Analyze
	   "png"	PNG
	   "bin"	Raw binary
	   // not supported "ascii"	Raw ascii
	   */
class XMedConIO{
 public:
  int sx,sy,sz,bit,type;
  double px,py,pz;
  double intercept;
  std::string imagedata; 
  int CMAX,SMAX;
  //unsigned long IMAX;

  //unsigned long long LMAX;
  
  char tmpbinaryheaderfilename[256];
  char tmpbinaryimagefilename[256];
  char tmpouttxt[256];
  int numberofimages;
  /*
    XMedCon Formats
    -c, --convert <format> ...
    Notation	Format
    "acr"	Acr/Nema 2.0
    "anlz"	Analyze (SPM)
    "conc"	Concorde/microPET
    "dicom"	DICOM 3.0 (NM)
    "ecat6" or "ecat"	CTI ECAT 6.4
    "ecat7"	CTI ECAT 7.2
    "gif"	Gif89a // unsupported for RIKEN's system
    "intf"	InterFile 3.3
    "inw"	INW (RUG)
    "nifti"	NIfTI Analyze
    "png"	PNG
    "bin"	Raw binary // unsupported for RIKEN's system
    "ascii"	Raw ascii // unsupported for RIKEN's system
    
  */
  /* pixel types                             */
  //#define BIT1      1    /*  1-bit           */unsupported
  //#define BIT8_S    2    /*  8-bit signed    */
  //#define BIT8_U    3    /*  8-bit unsigned  */
  //#define BIT16_S   4    /* 16-bit signed    */
  //#define BIT16_U   5    /* 16-bit unsigned  */
  //#define BIT32_S   6    /* 32-bit signed    */
  //#define BIT32_U   7    /* 32-bit unsigned  */
  //#define BIT64_S   8    /* 64-bit signed    */
  //#define BIT64_U   9    /* 64-bit unsigned  */
  //#define FLT32    10    /* 32-bit float     */
  //#define FLT64    11    /* 64-bit double    */
  //#define ASCII    12    /* ascii            */ unsupported
  //#define VAXFL32  13    /* 32-bit vaxfloat  */
  //#define COLRGB   20    /* RGB triplets     */
  
  XMedConIO(){
    sx=sy=-1;sz=1;numberofimages=1;
    px=py=pz=1.0;
	intercept = 0.0;
    bit=-1;  type=-1;
    CMAX = 2*SCHAR_MAX+1;
    SMAX = 2*SHRT_MAX+1;
    //IMAX = 2*INT_MAX+1;
    //LMAX = 2*LONG_MAX+1;
    sprintf(tmpbinaryheaderfilename,"predef_tmp.rpi");
    sprintf(tmpbinaryimagefilename,"predef_tmp.bin");
    sprintf(tmpouttxt,"tmp.medcon.txt");
  }
  virtual ~XMedConIO(){}

	/***************************************************/
	//INPUT from medcon supported file formats: header and binary image  
  	bool setFileGetHeaderInfoWin(const char *filename)
	{
		if(filename==NULL)return false;
		std::string fname = std::string(filename);
		int filenamesize = (int) fname.size();
		if(filenamesize<=0)return false;

		fflush(stdout);
		char *cmd = new char[(1024+filenamesize)];
		//sprintf(cmd,"medcon -n -f %s -c - bin 2> medcon_tmp_input.txt 1> medcon_tmp_data.bin",filename);
		//system(cmd);

		sprintf(cmd,"%s -n -f %s -c - bin 2> medcon_tmp_input.txt",MEDCON, filename);

		FILE *dfp = _popen(cmd,"rb");//for Windows
    
		if(dfp!=NULL){ 
			while(!feof(dfp)){
				char tmpc = fgetc(dfp);
	
				//for Windows
				if(tmpc!='\x0A') imagedata+=tmpc;
			}
			fclose(dfp);	
		}else{
			fprintf(stderr, "error in setFileGetHeaderInfoWin: NULL\n");
			fprintf(stderr, "cmd = %s\n",cmd);
		}

		// FILE *hfp = fopen("medcon_tmp_input.txt","r"); for UNIX
		FILE *hfp = fopen("medcon_tmp_input.txt","rb");
		std::string header;
		while(!feof(hfp)){
			//sz_b = fread(buffer_b,sizeof(buffer_b),BUFSIZ,hfp);
			//if(sz_b<=0||sz_b==EOF)break;
			//for(int i=0;i<sz_b;i++)header+=buffer_b[i];
			char tmpc = fgetc(dfp);
			//for Windows
			if(tmpc!='\x0A') header+=tmpc;
		}
    /*
      while(1){
      sz_b = fread(buffer_b,sizeof(buffer_b),BUFSIZ,hfp);
      if(sz_b<=0||sz_b==EOF)break;
      for(int i=0;i<sz_b;i++)header+=buffer_b[i];
      
    }
    */
		fclose(hfp);
    /*
    FILE *dfp = fopen("medcon_tmp_data.bin","rb");
     while(!feof(dfp)){
       char tmpc = fgetc(dfp);
       //if(tmpc!='\x0A')
       if(tmpc!=EOF)
       imagedata+=tmpc;
     }
    fclose(dfp);
    */
    //Data Size(Byte): 524288
    //printf("header size = %d\n",header.size());
    //printf("imagedata size = %d\n",imagedata.size());
    //std::cout << header << std::endl;
		delete [] cmd;
		std::istringstream headerstream(header);
		std::string tmp1,tmp2;
      
		while(!headerstream.fail()){
			headerstream >> tmp1;
			if(tmp1=="mwidth"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream >> sx;
				//for debug std::cout << "width = "<< sx << std::endl;
			}else if(tmp1=="mheight"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream  >> sy;
				//for debug std::cout << "height = "<< sy << std::endl;
			}else if(tmp1=="number"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream >> sz;
				//for debug std::cout << "number = "<< tmp2 << " " << sz <<std::endl;
			}else if(tmp1=="pixdim[1]"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream >> px;
				//for debug std::cout << "px = "<< px << std::endl;
			}else if(tmp1=="pixdim[2]"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream >> py;
				//for debug std::cout << "py = "<< py << std::endl;
			}else if(tmp1=="pixdim[3]"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream >> pz;
				//for debug std::cout << "pz = "<< pz << std::endl;	
			}else if(tmp1=="bits"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream >> bit;
				//for debugstd::cout << "bit = "<< bit << std::endl;
			}else if(tmp1=="type"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream >> type;
				//for debug  std::cout << "type = "<< type << std::endl;
			}else if(tmp1=="intercept"){
				headerstream >> tmp2;
				if(tmp2==":") headerstream >> intercept;
			}
	}
      /* for debug print header */
      //std::cout << header << std::endl;
      
      

  
    // is data there ?
    if(imagedata.size()<=0)return false;
    
    // check format consistency
    if(type<=1||type==12)return false;
    if(!(bit==8||bit==16||bit==24||bit==32||bit==64))return false;
	if(type==2||type==3)if(bit!=8)return false;
    if(type==4||type==5)if(bit!=16)return false;
    if(type==20)if(bit!=24)return false;
    if(type==6||type==7||type==10||type==13)if(bit!=32)return false;
    if(type==8||type==9||type==11)if(bit!=64)return false;
    return true;
  }

  
  /***************************************************/
  //setFileGetHeaderInfo should be called before call ReadRawBinaryGray
  //single,
  bool ReadRawBinaryGray(double **data){
    if(sx<=0||sy<=0||px<=0.0||py<=0.0||pz<=0.0)return false;
    	int i,j,posi=0;
    if(bit==8){
      if(type==2){//char
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    data[i][j] = ((double)(imagedata[posi]));
	    posi++;
	  }
      }else if(type==3){//unsigned char
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    data[i][j] = ((double)(uchar2int(imagedata[posi])));
	    posi++;
	  }
	
      }else{
	return false;
      }
      
    }else if(bit==16){
      if(type==4){//short
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 2*(i*sx+j);
	    data[i][j] = ((double)(getShortFromChar(posi)));
	  }
      }else if(type==5){//unsigned short
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 2*(i*sx+j);
	    data[i][j] = ((double)(getUShortFromChar(posi)));
	    
	  }
      }else{
	return false;
      }
    }else if(bit==32){
      if(type==10||type==13){//float
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 4*(i*sx+j);
	    data[i][j] = ((double)(getFloatFromChar(posi)));
	  }
	
      }else if(type==6){//int
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 4*(i*sx+j);
	    data[i][j] = ((double)(getIntFromChar(posi)));
	    
	  }
      }else if(type==7){//unsigned int
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 4*(i*sx+j);
	    data[i][j] = ((double)(getUIntFromChar(posi)));
	    
	  }
	
      }else if(bit==64){
	if(type==8){//long  not correctly supported for 32bit machine
	  for(i=0;i<sy;i++)
	    for(j=0;j<sx;j++){
	      posi = 8*(i*sx+j);
	      data[i][j] = ((double)(getLongFromChar(posi)));
	    }
	  
	}else if(type==9){//unsigned long not correctly supported for 32bit machine
	  for(i=0;i<sy;i++)
	    for(j=0;j<sx;j++){
	      posi = 8*(i*sx+j);
	      data[i][j] = ((double)(getULongFromChar(posi)));
	    }
	  
	  
	}else if(type==11){//double
	  for(i=0;i<sy;i++)
	    for(j=0;j<sx;j++){
	      posi = 8*(i*sx+j);
	      data[i][j] = ((double)(getDoubleFromChar(posi)));
	    }
	  
	}else{
	  return false;
	}
	
      }else{
	return false;
      }
      
    }
    
    return true;
  }

  
  //multi,
  bool ReadRawBinaryGray(double ***data){
    if(sx<=0||sy<=0||sz<=0||px<=0.0||py<=0.0||pz<=0.0)return false;
    int i,j,kz,posi=0;
    if(bit==8){
      if(type==2){//char
	for(kz=0;kz<sz;kz++)
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    data[kz][i][j] = ((double)(imagedata[posi]));
	    posi++;
	  }
      }else if(type==3){//unsigned char
		for(kz=0;kz<sz;kz++)
		  for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    data[kz][i][j] = ((double)(uchar2int(imagedata[posi])));
	    posi++;
	  }
	
      }else{
	return false;
      }
      
    }else if(bit==16){
      if(type==4){//short
		for(kz=0;kz<sz;kz++)for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 2*(kz*sx*sy+i*sx+j);
	    data[kz][i][j] = ((double)(getShortFromChar(posi)));
	  }
      }else if(type==5){//unsigned short
		for(kz=0;kz<sz;kz++)for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 2*(kz*sx*sy+i*sx+j);
	    data[kz][i][j] = ((double)(getUShortFromChar(posi)));
	    
	  }
      }else{
	return false;
      }
    }else if(bit==32){
      if(type==10||type==13){//float
	for(kz=0;kz<sz;kz++)for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 4*(kz*sx*sy+i*sx+j);
	    data[kz][i][j] = ((double)(getFloatFromChar(posi)));
	  }
	
      }else if(type==6){//int
	for(kz=0;kz<sz;kz++)for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 4*(kz*sx*sy+i*sx+j);
	    data[kz][i][j] = ((double)(getIntFromChar(posi)));
	    
	  }
      }else if(type==7){//unsigned int
	for(kz=0;kz<sz;kz++)for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    posi = 4*(kz*sx*sy+i*sx+j);
	    data[kz][i][j] = ((double)(getUIntFromChar(posi)));
	    
	  }
	
      }else if(bit==64){
	if(type==8){//long  not correctly supported for 32bit machine
	  for(kz=0;kz<sz;kz++)for(i=0;i<sy;i++)
	    for(j=0;j<sx;j++){
	      posi = 8*(kz*sx*sy+i*sx+j);
	      data[kz][i][j] = ((double)(getLongFromChar(posi)));
	    }
	  
	}else if(type==9){//unsigned long not correctly supported for 32bit machine
	  for(kz=0;kz<sz;kz++)for(i=0;i<sy;i++)
	    for(j=0;j<sx;j++){
	      posi = 8*(kz*sx*sy+i*sx+j);
	      data[kz][i][j] = ((double)(getULongFromChar(posi)));
	    }
	  
	  
	}else if(type==11){//double
	  for(kz=0;kz<sz;kz++)for(i=0;i<sy;i++)
	    for(j=0;j<sx;j++){
	      posi = 8*(kz*sx*sy+i*sx+j);
	      data[kz][i][j] = ((double)(getDoubleFromChar(posi)));
	    }
	  
	}else{
	  return false;
	}
	
      }else{
	return false;
      }
      
    }
    
    return true;
  }
  /***************************************************/
  //setFileGetHeaderInfo should be called before call ReadRawBinaryColor
  //Single
  bool ReadRawBinaryColor(double **dataR,double **dataG,double **dataB){
    if(sx<=0||sy<=0||px<=0.0||py<=0.0||pz<=0.0)return false;
    int i,j,posi=0;
    if(bit==24){
      for(i=0;i<sy;i++)
	for(j=0;j<sx;j++){
	  
	  dataR[i][j] = ((double)(uchar2int(imagedata[3*posi])));
	  dataG[i][j] = ((double)(uchar2int(imagedata[3*posi+1])));
	  dataB[i][j] = ((double)(uchar2int(imagedata[3*posi+2])));
	  //if(dataR[i][j]<0||dataG[i][j]<0||dataB[i][j]<0)printf("%lf %lf %lf\n",dataR[i][j],dataG[i][j],dataB[i][j]);
	  
	  posi++;
	  
	}
      
    }else{
      return false;
    }
    
    return true;
  }
   //Multi
  bool ReadRawBinaryColor(double ***dataR,double ***dataG,double ***dataB){
    if(sx<=0||sy<=0||sz<=0||px<=0.0||py<=0.0||pz<=0.0)return false;
    int i,j,kz,posi=0;
    if(bit==24){
      for(kz=0;kz<sz;kz++)
      for(i=0;i<sy;i++)
	for(j=0;j<sx;j++){
	  
	  dataR[kz][i][j] = ((double)(uchar2int(imagedata[3*posi])));
	  dataG[kz][i][j] = ((double)(uchar2int(imagedata[3*posi+1])));
	  dataB[kz][i][j] = ((double)(uchar2int(imagedata[3*posi+2])));
	  //if(dataR[i][j]<0||dataG[i][j]<0||dataB[i][j]<0)printf("%lf %lf %lf\n",dataR[i][j],dataG[i][j],dataB[i][j]);
	  
	  posi++;
	  
	}
      
    }else{
      return false;
    }
    
    return true;
  }

   /**********************************************************/
  // Print Header for debug.
	void PrintHeaderInfo(char* filename)
	{
		fprintf(stderr, "File Name: %s\n",filename);
		fprintf(stderr, "Number of Images in the File %d\n",sz);
		fprintf(stderr, "Bits: %d\n",bit);
		fprintf(stderr, "type: %d\n",type);
		fprintf(stderr, "Size(X,Y): (%d,%d)\n",sx,sy);
		fprintf(stderr, "Pitch(Z:Y:Z): (%lf:%lf:%lf)\n",px,py,pz);
		fprintf(stderr, "Data Size(Byte): %d\n",imagedata.size());
	}
    //Print image to PGM for debug.
     void PrintPGM(char* filename,double **data){
	int maxval=0;
	Shift(data);
	maxval = ((int)(data[0][0]));
	int i,j;
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    if(maxval<data[i][j])maxval=((int)(data[i][j]));
	  }
	FILE *out = fopen(filename,"w");
	fprintf(out,"P2\n");
	fprintf(out,"%d %d\n",sx,sy);
	fprintf(out,"%d\n",maxval);
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    fprintf(out,"%d ",((int)(data[i][j])));
	  }
	fprintf(out,"\n");
	fclose(out);
     }
    //Print image to PPM for debug.
     void PrintPPM(char * filename,double **dR,double **dG,double **dB){
	int maxval=0;
	Shift(dR,dG,dB);
	maxval = ((int)(dR[0][0]));
	int i,j;
	for(i=0;i<sy;i++)
	    for(j=0;j<sx;j++){
		if(maxval<dR[i][j])maxval=((int)(dR[i][j]));
		if(maxval<dG[i][j])maxval=((int)(dG[i][j]));
		if(maxval<dB[i][j])maxval=((int)(dB[i][j]));
		
	    }
	FILE *out = fopen(filename,"w");
	fprintf(out,"P3\n");
	fprintf(out,"%d %d\n",sx,sy);
	fprintf(out,"%d\n",maxval);
	for(i=0;i<sy;i++)
	  for(j=0;j<sx;j++){
	    fprintf(out,"%d %d %d ",((int)(dR[i][j])),((int)(dG[i][j])),((int)(dB[i][j])));
	  }
	fprintf(out,"\n");
	fclose(out);
     }
     //Shift Intensity to positive value.
     void Shift(double **data){
       double min = getMinVal(data);
       if(min<0.0){
	 int i,j;
	 for(i=0;i<sy;i++)
	   for(j=0;j<sx;j++){
	     data[i][j]-=min;
	   }
       }
     }
     
     void Shift(double **dR,double **dG,double **dB){
       
       double minR = getMinVal(dR);
       double minG = getMinVal(dG);
       double minB = getMinVal(dB);
       if(minR<0.0||minG<0.0||minB<0.0){
	 double min=minR;
	 if(min>minG)min=minG;
	 if(min>minB)min=minB;
	 int i,j;
	 for(i=0;i<sy;i++)
	   for(j=0;j<sx;j++){
	     dR[i][j]-=min;
	     dG[i][j]-=min;
	     dB[i][j]-=min;
	     
	   }
	 
       }
       
     }
     
    
    double getMinVal(double **data){
      double min = data[0][0];
      int i,j;
      for(i=0;i<sy;i++)
	for(j=0;j<sx;j++){
	  if(data[i][j]<min)min=data[i][j];
	}
      return min;
    }
    

  /*****************************
Binary Read Functions: imagedata shoud not be empty() 
  ************************/
 private:
    int uchar2int(char val){
      int ival=((int)(val));
      if(val<0){
	ival += CMAX;
      }
      return ival;
    }
    
  short getShortFromChar(int posi){
    short sval;
    char *myint = (char *)(&sval);
    myint[0] = imagedata[posi];
    myint[1] = imagedata[posi+1];
    return sval;
  }
  unsigned short getUShortFromChar(int posi){
    unsigned short sval;
    char *myint = (char *)(&sval);
    myint[0] = imagedata[posi];
    myint[1] = imagedata[posi+1];
    return sval;
  }
  int getIntFromChar(int posi){
    int sval;
    char *myint = (char *)(&sval);
    myint[0] = imagedata[posi];
    myint[1] = imagedata[posi+1];
    myint[2] = imagedata[posi+2];
    myint[3] = imagedata[posi+3];
    return sval;
  }
  unsigned int getUIntFromChar(int posi){
    unsigned int sval;
    char *myint = (char *)(&sval);
    myint[0] = imagedata[posi];
    myint[1] = imagedata[posi+1];
    myint[2] = imagedata[posi+2];
    myint[3] = imagedata[posi+3];
    return sval;
  }
  float getFloatFromChar(int posi){
    float sval;
    char *myint = (char *)(&sval);
    myint[0] = imagedata[posi];
    myint[1] = imagedata[posi+1];
    myint[2] = imagedata[posi+2];
    myint[3] = imagedata[posi+3];
    return sval;
  }
  
  long getLongFromChar(int posi){
    long sval;
    char *myint = (char *)(&sval);
    int i;
    for(i=0;i<8;i++)myint[i] = imagedata[posi+i];
    
    
    return sval;
  }
  unsigned long getULongFromChar(int posi){
    unsigned long sval;
    char *myint = (char *)(&sval);
    int i;for(i=0;i<8;i++)myint[i] = imagedata[posi+i];
    
    return sval;
  }
  double getDoubleFromChar(int posi){
    double sval;
    char *myint = (char *)(&sval);
    int i;
    for(i=0;i<8;i++)myint[i] = imagedata[posi+i];
    
    
    return sval;
  }
  /**********************
Binary Write Functions

  *************************/
   public:
  //single 
  bool SetHeader(int dtype,int dbit,int dx,int dy,double dpx,double dpy,double dpz){
       //true: supported
       //false: unsupported
	
       type = dtype;if(type<=1||type==12)return false;
       bit = dbit;if(!(bit==8||bit==16||bit==24||bit==32||bit==64))return false;
       if(type==2||type==3)if(bit!=8)return false;
       if(type==4||type==5)if(bit!=16)return false;
       if(type==20)if(bit!=24)return false;
       if(type==6||type==7||type==10||type==13)if(bit!=32)return false;
       if(type==8||type==9||type==11)if(bit!=64)return false;
       
       sx = dx;if(sx<=0)return false;
       sy = dy;if(sy<=0)return false;
       px = dpx;if(px<=0.0)return false;
       py = dpy;if(py<=0.0)return false;
       pz = dpz;if(pz<=0.0)return false;
       //Currently Raw Binary output to XMedCon is only supported for double or float, so let us return false where...
       
       if(type==11||type==10)return true;
       
       
       return false;
       
     }
  //multiple
  bool SetHeader(int dtype,int dbit,int dx,int dy,int dz,double dpx,double dpy,double dpz){
       //true: supported
       //false: unsupported
	
       type = dtype;if(type<=1||type==12)return false;
       bit = dbit;if(!(bit==8||bit==16||bit==24||bit==32||bit==64))return false;
       if(type==2||type==3)if(bit!=8)return false;
       if(type==4||type==5)if(bit!=16)return false;
       if(type==20)if(bit!=24)return false;
       if(type==6||type==7||type==10||type==13)if(bit!=32)return false;
       if(type==8||type==9||type==11)if(bit!=64)return false;
       
       sx = dx;if(sx<=0)return false;
       sy = dy;if(sy<=0)return false;
       sz = dz;if(sz<=0)return false;
       px = dpx;if(px<=0.0)return false;
       py = dpy;if(py<=0.0)return false;
       pz = dpz;if(pz<=0.0)return false;
       //Currently Raw Binary output to XMedCon is only supported for double or float, so let us return false where...
       
       if(type==11||type==10)return true;
       
       
       return false;
       
     }


     /**************************
Convert Formats Function

    ***********************/
 public:
  //single
     bool Save(double **data,char *outfilename,char *ext){
	
       if(!(strcmp(ext,"acr")==0||
	    strcmp(ext,"anlz")==0||
	    strcmp(ext,"conc")==0||
	    strcmp(ext,"dicom")==0||
	    strcmp(ext,"ecat6")==0||
	    strcmp(ext,"ecat7")==0||
	    strcmp(ext,"intf")==0||
	    strcmp(ext,"inw")==0||
	    strcmp(ext,"nifti")==0||
	    strcmp(ext,"png")==0)){
	 printf("Unsupported Formats !!\n");
	 return false;
       }
       /* the String "ext" should be one of
	  XMedCon Formats
	  "acr"	Acr/Nema 2.0
	  "anlz"	Analyze (SPM)
	  "conc"	Concorde/microPET
	  "dicom"	DICOM 3.0 (NM)
	  "ecat6" or "ecat"	CTI ECAT 6.4
	  "ecat7"	CTI ECAT 7.2
	  // not supported "gif"	Gif89a
	   "intf"	InterFile 3.3
	   "inw"	INW (RUG)
	   "nifti"	NIfTI Analyze
	   "png"	PNG
	   "bin"	Raw binary
	   // not supported "ascii"	Raw ascii
	   */
       
       if(!PrintHeader(tmpbinaryheaderfilename)){
	 printf("PrintHeader: error !\n");
	 return false;
       }
       if(!PrintRawBinaryGray(tmpbinaryimagefilename,data)){
	 printf("PrintRawBinaryGray: error !\n");
	 return false;
       }
	//difference between in Java and this is -s --silence and > dump message to tmpouttxt
	
   
	   
	int filenamesizes = (int)( strlen(tmpbinaryimagefilename)+strlen(outfilename)+strlen(tmpbinaryheaderfilename)+strlen(tmpouttxt) );
	char *cmd = new char[1024+filenamesizes];
	sprintf(cmd,"%s -s -n -w -noprefix -f %s -o %s -c %s -i < %s > %s",MEDCON,tmpbinaryimagefilename,outfilename,ext,tmpbinaryheaderfilename,tmpouttxt);
	
	
	system(cmd);
	delete [] cmd;
	char *cmd2 = new char[256+strlen(tmpbinaryheaderfilename)];
	sprintf(cmd2,"rm -rf %s",tmpbinaryheaderfilename);
	system(cmd2);
	delete [] cmd2;
	char *cmd3 = new char[256+strlen(tmpbinaryimagefilename)];
	sprintf(cmd3,"rm -rf %s",tmpbinaryimagefilename);
	system(cmd3);
	delete [] cmd3;
	char *cmd4 = new char[256+strlen(tmpouttxt)];
	sprintf(cmd4,"rm -rf %s",tmpouttxt);
	system(cmd4);
	delete [] cmd4;
	
	return true;
     }
     //multiple
     bool Save(double ***data,char *outfilename,char *ext){
	
       if(!(strcmp(ext,"acr")==0||
	    strcmp(ext,"anlz")==0||
	    strcmp(ext,"conc")==0||
	    strcmp(ext,"dicom")==0||
	    strcmp(ext,"ecat6")==0||
	    strcmp(ext,"ecat7")==0||
	    strcmp(ext,"intf")==0||
	    strcmp(ext,"inw")==0||
	    strcmp(ext,"nifti")==0||
	    strcmp(ext,"png")==0)){
	 printf("Unsupported Formats !!\n");
	 return false;
       }
       /* the String "ext" should be one of
	  XMedCon Formats
	  "acr"	Acr/Nema 2.0
	  "anlz"	Analyze (SPM)
	  "conc"	Concorde/microPET
	  "dicom"	DICOM 3.0 (NM)
	  "ecat6" or "ecat"	CTI ECAT 6.4
	  "ecat7"	CTI ECAT 7.2
	  // not supported "gif"	Gif89a
	   "intf"	InterFile 3.3
	   "inw"	INW (RUG)
	   "nifti"	NIfTI Analyze
	   "png"	PNG
	   "bin"	Raw binary
	   // not supported "ascii"	Raw ascii
	   */
       
       if(!PrintHeader(tmpbinaryheaderfilename)){
	 printf("PrintHeader: error !\n");
	 return false;
       }
       if(!PrintRawBinaryGray(tmpbinaryimagefilename,data)){
	 printf("PrintRawBinaryGray: error !\n");
	 return false;
       }
	//difference between in Java and this is -s --silence and > dump message to tmpouttxt
	
       int  filenamesizes = (int)(strlen(tmpbinaryimagefilename)+strlen(outfilename)+strlen(tmpbinaryheaderfilename)+strlen(tmpouttxt) );
	char *cmd = new char[1024+filenamesizes];
	sprintf(cmd,"%s -s -n -w -noprefix -f %s -o %s -c %s -i < %s > %s",MEDCON,tmpbinaryimagefilename,outfilename,ext,tmpbinaryheaderfilename,tmpouttxt);
	
	
	system(cmd);
	delete [] cmd;
	char *cmd2 = new char[256+strlen(tmpbinaryheaderfilename)];
	sprintf(cmd2,"rm -rf %s",tmpbinaryheaderfilename);
	system(cmd2);
	delete [] cmd2;
	char *cmd3 = new char[256+strlen(tmpbinaryimagefilename)];
	sprintf(cmd3,"rm -rf %s",tmpbinaryimagefilename);
	system(cmd3);
	delete [] cmd3;
	char *cmd4 = new char[256+strlen(tmpouttxt)];
	sprintf(cmd4,"rm -rf %s",tmpouttxt);
	system(cmd4);
	delete [] cmd4;
	
	return true;
     }

     //color, single
     bool Save(double **dataR,double **dataG,double **dataB,char *outfilename,char *ext){
	
       if(!(strcmp(ext,"acr")==0||
	    strcmp(ext,"anlz")==0||
	    strcmp(ext,"conc")==0||
	    strcmp(ext,"dicom")==0||
	    strcmp(ext,"ecat6")==0||
	    strcmp(ext,"ecat7")==0||
	    strcmp(ext,"intf")==0||
	    strcmp(ext,"inw")==0||
	    strcmp(ext,"nifti")==0||
	    strcmp(ext,"png")==0)){
	 printf("Unsupported Formats !!\n");
	 return false;
       }
       type=20;
       bit = 24;
       
       /* the String "ext" should be one of
	  XMedCon Formats
	  "acr"	Acr/Nema 2.0
	  "anlz"	Analyze (SPM)
	  "conc"	Concorde/microPET
	  "dicom"	DICOM 3.0 (NM)
	  "ecat6" or "ecat"	CTI ECAT 6.4
	  "ecat7"	CTI ECAT 7.2
	  // not supported "gif"	Gif89a
	   "intf"	InterFile 3.3
	   "inw"	INW (RUG)
	   "nifti"	NIfTI Analyze
	   "png"	PNG
	   "bin"	Raw binary
	   // not supported "ascii"	Raw ascii
	   */
       
       if(!PrintHeader(tmpbinaryheaderfilename)){
	 printf("PrintHeader: error !\n");
	 return false;
       }
       if(!PrintRawBinaryColor(tmpbinaryimagefilename,dataR,dataG,dataB)){
	 printf("PrintRawBinaryGray: error !\n");
	 return false;
       }
	//difference between in Java and this is -s --silence and > dump message to tmpouttxt
	
       int filenamesizes = (int)(strlen(tmpbinaryimagefilename)+strlen(outfilename)+strlen(tmpbinaryheaderfilename)+strlen(tmpouttxt));
	char *cmd = new char[1024+filenamesizes];
	sprintf(cmd,"%s -s -n -w -noprefix -f %s -o %s -c %s -i < %s > %s",MEDCON,tmpbinaryimagefilename,outfilename,ext,tmpbinaryheaderfilename,tmpouttxt);
	
	
	system(cmd);
	delete [] cmd;
	char *cmd2 = new char[256+strlen(tmpbinaryheaderfilename)];
	sprintf(cmd2,"rm -rf %s",tmpbinaryheaderfilename);
	system(cmd2);
	delete [] cmd2;
	char *cmd3 = new char[256+strlen(tmpbinaryimagefilename)];
	sprintf(cmd3,"rm -rf %s",tmpbinaryimagefilename);
	system(cmd3);
	delete [] cmd3;
	char *cmd4 = new char[256+strlen(tmpouttxt)];
	sprintf(cmd4,"rm -rf %s",tmpouttxt);
	system(cmd4);
	delete [] cmd4;
	
	return true;
     }

      //color, multiple
     bool Save(double ***dataR,double ***dataG,double ***dataB,char *outfilename,char *ext){
	
       if(!(strcmp(ext,"acr")==0||
	    strcmp(ext,"anlz")==0||
	    strcmp(ext,"conc")==0||
	    strcmp(ext,"dicom")==0||
	    strcmp(ext,"ecat6")==0||
	    strcmp(ext,"ecat7")==0||
	    strcmp(ext,"intf")==0||
	    strcmp(ext,"inw")==0||
	    strcmp(ext,"nifti")==0||
	    strcmp(ext,"png")==0)){
	 printf("Unsupported Formats !!\n");
	 return false;
       }
       type=20;
       bit = 24;
       
       /* the String "ext" should be one of
	  XMedCon Formats
	  "acr"	Acr/Nema 2.0
	  "anlz"	Analyze (SPM)
	  "conc"	Concorde/microPET
	  "dicom"	DICOM 3.0 (NM)
	  "ecat6" or "ecat"	CTI ECAT 6.4
	  "ecat7"	CTI ECAT 7.2
	  // not supported "gif"	Gif89a
	   "intf"	InterFile 3.3
	   "inw"	INW (RUG)
	   "nifti"	NIfTI Analyze
	   "png"	PNG
	   "bin"	Raw binary
	   // not supported "ascii"	Raw ascii
	   */
       
       if(!PrintHeader(tmpbinaryheaderfilename)){
	 printf("PrintHeader: error !\n");
	 return false;
       }
       if(!PrintRawBinaryColor(tmpbinaryimagefilename,dataR,dataG,dataB)){
	 printf("PrintRawBinaryGray: error !\n");
	 return false;
       }
	//difference between in Java and this is -s --silence and > dump message to tmpouttxt
	
       int filenamesizes = (int)(strlen(tmpbinaryimagefilename)+strlen(outfilename)+strlen(tmpbinaryheaderfilename)+strlen(tmpouttxt));
	char *cmd = new char[1024+filenamesizes];
	sprintf(cmd,"%s -s -n -w -noprefix -f %s -o %s -c %s -i < %s > %s",MEDCON,tmpbinaryimagefilename,outfilename,ext,tmpbinaryheaderfilename,tmpouttxt);
	
	
	system(cmd);
	delete [] cmd;
	char *cmd2 = new char[256+strlen(tmpbinaryheaderfilename)];
	sprintf(cmd2,"rm -rf %s",tmpbinaryheaderfilename);
	system(cmd2);
	delete [] cmd2;
	char *cmd3 = new char[256+strlen(tmpbinaryimagefilename)];
	sprintf(cmd3,"rm -rf %s",tmpbinaryimagefilename);
	system(cmd3);
	delete [] cmd3;
	char *cmd4 = new char[256+strlen(tmpouttxt)];
	sprintf(cmd4,"rm -rf %s",tmpouttxt);
	system(cmd4);
	delete [] cmd4;
	
	return true;
     }



 private:
     bool PrintHeader(char* filename){
       FILE *out = fopen(filename,"w");	
       if(out==NULL)return false;
       
       fprintf(out,"# RPI v0.1 - BEGIN #\n");
       fprintf(out,"#\n");
       fprintf(out,"# Total number of images?\n");
       fprintf(out,"%d\n",sz);
       fprintf(out,"# General header offset (bytes)?\n");
       fprintf(out,"0\n");
       fprintf(out,"# Image   header offset (bytes)?\n");
       fprintf(out,"0\n");
       fprintf(out,"# Repeated image header?\n");
       fprintf(out,"no\n");
       fprintf(out,"# Swap pixel bytes?\n");
       fprintf(out,"no\n");// yes in Java version of XMedConIO.java but no in C/C++ version
       fprintf(out,"# Identical images?\n");
       fprintf(out,"yes\n");
       fprintf(out,"# Absolute offset in bytes?\n");
       fprintf(out,"0\n");
       fprintf(out,"# Image columns?\n");
       fprintf(out,"%d\n",sx);
       fprintf(out,"# Image rows?\n");
       fprintf(out,"%d\n",sy);
       fprintf(out,"# Pixel data type?\n");
       fprintf(out,"%d\n",type);
       fprintf(out,"# Redo input?\n");
       fprintf(out,"no\n");
       fprintf(out,"# RPI v0.1 - END #\n");
       
       fclose(out);
	return true;
    }
     //single
     bool PrintRawBinaryGray(char *filename,double **data){
       int i;
       FILE *out = fopen(filename,"wb");
       if(out==NULL)return false;
       if(type==10){//float
	 size_t len = sizeof(float);
	 
	 
	 for(i=0;i<sy;i++)
	   fwrite(data[i],len,sx,out);
	 
	 
       }else if(type==11){//double
	 
	 size_t len = sizeof(double);
	 	 
	 for(i=0;i<sy;i++)
	   fwrite(data[i],len,sx,out); 
	 
       }else{
	 fclose(out);
	 return false;
       }
       fclose(out);
       return true;
     }
     //multiple
     bool PrintRawBinaryGray(char *filename,double ***data){
       int i,kz;
       FILE *out = fopen(filename,"wb");
       if(out==NULL)return false;
       if(type==10){//float
	 size_t len = sizeof(float);
	 
	 for(kz=0;kz<sz;kz++)
	 for(i=0;i<sy;i++)
	   fwrite(data[kz][i],len,sx,out);
	 
	 
       }else if(type==11){//double
	 
	 size_t len = sizeof(double);
	 for(kz=0;kz<sz;kz++)	 
	 for(i=0;i<sy;i++)
	   fwrite(data[kz][i],len,sx,out); 
	 
       }else{
	 fclose(out);
	 return false;
       }
       fclose(out);
       return true;
     }
     //single
     bool PrintRawBinaryColor(char *filename,double **dR,double **dG,double **dB){
	int i,j;
       FILE *out = fopen(filename,"wb");
       if(out==NULL)return false;
       if(type==20&&bit==24){//RGB color
	 size_t len = sizeof(unsigned char);
	 
	 unsigned char val=0;
	 for(i=0;i<sy;i++)
	   for(j=0;j<sx;j++){
	     val = ((unsigned char)(dR[i][j]));
	     fwrite(&val,len,1,out);
	     val = ((unsigned char)(dG[i][j]));
	     fwrite(&val,len,1,out);
	     val = ((unsigned char)(dB[i][j]));
	     fwrite(&val,len,1,out);
	     
	   }
       }else{
	 fclose(out);
	 return false;
       }
       fclose(out);
       return true;
     }
      //multiple
     bool PrintRawBinaryColor(char *filename,double ***dR,double ***dG,double ***dB){
       int i,j,kz;
       FILE *out = fopen(filename,"wb");
       if(out==NULL)return false;
       if(type==20&&bit==24){//RGB color
	 size_t len = sizeof(unsigned char);
	 
	 unsigned char val=0;
	 for(kz=0;kz<sz;kz++)
	   for(i=0;i<sy;i++)
	     for(j=0;j<sx;j++){
	       val = ((unsigned char)(dR[kz][i][j]));
	       fwrite(&val,len,1,out);
	       val = ((unsigned char)(dG[kz][i][j]));
	       fwrite(&val,len,1,out);
	       val = ((unsigned char)(dB[kz][i][j]));
	       fwrite(&val,len,1,out);
	     
	   }
       }else{
	 fclose(out);
	 return false;
       }
       fclose(out);
       return true;
     }
     
};

#undef W
#undef R
#undef MEDCON