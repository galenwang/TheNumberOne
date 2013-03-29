#include "stdafx.h"
#include "travelmap.h"

#include "tiffio.h"
#include "geotiffio.h"

//#define  SUPPORT_JPEG_TEXTURE
#ifdef SUPPORT_JPEG_TEXTURE
#include "studio.h"
#include "MagickWand.h"
#endif

#include "common.h"
#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osgDB/FileNameUtils>
#include <osgDB/FileUtils>
#include <osgDB/WriteFile>
#include <osg/CoordinateSystemNode>
#include <osg/Node>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/ShapeDrawable>
#include <osg/LineSegment>
#include <osgUtil/IntersectVisitor>
#include "sgstream.hxx"
#include "interface.h"
#include "TerrainProcessingUtils.h"
#include "CorrectionSpaceMesh.h"

#pragma warning( disable: 4309 4018)

//#define NOMERGE
#define TIFFTAG_GEOPIXELSCALE       33550
#define TIFFTAG_INTERGRAPH_MATRIX   33920
#define TIFFTAG_GEOTIEPOINTS        33922
#define TIFFTAG_GEOTRANSMATRIX      34264
#define TIFFTAG_GEOKEYDIRECTORY     34735
#define TIFFTAG_GEODOUBLEPARAMS     34736
#define TIFFTAG_GEOASCIIPARAMS      34737
#define TIFFTAG_CFAREPEATPATTERNDIM	33421	/* dimensions of CFA pattern */
#define TIFFTAG_CFAPATTERN			33422	/* color filter array pattern */
#define TEMP_TEXTFILE				"dest"
#define TEMP_TEXTGEOFILE			"T"
#define TEMP_HGTFILE				"hgt"
#define TEMP_HGTGEOFILE			 	"D"
#define TEMP_EXTENDNAME				"tif"
#define USER_HIGHRESOLUTION			"User"
#define SYSTEM_LOWRESOLUTION		"System"
#define LOG_FILE					"TerrainCreator.log"

/* 定义并记录当前状态 */
#define STATUS_UNWORK				0
#define STATUS_SUCCESS				1
#define STATUS_HGTBAD				2
#define STATUS_TEXTBAD				3
#define STATUS_HGTGEOBAD			4
#define STATUS_TEXTGEOBAD			5

/* 当前是纹理还是高程 */
#define STATUS_TEXTURE				1
#define STATUS_HGT					2
#define MBUFFSIZE 					0x100000

/* 基本字符串空间定义，防止缓冲区溢出。*/
#define STRING_STANDARD_LEN	128

using std::string;

/************************************************
 *
 *    引用外部过程及变量
 *
 ************************************************/
#ifdef SUPPORT_JPEG_TEXTURE
extern "C"/*MagickExport*/{
	void  MagickComponentTerminus(void);
	void MagickCoreGenesis(const char *,const MagickBooleanType);
	void MagickCoreTerminus(void);
	MagickBooleanType  MagickCommandGenesis(ImageInfo *,MagickCommand,int,char **,char **, ExceptionInfo *);
	MagickBooleanType  ConvertImageCommand(ImageInfo *,int,char **,char **,ExceptionInfo *);
}
#endif
extern int LoadShpAndDrawIt(char *shpFn);
extern int LoadShpAndCalculateMM(char *shpFn);
extern int CurrentTerrainIsCut;
extern void MessageToDialog01(char *msg);
extern void MessageToDialog02(char *msg);
extern void UpdateDialog01();
extern int UpdateIVEByBTG(osg::ref_ptr<osg::Group> btg, char *fnIve);
extern MaxMinLonLat currLayer;
extern std::string cultures[];
extern int culture_cc;
extern std::vector<TerPoint> curPnt;
extern std::vector<TerTextcoord> curText;
extern std::vector<TriType> curTri;
extern double dXpart, dYpart, dZpart;
extern float  fRadiusofBS;
extern int iCreationTime, totalNormal;
extern unsigned short iVersion, iMagicNum, iNumOfToplvlObject; 
extern double gbCenterX, gbCenterY, gbCenterZ;
extern LimitValue currMes;
extern bool foundAround;
extern int ReadFltAndCreateLightsDataStru(char *fname);
extern std::vector<SimPoint> CultPoints;
extern int AddLODIntoLightPoints();
extern int CreateLightsDataStructFromCultPoints();
extern std::vector<SimPoint> TreePoints;
extern int AddLODIntoTreePoints();	
extern int CreateTreesDataStructFromTreePoints();
extern void AnalysisBTG(char *btgName);
extern mt random_seed;
extern int currShpType;
extern int CreateSkirtDataStruFromSkirtGroup(osg::ref_ptr<osg::Group> gp);
extern int CreateRoadDataStruAndBtg(osg::ref_ptr<osg::Group> gp);
extern int CreateLakeDataStruAndBtg(osg::ref_ptr<osg::Group> gp);
extern int CreateNewStyleLakeDataStruAndBtg(osg::ref_ptr<osg::Group> gp);
extern int splitPathName(std::string matname);
extern std::vector<std::string> materialParts;
extern int ReadFltAndCreateAptLight(char *fname);

/************************************************
 *
 *    过程说明置顶
 *
 ************************************************/
double GetElevation(double lon, double lat, double pre);
int ReadALineToArray(FILE *fi);
int TravelDirectoryToReadAc(char *dPath);
int TravelDirectoryToReadShp(char *dPath);
int TravelDirectoryToReadFlt(char *dPath);
int CreatePerDegreeTerrain(int noVpb);
int MakeSTGFile(Bucket *Bt, char *main, int type);
int ReadTerrainIveAndCut(int lon, int lat);
int DoPartVpbAndWriteStg(int idx);
int DoHighVpbAndWriteStg(int idx);
int WriteAllStgs(char *ctileID, int type);
int CreateAirport( int argc, char **argv );
void ChangeCoordSystem(osg::Geode* gd, const osg::Matrix mt);
int MoveHigh50ToLow();
int ReadShpAndMakeBTG(int idx);
void MessageToDialog(char *msg);
int XchgModelBtgToDestination();
int ReadAptOSGAndOutputLineLoopIve();
int DoAllGA();
int DoAllFlt();
int CheckIfHas50(int iW, int iS);
void SetCurrHTI(HighTerrainInfo *curr);
int CreateUserSceneryDirectory();
void AddCorrectionSpaceMesh(osg::ref_ptr<osg::Group> gp);
int DeleteAllPrevCultureFiles(int lon, int lat);
int RemoveAllTileNodes();
int openhgt(TGhgt *hgt, char *file_name);
double GetElevationByLonLat(double lon, double lat);
int WriteHTxtInfo2Config(FILE *fo);
int WriteHHgtInfo2Config(FILE *fo);
int WriteHighInfo2Config(FILE *fo);
int WriteShpInfo2Config(FILE *fo);
int WriteBldInfo2Config(FILE *fo);

/************************************************
 *
 *    全局变量
 *
 ************************************************/
int ifHasIslandInCurrentTerrain;
int ifFoundBasePath;
float  mRoadWidthInShp = 50.0;
HANDLE MsgMutex;
int MessageTo;
std::vector<MapNode>  totalnode;
std::vector<HighNode> highTxtnode;
std::vector<HighNode> highHgtnode;
std::vector<CultNode> totalcult;
std::vector<LocateBuilding> totalUserBld;
std::vector<HighNode> lowTxtnode;
std::vector<CultNode> totalShp;
bool rebuildTerrainForGeneapt;
bool CurrentCreateGeneralAirport;
bool CurrentLampType;
TerrainInput currTi;
GeneralAirport currGA;
int nodeID;
int hTxtNum;
int hHgtNum;
int cultID;
int bldNum;
int lTxtNum;
int shpID;
int isCurrTxtOrHgt;			//txt = 1; hgt = 2;
int currHighLevel;			//当前的纹理级别，根据目录名称来确定
int updateTTNode;
char cf[64][STRING_STANDARD_LEN];
FILE *fi, *fo, *fLog;
gzFile gzfo;
MapNode mergedTif;						//整合后的图片相关信息
MapNode mergedHgt;						//整合后的高程相关信息

/* 全局路径 */
char hgtroot[_MAX_PATH];
char iveroot[_MAX_PATH];
char objroot[_MAX_PATH];
char tempDirectory[_MAX_PATH];
char workDirectory[_MAX_PATH];
char backDirectory[_MAX_PATH];
char currDirectory[_MAX_PATH];
char findJpgDirectory[_MAX_PATH];
char configFile[_MAX_PATH];
char userConfig[_MAX_PATH];
char temTextFile[_MAX_PATH];
char temTextGeoFile[_MAX_PATH];
char temHgtFile[_MAX_PATH];
char temHgtGeoFile[_MAX_PATH];
char temLogFile[_MAX_PATH];
char baseTextDir[_MAX_PATH];
char baseHgtDir[_MAX_PATH];
char fsData[_MAX_PATH];
char debugFile[_MAX_PATH];

#define MAXFINDNODES 25
std::vector<MapNode> findResult;
int resultID, sortrows;
int sortlist[10][10], sortcolums[10];
unsigned char *newPicPtr = NULL;

#define MAXPICTURES 36
typedef struct _TILEPICS_
{
	int 	tile;
	double	W;
	double	S;
	double	E;
	double	N;
	int		picNum;
	int		pics[MAXPICTURES];
}TilePics;

std::vector<TilePics> currTP;
int currTPID;
int allPics[MAXPICTURES];
int currPics;

/* 获取高程的全局变量*/
// Create root node of the elevation calculation.
#define MAXTOPENEDILES 64
osg::ref_ptr<osg::Group> rootNode;
typedef struct _OPENTILES_
{
	osg::ref_ptr<osg::Node> node;
	int tile;
}OpenTiles;
std::vector<OpenTiles> openedTiles;
TGhgt currHgt;

/* 批量单ive制作地形所需变量及数据结构 */

int TextDirIdx;
std::vector<TextDirInfo> allDirs;
PerZone *currPZ;

int needUpdate;

std::vector<string> dirsNoCultDone;
int  numDNCD;
std::vector<CultNode> findResultAC; 
int  numFRAC;
std::vector<LocateBuilding> findResultFlt;
int numFRFLT;
std::vector<CultNode> fRtoCutIsland; 

/* 调整地景点高度，使得地景与所加机场平滑过渡 */
std::vector< osg::ref_ptr<CorrectionSpaceMesh> > csmList;
	
/* 读取colors.txt，保存在这里 */
std::vector<ColorName> colors;

/* flt机场坐标变换全局变量 */
osg::Matrix coordMatrix;
osg::Matrix rotateMatrix;

/* flt机场轮廓所占区域的极值 */
MaxMinLonLat currFlt;

/* 细节纹理收集 */
std::vector<std::string> allDetailedTxt;
	
/* Shp文件管理 */
std::vector<ShpManagement> allShps;
	
/* 初始化检查项目 */
InitializationCheck initCheck;	

/* 管理所有的文化信息制作项目 */
std::vector<HighCultureInfo> allHCI;

/* 给路面湖面添加裙边 */
osg::ref_ptr<osg::Group> allLakeSkirt;
osg::ref_ptr<osg::Group> allRoadSkirt;
osg::ref_ptr<osg::Group> allDynamicGrass;
osg::ref_ptr<osg::Group> allIslands;
osg::ref_ptr<osg::Group> allLakeInCurrLonLat;


/**************************************************
 *
 *   把 currTi.lod 分解为 3 级 LOD 字符串
 *
 **************************************************/
string getCurrLODString()
{
	string lodstr;
	int iLod = currTi.lod;
	if(iLod < 10)
	{
		iLod += 0x30;
		lodstr = iLod;
		lodstr += " ";
		lodstr += iLod;
		lodstr += " ";
		lodstr += iLod;
	}
	else if ((iLod >= 10) && (iLod < 100))
	{
		int iLodL0 = (iLod / 10) + 0x30;
		int iLodL1 = (iLod % 10) + 0x30;
		lodstr = iLodL0;
		lodstr += " ";
		lodstr += iLodL1;
		lodstr += " ";
		lodstr += iLodL1;
	}else if(iLod >= 100)
	{
		int iLodL0 = (iLod / 100) + 0x30;
		int iLodL1 = (iLod / 10 - (iLodL0-0x30) * 10) + 0x30;
		int iLodL2 = (iLod % 10) + 0x30;
		lodstr = iLodL0;
		lodstr += " ";
		lodstr += iLodL1;
		lodstr += " ";
		lodstr += iLodL2;
	}
	return lodstr;
}

int getCurrLodInt(const char *lodstr)
{
	int j0, j1, d0, j;
	char OneLine[128];
	for(j=0; j<strlen(lodstr);j++)	OneLine[j] = lodstr[j];
	OneLine[j] = 0;
	memset(cf, 0, sizeof(cf));

	j0=j1=d0=0;
	for(unsigned int i=0;i<(int)strlen(OneLine);i++)
	{
		if(!((OneLine[i]==' ')||(OneLine[i]==0x0a)||(OneLine[i]==0x09)))  cf[j0][j1++] = OneLine[i];
		if(((OneLine[i]==' ')||(OneLine[i]==9))&&((OneLine[i-1]!=' ')&&(OneLine[i-1]!=9))&&(i!=0)) {j0++; j1=0;}
	}
	
	unsigned int iLodL0 = 7, iLodL1, iLodL2;
	if(strlen(cf[0]) > 0)	iLodL0 = atoi(cf[0]);
	if(strlen(cf[1]) > 0)	iLodL1 = atoi(cf[1]);	else	iLodL1 = iLodL0;
	if(strlen(cf[2]) > 0)	iLodL2 = atoi(cf[2]);	else	iLodL2 = iLodL1;
	return iLodL0*100 + iLodL1*10 + iLodL2;
}

int getLodLvl0()
{
	int iLod = currTi.lod;
	return iLod / 100;
}

int getLodLvl1()
{
	int iLod = currTi.lod;
	int iLodL0 = iLod / 100;
	return iLod / 10 - iLodL0 * 10;
}

int getLodLvl2()
{
	int iLod = currTi.lod;
	return iLod % 10;
}

#ifdef	USE_HCI
/**************************************************
 *
 *   根据tileID号，从高清地景制作列表中找到相关项目，赋值allHCI
 *
 **************************************************/
void InsertHCI(int tileID)
{
	/* 根据tileID，从allHCI中查找，看是否已经添加过。*/
	for(unsigned int i = 0; i<allHCI.size(); i++)
	{
		if(allHCI[i].tileID == tileID)	return;
	}
	
	/* 从allHTI中查找同tileID的项目。*/
	int iFlag = 0;
	HighTerrainInfo tmpHTI;
	for(unsigned int i = 0; i<allHTI.size(); i++)
	{
		if(allHTI[i].tileID == tileID)
		{
			iFlag = 1;
			strcpy(tmpHTI.texFname, allHTI[i].texFname);  
			strcpy(tmpHTI.hgtFname, allHTI[i].hgtFname);  
			tmpHTI.startLongi = allHTI[i].startLongi	;
			tmpHTI.startLati  = allHTI[i].startLati 	;
			tmpHTI.offsetLongi= allHTI[i].offsetLongi	;
			tmpHTI.offsetLati = allHTI[i].offsetLati 	;
			tmpHTI.highLevel  = allHTI[i].highLevel 	;
			tmpHTI.tileID     = allHTI[i].tileID    	;
			tmpHTI.ID         = allHTI[i].ID        	;
			tmpHTI.isOK       = allHTI[i].isOK      	;
			break;
		}
	}
	
	if(iFlag)
	{
		HighCultureInfo hci;
		strcpy(hci.texFname, tmpHTI.texFname);
		hci.startLongi = tmpHTI.startLongi	;
		hci.startLati  = tmpHTI.startLati 	;
		hci.offsetLongi= tmpHTI.offsetLongi	;
		hci.offsetLati = tmpHTI.offsetLati 	;
		hci.tileID     = tmpHTI.tileID    	;
		hci.ID         = allHCI.size()     	;
		hci.isOK       = 1      			;
		allHCI.push_back(hci);
	}
}

/* 根据tileID号，检查当前文化信息是否做过。*/
int CheckHCI(int tileID)
{
	/* 根据tileID，从allHCI中查找，看是否已经添加过。*/
	int iFlag = 0;
	for(unsigned int i = 0; i<allHCI.size(); i++)
	{
		if(allHCI[i].tileID == tileID)
		{
			if(allHCI[i].isOK == 1)	return 1;
		}
	}
	return iFlag;
}
#endif	//USE_HCI

/* 写配置，把用户高分辨率文化信息制作情况写入配置文件。 */
int WriteHighCultureInfo2Config(FILE *fo)
{
	unsigned int i;
	
	fprintf(fo, "\n%d",  allHCI.size());
	//再写入Map表
	for(i=0; i<allHCI.size(); i++)
	{
		fprintf(fo, "\n%s %11.6f %11.6f %11.6f %11.6f %9d %2d", allHCI[i].texFname, allHCI[i].startLongi, allHCI[i].startLati, allHCI[i].offsetLongi, allHCI[i].offsetLati, allHCI[i].tileID, allHCI[i].isOK);
	}
	return 1;
}

/* 读配置，从配置文件中读取用户高分辨率文化信息制作情况。 */
int ReadHighCultureInfoFromConfig(FILE *fi)
{
	unsigned int i;
	
	ReadALineToArray(fi);	
	int numHCI = atoi(cf[0]);	allHCI.clear();
	for(i=0; i<numHCI; i++)
	{
		HighCultureInfo hci;
		ReadALineToArray(fi);
		strcpy(hci.texFname, cf[0]);
		hci.startLongi = atof(cf[1]);
		hci.startLati  = atof(cf[2]);
		hci.offsetLongi= atof(cf[3]);
		hci.offsetLati = atof(cf[4]);
		hci.tileID 	   = atoi(cf[5]);	
		hci.isOK	   = atoi(cf[6]);
		hci.ID		   = i;
		allHCI.push_back(hci);
	}
	return 1;
}

/**************************************************
 *
 *   将地景配置文件写入Log
 *
 **************************************************/
void WriteTerrainCfgIntoLogFile()
{
	/* 写用户配置数据结构进入Log */
	fprintf(fLog, "\n当前地景生成信息：\n");
	fprintf(fLog, "当前地景生成LOD级别: %d\n", currTi.lod );
	WriteHTxtInfo2Config(fLog);
	WriteHHgtInfo2Config(fLog);
	WriteHighInfo2Config(fLog);
	fprintf(fLog, "\n");
}

void WriteCultureCfgIntoLogFile()
{
	/* 写用户配置数据结构进入Log */
	fprintf(fLog, "\n");
	fprintf(fLog, "\n当前文化信息生成信息：\n");
	WriteShpInfo2Config(fLog);
	WriteBldInfo2Config(fLog);
	WriteHighCultureInfo2Config(fLog);
	fprintf(fLog, "\n");
}

/************************************************
 *
 *    获取文件长度
 *
 ************************************************/

fpos_t GetFileSize(char *fn)
{
	FILE *ft;
	fpos_t iFileSize;
	
		
	ft = fopen(fn, "rb");
	if(ft == NULL)	return 0;

	/*获取文件长度*/
	fseek(ft, 0, SEEK_END);
	fgetpos(ft, &iFileSize);
	fclose(ft);
	
	return iFileSize;
}

/************************************************
 *
 *    把指定经纬度目录下的所有btg文件输出为txt文件，以供分析
 *
 ************************************************/
void ConvertBTGToTXT(int lon, int lat)
{
	Bucket Bt;
	char basePath[STRING_STANDARD_LEN];
	char fullPath[_MAX_PATH];
	set_bucket(&Bt, lon, lat);
	gen_base_path(&Bt, basePath);
	sprintf(fullPath, "%s\\%s", iveroot, basePath);		
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
	
		/* 在这个路径下查找 *.btg.gz */
		io = chdir(fullPath);
		
		dio = _findfirst("*.btg.gz", &sdff);
		if((!dio)||(dio==-1)) goto CBTT_Continue_1;
		dfio = 0;
	
		while(!dfio)
		{
			/* 把这个文化信息btg读入到一个 group里面 */				
			AnalysisBTG(sdff.name);	

			/* 寻找下一个 */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}
CBTT_Continue_1:
	return;
}


/************************************************
 *
 *    从高程文件中读取当前点的高程数据
 *
 ************************************************/

// Return the elevation of the closest non-void grid point to lon, lat
double closest_nonvoid_elev(TGhgt *a, double lon, double lat ) 
{
    int row, col;
    double mindist = 99999999999.9;
    double minelev = -9999.0;
    Geodesy p0, p1;
    Carton  cp0, cp1;

    p0._lon = lon; p0._lat = lat; p0._elevation = 0.0;
    GeodToCart(&p0, &cp0);

    for ( row = 0; row < a->rows; row++ ) {
        for ( col = 0; col < a->cols; col++ ) {
            double dist, elev;

            p1._lon = a->originx + col * a->col_step;
            p1._lat = a->originy + row * a->row_step;
            p1._elevation = 0.0;
            GeodToCart(&p1, &cp1);	

            dist = distance3D( &cp0, &cp1 );
            elev = a->data[col][row];
            if ( dist < mindist && elev > -9000 ) {
                mindist = dist;
                minelev = elev;
                // cout << "dist = " << mindist;
                // cout << "  elev = " << elev << endl;
            }
        }
    }

    if ( minelev > -9999.0 ) {
        return minelev;
    } else {
        return 0.0;
    }
}


// return the current altitude based on grid data.  We should rewrite
// this to interpolate exact values, but for now this is good enough
int xindex, yindex;
double altitude_from_grid(TGhgt *a, double lon, double lat ) 
{
    // we expect incoming (lon,lat) to be in arcsec for now
    double xlocal, ylocal, dx, dy, zA, zB, elev;
    int x1, x2, x3, y1, y2, y3;
    float z1, z2, z3;
 
    /* determine if we are in the lower triangle or the upper triangle 
       ______
       |   /|
       |  / |
       | /  |
       |/   |
       ------

       then calculate our end points
     */

    xlocal = (lon - a->originx) / a->col_step;
    ylocal = (lat - a->originy) / a->row_step;

    xindex = (int)(xlocal);
    yindex = (int)(ylocal);

    // printf("xindex = %d  yindex = %d\n", xindex, yindex);
    if ( xindex + 1 == a->cols ) {
      xindex--;
    }

    if ( yindex + 1 == a->rows ) {
      yindex--;
    }

    if ( (xindex < 0) || (xindex + 1 >= a->cols) ||
      (yindex < 0) || (yindex + 1 >= a->rows) ) {
      // cout << "WARNING: Attempt to interpolate value outside of array!!!" 
      //      << endl;
      return 0.0;
    }

    dx = xlocal - xindex;
    dy = ylocal - yindex;

    if ( dx > dy ) {
      // lower triangle
      // printf("  Lower triangle\n");

      x1 = xindex; 
      y1 = yindex; 
      z1 = a->data[x1][y1];

      x2 = xindex + 1; 
      y2 = yindex; 
      z2 = a->data[x2][y2];
				  
      x3 = xindex + 1; 
      y3 = yindex + 1; 
      z3 = a->data[x3][y3];

      if ( z1 < -9000 || z2 < -9000 || z3 < -9000 ) {
          // don't interpolate off a void
          return closest_nonvoid_elev(a, lon, lat );
      }

      zA = dx * (z2 - z1) + z1;
      zB = dx * (z3 - z1) + z1;
	
      if ( dx > SG_EPSILON ) {
          elev = dy * (zB - zA) / dx + zA;
      } else {
          elev = zA;
      }
    } else {
      // upper triangle
      // printf("  Upper triangle\n");
      
      x1 = xindex; 
      y1 = yindex; 
      z1 = a->data[x1][y1];
      
      x2 = xindex; 
      y2 = yindex + 1; 
      z2 = a->data[x2][y2];
		  		  
      x3 = xindex + 1; 
      y3 = yindex + 1; 
      z3 = a->data[x3][y3];
      
      if ( z1 < -9000 || z2 < -9000 || z3 < -9000 ) {
          // don't interpolate off a void
          return closest_nonvoid_elev(a, lon, lat );
      }
      
      zA = dy * (z2 - z1) + z1;
      zB = dy * (z3 - z1) + z1;
	    
      if ( dy > SG_EPSILON ) {
	        elev = dx * (zB - zA) / dy    + zA;
      } else {
	        elev = zA;
      }
    }
   
   return elev;
}


double GetElevationFromHgt(double lon, double lat)
{
	int iCurrLong, iCurrLat;
	char hgtFilename[_MAX_PATH], hgtPathname[_MAX_PATH];
	char fullPath[_MAX_PATH];
	double clon, clat, cElev;
	
	/* 获得当前位置高程文件路径名和文件名 */
	iCurrLong = (int) lon;
	iCurrLat  = (int) lat;

	sprintf(hgtFilename, "N%02dE%03d.hgt", iCurrLat, iCurrLong);
	sprintf(hgtPathname, "%s\\n%02d", hgtroot, iCurrLat);
	
	/* 打开并读取高程文件 */
	TGhgt * curHgt = &currHgt;
	if((curHgt->curLong != iCurrLong) || (curHgt->curLati != iCurrLat))
	{
		sprintf(fullPath, "%s\\%s", hgtPathname, hgtFilename);
		if(!openhgt(curHgt, fullPath))
		{
			return 1.0;
		}
	}
	
	/* 获取高程数据 */
	clon = lon * 3600.0; clat = lat * 3600.0;
	cElev = altitude_from_grid(curHgt, clon, clat);
	
	return cElev;
}
	


/************************************************
 *
 *    从一个文件中读取colors所有名称信息
 *
 ************************************************/
void ReadColorsFromFile()
{
	FILE *fc;
	char fname[STRING_STANDARD_LEN * 2];
	
	sprintf(fname, "%s\\colors.txt", workDirectory);
	fc = fopen(fname, "rt");
	if(fc == NULL) return;
	
	while(!feof(fc))
	{
		ColorName cn;

		ReadALineToArray(fc);
		cn.name = cf[0];
		cn.color.set(atof(cf[1]), atof(cf[2]), atof(cf[3]), atof(cf[4]));
		colors.push_back(cn);
	}
	fclose(fc);
}

/* 根据一个color的颜色值，在color表中查找它的位置 */
int FindIdxByColor(osg::Vec4 color)
{
	int retVal = -1;
	
	for(int i = 0; i<colors.size(); i++)
	{
		if(	(fabs(colors[i].color.x() - color.x()) < 0.0001) &&
			(fabs(colors[i].color.y() - color.y()) < 0.0001) &&
			(fabs(colors[i].color.z() - color.z()) < 0.0001) &&
			(fabs(colors[i].color.w() - color.w()) < 0.0001) )
			{	retVal = i; break;	}
	}
	return retVal;
}


/************************************************
 *
 *    从一个totalShp里面，根据不同的整数经纬度分出一个或多个findResultAC
 *    用于处理文化信息跨越经纬度的情况。
 *
 ************************************************/
int SplitByLonLatShp(int idx)
{
	int minLon, maxLon, minLat, maxLat, i, j;

	/* 先找出findResultAC所占据的最大最小经纬度整数值 */
	minLon = (int)totalShp[idx].minLong;
	maxLon = (int)totalShp[idx].maxLong;
	minLat = (int)totalShp[idx].minLati;
	maxLat = (int)totalShp[idx].maxLati;
	
	/* 做双重循环，把所有占据情况都计入findAllItemsAC */
	for(i=minLon; i<=maxLon; i++)
	{
		for(j=minLat; j<=maxLat; j++)
		{
			CultNode cn;
			strcpy(cn.fileName, totalShp[idx].fileName);
			strcpy(cn.pathName, totalShp[idx].pathName);
			if((float)i		< totalShp[idx].minLong)	cn.minLong = totalShp[idx].minLong; else cn.minLong = (double)(i+0.0001);
			if((float)(i+1) > totalShp[idx].maxLong)	cn.maxLong = totalShp[idx].maxLong; else cn.maxLong = (double)(i+0.9999);
			if((float)j		< totalShp[idx].minLati)	cn.minLati = totalShp[idx].minLati; else cn.minLati = (double)(j+0.0001);
			if((float)(j+1) > totalShp[idx].maxLati)	cn.maxLati = totalShp[idx].maxLati; else cn.maxLati = (double)(j+0.9999);
			
			/* 检查一下，当前块是否值得做 。*/
			if(fabs(cn.minLong - cn.maxLong) < 0.001) continue;
			if(fabs(cn.minLati - cn.maxLati) < 0.001) continue;
			
			/* 如果值得做，才增加需要制作的项目 。*/
			strcpy(cn.destFile, totalShp[idx].destFile);
			strcpy(cn.destPath, totalShp[idx].destPath);
			cn.destLongi    = i; 
			cn.destLati     = j;
			cn.destTileID   = totalShp[idx].destTileID ;
			cn.ID			= totalShp[idx].ID;
			cn.isOK			= totalShp[idx].isOK;
			findResultAC.push_back(cn);		numFRAC++;
		}
	}
	return 1;	
}

/************************************************
 *
 *    从一个totalcult里面，根据不同的整数经纬度分出一个或多个findResultAC
 *    用于处理文化信息跨越经纬度的情况。
 *
 ************************************************/
int SplitByLonLatAC(int idx)
{
	int minLon, maxLon, minLat, maxLat, i, j;

	/* 先找出findResultAC所占据的最大最小经纬度整数值 */
	minLon = (int)totalcult[idx].minLong;
	maxLon = (int)totalcult[idx].maxLong;
	minLat = (int)totalcult[idx].minLati;
	maxLat = (int)totalcult[idx].maxLati;
	
	/* 做双重循环，把所有占据情况都计入findAllItemsAC */
	for(i=minLon; i<=maxLon; i++)
	{
		for(j=minLat; j<=maxLat; j++)
		{
			CultNode cn;
			strcpy(cn.fileName, totalcult[idx].fileName);
			strcpy(cn.pathName, totalcult[idx].pathName);
			if((float)i		< totalcult[idx].minLong)	cn.minLong = totalcult[idx].minLong; else cn.minLong = (double)(i+0.0001);
			if((float)(i+1) > totalcult[idx].maxLong)	cn.maxLong = totalcult[idx].maxLong; else cn.maxLong = (double)(i+0.9999);
			if((float)j		< totalcult[idx].minLati)	cn.minLati = totalcult[idx].minLati; else cn.minLati = (double)(j+0.0001);
			if((float)(j+1)	> totalcult[idx].maxLati)	cn.maxLati = totalcult[idx].maxLati; else cn.maxLati = (double)(j+0.9999);
			
			/* 检查一下，当前块是否值得做 。*/
			if(fabs(cn.minLong - cn.maxLong) < 0.001) continue;
			if(fabs(cn.minLati - cn.maxLati) < 0.001) continue;
			
			/* 如果值得做，才增加需要制作的项目 。*/
			strcpy(cn.destFile, totalcult[idx].destFile);
			strcpy(cn.destPath, totalcult[idx].destPath);
			cn.destLongi    = i; 
			cn.destLati     = j;
			cn.destTileID   = totalcult[idx].destTileID ;
			cn.ID			= totalcult[idx].ID;
			cn.isOK			= totalcult[idx].isOK;
			findResultAC.push_back(cn);		numFRAC++;
		}
	}
	return 1;	
}

/************************************************
 *
 *    设置机场代码，根据findResultFlt前面已经设过的代码情况，设置当前代码，不能和前面的重复了。
 *
 ************************************************/
int SetCurrFltAPCode(int idx)
{
	int iFlag = 0, iCode;
	if(idx == 0) { findResultFlt[idx].iCode = 'A'; return 1; }
	iCode = 'A';
	while(1)
	{
		for(int i = 0; i < numFRFLT; i++)
		{
			if(iCode == findResultFlt[i].iCode)
			{	iFlag = 1; break; }		
		}
		if(!iFlag) break;
		if( iFlag){	iCode++; iFlag = 0;	}
	}
	findResultFlt[idx].iCode = iCode;
	return 1;
}

/************************************************
 *
 *    根据tile在findResultAC里面查找，是否有相同的tile存在，直到确定一个独一无二的。
 *
 ************************************************/
int findUniqueTile(int tile)
{
	int iFlag = 0, iTile;
	iTile = tile;
	while(1)
	{
		for(int i = 0; i < numFRAC; i++)
		{
			if(iTile == findResultAC[i].destTileID)
			{	iFlag = 1; break; }		
		}
		if(!iFlag) break;
		if( iFlag){	iTile++; iFlag = 0;	}
	}
	return iTile;
}

/************************************************
 *
 *    获取环境变量Path，用来添加软件自己的环境变量
 *
 ************************************************/
static void CreateEnvPath()
{
	std::vector<string> paths;
	char onepath[STRING_STANDARD_LEN];
	int i, j = 0;
	
	/* 获取当前环境变量PATH，放在paths列表里面 */
	paths.clear();
    char *envp = ::getenv( "PATH" );
    if( envp != NULL ) 
    {
    	for(i=0; i<strlen(envp); i++)
    	{
    		if(envp[i] == ';')	
    		{ 
    			onepath[j] = 0; j = 0;
    			string path = onepath;
    			paths.push_back(path);
    		}else{
    			onepath[j++] = envp[i];
    		}
    	}
		onepath[j] = 0; j = 0;
		string path = onepath;
		paths.push_back(path);    	
    }
    
    /* paths列表与当前工作目录逐一比较，如果有相同的就退出，如果没有相同的就继续。 */
    string currwork = workDirectory;
    for(i=0; i<paths.size(); i++)
    {
    	if(paths[i] == currwork)
    	{
    		return;
    	}
    }
    
    /* 将当前工作目录加入到列表中，构造新的path，写入环境变量里面 */
	char newpath[2048];
	sprintf(newpath, "set path=%s;%s", envp, workDirectory);
	::putenv(newpath);
}

/************************************************
 *
 *    备份或者恢复制作文化信息的地景库文件
 *
 ************************************************/

/* 复制文件，把src复制到dst里面 */
int CopyFile(char *src, char *dst)
{
	FILE *pInput;
	FILE *pOutput;
	fpos_t iFileSize, iDataLen;

	/* 打开文件 */	
	if((pInput=fopen(src,"rb"))==NULL||						
		(pOutput=fopen(dst,"wb"))==NULL)
	{
		return 0;
	}

	/*获取文件长度*/
	fseek(pInput, 0, SEEK_END);
	fgetpos(pInput, &iFileSize);
	rewind(pInput);

	/*读写文件体*/
	unsigned char * pcInputLine= new unsigned char[MBUFFSIZE];
	iDataLen = iFileSize;
	while(iDataLen > MBUFFSIZE)
	{
		fread(pcInputLine, 1, MBUFFSIZE, pInput);
		fwrite(pcInputLine, 1, MBUFFSIZE, pOutput);
		iDataLen -= MBUFFSIZE;
	}
	fread(pcInputLine, 1, iDataLen, pInput);
	fwrite(pcInputLine, 1, iDataLen, pOutput);
	
	/* 关闭文件 */
	delete pcInputLine;
	fclose(pInput);
	fclose(pOutput);
	return 1;
}

/* 把当前文件从当前目录复制到Temp目录下，并根据文件名创建和复制子目录下的相关文件 */
int BackUpFileToTemp(char *fn, char *fnpath)
{
	char fullSubPath[_MAX_PATH], currDir[_MAX_PATH], backSubDir[_MAX_PATH];
	char dstFile[_MAX_PATH];
	char subDir[_MAX_PATH], currRoot[_MAX_PATH];
	struct _finddata_t iveSubf;

	/* 首先备份当前根路径, 然后进入当前参数的路径 */
	if(strlen(fnpath) != 0)
	{
		_getcwd(currRoot, _MAX_PATH);
		int ioCur = chdir(fnpath);
		if(ioCur) 
		{
			char msg[STRING_STANDARD_LEN];
			sprintf(msg, "不能进入目标目录 -- %s", fnpath);
			MessageToDialog(msg);
			return 0;
		}
	}
	
	/* 复制顶层 */
	sprintf(dstFile, "%s\\%s", backDirectory, fn);
	
	/* 检查一下备份文件在不在，如果已经备份过了，就不用再备份了 --- 暂不需要 */
	if(!CopyFile(fn, dstFile))	return 0;
	
	/* 以下是处理PageLOD的底层的各个ive文件，所有底层文件存放在 XXXX_root_L0_X0_Y0 目录下 */
	strcpy(subDir, fn);
	int len = strlen(subDir);
	subDir[len - 4] = 0;
	sprintf(fullSubPath, "%s_root_L0_X0_Y0", subDir);

	/* 在备份路径下创建子目录 */
	_getcwd(currDir, _MAX_PATH);
	int io = chdir(backDirectory);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", backDirectory);
		MessageToDialog(msg);
		return 0;
	}
	
	io = chdir(fullSubPath);
	if(io) 
	{
		int dio = mkdir(fullSubPath);
		if(dio)
		{
			MessageToDialog("不能创建地景子目录。\n");
			return 0;
		}
	}
	sprintf(backSubDir, "%s\\%s", backDirectory, fullSubPath);

	/* 返回源数据目录并进入到子目录中 */
	io = chdir(currDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", currDir);
		MessageToDialog(msg);
		return 0;
	}

	io = chdir(fullSubPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", fullSubPath);
		MessageToDialog(msg);
		return 0;
	}

	/* 在子目录中查找ive文件 */
	int subffio = _findfirst("*.ive", &iveSubf);	
	if((!subffio)||(subffio == -1))
	{
		MessageToDialog("找不到ive文件。");
		return 0;
	}

	int subfio = 0, subidx = 0;
	while(!subfio)
	{
		subidx++;
		sprintf(dstFile, "%s\\%s", backSubDir, iveSubf.name);
		
		/* 调用复制函数 */
		if(!CopyFile(iveSubf.name, dstFile))	return 0;

		/* 查找下一个符合条件的文件 */
		subfio = _findnext(subffio, &iveSubf);
	}
	_findclose(subffio);

	/* 退出前返回源数据目录 */
	io = chdir(currDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", currDir);
		MessageToDialog(msg);
		return 0;
	}
	
	/* 所有工作完成之后，返回根路径 */
	if(strlen(fnpath) != 0)
	{
		int ioCur = chdir(currRoot);
		if(ioCur) 
		{
			char msg[STRING_STANDARD_LEN];
			sprintf(msg, "不能进入目标目录 -- %s", currRoot);
			MessageToDialog(msg);
			return 0;
		}
	}
	
	return 1;
}

/* 从Temp目录恢复文件到当前目录下，并根据文件名恢复子目录下的相关文件 */
int RestoreFileFromTemp(char *fn)
{
	char fullSubPath[_MAX_PATH], currDir[_MAX_PATH], backSubDir[_MAX_PATH];
	char dstFile[_MAX_PATH];
	char subDir[24];
	struct _finddata_t iveSubf;

	sprintf(dstFile, "%s\\%s", backDirectory, fn);

	/* 先比较顶层的两个文件长度,如果长度相等,说明该组文件没变化,那就没必要再复制了.*/
	fpos_t srcsize, tgtsize;
	if((srcsize = GetFileSize(fn)) == 0) return 0;
	if((tgtsize = GetFileSize(dstFile)) == 0) return 2;
	if(srcsize == tgtsize) return 0;

	/* 复制顶层 */
	if(!CopyFile(dstFile, fn))	return 2;
	
	/* 以下是处理PageLOD的底层的各个ive文件，所有底层文件存放在 XXXX_root_L9_X0_Y0 目录下 */
	strcpy(subDir, fn);
	int len = strlen(subDir);
	subDir[len - 4] = 0;
	sprintf(fullSubPath, "%s_root_L0_X0_Y0", subDir);
	sprintf(backSubDir, "%s\\%s", backDirectory, fullSubPath);

	/* 在备份路径下创建子目录 */
	_getcwd(currDir, _MAX_PATH);
	int io = chdir(backSubDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", backSubDir);
		MessageToDialog(msg);
		return 0;
	}
	
	/* 返回源数据目录并进入到子目录中 */
	io = chdir(currDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", currDir);
		MessageToDialog(msg);
		return 0;
	}

	io = chdir(fullSubPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", fullSubPath);
		MessageToDialog(msg);
		return 0;
	}

	/* 在子目录中查找ive文件 */
	int subffio = _findfirst("*.ive", &iveSubf);	
	if((!subffio)||(subffio == -1))
	{
		MessageToDialog("找不到ive文件。");
		return 0;
	}

	int subfio = 0, subidx = 0;
	while(!subfio)
	{
		subidx++;
		sprintf(dstFile, "%s\\%s", backSubDir, iveSubf.name);
		
		/* 调用复制函数 */
		if(!CopyFile(dstFile, iveSubf.name))	return 0;

		/* 查找下一个符合条件的文件 */
		subfio = _findnext(subffio, &iveSubf);
	}
	_findclose(subffio);

	/* 退出前返回源数据目录 */
	io = chdir(currDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", currDir);
		MessageToDialog(msg);
		return 0;
	}
	return 1;
}

/* 根据文件长度来决定备份还是恢复，文件长度在FILELENMIN和FILELENMAX之间，则备份，否则恢复。 */
#define FILELENMAX 717300
#define FILELENMIN 717200

/* 在恢复地景文件之前，先检查一下备份过的ive根文件长度，
   如果备份文件不存在，说明还没有被备份。
   如果本文件长度没变化并且有备份文件，说明文件没有被改动。则不需要做备份或者恢复操作
   文件长度没变化则返回 0， 不存在和有变化则返回 1；
  */
int CheckFileLengthBeforeRestore(int tID)
{
	char backIveFile[_MAX_PATH];
	char fullIveFile[_MAX_PATH];
	fpos_t iFileSizeSrc, iFileSizeBkp;
	int iFlag = 1;

	/* 读取当前的一个ive文件。 */
	sprintf(fullIveFile, "%d.ive", tID);
	iFileSizeSrc = GetFileSize(fullIveFile);

	/* 读取备份过的一个ive文件。 */
	sprintf(backIveFile, "%s\\%d.ive",backDirectory, tID);
	iFileSizeBkp = GetFileSize(backIveFile);
	if(iFileSizeBkp == 0)	return iFlag;

	/* 比较两个文件的长度 */
	if(iFileSizeSrc == iFileSizeBkp) iFlag = 0;

	return iFlag;
}


/************************************************
 *
 *    读入一个tile.btg.gz，然后只把里面的cultures部分写入一个新的Btg里面，文件名 A+tile.btg.gz
 *
 ************************************************/
int CreateCultureBtgFromSrc(int idx)
{
	char fn[STRING_STANDARD_LEN], afn[STRING_STANDARD_LEN], msg[_MAX_PATH];

	char fullPathFile[_MAX_PATH];
	sprintf(fullPathFile, "%s\\%s", findResultAC[idx].destPath, findResultAC[idx].destFile);
	strcpy(fn, findResultAC[idx].destFile);

	/* 先把tile.btg.gz写入各个stg中 */
	if(WriteDataStructToBTG(fullPathFile))
	{
		if(!WriteAllStgs(fn, 1)) return 0;
	}else{
		sprintf(msg, "打开文件失败 --- %s", fullPathFile);
		MessageToDialog(msg);
		return 0;
	}

	/* 根据tile.btg.gz，去掉非文化信息部分，存成A+tile.btg.gz*/
	sprintf(fullPathFile, "%s\\A%s", findResultAC[idx].destPath, findResultAC[idx].destFile);
	sprintf(afn, "A%s", fn);
	if(WriteDataStructToBTG_CulturePlusIsland(fullPathFile))
	{
		/* 把A+tile.btg.gz也写入各个stg中 */
		if(!WriteAllStgs(afn, 1)) return 0;
	}
	return 1;
}

int DeleteIslandFromCultureBTG()
{
	char fullPathFile[_MAX_PATH];
	int i;
	
	for(i=0; i<fRtoCutIsland.size(); i++)
	{
		if(fRtoCutIsland[i].isOK = 2)
		{
			fRtoCutIsland[i].isOK = 1;
			
			/* 根据tile.btg.gz，去掉非文化信息部分，存成A+tile.btg.gz*/
			sprintf(fullPathFile, "%s\\A%s", fRtoCutIsland[i].destPath, fRtoCutIsland[i].destFile);  

			/* 先读入这个文件到内存. */
			if(!ReadBTGToDataStruct(fullPathFile)) continue;

			/* 再删除掉这个文件， */
			DeleteFile(fullPathFile);

			/* 然后再根据内存中的数据重新生成*/
			if(!WriteDataStructToBTG_OnlyCulture(fullPathFile)) continue;
		}
	}
	fRtoCutIsland.clear();

	return 1;
}

/************************************************
 *
 *    管理文化信息挖面次序
 *
 ************************************************/
typedef struct _CUTLIST_
{
	int lon;
	int lat;
}CutList;
std::vector<CutList> allNeedCut;			//统一挖洞的列表，每一项表示在一个经纬度区间内挖洞。

/* 把当前的经纬度组合插入到CutList中，如果里面有相同的值，则不再重复插入 */
int InsertItemIntoCutList(int lon, int lat)
{
	int iFlag = 0;
	for(std::vector<CutList>::iterator p = allNeedCut.begin(); p != allNeedCut.end(); p++)
	{
		if((p->lon == lon)&&(p->lat == lat))
		{
			iFlag = 1; break;
		}
	}
	if(!iFlag)
	{
		CutList cl; cl.lon = lon; cl.lat = lat; allNeedCut.push_back(cl);
	}
	return 1;
}

/* 根据CutList的内容，对地景块进行挖洞操作 */
int DoTerrainCut()
{
	char msg[STRING_STANDARD_LEN];

	for(std::vector<CutList>::iterator p = allNeedCut.begin(); p != allNeedCut.end(); p++)
	{
		/* 先把当前经纬度目录下的btg输出为txt */
		ConvertBTGToTXT(p->lon, p->lat);
		
		/* 将当前块（1个经度x1个纬度所决定的空间）内的所有的文化信息Btg文件，然后加上通用机场和用户机场的轮廓ive文件，共同组成一个组Group，根据这个Group对地景进行挖洞。 */
		if(!ReadTerrainIveAndCut(p->lon, p->lat))
		{   
			sprintf(msg, "转换IVE失败。区域：经度 %d 纬度 %d", p->lon, p->lat);
			MessageToDialog(msg);
		}
	}
	allNeedCut.clear();
	return 1;
}

/************************************************
 *
 *    管理用户自定义三维物体列表
 *
 ************************************************/
/* 从totalUserBld库里面，找出符合dirsNoCultDone里面一条terPname的所有记录，包括做好的和没做好的 */
int FindAllFltItemsWithSameDir(int idx)
{
	int i;
	
	if(idx >= numDNCD) return 0;
	findResultFlt.clear();	numFRFLT = 0;
	for(i=0; i<bldNum; i++)
	{
		if(!strcmp(dirsNoCultDone[idx].c_str(), totalUserBld[i].terPname))
		{
			LocateBuilding lb;
			strcpy(lb.fltPname, totalUserBld[i].fltPname);
			strcpy(lb.fltFname, totalUserBld[i].fltFname);
			strcpy(lb.locFname, totalUserBld[i].locFname);
			strcpy(lb.terPname, totalUserBld[i].terPname);
			lb.fLongi       = totalUserBld[i].fLongi     ;
			lb.fLati        = totalUserBld[i].fLati      ; 
			lb.fElevation   = totalUserBld[i].fElevation ;
			lb.fHeading     = totalUserBld[i].fHeading   ;
			lb.iCode        = totalUserBld[i].iCode      ;
			lb.ID           = totalUserBld[i].ID         ;
			lb.isOK         = totalUserBld[i].isOK       ;
			findResultFlt.push_back(lb);	numFRFLT++;
		}
	}
	return numFRFLT;
}

/* 找出所有未定位用户三维物体的地景目录集合 */
int FindOutDirsNoBldDone()
{
	int i, j, iFlag;
	
	dirsNoCultDone.clear(); numDNCD = 0;
	for(i=0; i<bldNum; i++)
	{
		if((totalUserBld[i].isOK == 0) && (strcmp("NULL", totalUserBld[i].terPname)))
		{
			iFlag = 0;
			for(j=0; j<numDNCD; j++)
			{
				if(!strcmp(dirsNoCultDone[j].c_str(), totalUserBld[i].terPname))
				{
					iFlag = 1; break;
				}
			}
			if(!iFlag)
			{
				dirsNoCultDone.push_back(totalUserBld[i].terPname); numDNCD++;
			}
		}
	}
	return numDNCD;
}

/* 已经定位完成的三维物体更新到所有三维物体列表中 */
int UpdateTotalUserBldWithNewFindResult()
{
	int i, idx;
	for(i=0;i<numFRFLT;i++)
	{
		idx = findResultFlt[i].ID;
		strcpy(totalUserBld[idx].fltPname, findResultFlt[i].fltPname);
		strcpy(totalUserBld[idx].fltFname, findResultFlt[i].fltFname);
		strcpy(totalUserBld[idx].locFname, findResultFlt[i].locFname);
		strcpy(totalUserBld[idx].terPname, findResultFlt[i].terPname);
		totalUserBld[idx].fLongi       = findResultFlt[i].fLongi     ;
		totalUserBld[idx].fLati        = findResultFlt[i].fLati      ; 
		totalUserBld[idx].fElevation   = findResultFlt[i].fElevation ;
		totalUserBld[idx].fHeading     = findResultFlt[i].fHeading   ;
		totalUserBld[idx].iCode        = findResultFlt[i].iCode      ;
		totalUserBld[idx].ID           = findResultFlt[i].ID         ;
		totalUserBld[idx].isOK         = findResultFlt[i].isOK       ;
	}
	return 1;
}

/* 写配置，把所有三维物体列表写入配置文件 */
int WriteBldInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  bldNum);
	//再写入表
	for(i=0; i<bldNum; i++)
	{
		fprintf(fo, "\n%s %s %s %s %13.8f %13.8f %13.8f %13.8f %d %5d %1d", totalUserBld[i].fltPname, totalUserBld[i].fltFname, totalUserBld[i].locFname, totalUserBld[i].terPname, totalUserBld[i].fLongi, totalUserBld[i].fLati, totalUserBld[i].fElevation, totalUserBld[i].fHeading, totalUserBld[i].iCode, totalUserBld[i].ID, totalUserBld[i].isOK);
	}

	return 1;
}

/* 读配置，从配置文件中读取所有三维物体列表 */
int ReadBldInfoFromConfig(FILE *fi)
{
	int i, blen;
	
	blen = strlen(baseTextDir);
	ReadALineToArray(fi);	bldNum = atoi(cf[0]);	totalUserBld.clear();
	for(i=0; i<bldNum; i++)
	{
		LocateBuilding lb;
		ReadALineToArray(fi);
		strcpy(lb.fltPname, cf[0]);
		strcpy(lb.fltFname, cf[1]);
		strcpy(lb.locFname, cf[2]);
		strcpy(lb.terPname, cf[3]);
		lb.fLongi       = atof(cf[4]);
		lb.fLati        = atof(cf[5]); 
		lb.fElevation   = atof(cf[6]);
		lb.fHeading     = atof(cf[7]);
		lb.iCode        = atof(cf[8]);
		lb.ID           = i;
		lb.isOK         = atof(cf[10]);
		totalUserBld.push_back(lb);
	}
	return 1;
}


/************************************************
 *
 *    增加按需生成的功能
 *
 ************************************************/

typedef struct _TERRAINONDEMAND_
{
	float maxLongi;
	float maxLati;
	float minLongi;
	float minLati;
}TerrainOnDemand;

TerrainOnDemand inputParam;

/* 初始化按需生成输入参数*/
int InitOnDemand()
{
	inputParam.maxLongi = 0.0;
	inputParam.maxLati  = 0.0; 
	inputParam.minLongi = 0.0;
	inputParam.minLati  = 0.0; 
	return 1;
}

/* 设置按需生成的区域 */
int SetOnDemand()
{
	inputParam.maxLongi = currTi.maxLongi;
	inputParam.maxLati  = currTi.maxLati ; 
	inputParam.minLongi = currTi.minLongi;
	inputParam.minLati  = currTi.minLati ;

	int i, j, lon, lat;
	int maxlon, maxlat, minlon, minlat;
	maxlon = (int)inputParam.maxLongi; maxlat = (int)inputParam.maxLati; minlon = (int)inputParam.minLongi; minlat = (int)inputParam.minLati;
	for(i=0; i<TextDirIdx; i++)
	{
		for(j=0; j<4; j++)
		{
			lon = allDirs[i].zone[j].lon; lat = allDirs[i].zone[j].lat;
			if((lon >= minlon) && (lon <= maxlon) && (lat >= minlat) && (lat <= maxlat))
			{
				allDirs[i].zone[j].isOK = STATUS_UNWORK;
			}
		}
	}
	 
	return 1;
}

/* 判断当前所作内容是否在按需生成范围内 */
int CheckOnDemand()
{
	if(	(inputParam.maxLongi == 0.0) &&
		(inputParam.maxLati  == 0.0) &&
		(inputParam.minLongi == 0.0) &&
		(inputParam.minLati  == 0.0) )
		return 1;
	
	if(	(inputParam.maxLongi >= currTi.maxLongi) &&
		(inputParam.maxLati  >= currTi.maxLati ) &&
		(inputParam.minLongi <= currTi.minLongi) &&
		(inputParam.minLati  <= currTi.minLati ) )
		return 1;
	else
		return 0;
}      

/************************************************
 *
 *  从磁盘里面找到特定的目录名。
 *  目前视景系统运行需要的特定目录如下：
 *  BaseTexture: 里面包含地景纹理源数据，基本源数据分为两种：一种是JPEG+Map数据，一种是GeoTiff数据，在User\50m目录下。
 *               用户高分辨率地景放在下面的User子目录下。用户应将数据复制到 User\1m  User\5m  User\10m 目录下
 *  BaseHgt:     里面包含高程源数据。用户的高分辨率高程数据应复制到 User 下面。
 *
 ************************************************/
char PathFound[_MAX_PATH];
int  PathIsFound;

/* 在子目录里面查找pathname名字。*/
int TravelSubDirectoryForPathname(char *dir,  char *pathname, int level)
{
	int dio, dfio, dirIo;
	struct _finddata_t sdff;

	//如果已经找到了，就不需要再找了
	if(PathIsFound) return 1;

	//查找当前路径下的子目录；
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//如果不是子目录也不是文件
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR_PATHNAME;

		//如果是子目录，进入当前子目录
		if(sdff.attrib == _A_SUBDIR)
		{
			//如果找到了。
			if((!strcmp(pathname, sdff.name)))
			{
				if(!strcmp(pathname, sdff.name))
				{
					_getcwd(path, _MAX_PATH);
					if(path[strlen(path) - 1] == '\\')
						sprintf(PathFound, "%s%s", path, sdff.name);
					else
						sprintf(PathFound, "%s\\%s", path, sdff.name);
					PathIsFound = 1;
					return 1;
				}
			}else{	
				//进入子目录
				if(level >= 4)
				{

				}
				else
				{
					dirIo = chdir(sdff.name);
					if(!dirIo)
					{
						_getcwd(path, _MAX_PATH);
				
						//遍历子目录
						level++;
						TravelSubDirectoryForPathname(path, pathname, level);

						//退回到上一层目录
						dirIo = chdir("..");
						level--;
					}
				}
			}
		}

GET_NEXT_SUBDIR_PATHNAME:
		//寻找下一个
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* 从磁盘上查找特定目录 */
int FindSpecDirFromDisk(char *dirname)
{
	int i, drv, dio;
	char path[_MAX_PATH];
	
	memset(PathFound, 0, sizeof(PathFound)); PathIsFound = 0;
	//先进入到输出根目录
	drv = 'C' - 'A' + 1;
	for(i=0; i<8; i++)
	{
		//如果已经找到了，就不需要再找了
		if(PathIsFound) break;
		
		dio = _chdrive(drv);
		if(!dio)
		{
			sprintf(path, "%c:\\", drv + 'A' - 1);
			dio = chdir(path);
			
			//遍历子目录
			sprintf(path, "%c:", drv + 'A' - 1);
			TravelSubDirectoryForPathname(path, dirname, 1);
		}
		drv++;
	}
	
	return 1;
}

/* 在工作目录下读取目录名称文件，如果读不到，就从磁盘中查找，如果查找到，就写这个目录名称文件
   如果查找不到，返回0。
 */
int FindBasePathName()
{
	FILE *fpath;
	int io, isValid = 1;

	/* 首先要确定在工作目录下*/
	io = chdir(workDirectory);
	if(io)  return 0;
	
	/* 读取PathName文件 */
	fpath = fopen("secret.pth", "rt");
	if(fpath != NULL )
	{
		fgets(baseTextDir, _MAX_PATH, fpath);
		fgets(baseHgtDir,  _MAX_PATH, fpath);
		fgets(fsData,  _MAX_PATH, fpath);
		fclose(fpath);
		
		/* 检验这个路径是否有效 */
		int len = strlen(baseTextDir);
		if(baseTextDir[len - 1] == 0xa) baseTextDir[len - 1] = 0;
		len = strlen(baseHgtDir);
		if(baseHgtDir[len - 1] == 0xa) baseHgtDir[len - 1] = 0;
		len = strlen(fsData);
		if(fsData[len - 1] == 0xa) fsData[len - 1] = 0;
		io = chdir(baseTextDir);
		if(io)  isValid = 0;
		io = chdir(baseHgtDir);
		if(io)  isValid = 0;
		io = chdir(fsData);
		if(io)  isValid = 0;
		if(isValid)	return 1;
	}
	
	/* 如果没找到这个文件 */
	int iFlag = 1;
	memset(baseTextDir, 0, sizeof(baseTextDir)); memset(baseHgtDir, 0, sizeof(baseHgtDir));	memset(fsData, 0, sizeof(fsData));
	FindSpecDirFromDisk("BaseTexture");
	if(PathIsFound) {	strcpy(baseTextDir, PathFound);}	else {	iFlag *= 0;	initCheck.pathOfTexture = false;	}
	FindSpecDirFromDisk("BaseHgt");
	if(PathIsFound) {	strcpy(baseHgtDir,  PathFound);}	else {	iFlag *= 0;	initCheck.pathOfHgt = false;	}
	FindSpecDirFromDisk("FG");
	if(PathIsFound) {	strcpy(fsData,		PathFound);}	else {	iFlag *= 0;	initCheck.pathOfTarget = false;	}
	ifFoundBasePath = iFlag;
	
	/* 如果没能找到上述三个目录中任意一个，则退出。*/
	if(!iFlag)	return 0;
	
	/* 首先要确定回到工作目录下*/
	io = chdir(workDirectory);
	if(io)  return 0;

	/* 如果两个目录都找到了，写文件返回*/
	int textlen = strlen(baseTextDir);
	int hgtlen  = strlen(baseHgtDir);
	int fslen   = strlen(fsData);
	if(textlen && hgtlen && fslen)
	{
		fpath = fopen("secret.pth", "wt");
		if(fpath != NULL )
		{
			fprintf(fpath, "%s\n", baseTextDir);
			fprintf(fpath, "%s\n", baseHgtDir);
			fprintf(fpath, "%s\n", fsData);
			fclose(fpath);
			return 1;
		}	
	}
	
	return 0;
}


/************************************************
 *
 *    找出所有未制作文化信息的shp文件的目录集合
 *
 ************************************************/
int FindOutDirsNoShpDone()
{
	int i, j, iFlag;
	
	dirsNoCultDone.clear(); numDNCD = 0;
	for(i=0; i<shpID; i++)
	{
		if(totalShp[i].isOK == 0)
		{
			iFlag = 0;
			for(j=0; j<numDNCD; j++)
			{
				if(!strcmp(dirsNoCultDone[j].c_str(), totalShp[i].pathName))
				{
					iFlag = 1; break;
				}
			}
			if(!iFlag)
			{
				dirsNoCultDone.push_back(totalShp[i].pathName); numDNCD++;
			}
		}
	}
	return numDNCD;
}

/* 从totalShp库里面，找出符合dirsNoCultDone里面一条pathName的所有记录，包括做好的和没做好的 */
int FindAllShpItemsWithSameDir(int idx)
{
	int i;
	
	if(idx >= numDNCD) return 0;
	findResultAC.clear();	numFRAC = 0;
	for(i=0; i<shpID; i++)
	{
		if(!strcmp(dirsNoCultDone[idx].c_str(), totalShp[i].pathName))
		{
			/* 如果一个shp文化信息跨越多个经纬度，则按每个经纬度细分一下。*/
			SplitByLonLatShp(i);
		}
	}
	return numFRAC;
}

/* 已经制作完成的shp文化信息更新到所有shp列表中 */
int UpdatetotalShpWithNewFindResult()
{
	int i, idx;
	for(i=0;i<numFRAC;i++)
	{
		idx = findResultAC[i].ID;
		strcpy(totalShp[idx].fileName, findResultAC[i].fileName);
		strcpy(totalShp[idx].pathName, findResultAC[i].pathName);
		totalShp[idx].minLong			= findResultAC[i].minLong;
		totalShp[idx].maxLong			= findResultAC[i].maxLong;
		totalShp[idx].minLati			= findResultAC[i].minLati;	
		totalShp[idx].maxLati			= findResultAC[i].maxLati;
		strcpy(totalShp[idx].destFile, findResultAC[i].destFile);
		strcpy(totalShp[idx].destPath, findResultAC[i].destPath);
		totalShp[idx].destLongi     	= findResultAC[i].destLongi  ; 
		totalShp[idx].destLati      	= findResultAC[i].destLati   ;
		totalShp[idx].destTileID    	= findResultAC[i].destTileID ;
		totalShp[idx].ID				= findResultAC[i].ID;
		totalShp[idx].isOK				= 1;
		
		/* 如果isOK == 2，说明要清除一个岛屿，添加一个fRtoCutIsland */
		if(findResultAC[i].isOK == 2)
		{
			CultNode cn;
			strcpy(cn.fileName,  findResultAC[i].fileName);     
			strcpy(cn.pathName,  findResultAC[i].pathName);     
			cn.minLong =    	 findResultAC[i].minLong;     
			cn.maxLong =    	 findResultAC[i].maxLong;     
			cn.minLati =    	 findResultAC[i].minLati;	    
			cn.maxLati =    	 findResultAC[i].maxLati;     
			strcpy(cn.destFile,  findResultAC[i].destFile);     
			strcpy(cn.destPath,  findResultAC[i].destPath);     
			cn.destLongi  = 	 findResultAC[i].destLongi  ; 
			cn.destLati   = 	 findResultAC[i].destLati   ; 
			cn.destTileID = 	 findResultAC[i].destTileID ; 
			cn.ID   =			 findResultAC[i].ID;          
			cn.isOK =			 findResultAC[i].isOK;          
			fRtoCutIsland.push_back(cn);
		}
	}
	return 1;
}


/************************************************
 *
 *    找出所有未制作文化信息的ac文件的目录集合
 *
 ************************************************/
int FindOutDirsNoCultDone()
{
	int i, j, iFlag;
	
	dirsNoCultDone.clear(); numDNCD = 0;
	for(i=0; i<cultID; i++)
	{
		if(totalcult[i].isOK == 0)
		{
			iFlag = 0;
			for(j=0; j<numDNCD; j++)
			{
				if(!strcmp(dirsNoCultDone[j].c_str(), totalcult[i].pathName))
				{
					iFlag = 1; break;
				}
			}
			if(!iFlag)
			{
				dirsNoCultDone.push_back(totalcult[i].pathName); numDNCD++;
			}
		}
	}
	return numDNCD;
}

/* 从totalcult库里面，找出符合dirsNoCultDone里面一条pathName的所有记录，包括做好的和没做好的 */
int FindAllAcItemsWithSameDir(int idx)
{
	int i;
	
	if(idx >= numDNCD) return 0;
	findResultAC.clear();	numFRAC = 0;
	for(i=0; i<cultID; i++)
	{
		if(!strcmp(dirsNoCultDone[idx].c_str(), totalcult[i].pathName))
		{
			/* 如果一个AC文化信息跨越多个经纬度，则按每个经纬度细分一下。*/
			SplitByLonLatAC(i);
		}
	}
	return numFRAC;
}

/* 已经制作完成的ac文化更新到所有ac列表中 */
int UpdateTotalcultWithNewFindResult()
{
	int i, idx;
	for(i=0;i<numFRAC;i++)
	{
		idx = findResultAC[i].ID;
		strcpy(totalcult[idx].fileName, findResultAC[i].fileName);
		strcpy(totalcult[idx].pathName, findResultAC[i].pathName);
		totalcult[idx].minLong			= findResultAC[i].minLong;
		totalcult[idx].maxLong			= findResultAC[i].maxLong;
		totalcult[idx].minLati			= findResultAC[i].minLati;	
		totalcult[idx].maxLati			= findResultAC[i].maxLati;
		strcpy(totalcult[idx].destFile, findResultAC[i].destFile);
		strcpy(totalcult[idx].destPath, findResultAC[i].destPath);
		totalcult[idx].destLongi     	= findResultAC[i].destLongi  ; 
		totalcult[idx].destLati      	= findResultAC[i].destLati   ;
		totalcult[idx].destTileID    	= findResultAC[i].destTileID ;
		totalcult[idx].ID				= findResultAC[i].ID;
		totalcult[idx].isOK				= 1;

		/* 如果isOK == 2，说明要清除一个岛屿，添加一个fRtoCutIsland */
		if(findResultAC[i].isOK == 2)
		{
			CultNode cn;
			strcpy(cn.fileName,  findResultAC[i].fileName);     
			strcpy(cn.pathName,  findResultAC[i].pathName);     
			cn.minLong =    	 findResultAC[i].minLong;     
			cn.maxLong =    	 findResultAC[i].maxLong;     
			cn.minLati =    	 findResultAC[i].minLati;	    
			cn.maxLati =    	 findResultAC[i].maxLati;     
			strcpy(cn.destFile,  findResultAC[i].destFile);     
			strcpy(cn.destPath,  findResultAC[i].destPath);     
			cn.destLongi  = 	 findResultAC[i].destLongi  ; 
			cn.destLati   = 	 findResultAC[i].destLati   ; 
			cn.destTileID = 	 findResultAC[i].destTileID ; 
			cn.ID   =			 findResultAC[i].ID;          
			cn.isOK =			 findResultAC[i].isOK;          
			fRtoCutIsland.push_back(cn);
		}
	}
	return 1;
}

/************************************************
 *
 *    根据全局变量决定将消息输出到哪里的对话框
 *
 ************************************************/
void MessageToDialog(char *msg)
{
	bool fDone = false;
	int dw;
	
	while(!fDone)
	{
		dw = WaitForSingleObject(MsgMutex, INFINITE);
		if(dw == WAIT_OBJECT_0)
		{
			if(MessageTo==0)
				MessageToDialog01(msg);
			if(MessageTo==1)
				MessageToDialog02(msg);
			Sleep(50);
			ReleaseMutex(MsgMutex);
			fDone = true;
		}
	}

#define WRITELOG
#ifdef WRITELOG	
	/* 同时要写Log文件。*/
	/* 消息前面加上时间头 */
	struct tm *local;  
	time_t t;  
	t = time(NULL);  
	local = localtime(&t); 
	
	/* 拼装完整的消息 */
	char fullmsg[128];
	sprintf(fullmsg, "%02d:%02d:%02d: %s", local->tm_hour, local->tm_min, local->tm_sec, msg);
	
	/* 写入记录文件。*/
	fprintf(fLog, "%s\n", fullmsg);
#endif
}


/************************************************
 *
 *    处理每一个目录名称
 *
 ************************************************/
int CheckTextureDirName(char *dirname, char *currPath)
{
	int i, slen, ilat, ilon, jlat, jlon, flag;
	char clat[2][5], clon[2][5];
	int  lat[2], lon[2], minlat, minlon;
	
	/* 从字符串里面读取四个数据 */
	slen = strlen(dirname);  if(slen < 6) return 1;
	jlat = jlon = 0; ilat = ilon = flag = -1;
	memset(clat, 0, sizeof(clat)); memset(clon, 0, sizeof(clon));
	for(i=0; i<slen; i++)
	{
		if((flag == 1)&&(dirname[i]>=0x30)&&(dirname[i]<=0x39)) clat[ilat][jlat++] = dirname[i];
		if((flag == 0)&&(dirname[i]>=0x30)&&(dirname[i]<=0x39)) clon[ilon][jlon++] = dirname[i];
		if((dirname[i]=='N')||(dirname[i]=='n')||(dirname[i]=='S')||(dirname[i]=='s')) {flag = 1;  jlat = 0; ilat++; }
		if((dirname[i]=='W')||(dirname[i]=='w')||(dirname[i]=='E')||(dirname[i]=='e')) {flag = 0;  jlon = 0; ilon++; }
	}
	
	/* 判断当前路径是否合法 */
	if((strlen(clon[0])==0)||(strlen(clat[0])==0)||flag==-1)
		return 1;
	
	/* 转换成整形, 取小值 */
	lat[0] = atoi(clat[0]);  lat[1] = atoi(clat[1]); lon[0] = atoi(clon[0]); lon[1] = atoi(clon[1]);
	if(lat[0]<lat[1]) minlat = lat[0]; else minlat = lat[1]; if(lon[0]<lon[1]) minlon = lon[0]; else minlon = lon[1]; 
	
	/* 填写数据结构 */
	TextDirInfo tdi;
	sprintf(tdi.dirName, "%s\\%s", currPath, dirname);
	tdi.zone[0].isOK = STATUS_UNWORK;	tdi.zone[0].lat = minlat;	tdi.zone[0].lon = minlon;
	tdi.zone[1].isOK = STATUS_UNWORK;	tdi.zone[1].lat = minlat;	tdi.zone[1].lon = minlon+1;
	tdi.zone[2].isOK = STATUS_UNWORK;	tdi.zone[2].lat = minlat+1;	tdi.zone[2].lon = minlon;
	tdi.zone[3].isOK = STATUS_UNWORK;	tdi.zone[3].lat = minlat+1;	tdi.zone[3].lon = minlon+1;
	allDirs.push_back(tdi);		TextDirIdx++;
	
	return 1;
}

/*  检查当前块是否已经做完，做完返回0， 没做设置currTi, 返回1 */
int SetDirInfo2CurrTi(int o, int i)
{
	if(allDirs[o].zone[i].isOK) return 0;
	currTi.minLongi = (float)(allDirs[o].zone[i].lon);
	currTi.minLati  = (float)(allDirs[o].zone[i].lat);
	currTi.maxLongi = (float)(allDirs[o].zone[i].lon + 1);
	currTi.maxLati  = (float)(allDirs[o].zone[i].lat + 1);
	currTi.lod      = 777;
	currTi.scale    = 1.0;
	strcpy(findJpgDirectory, allDirs[o].dirName);
	currPZ = &allDirs[o].zone[i];
	
	return 1;
}

/* 写配置，将目录信息写入配置文件中 */
int WriteDirInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  TextDirIdx);
	for(i=0; i<TextDirIdx; i++)
	{
		fprintf(fo, "\n%s %3d %3d %d", allDirs[i].dirName, allDirs[i].zone[0].lon, allDirs[i].zone[0].lat, allDirs[i].zone[0].isOK);
		fprintf(fo, "\n%s %3d %3d %d", allDirs[i].dirName, allDirs[i].zone[1].lon, allDirs[i].zone[1].lat, allDirs[i].zone[1].isOK);
		fprintf(fo, "\n%s %3d %3d %d", allDirs[i].dirName, allDirs[i].zone[2].lon, allDirs[i].zone[2].lat, allDirs[i].zone[2].isOK);
		fprintf(fo, "\n%s %3d %3d %d", allDirs[i].dirName, allDirs[i].zone[3].lon, allDirs[i].zone[3].lat, allDirs[i].zone[3].isOK);
	}
	return 1;
}

/* 读配置，从配置文件中读取所有目录信息，用于图片拼接方式的地景制作 */
int ReadDirInfoFromConfig(FILE *fi)
{
	int i, j, k, blen;
	char currPath[STRING_STANDARD_LEN];
	
	blen = strlen(baseTextDir);	
	ReadALineToArray(fi);	TextDirIdx = atoi(cf[0]);	allDirs.clear();
	for(i=0; i<TextDirIdx; i++)
	{
		ReadALineToArray(fi);
		
		/*如果目录名前半部分与baseTextDir不同，则用baseTextDir代替 */
		for(j=0; j<blen; j++) currPath[j] = cf[0][j]; currPath[j] = 0;
		if(strcmp(currPath, baseTextDir)) 	
		{
			j=1; while( !((cf[0][j-1]=='\\')&&(cf[0][j]=='n')) ) j++;
			blen = strlen(cf[0]) - j;
			for(k=0; k<blen; k++) currPath[k] = cf[0][j++]; currPath[k] = 0;
			sprintf(cf[0], "%s\\%s", baseTextDir, currPath);
		}
		TextDirInfo tdi;
		strcpy(tdi.dirName, cf[0]); 
							  tdi.zone[0].lon = atoi(cf[1]);  tdi.zone[0].lat = atoi(cf[2]);  tdi.zone[0].isOK = atoi(cf[3]); 
		ReadALineToArray(fi); tdi.zone[1].lon = atoi(cf[1]);  tdi.zone[1].lat = atoi(cf[2]);  tdi.zone[1].isOK = atoi(cf[3]); 
		ReadALineToArray(fi); tdi.zone[2].lon = atoi(cf[1]);  tdi.zone[2].lat = atoi(cf[2]);  tdi.zone[2].isOK = atoi(cf[3]); 
		ReadALineToArray(fi); tdi.zone[3].lon = atoi(cf[1]);  tdi.zone[3].lat = atoi(cf[2]);  tdi.zone[3].isOK = atoi(cf[3]); 
		allDirs.push_back(tdi);
	}
	return 1;
}

/************************************************
 *
 *    读写Map信息。与配置文件交换
 *
 ************************************************/
/* 写配置，将Map信息写入配置文件中 */
int WriteMapInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  nodeID);
	//再写入Map表
	for(i=0; i<nodeID; i++)
	{
		fprintf(fo, "\n%11.6f %11.6f %11.6f %11.6f %5d %5d %s %s", totalnode[i].minLong, totalnode[i].maxLong, totalnode[i].minLati, totalnode[i].maxLati, totalnode[i].image_width, totalnode[i].image_height, totalnode[i].fileName, totalnode[i].pathName);
	}
	return 1;
}

/* 读配置，从配置文件中读取所有Map信息，用于图片拼接方式的地景制作 */
int ReadMapInfoFromConfig(FILE *fi)
{
	int i, j, k, blen;
	char currPath[STRING_STANDARD_LEN];
	
	blen = strlen(baseTextDir);
	while(nodeID == 0)	{ReadALineToArray(fi);	nodeID = atoi(cf[0]);}
	for(i=0; i<nodeID; i++)
	{
		MapNode mn;
		ReadALineToArray(fi);
		mn.minLong = atof(cf[0]);
		mn.maxLong = atof(cf[1]);
		mn.minLati = atof(cf[2]);	
		mn.maxLati = atof(cf[3]);
		mn.image_width =  atoi(cf[4]); 
		mn.image_height = atoi(cf[5]);
		strcpy(mn.fileName, cf[6]);
		
		/*如果目录名前半部分与baseTextDir不同，则用baseTextDir代替 */
		for(j=0; j<blen; j++) currPath[j] = cf[7][j]; currPath[j] = 0;
		if(strcmp(currPath, baseTextDir)) 	
		{
			j=1; while( !((cf[7][j-1]=='\\')&&(cf[7][j]=='n')) ) j++;
			blen = strlen(cf[7]) - j;
			for(k=0; k<blen; k++) currPath[k] = cf[7][j++]; currPath[k] = 0;
			sprintf(cf[7], "%s\\%s", baseTextDir, currPath);
		}
		strcpy(mn.pathName, cf[7]);
		mn.ID = i;
		totalnode.push_back(mn);
	}
	return 1;
}

/************************************************
 *
 *    读写HighTxt信息。与配置文件交换
 *
 ************************************************/
/* 写配置，将用户高分辨率地景纹理信息写入配置文件中 */
int WriteHTxtInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  hTxtNum);
	//再写入HighTxt表
	for(i=0; i<hTxtNum; i++)
	{
		fprintf(fo, "\n%11.6f %11.6f %11.6f %11.6f %5d %5d %1d %s %s", highTxtnode[i].minLong, highTxtnode[i].maxLong, highTxtnode[i].minLati, highTxtnode[i].maxLati, highTxtnode[i].image_width, highTxtnode[i].image_height, highTxtnode[i].highLevel, highTxtnode[i].fileName, highTxtnode[i].pathName);
	}
	return 1;
}

/* 读配置，从配置文件中读取所有用户高分辨率地景纹理信息，用于Tif方式的地景制作 */
int ReadHTxtInfoFromConfig(FILE *fi)
{
	int i, j, k, blen;
	char currPath[STRING_STANDARD_LEN];
	
	blen = strlen(baseTextDir);
	ReadALineToArray(fi);	hTxtNum = atoi(cf[0]);	highTxtnode.clear();
	for(i=0; i<hTxtNum; i++)
	{
		HighNode hn;
		ReadALineToArray(fi);
		hn.minLong = atof(cf[0]);
		hn.maxLong = atof(cf[1]);
		hn.minLati = atof(cf[2]);	
		hn.maxLati = atof(cf[3]);
		hn.image_width =  atoi(cf[4]); 
		hn.image_height = atoi(cf[5]);
		hn.highLevel = atoi(cf[6]);
		strcpy(hn.fileName, cf[7]);

		/*如果目录名前半部分与baseTextDir不同，则用baseTextDir代替 */
		for(j=0; j<blen; j++) currPath[j] = cf[8][j]; currPath[j] = 0;
		if(strcmp(currPath, baseTextDir)) 	
		{
			j=1; while( !((cf[8][j-1]=='\\')&&(cf[8][j]=='U')) ) j++;
			blen = strlen(cf[8]) - j;
			for(k=0; k<blen; k++) currPath[k] = cf[8][j++]; currPath[k] = 0;
			sprintf(cf[8], "%s\\%s", baseTextDir, currPath);
		}
		strcpy(hn.pathName, cf[8]);
		hn.ID = i;
		highTxtnode.push_back(hn);
	}
	return 1;
}


/************************************************
 *
 *    读写LowTxt信息。与配置文件交换
 *
 ************************************************/
/* 写配置，将系统缺省低分辨率地景纹理信息写入配置文件中 */
int WriteLTxtInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  lTxtNum);
	//再写入LowTxt表
	for(i=0; i<lTxtNum; i++)
	{
		fprintf(fo, "\n%11.6f %11.6f %11.6f %11.6f %5d %5d %1d %s %s", lowTxtnode[i].minLong, lowTxtnode[i].maxLong, lowTxtnode[i].minLati, lowTxtnode[i].maxLati, lowTxtnode[i].image_width, lowTxtnode[i].image_height, lowTxtnode[i].highLevel, lowTxtnode[i].fileName, lowTxtnode[i].pathName);
	}
	return 1;
}

/* 读配置，从配置文件中读取所有系统缺省分辨率地景纹理信息，用于Tif方式的地景制作 */
int ReadLTxtInfoFromConfig(FILE *fi)
{
	int i, j, k, blen;
	char currPath[STRING_STANDARD_LEN];
	
	blen = strlen(baseTextDir);
	ReadALineToArray(fi);	lTxtNum = atoi(cf[0]);	lowTxtnode.clear();
	for(i=0; i<lTxtNum; i++)
	{
		HighNode hn;
		ReadALineToArray(fi);
		hn.minLong = atof(cf[0]);
		hn.maxLong = atof(cf[1]);
		hn.minLati = atof(cf[2]);	
		hn.maxLati = atof(cf[3]);
		hn.image_width =  atoi(cf[4]); 
		hn.image_height = atoi(cf[5]);
		hn.highLevel = atoi(cf[6]);
		strcpy(hn.fileName, cf[7]);

		/*如果目录名前半部分与baseTextDir不同，则用baseTextDir代替 */
		for(j=0; j<blen; j++) currPath[j] = cf[8][j]; currPath[j] = 0;
		if(strcmp(currPath, baseTextDir)) 	
		{
			j=1; while( !((cf[8][j-1]=='\\')&&(cf[8][j]=='U')) ) j++;
			blen = strlen(cf[8]) - j;
			for(k=0; k<blen; k++) currPath[k] = cf[8][j++]; currPath[k] = 0;
			sprintf(cf[8], "%s\\%s", baseTextDir, currPath);
		}
		strcpy(hn.pathName, cf[8]);
		hn.ID = i;
		lowTxtnode.push_back(hn);
	}
	return 1;
}


/************************************************
 *
 *    读写HighHgt信息。与配置文件交换
 *
 ************************************************/
/* 写配置，将用户高分辨率地景高程信息写入配置文件中 */
int WriteHHgtInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  hHgtNum);
	//再写入HighTxt表
	for(i=0; i<hHgtNum; i++)
	{
		fprintf(fo, "\n%11.6f %11.6f %11.6f %11.6f %5d %5d %1d %s %s", highHgtnode[i].minLong, highHgtnode[i].maxLong, highHgtnode[i].minLati, highHgtnode[i].maxLati, highHgtnode[i].image_width, highHgtnode[i].image_height, highHgtnode[i].highLevel, highHgtnode[i].fileName, highHgtnode[i].pathName);
	}
	return 1;
}

/* 读配置，从配置文件中读取所有用户高分辨率地景高程信息，用于Tif方式的地景制作 */
int ReadHHgtInfoFromConfig(FILE *fi)
{
	int i, j, k, blen;
	char currPath[STRING_STANDARD_LEN];
	
	blen = strlen(baseHgtDir);
	ReadALineToArray(fi);	hHgtNum = atoi(cf[0]);	highHgtnode.clear();
	for(i=0; i<hHgtNum; i++)
	{
		HighNode hn;
		ReadALineToArray(fi);
		hn.minLong = atof(cf[0]);
		hn.maxLong = atof(cf[1]);
		hn.minLati = atof(cf[2]);	
		hn.maxLati = atof(cf[3]);
		hn.image_width =  atoi(cf[4]); 
		hn.image_height = atoi(cf[5]);
		hn.highLevel = atoi(cf[6]);
		strcpy(hn.fileName, cf[7]);

		/*如果目录名前半部分与baseTextDir不同，则用baseTextDir代替 */
		for(j=0; j<blen; j++) currPath[j] = cf[8][j]; currPath[j] = 0;
		if(strcmp(currPath, baseHgtDir)) 	
		{
			j=1; while( !((cf[8][j-1]=='\\')&&(cf[8][j]=='U')) ) j++;
			blen = strlen(cf[8]) - j;
			for(k=0; k<blen; k++) currPath[k] = cf[8][j++]; currPath[k] = 0;
			sprintf(cf[8], "%s\\%s", baseHgtDir, currPath);
		}
		strcpy(hn.pathName, cf[8]);
		hn.ID = i;
		highHgtnode.push_back(hn);
	}
	return 1;
}

/************************************************
 *
 *    读写Shp文化初始信息。与配置文件交换
 *
 ************************************************/
/* 写配置，将用户shp格式文化信息写入配置文件中 */
int WriteShpInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  shpID);
	//再写入Map表
	for(i=0; i<shpID; i++)
	{
		fprintf(fo, "\n%s %s %11.6f %11.6f %11.6f %11.6f %s %s %3d %2d %d %d %d", totalShp[i].fileName, totalShp[i].pathName, totalShp[i].minLong, totalShp[i].maxLong, totalShp[i].minLati, totalShp[i].maxLati, totalShp[i].destFile, totalShp[i].destPath, totalShp[i].destLongi, totalShp[i].destLati, totalShp[i].destTileID, totalShp[i].ID, totalShp[i].isOK );
	}
	return 1;
}

/* 读配置，从配置文件中读取所有用户shp格式信息。 */
int ReadShpInfoFromConfig(FILE *fi)
{
	int i;
	
	ReadALineToArray(fi);	shpID = atoi(cf[0]);	totalShp.clear();
	for(i=0; i<shpID; i++)
	{
		CultNode cn;
		ReadALineToArray(fi);
		strcpy(cn.fileName,  cf[0]);
		strcpy(cn.pathName,  cf[1]);
		cn.minLong =    atof(cf[2]);
		cn.maxLong =    atof(cf[3]);
		cn.minLati =    atof(cf[4]);	
		cn.maxLati =    atof(cf[5]);
		strcpy(cn.destFile,  cf[6]);
		strcpy(cn.destPath,  cf[7]);
		cn.destLongi  = atoi(cf[8]); 
		cn.destLati   = atoi(cf[9]);
		cn.destTileID = atoi(cf[10]);
		cn.ID 		= i;
		cn.isOK 		= atoi(cf[12]);
		totalShp.push_back(cn);
	}
	return 1;
}



/************************************************
 *
 *    读写Ac文化初始信息。与配置文件交换
 *
 ************************************************/
/* 写配置，将用户ac格式文化信息写入配置文件中 */
int WriteCultInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  cultID);
	//再写入Map表
	for(i=0; i<cultID; i++)
	{
		fprintf(fo, "\n%s %s %11.6f %11.6f %11.6f %11.6f %s %s %3d %2d %d %d %d", totalcult[i].fileName, totalcult[i].pathName, totalcult[i].minLong, totalcult[i].maxLong, totalcult[i].minLati, totalcult[i].maxLati, totalcult[i].destFile, totalcult[i].destPath, totalcult[i].destLongi, totalcult[i].destLati, totalcult[i].destTileID, totalcult[i].ID, totalcult[i].isOK );
	}
	return 1;
}

/* 读配置，从配置文件中读取所有用户ac格式信息。 */
int ReadCultInfoFromConfig(FILE *fi)
{
	int i;
	
	ReadALineToArray(fi);	cultID = atoi(cf[0]);	totalcult.clear();
	for(i=0; i<cultID; i++)
	{
		CultNode cn;
		ReadALineToArray(fi);
		strcpy(cn.fileName,  cf[0]);
		strcpy(cn.pathName,  cf[1]);
		cn.minLong =    atof(cf[2]);
		cn.maxLong =    atof(cf[3]);
		cn.minLati =    atof(cf[4]);
		cn.maxLati =    atof(cf[5]);
		strcpy(cn.destFile,  cf[6]);
		strcpy(cn.destPath,  cf[7]);
		cn.destLongi  = atoi(cf[8]); 
		cn.destLati   = atoi(cf[9]);
		cn.destTileID = atoi(cf[10]);
		cn.ID 		= i;
		cn.isOK 		= atoi(cf[12]);
		totalcult.push_back(cn);
	}
	return 1;
}

/************************************************
 *
 *    管理高精度地景块制作列表
 *
 ************************************************/
std::vector<HighTerrainInfo> allHTI;
unsigned int numHTI;
std::vector<HighTerrainInfo> lowHTI;
unsigned int numLowHTI;
std::vector<PerZone>  ZoneWithHTI;
int numZWHTI;
std::vector<HighTerrainInfo> findResultHTI;
int numFRHTI;

typedef struct _TERRAINPOINT_
{
	double px;
	double py;				//实际直角坐标系坐标   需要计算转换过来的
	double pz;				
	double lon;
	double lat;				//经纬度坐标，需要计算得来
	double elev;
	double tU;
	double tV;				//纹理坐标   直接可读取的
}TerrainPoint;

std::vector<TerrainOnDemand> doTerrain;
int numDT;


/**************************************************
 *
 *   判断当前这个tileID在地景库中是否有地景存在，存在返回1， 否则返回0
 *
 **************************************************/

int IfTileIDExistInLib(HighTerrainInfo *hti)
{
	char basePath[STRING_STANDARD_LEN];
	char fullPath[_MAX_PATH];
	float W, S, offLon, offLat;
	Bucket b;
	FILE *bin;
	int iFlag = 0;

	//计算相关路径信息
	W = hti->startLongi;
	S = hti->startLati;
	offLon = hti->offsetLongi;
	offLat = hti->offsetLati; 
	set_bucket(&b, W+0.01, S+0.01);

	/* 拼装输出路径及文件 */
	gen_base_path(&b, basePath);
	sprintf(fullPath, "%s\\%s\\%d.ive", iveroot, basePath, hti->tileID);
	
	/* 看看是否能够打开文件 */
	bin = fopen(fullPath, "rb");
	if(bin != NULL)
	{
		iFlag = 1;	fclose(bin);
	}
	return iFlag;
}


/* 在allHTI搜索出共占用了多少平方度  */
int CheckZoneInAllHti()
{
	int i, j, ilon, ilat, iFlag;
	numZWHTI = 0;	ZoneWithHTI.clear();
	for(i=0; i<numHTI; i++)
	{
		iFlag = 0;
		ilon = (int)allHTI[i].startLongi;
		ilat = (int)allHTI[i].startLati;
		for(j=0; j<numZWHTI; j++)
		{
			if ((ZoneWithHTI[j].lon == ilon) && (ZoneWithHTI[j].lat == ilat)) 
			{
				iFlag = 1; break;
			}
		}
		if(!iFlag)
		{
			PerZone pz;	pz.lon = ilon; pz.lat = ilat; ZoneWithHTI.push_back(pz); numZWHTI++;
		}
	}
	return 1;
}

/* 在allHTI里面找出在同一个Zone(平方度)中的所有项目，存放在findResultHTI里面 */
int FindAllHTIInSameZone(int idx)
{
	int i, ilon, ilat;
	findResultHTI.clear(); numFRHTI=0;
	for(i=0; i<numHTI; i++)
	{
		ilon = (int)allHTI[i].startLongi;
		ilat = (int)allHTI[i].startLati;
		if((ZoneWithHTI[idx].lon == ilon) && (ZoneWithHTI[idx].lat == ilat)) 
		{
			HighTerrainInfo hti;
			strcpy(hti.texFname, allHTI[i].texFname);  
			strcpy(hti.hgtFname, allHTI[i].hgtFname);  
			hti.startLongi = allHTI[i].startLongi	;
			hti.startLati  = allHTI[i].startLati 	;
			hti.offsetLongi= allHTI[i].offsetLongi	;
			hti.offsetLati = allHTI[i].offsetLati 	;
			hti.highLevel  = allHTI[i].highLevel 	;
			hti.tileID     = allHTI[i].tileID    	;
			hti.ID         = allHTI[i].ID        	;
			hti.isOK       = allHTI[i].isOK      	;
			findResultHTI.push_back(hti);	numFRHTI++;
		}
	}	
	return 1;
}


/* 检查某一区块的高精度地景是否已经全部生成完毕 */
int CheckHTinZoneIsDone()
{
	int i, iFlag = 1;
	for(i=0; i<numFRHTI; i++)
	{
		if(findResultHTI[i].isOK == 0)
		{
			iFlag = 0; break;
		}
	}
	return iFlag;
}

/******************************************************************************************************
 *  这里插入计算高精度地景的新算法，支持更加复杂的情况 ---- 2012.01.20
 *****************************************************************************************************/
/* 把findResultHTI数组里面的所有经纬度信息，分别进行从小到大的排序 */
double LongiSortList[32], LatiSortList[32];
int numSLon, numSLat;
#define DELTA_THREAHOLD 0.0001
int SortTheFRHTI(TerrainOnDemand *curr)
{
	int i, j;
	
	/* 先把边界信息加进去 */
	memset(LongiSortList, 0, sizeof(LongiSortList)); memset(LatiSortList, 0, sizeof(LatiSortList));
	LongiSortList[0] = curr->minLongi; LongiSortList[1] = curr->maxLongi;  numSLon = 2;
	LatiSortList[0]  = curr->minLati;  LatiSortList[1]  = curr->maxLati;   numSLat = 2;
	
	/* 再把findResultHTI里面的最大最小经纬度值都加进去，如果有重复的，则忽略不计 */
	for(i=0; i<numFRHTI; i++)
	{
		double minFrLon, minFrLat, maxFrLon, maxFrLat;
		int iFlag = 0;
		minFrLon = findResultHTI[i].startLongi; maxFrLon = findResultHTI[i].startLongi + findResultHTI[i].offsetLongi;
		minFrLat = findResultHTI[i].startLati;  maxFrLat = findResultHTI[i].startLati  + findResultHTI[i].offsetLati;

		/* 如果高精度块比当前精度块大，则限制在当前精度块内 */
		if(minFrLon < curr->minLongi) minFrLon = curr->minLongi;
		if(minFrLat < curr->minLati ) minFrLat = curr->minLati ;
		if(maxFrLon > curr->maxLongi) maxFrLon = curr->maxLongi;
		if(maxFrLat > curr->maxLati ) maxFrLat = curr->maxLati ;
		
		/* 检测 minFrLon 不重复后，插入数组内 */
		iFlag = 0;
		for(j=0; j<numSLon; j++)
		{
			double delta = fabs(LongiSortList[j] - minFrLon);
			if(delta < DELTA_THREAHOLD)	{ iFlag = 1; break; }
		}
		if(!iFlag) { LongiSortList[numSLon++] = minFrLon; }

		/* 检测 maxFrlon 不重复后，插入数组内 */
		iFlag = 0;
		for(j=0; j<numSLon; j++)
		{
			double delta = fabs(LongiSortList[j] - maxFrLon);
			if(delta < DELTA_THREAHOLD)	{ iFlag = 1; break; }
		}
		if(!iFlag) { LongiSortList[numSLon++] = maxFrLon; }

		/* 检测 minFrLat 不重复后，插入数组内 */
		iFlag = 0;
		for(j=0; j<numSLat; j++)
		{
			double delta = fabs(LatiSortList[j] - minFrLat);
			if(delta < DELTA_THREAHOLD)	{ iFlag = 1; break; }
		}
		if(!iFlag) { LatiSortList[numSLat++] = minFrLat; }

		/* 检测 maxFrLat 不重复后，插入数组内 */
		iFlag = 0;
		for(j=0; j<numSLat; j++)
		{
			double delta = fabs(LatiSortList[j] - maxFrLat);
			if(delta < DELTA_THREAHOLD)	{ iFlag = 1; break; }
		}
		if(!iFlag) { LatiSortList[numSLat++] = maxFrLat; }
	}
	
	/* 接下来对这两个数组 LongiSortList,LatiSortList 进行从低到高的排序 */
	for(i=0; i<numSLon - 1; i++)
	{
		for(j=i+1; j<numSLon; j++)
		{
			double Ex;
			if(LongiSortList[i] > LongiSortList[j]) {Ex = LongiSortList[i]; LongiSortList[i] = LongiSortList[j]; LongiSortList[j] = Ex;}
		}
	}

	for(i=0; i<numSLat - 1; i++)
	{
		for(j=i+1; j<numSLat; j++)
		{
			double Ex;
			if(LatiSortList[i] > LatiSortList[j]) {Ex = LatiSortList[i]; LatiSortList[i] = LatiSortList[j]; LatiSortList[j] = Ex;}
		}
	}
	return 1;	
}

/* 定义一个数据结构，记录二维数组每一个交叉点（格点）的信息 */
typedef struct _PERGRIDINFO_
{
	int 	RowIndex;
	int		ColIndex;			//在二维数组中的位置；
	double	MinLongi;
	double	MaxLongi;
	double	MinLati;
	double	MaxLati;			//本格所包含的区域
	double	MidLongi;
	double	MidLati;			//本格的中心点经纬度
	int		Info;				//本格的状态	0：未被更高精度覆盖； 1：已被更高精度覆盖
}PerGridInfo;
PerGridInfo AllGrids[32][32];

/* 将排序好的LongiSortList，LatiSortList组成二维网格，并检查记录每个网格点的状态 */
int MakeGridFindInfo()
{
	int i, j, k, iFlag;
	memset(AllGrids, 0, sizeof(AllGrids));
	
	/* 先对网格进行赋值 行为经度，列为纬度 */
	for(i=0; i<numSLat-1; i++)
	{
		for(j=0; j<numSLon-1; j++)
		{
			AllGrids[i][j].RowIndex = i;
			AllGrids[i][j].ColIndex = j;
			AllGrids[i][j].MinLongi = LongiSortList[j];
			AllGrids[i][j].MaxLongi = LongiSortList[j+1];
			AllGrids[i][j].MinLati  = LatiSortList [i];
			AllGrids[i][j].MaxLati  = LatiSortList [i+1];
			AllGrids[i][j].MidLongi = (AllGrids[i][j].MinLongi + AllGrids[i][j].MaxLongi) / 2.0;
			AllGrids[i][j].MidLati  = (AllGrids[i][j].MinLati  + AllGrids[i][j].MaxLati ) / 2.0;
			
			/* 根据中点判断被占用*/
			iFlag = 0;
			for(k=0; k<numFRHTI; k++)
			{
				double minFrLon, minFrLat, maxFrLon, maxFrLat;
				minFrLon = findResultHTI[k].startLongi; maxFrLon = findResultHTI[k].startLongi + findResultHTI[k].offsetLongi;
				minFrLat = findResultHTI[k].startLati;  maxFrLat = findResultHTI[k].startLati  + findResultHTI[k].offsetLati;
				if(	(AllGrids[i][j].MidLongi > minFrLon) &&
					(AllGrids[i][j].MidLongi < maxFrLon) &&
					(AllGrids[i][j].MidLati  > minFrLat) &&
					(AllGrids[i][j].MidLati  < maxFrLat) )
				{
					iFlag = 1; break;
				}
			}
			if(iFlag)	{AllGrids[i][j].Info = 1;}
		}
	}
	return 1;
}

/* 根据网格里面的信息，来收集低精度地景块所需制作的doTerrain */
typedef struct _RECORDGRIDINFO_
{
	int 	Row;
	int		Col;			//在二维数组中的位置；
}RecordGridInfo;

/* 根据二维网格内各个块的分布情况，（其中Info>0是被高分辨率地景覆盖，Info=0没有被覆盖），按照最大方块原则，制作低分辨率地景。*/
int CalcDTnByGrids(HighTerrainInfo *curr)
{
	int i, j, k;
	RecordGridInfo RecGrid[32][32];
	int numGRow, numGCol;
	
	while(1)
	{
		int StartRow, StartCol, iFlag, NextRow, NextCol, Lines, bkFlag ;
		
		/* 首先查找当前最左上角的 Info = 0 的块，记录开始行StartRow和开始列StartCol */
		iFlag = 0;
		for(i=0; i<numSLat-1; i++)
		{
			bkFlag = 0;
			for(j=0; j<numSLon-1; j++)
			{
				if(AllGrids[i][j].Info == 0){StartRow = i; StartCol = j; iFlag = 1; bkFlag = 1; break;}
			}
			if(bkFlag) break;
		}
		
		/* 如果找不到Info == 0 的块，则表示已经做完，退出循环 */
		if(!iFlag) break;
		
		/* 根据开始行StartRow和开始列StartCol，向下查找，找到第一个 Info != 0的块，或者尽头，记录在 RecGrid。*/
		memset(RecGrid, -1, sizeof(RecGrid)); RecGrid[0][0].Row = StartRow; RecGrid[0][0].Col = StartCol; 	numGRow = numGCol = 1;
		NextRow = StartRow + 1; Lines = 1;
		AllGrids[StartRow][StartCol].Info = 2;
		while((AllGrids[NextRow][StartCol].Info == 0) && (NextRow < (numSLat-1)))
		{
			AllGrids[NextRow][StartCol].Info = 2;
			RecGrid[numGRow][0].Row = NextRow; RecGrid[numGRow][0].Col = StartCol; NextRow++; Lines++; numGRow++;
		}
		
		/* 再从 StartCol+1 开始向右查找，一列一列地查，列的所有行元素Info都为0是正确的列 */
		NextCol = StartCol + 1;
		while(NextCol < (numSLon-1))
		{
			int sFlag = 0;
			for(k=0; k<Lines; k++) 	{	sFlag += AllGrids[StartRow + k][NextCol].Info;	}
			if(sFlag == 0)
			{
				for(k=0; k<Lines; k++)
				{
					AllGrids[StartRow + k][NextCol].Info = 2;
					RecGrid[k][numGCol].Row = StartRow + k; RecGrid[k][numGCol].Col = NextCol;
				}
				NextCol++; numGCol++;
				if(NextCol >= (numSLon-1)) break;
			}else{
				break;
			}
		}
		
		/* 给DoTerrain数组赋值 */
		if(curr)
		{
			TerrainOnDemand tod;
			tod.minLongi = AllGrids[RecGrid[0][0].Row][RecGrid[0][0].Col].MinLongi ; 
			tod.minLati  = AllGrids[RecGrid[0][0].Row][RecGrid[0][0].Col].MinLati  ;
			tod.maxLongi = AllGrids[RecGrid[numGRow-1][numGCol-1].Row][RecGrid[numGRow-1][numGCol-1].Col].MaxLongi ; 
			tod.maxLati  = AllGrids[RecGrid[numGRow-1][numGCol-1].Row][RecGrid[numGRow-1][numGCol-1].Col].MaxLati  ;
	
			HighTerrainInfo tmp;
			tmp.startLongi = tod.minLongi;
			tmp.startLati  = tod.minLati;
			tmp.offsetLongi= tod.maxLongi - tod.minLongi;
			tmp.offsetLati = tod.maxLati  - tod.minLati ;
			tmp.highLevel  = curr->highLevel;
			strcpy(tmp.texFname, curr->texFname);
			strcpy(tmp.hgtFname, curr->hgtFname);
			tmp.tileID     = curr->tileID; 
			tmp.isOK       = 0;
			SetCurrHTI(&tmp);
		}else{
			TerrainOnDemand tod;
			tod.minLongi = AllGrids[RecGrid[0][0].Row][RecGrid[0][0].Col].MinLongi ; 
			tod.minLati  = AllGrids[RecGrid[0][0].Row][RecGrid[0][0].Col].MinLati  ;
			tod.maxLongi = AllGrids[RecGrid[numGRow-1][numGCol-1].Row][RecGrid[numGRow-1][numGCol-1].Col].MaxLongi ; 
			tod.maxLati  = AllGrids[RecGrid[numGRow-1][numGCol-1].Row][RecGrid[numGRow-1][numGCol-1].Col].MaxLati  ;
			doTerrain.push_back(tod);	numDT++;
		}
	}
	return 1;
}

/* 做总体调用 */
/* 检查1个经纬度区间内，没有被高精度地景覆盖的部分 */
int CheckNoHighInZoneEx()
{
	TerrainOnDemand	tod;	
	tod.minLongi  = currTi.minLongi; tod.maxLongi  = currTi.maxLongi;
	tod.minLati   = currTi.minLati;  tod.maxLati   = currTi.maxLati;
	
	doTerrain.clear();	numDT = 0;
	if(!SortTheFRHTI(&tod))	return 0;
	if(!MakeGridFindInfo()) return 0;
	if(!CalcDTnByGrids(NULL)) return 0;

	return 1;
}

/* 计算一个高分辨率纹理块被另一个更高分辨率的纹理块全覆盖 */
int CalcuHighOverlayByHigherEx(HighTerrainInfo *curr)
{
	TerrainOnDemand	tod;	
	tod.minLongi  = curr->startLongi; tod.maxLongi  = curr->startLongi + curr->offsetLongi ;
	tod.minLati   = curr->startLati;  tod.maxLati   = curr->startLati  + curr->offsetLati  ;
	
	if(!SortTheFRHTI(&tod))	return 0;
	if(!MakeGridFindInfo()) return 0;
	if(!CalcDTnByGrids(curr)) return 0;

	return 1;
}


/******************************************************************************************************
 *  新算法插入完毕 ---- 2012.01.21
 *****************************************************************************************************/

//根据收集到的doTerrain，做低精度地景块
int DoLowRTerrain()
{
	int i;
	char msg[STRING_STANDARD_LEN];
	
	for(i=0; i<numDT; i++)
	{
		sprintf(msg, "制作第 %d 块低精度区块。", i+1);
		MessageToDialog(msg);

		if(!DoPartVpbAndWriteStg(i))
		{
			MessageToDialog("执行VPB失败\n");
			return 0;			
		}
	}
	return 1;
}

/* 写回到记录总集里面 */
int WriteBackToallHTI()
{
	int i;
	for(i=0; i<numFRHTI; i++)
	{
		int j = findResultHTI[i].ID;
		
		strcpy(allHTI[j].texFname, findResultHTI[i].texFname);  
		strcpy(allHTI[j].hgtFname, findResultHTI[i].hgtFname);  
		allHTI[j].startLongi = findResultHTI[i].startLongi	;
		allHTI[j].startLati  = findResultHTI[i].startLati 	;
		allHTI[j].offsetLongi= findResultHTI[i].offsetLongi	;
		allHTI[j].offsetLati = findResultHTI[i].offsetLati 	;
		allHTI[j].highLevel  = findResultHTI[i].highLevel 	;
		allHTI[j].tileID     = findResultHTI[i].tileID    	;
		allHTI[j].ID         = findResultHTI[i].ID        	;
		allHTI[j].isOK       = findResultHTI[i].isOK      	;		
	}
	return 1;
}

/* 根据经度和纬度的整数值，查找findJpgDirectory */
int SearchDirWithLonLat(int lon, int lat)
{
	int i, j, retval = 0;
	for(i=0; i<TextDirIdx; i++)
	{
		for(j=0; j<4; j++)
		{
			if((lon == allDirs[i].zone[j].lon) && (lat == allDirs[i].zone[j].lat))
			{
				strcpy(findJpgDirectory, allDirs[i].dirName); retval = 1; break;
			}
		}
	}
	return retval;
}


/* 制作高分辨率地景的流程 */
int DoHightResolutionTerrain()
{
	int n, j;
	PerZone mPZ;
	char msg[STRING_STANDARD_LEN];
	
	/* 首先检查所有高分辨率纹理总共占据了多少经度x纬度（平方度）*/
	CheckZoneInAllHti();
	
	/* 对每个被占用的平方度做循环*/
	for(n=0; n<numZWHTI; n++)
	{
		int lowEnable = 1;		//本经纬度块是否能做低分辨率地景数据？ 2012-11-1

		/* 先收集本平方度内所有的高精度区块，结果放到 findResultHTI里面。 */
		FindAllHTIInSameZone(n);
		mPZ.lon = ZoneWithHTI[n].lon; mPZ.lat = ZoneWithHTI[n].lat; mPZ.isOK = 0;
		
		/* 先制作这个平方度内的基本低分辨率地景 *//* 生成基础地形 */
		sprintf(msg, "开始制作经度%d 纬度%d 区域的基本地景.", mPZ.lon, mPZ.lat);
		MessageToDialog(msg);
		currTi.minLongi = (float)(mPZ.lon);
		currTi.minLati  = (float)(mPZ.lat);
		currTi.maxLongi = (float)(mPZ.lon + 1);
		currTi.maxLati  = (float)(mPZ.lat + 1);
		currPZ = &mPZ;

		/* 首先要拼出低精度地景的高程数据和纹理数据 */
		int iW = (int)currTi.minLongi, iS = (int)currTi.minLati;
		if(CheckIfHas50(iW, iS) != -1)
		{
			if(!CreatePerDegreeTerrain(2))
			{
				sprintf(msg, "制作经度%d 纬度%d 区域的基本地景出错.", mPZ.lon, mPZ.lat);
				MessageToDialog(msg);
				lowEnable = 0;
			}
		}else{
			if(!CreatePerDegreeTerrain(1))
			{
				sprintf(msg, "制作经度%d 纬度%d 区域的基本地景出错.", mPZ.lon, mPZ.lat);
				MessageToDialog(msg);
				lowEnable = 0;
			}
		}

		/* 把本区块没做过高精度地景的tile全部做完 */
		for(j=0;j<numFRHTI;j++)
		{
			/* 略过已经做完的 */
			if(findResultHTI[j].isOK) continue;
			
			/* 开始制作一个tile的高精度地景 */
			sprintf(msg, "开始制作Tile %d 区域的高精度地景.", findResultHTI[j].tileID);
			MessageToDialog(msg);

			/* 装配vpb批处理命令。*/
			MessageToDialog("调用VPB生成地景");

			/* 装配vpb批处理命令。并调用vpb，制作当前块高精度地景*/
			if(!DoHighVpbAndWriteStg(j))
			{
				MessageToDialog("执行VPB失败\n");
				continue;		//return 0;
			}
			
			/* 记录这个块已经做好 */
			findResultHTI[j].isOK = 1;
		}

		/* 计算需要生成的低精度区块系列 doTerrain[]*/
		if(!CheckNoHighInZoneEx()) return 0;
		
		/* 如果numDT为0，说明该地块全部覆盖了高精度地景。不需要再做低精度地景。 */
		if((numDT != 0)&&(lowEnable))
		{
			/* 分别制作各个低精度区块 */
			if(!DoLowRTerrain())
			{   
				sprintf(msg, "制作低精度区块--经度 %d 纬度 %d 失败。", mPZ.lon, mPZ.lat);
				MessageToDialog(msg);
				return 0;
			}
		}

		/* 做完了以后，还要把做好的结果写回到allHTI */
		WriteBackToallHTI();

		/* 把allHTI写回配置文件 */
		CheckConfig(1);
		
		sprintf(msg, "经度 %d 纬度 %d 区块内高精度地景制作完成。", mPZ.lon, mPZ.lat);
		MessageToDialog(msg);
	}
	return 1;
}

/************************************************
 *
 *    填充高精度制作系列allHTI
 *
 ************************************************/

/* 设置当前HTI, 如果先前设过了，则跳过 */
void SetCurrHTI(HighTerrainInfo *curr)
{
	int i, iFlag, iNum;
	
	/* 根据经纬度来检查是否之前的项目有重复。*/
	iFlag = 0;
	for(i=0;i<numHTI;i++)
	{
		if(	(fabs(allHTI[i].startLongi  - curr->startLongi ) < 0.001) &&
			(fabs(allHTI[i].startLati   - curr->startLati  ) < 0.001) && 
			(fabs(allHTI[i].offsetLongi - curr->offsetLongi) < 0.001) && 
			(fabs(allHTI[i].offsetLati  - curr->offsetLati ) < 0.001)  )
		{ iFlag = 1; break; }
	}
	
	/* 如果没有重复，将新项目添加到列表最后一项。*/
	if(!iFlag)
	{
		HighTerrainInfo hti;

		/* 保证tileID不重复 */
		iNum = 0;
		for(i=0;i<numHTI;i++)
		{
			int iTmp = (int)(allHTI[i].tileID / 100);
			if(iTmp == curr->tileID)	iNum++;
		}
		
		hti.startLongi  = curr->startLongi ; 
		hti.startLati   = curr->startLati  ; 
		hti.offsetLongi = curr->offsetLongi; 
		hti.offsetLati  = curr->offsetLati ; 
		hti.highLevel   = curr->highLevel  ; 
		hti.tileID      = (curr->tileID * 100) + iNum; 	/* 赋予当前级别地块一个不重复的ID号 */
		strcpy(hti.texFname, curr->texFname);
		strcpy(hti.hgtFname, curr->hgtFname);
		hti.ID          = numHTI;
		hti.isOK		= IfTileIDExistInLib(&hti);
		allHTI.push_back(hti);	numHTI++;
	}
}

/* 设置比较HTI */
void SetHigherHTIEx(HighTerrainInfo *curr)
{
	HighTerrainInfo hti;
	hti.startLongi  = curr->startLongi ; 
	hti.startLati   = curr->startLati  ; 
	hti.offsetLongi = curr->offsetLongi; 
	hti.offsetLati  = curr->offsetLati ; 
	hti.highLevel   = curr->highLevel  ; 
	hti.tileID      = curr->tileID     ; 
	strcpy(hti.texFname, curr->texFname);
	strcpy(hti.hgtFname, curr->hgtFname);
	hti.ID          = curr->ID	     ; 
	hti.isOK		= curr->isOK;
	findResultHTI.push_back(hti);	numFRHTI++;
}

/* 判断当前的高分辨率块，是否包含更高一层分辨率的子块 */
int DetectHTContainHigherHT(HighTerrainInfo *curr, int higherLevel)
{
	int i;

	/* 判断更高一级分辨率的块 与 当前分辨率的块 是否有相交，如果有就 SetHigherHTI */
	float curW = curr->startLongi, curE = curr->startLongi + curr->offsetLongi, curS = curr->startLati, curN = curr->startLati + curr->offsetLati;
	findResultHTI.clear(); 	numFRHTI=0;
	
	for(i=0; i<numHTI; i++)
	{
		if(allHTI[i].highLevel  <= higherLevel)
		{
			float tmpW = allHTI[i].startLongi, tmpE = allHTI[i].startLongi + allHTI[i].offsetLongi, tmpS = allHTI[i].startLati, tmpN = allHTI[i].startLati + allHTI[i].offsetLati;
			
			/* 高精度地景数据与当前精度地景数据有交叉，或者当前精度地景数据完全覆盖高精度地景数据 */
			if(	(((tmpW > curW) && (tmpW < curE))	|| ((tmpE > curW) && (tmpE < curE))	|| ((curW > tmpW) && (curW < tmpE))	|| ((curE > tmpW) && (curE < tmpE))) &&
				  (((tmpS > curS) && (tmpS < curN))	|| ((tmpN > curS) && (tmpN < curN)) || ((curS > tmpS) && (curS < tmpN))	|| ((curN > tmpS) && (curN < tmpN))) )
			{
				SetHigherHTIEx(&allHTI[i]);
			}
			
			/* 高精度地景数据完全覆盖了当前精度地景数据，这块地形就不做，直接退出 */
			else if( (tmpW <= curW) && (tmpE >= curE) && (tmpS <= curS) && (tmpN >= curN) )
			{
				return 1;
			}
			
			/* 高精度地景数据与当前精度地景数据没有重合 */
			else
			{
			}
			
		}
	}
	
	/* 如果有相覆盖的情况，则计算如何分割当前块curr */
	if(numFRHTI==0)
	{
		/* 没有查找到覆盖情况，退出 */
		return 0;
	}else{
		/* 有查找到覆盖情况 计算覆盖方法 */
		if(!CalcuHighOverlayByHigherEx(curr))	return 0;			
	}
	return 1;
}



/* 计算不同精度纹理数据叠加 */
/*  
   纹理叠加的情况相当复杂。
   按覆盖面来分类：有全覆盖的，有半覆盖的；
   按精度级别来分：有二层覆盖的，有三层覆盖的；
   按区间内覆盖数量分：有单数覆盖的，有多数覆盖的；
   现在新算法原则上支持所有的情况。
*/
int CalculateTextOverlay(HighTerrainInfo *curr)
{
	/* 如果当前已经是最高精度的了，不用计算，直接赋值退出 */
	if(curr->highLevel == 1)
	{
		SetCurrHTI(curr);
		return 1;
	}
	
	/* 查找高于当前精度的，并且在同一经纬度区间内的所有HTI ，先算 5米精度 覆盖 1米精度 */
	if(curr->highLevel == 5)
	{
		if(!DetectHTContainHigherHT(curr, 1))
		{
			SetCurrHTI(curr);
			return 1;
		}
	}

	/* 算 10 米精度 覆盖 1米 和 5米精度*/
	if(curr->highLevel == 10)
	{
		if (!DetectHTContainHigherHT(curr, 5))
		{
			SetCurrHTI(curr);
			return 1;
		}
	}

	return 1;
}


/* 根据整数经度和纬度信息，填充allHTI，计算出tileIDs */
int CalculateHTIFromLonLat(int i)
{
	float minLongi;
	float maxLongi;
	float minLati ;
	float maxLati ;
	unsigned int j;
	int lon, lat, lonp1, latp1;
	char fnPath[_MAX_PATH], tileStr[32];
	HighTerrainInfo currHTI;

	sprintf(fnPath, "%s\\%s", highTxtnode[i].pathName, highTxtnode[i].fileName);
	lon = (int)highTxtnode[i].minLong;
	while(lon <= (int)highTxtnode[i].maxLong)
	{
		lat = (int)highTxtnode[i].minLati;
		while(lat <= (int)highTxtnode[i].maxLati)
		{
			/* 计算当前所占区域 */
			lonp1 = lon + 1;	latp1 = lat + 1;
			if(highTxtnode[i].minLong >= (float)lon  ) minLongi = highTxtnode[i].minLong; else minLongi = (float)lon ;
			if(highTxtnode[i].minLati >= (float)lat  ) minLati  = highTxtnode[i].minLati; else minLati  = (float)lat ;
			if(highTxtnode[i].maxLong <= (float)lonp1) maxLongi = highTxtnode[i].maxLong; else maxLongi = (float)lonp1;
			if(highTxtnode[i].maxLati <= (float)latp1) maxLati  = highTxtnode[i].maxLati; else maxLati  = (float)latp1;
			
			if(highTxtnode[i].highLevel < 10)
				sprintf(tileStr, "%d%d0%d", lon, lat, highTxtnode[i].highLevel);
			else
				sprintf(tileStr, "%d%d%d", lon, lat, highTxtnode[i].highLevel);
			currHTI.startLongi  = minLongi;
			currHTI.startLati   = minLati ;
			currHTI.offsetLongi = maxLongi - minLongi;
			currHTI.offsetLati  = maxLati  - minLati ;
			currHTI.highLevel   = highTxtnode[i].highLevel;
			currHTI.tileID      = atoi(tileStr);
			strcpy(currHTI.texFname, fnPath);

			/* 如果面积太小就不值得一做了 */
			/* 2012-11-1 面积太小也得做，否则就容易有缝，这句话要注释掉。*/
			//////////////////////////	if((currHTI.offsetLongi < 0.01) || (currHTI.offsetLati < 0.01)) goto CHFL_NextBlock;

			/* 从用户高程数据里面查找，是否有匹配的高程文件，如果没有匹配的用户高程文件，则使用系统高程文件 */
			for(j=0; j<hHgtNum; j++)
			{
				if(	(minLongi >= highHgtnode[j].minLong) &&
					(maxLongi <= highHgtnode[j].maxLong) &&
					(minLati  >= highHgtnode[j].minLati) &&
					(maxLati  <= highHgtnode[j].maxLati) ) 
					break;
			}
			if(j==hHgtNum)	sprintf(currHTI.hgtFname, "%s\\D%d%d.tif", tempDirectory, lon, lat);
			else			sprintf(currHTI.hgtFname, "%s\\%s", highHgtnode[j].pathName, highHgtnode[j].fileName);

			/* 计算不同精度纹理数据叠加 */
			CalculateTextOverlay(&currHTI);
//////////CHFL_NextBlock:
			lat++;
		}//while(lat)
		lon++;
	}//while(lon)

	return 1;
}

/* 根据totalnode里面纹理文件扩展名为tif的所有信息，填充allHTI表格*/
int ReadNodeListToCreateHTI()
{
	
	/*  如果是高程，就退出 */
	if(isCurrTxtOrHgt == STATUS_HGT) return 1;
	
	/* 检查highTxtnode里面，这些是高精度纹理数据 */
	/* 先整理1m精度的地景纹理数据 */
	for(int i=0;i<hTxtNum;i++)
	{
		if(highTxtnode[i].highLevel != 1) continue;

		/* 将这一块内所有的tile都计算出来 */
		CalculateHTIFromLonLat(i);
		
	}//for(i)

	/* 再整理5m精度的地景纹理数据 */
	for(int i=0;i<hTxtNum;i++)
	{
		if(highTxtnode[i].highLevel != 5) continue;

		/* 将这一块内所有的tile都计算出来 */
		CalculateHTIFromLonLat(i);
		
	}//for(i)

	/* 再整理10m精度的地景纹理数据 */
	for(int i=0;i<hTxtNum;i++)
	{
		if(highTxtnode[i].highLevel != 10) continue;

		/* 将这一块内所有的tile都计算出来 */
		CalculateHTIFromLonLat(i);
		
	}//for(i)

	return 1;	
}

/* 写配置，把用户高分辨率地景制作情况写入配置文件。 */
int WriteHighInfo2Config(FILE *fo)
{
	unsigned int i;
	
	fprintf(fo, "\n%d",  numHTI);
	//再写入Map表
	for(i=0; i<numHTI; i++)
	{
		fprintf(fo, "\n%s %s %11.6f %11.6f %11.6f %11.6f %d %8d %2d", allHTI[i].texFname, allHTI[i].hgtFname, allHTI[i].startLongi, allHTI[i].startLati, allHTI[i].offsetLongi, allHTI[i].offsetLati, allHTI[i].highLevel, allHTI[i].tileID, allHTI[i].isOK);
	}
	return 1;
}

/* 读配置，从配置文件中读取用户高分辨率地景制作情况。 */
int ReadHighInfoFromConfig(FILE *fi)
{
	unsigned int i;
	
	ReadALineToArray(fi);	numHTI = atoi(cf[0]);	allHTI.clear();
	for(i=0; i<numHTI; i++)
	{
		HighTerrainInfo hti;
		ReadALineToArray(fi);
		strcpy(hti.texFname, cf[0]);
		strcpy(hti.hgtFname, cf[1]);
		hti.startLongi = atof(cf[2]);
		hti.startLati  = atof(cf[3]);
		hti.offsetLongi= atof(cf[4]);
		hti.offsetLati = atof(cf[5]);
		hti.highLevel  = atoi(cf[6]);	
		hti.tileID 	   = atoi(cf[7]);	
		hti.isOK	   = atoi(cf[8]);
		hti.ID		   = i;
		allHTI.push_back(hti);
	}
	return 1;
}

/************************************************
 *
 *    根据分辨率为50m的Geotiff纹理数据做低精度地景
 *
 ************************************************/
/* 设置当前LowHTI, 如果先前设过了，则跳过 */
void SetCurrHTIToLowHTI(HighTerrainInfo *curr)
{
	int i, iFlag;
	
	iFlag = 0;
	for(i=0;i<numLowHTI;i++)
	{
		if(	(fabs(lowHTI[i].startLongi  - curr->startLongi ) < 0.001) &&
			(fabs(lowHTI[i].startLati   - curr->startLati  ) < 0.001) && 
			(fabs(lowHTI[i].offsetLongi - curr->offsetLongi) < 0.001) && 
			(fabs(lowHTI[i].offsetLati  - curr->offsetLati ) < 0.001)  )
		{ iFlag = 1; break; }
	}
	
	if(!iFlag)
	{
		HighTerrainInfo hti;
		hti.startLongi  = curr->startLongi ; 
		hti.startLati   = curr->startLati  ; 
		hti.offsetLongi = curr->offsetLongi; 
		hti.offsetLati  = curr->offsetLati ; 
		hti.highLevel   = curr->highLevel  ; 
		hti.tileID      = curr->tileID     ; 
		strcpy(hti.texFname, curr->texFname);
		strcpy(hti.hgtFname, curr->hgtFname);
		hti.ID          = numLowHTI;
		hti.isOK        = 0;
		lowHTI.push_back(hti);	numLowHTI++;
	}
}

/* 先校验50m精度的纹理数据，根据文件名来校验 */
int CheckTextureNameTif50()
{
	int n, i, slen, ilat, ilon, jlat, jlon, flag;
	char clat[2][5], clon[2][5], dirname[24];
	int  lat[2], lon[2], minlat, minlon, maxlon, maxlat;

	for(n=0; n<lTxtNum; n++)
	{
		if(lowTxtnode[n].highLevel != 50) continue;
			
		/* 先获取文件名，不带扩展名的。 */
		memset(dirname, 0, sizeof(dirname));
		for(i=0; i<strlen(lowTxtnode[n].fileName); i++) 
		{
			if(lowTxtnode[n].fileName[i] == '.') break;
			dirname[i] = lowTxtnode[n].fileName[i];
		}
		
		/* 从字符串里面读取四个数据 */
		slen = strlen(dirname);  if(slen < 6) return 1;
		jlat = jlon = 0; ilat = ilon = flag = -1;
		memset(clat, 0, sizeof(clat)); memset(clon, 0, sizeof(clon));
		for(i=0; i<slen; i++)
		{
			if((flag == 1)&&(dirname[i]>=0x30)&&(dirname[i]<=0x39)) clat[ilat][jlat++] = dirname[i];
			if((flag == 0)&&(dirname[i]>=0x30)&&(dirname[i]<=0x39)) clon[ilon][jlon++] = dirname[i];
			if((dirname[i]=='N')||(dirname[i]=='n')||(dirname[i]=='S')||(dirname[i]=='s')) {flag = 1;  jlat = 0; ilat++; }
			if((dirname[i]=='W')||(dirname[i]=='w')||(dirname[i]=='E')||(dirname[i]=='e')) {flag = 0;  jlon = 0; ilon++; }
		}
		
		/* 判断当前路径是否合法 */
		if((strlen(clon[0])==0)||(strlen(clat[0])==0)||flag==-1)
			return 1;
		
		/* 转换成整形, 取小值 */
		lat[0] = atoi(clat[0]);  lat[1] = atoi(clat[1]); lon[0] = atoi(clon[0]); lon[1] = atoi(clon[1]);
		if(lat[0]<lat[1]) minlat = lat[0]; else minlat = lat[1]; if(lon[0]<lon[1]) minlon = lon[0]; else minlon = lon[1]; 
		if(lat[0]<lat[1]) maxlat = lat[1]; else maxlat = lat[0]; if(lon[0]<lon[1]) maxlon = lon[1]; else maxlon = lon[0]; 
		
		/* 填写数据结构 */
		lowTxtnode[n].minLong = minlon;
		lowTxtnode[n].minLati = minlat;
		lowTxtnode[n].maxLong = maxlon;
		lowTxtnode[n].maxLati = maxlat;
	}

	return 1;
}

/* 检测是否存在50m精度的低分辨率数据 */
int CheckIfHas50(int iW, int iS)
{
	int ret = -1;
	for(int i=0; i<lTxtNum; i++)
	{
		if(	((int)lowTxtnode[i].minLong <= iW) &&
			((int)lowTxtnode[i].minLati <= iS) &&
			((int)lowTxtnode[i].maxLong >= (iW + 1)) &&
			((int)lowTxtnode[i].maxLati >= (iS + 1)) )
		{	ret = i; break; }
	}

	return ret;
}

/* 先搜集数据 */
int GetHTIListOfR50()
{
	double minLon, minLat;
	char fnPath[_MAX_PATH], tileStr[32];
	HighTerrainInfo currHTI;
	int  iMinLon, iMinLat;
	
	/* 根据名字校验一下，以免出错。 */
	CheckTextureNameTif50();
	
	for(int i=0; i<lTxtNum; i++)
	{
		if(lowTxtnode[i].highLevel == 50)
		{
			sprintf(fnPath, "%s\\%s", lowTxtnode[i].pathName, lowTxtnode[i].fileName);
			minLon = lowTxtnode[i].minLong;
			while(minLon < lowTxtnode[i].maxLong)
			{
				minLat = lowTxtnode[i].minLati;
				while(minLat < lowTxtnode[i].maxLati)
				{
					
					/* 制作当前HTI*/
					iMinLon = (int)(minLon + 0.00001);
					iMinLat = (int)(minLat + 0.00001);
					sprintf(tileStr, "%d%d", iMinLon, iMinLat);
					currHTI.startLongi  = minLon;
					currHTI.startLati   = minLat ;
					currHTI.offsetLongi = 1.0;
					currHTI.offsetLati  = 1.0;
					currHTI.highLevel   = lowTxtnode[i].highLevel;
					currHTI.tileID      = atoi(tileStr);
					strcpy(currHTI.texFname, fnPath);
					sprintf(currHTI.hgtFname, "%s\\D%d%d.tif", tempDirectory, iMinLon, iMinLat);
					
					/* 写入低分辨率地景制作情况列表当中*/
					SetCurrHTIToLowHTI(&currHTI);
					
					/* 下一个纬度 */					
					minLat += 1.0;
				}
				/* 下一个经度 */
				minLon += 1.0;
			}			
		}//if
	}//for
	
	return 1;
}

int WriteLowInfo2Config(FILE *fo)
{
	unsigned int i;
	
	fprintf(fo, "\n%d",  numLowHTI);
	//再写入Map表
	for(i=0; i<numLowHTI; i++)
	{
		fprintf(fo, "\n%s %s %11.6f %11.6f %11.6f %11.6f %d %8d %2d", lowHTI[i].texFname, lowHTI[i].hgtFname, lowHTI[i].startLongi, lowHTI[i].startLati, lowHTI[i].offsetLongi, lowHTI[i].offsetLati, lowHTI[i].highLevel, lowHTI[i].tileID, lowHTI[i].isOK);
	}
	return 1;
}

int ReadLowInfoFromConfig(FILE *fi)
{
	unsigned int i;
	
	ReadALineToArray(fi);	numLowHTI = atoi(cf[0]);	lowHTI.clear();
	for(i=0; i<numLowHTI; i++)
	{
		HighTerrainInfo hti;
		ReadALineToArray(fi);
		strcpy(hti.texFname, cf[0]);
		strcpy(hti.hgtFname, cf[1]);
		hti.startLongi = atof(cf[2]);
		hti.startLati  = atof(cf[3]);
		hti.offsetLongi= atof(cf[4]);
		hti.offsetLati = atof(cf[5]);
		hti.highLevel  = atoi(cf[6]);	
		hti.tileID 	   = atoi(cf[7]);	
		hti.isOK	   = atoi(cf[8]);
		hti.ID		   = i;
		lowHTI.push_back(hti);
	}
	return 1;
}

/* 制作低分辨率地景的流程 */
int DoLowResolutionTerrain()
{
	int n;
	PerZone mPZ;
	char msg[STRING_STANDARD_LEN];
	
	/* 对所有LowHTI做循环*/
	for(n=0; n<numLowHTI; n++)
	{
		/* 检查是否做过，做过的就跳过 */
		if(lowHTI[n].isOK) continue;
		
		/* 先收集本平方度内所有的高精度区块，结果放到 findResultHTI里面。 */
		mPZ.lon = lowHTI[n].startLongi; mPZ.lat = lowHTI[n].startLati; mPZ.isOK = 0;

		/* 先制作这个平方度内的基本低分辨率地景 *//* 生成基础地形 */
		sprintf(msg, "开始制作经度%d 纬度%d 区域的基本地景.", mPZ.lon, mPZ.lat);
		MessageToDialog(msg);
		currTi.minLongi = (float)(mPZ.lon);
		currTi.minLati  = (float)(mPZ.lat);
		currTi.maxLongi = (float)(mPZ.lon + 1);
		currTi.maxLati  = (float)(mPZ.lat + 1);
		currPZ = &mPZ;

		/* 首先要拼出低精度地景的高程数据 */
		if(!CreatePerDegreeTerrain(2))
		{
			sprintf(msg, "制作经度%d 纬度%d 区域的高程出错.", mPZ.lon, mPZ.lat);
			MessageToDialog(msg);
			lowHTI[n].isOK = 2;
			continue;
		}

		/* 拼装一个findResultHTI */
		strcpy(findResultHTI[0].texFname, lowHTI[n].texFname);  
		strcpy(findResultHTI[0].hgtFname, lowHTI[n].hgtFname);  
		findResultHTI[0].startLongi = lowHTI[n].startLongi	;
		findResultHTI[0].startLati  = lowHTI[n].startLati 	;
		findResultHTI[0].offsetLongi= lowHTI[n].offsetLongi	;
		findResultHTI[0].offsetLati = lowHTI[n].offsetLati 	;
		findResultHTI[0].highLevel  = lowHTI[n].highLevel 	;
		findResultHTI[0].tileID     = lowHTI[n].tileID    	;
		findResultHTI[0].ID         = lowHTI[n].ID        	;
		findResultHTI[0].isOK       = lowHTI[n].isOK      	;

		/* 装配vpb批处理命令。*/
		MessageToDialog("调用VPB生成地景");
		if(!DoHighVpbAndWriteStg(0))
		{
			MessageToDialog("执行VPB失败\n");
			lowHTI[n].isOK = 4;
			continue;
		}

		/* 记录这个块已经做好 */
		lowHTI[n].isOK = 1;

		/* 把allHTI写回配置文件 */
		CheckConfig(1);
	}
	return 1;
}

/*============================================================================================
 *  下面介绍了场景中任意一点高度的计算过程。
 *============================================================================================ */

/************************************************
 *
 *    管理openedTiles数组。这是为计算地景上某一点的高程所用。数组内保存着一个经纬度内的所有三维地景数据
 *	  计算高程的方法是：读入这一片地景的三维数据，然后在计算高程确定某一点(由经度、纬度决定)，做一个垂直线与地景做碰撞检测，
 *    检测返回的结果就是这一点的高程。
 *
 ************************************************/
 
 /* 插入一个地景文件到 高程计算组rootNode 里面 */
int InsertATile(int tile, osg::ref_ptr<osg::Node> node)
{
	OpenTiles ots;
	
	ots.tile = tile;
	ots.node = node;
	openedTiles.push_back(ots);
	
	rootNode->addChild(node);
	return 1;
}


/* 将当前地景目录下所有的根地景文件插入到 高程计算组rootNode 里面。*/
int InsertTilesInCurrDir(double lon, double lat)
{
	char fullpath[_MAX_PATH], head[24], lowPart[5], iveHead[48];
	int io, dio, dfio, dirIo, j, iio, iveio, iFlag = 0;
	struct _finddata_t sdff;
	struct _finddata_t iveff;

	/* 首先获取当前目录 */
	Bucket tile_bucket;
	set_bucket(&tile_bucket, lon, lat); 
	string root = iveroot;
	char base[32];
	gen_base_path(&tile_bucket, base);
	string base_path = root + "\\" + base; 
	sprintf(fullpath, "%s", base_path.c_str());
	int ilon = (int)lon;
	int ilat = (int)lat;
	int tileID = ilon*100 + ilat;
	
	/* 进入这个目录，然后搜索当前目录下扩展名为ive的所有文件。 */
	io = chdir(fullpath);
	if(io) return 0;
	
	/* 查找当前路径下的所有文件； */
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;
	
	/* 获得当前tileID */
	int currTileID = 0;
	if(openedTiles.size()>0) currTileID = openedTiles[0].tile;

	while(!dfio)
	{

		//如果不是子目录也不是文件
		if((sdff.name[0] == '.')) goto ITICD_GET_NEXT_SUBDIR;

		//如果是子目录，进入当前子目录
		if(sdff.attrib == _A_SUBDIR)
		{
			_getcwd(fullpath, _MAX_PATH);

			//进入子目录
			dirIo = chdir(sdff.name);
			if(dirIo == -1) goto ITICD_GET_NEXT_SUBDIR;
			
			for(j=0; j<strlen(sdff.name); j++)
			{ 
				if(sdff.name[j] != '_')	head[j] = sdff.name[j];
				else	break;
			}
			head[j++] = 0; lowPart[0] = sdff.name[j]; lowPart[1] = 0;
			tileID = atoi(head);

			/* 增加判断低精度地形的功能	*/
			if(strlen(head) < 8)
			{
				sprintf(head, "%s_%s", head, lowPart);
			}
			
			/* 比较当前目录头的名称是否和当前tileID一样，如果一样，说明这个目录无效，直接读取下一个目录 */
			if(currTileID == tileID)	goto ITICD_OUT_TO_UPPER_DIR;
			
			int level = 3;	//currTi.lod - 2;
			int isSucess = 1;
			while(1)
			{
				sprintf(iveHead, "%s_L%d*.ive", head, level);
				iio = _findfirst(iveHead, &iveff);
				if((!iio)||(iio==-1)) level--;
				else break;
				if(level < 0) {	isSucess = 0;	break;	}
			}

			/* 如果查不到任何文件 */
			if(isSucess == 0)	{	_findclose(iio);	break;	}
			
			/* 在将整个目录下的所有地景文件加入到组之前，先初始化管理数组 */
			RemoveAllTileNodes();
			
			iveio = 0;
			while(!iveio)
			{
				/* 打开文件 */
				osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(iveff.name);
				if(node == NULL) goto ITICD_Continue;

				/* 插入新打开的tile */
				InsertATile(tileID, node);
ITICD_Continue:
				/* 寻找下一个 */
				iveio = _findnext(iio, &iveff);
			}
			_findclose(iio);
			
			double elev = GetElevationByLonLat(lon, lat);
			if(elev != -1.0)
			{
				iFlag = 1;
			}
ITICD_OUT_TO_UPPER_DIR:			
			//退回到上一层目录
			dirIo = chdir("..");
			
			if(iFlag) break;
		}
ITICD_GET_NEXT_SUBDIR:

		if(iFlag) break;

		/* 寻找下一个 */
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);

	if(iFlag == 0)	return 0; else	return 1;
}

/* 从 openedTiles 高程计算管理数组中移除所有的数据。 */
int RemoveAllTileNodes()
{
	unsigned int i;
	
	for(i=0; i<openedTiles.size(); i++)	
	{
		if(openedTiles[i].node)
		{
			rootNode->removeChild(openedTiles[i].node);
			openedTiles[i].tile = 0;
		}
	}
	openedTiles.clear();
	
	return 1;
}

/* 根据一个点，来加载相应rootNode的地景， */
int LoadRootNodeByVertex(osg::Vec3 v)
{
	Carton ct;	ct.x = v[0]; ct.y = v[1]; ct.z = v[2];
	Geodesy gd;
	CartToGeod(&ct, &gd);
	
	/* 直接查找 经度纬度.ive文件。*/
	if(!InsertTilesInCurrDir(gd._lon, gd._lat))
		return 0;

	return 1;
}


/**************************************************
 *
 *   获取高度。获取地景上某一点（由经纬度决定）的高度，方法在openedTiles的地方有介绍。
 *
 **************************************************/
double GetElevation(double lon, double lat, double pre)
{
	double alt = 1.0; 

#if 0
	alt = GetElevationFromHgt(lon, lat);
	return alt;
#endif

	/* 直接查找 经度纬度.ive文件	，如果找不到这块地景，就返回pre，作为缺省高度。*/
	if(rootNode->getNumChildren() == 0)
	{
		if(!InsertTilesInCurrDir(lon, lat))
		{
			alt = GetElevationFromHgt(lon, lat);
			return alt;
		}
	}
	
	/* 创建相交线段，是一根在经度lon和纬度lat的垂直的线段。用这根线段与地景做相交。*/
	Geodesy stg; stg._lon = lon; stg._lat = lat; stg._elevation =  9999999.0;
	Carton  stc;
	GeodToCart(&stg, &stc);
	osg::Vec3d start; start.set(stc.x, stc.y, stc.z); 
	Geodesy eng; eng._lon = lon; eng._lat = lat; eng._elevation =  -9999999.0;
	Carton  enc;
	GeodToCart(&eng, &enc);
	osg::Vec3d end;   end.set(enc.x, enc.y, enc.z); 

	while(1)
	{
		/* 创建相交器  Create intersector  */
		osgUtil::IntersectVisitor intersectVisitor; 
		osg::LineSegment* lineSegment = new osg::LineSegment(start, end); 
		intersectVisitor.addLineSegment(lineSegment); 
		rootNode->accept(intersectVisitor); 
		
		/* 获取相交结果 Get intersections */
		bool hits = intersectVisitor.hits(); 
		if (hits) 
		{ 
			int nHits = intersectVisitor.getNumHits(lineSegment); 
			for (int i = 0; i < nHits; ++i) 
			{ 
				const osgUtil::Hit& hit = intersectVisitor.getHitList(lineSegment)[i]; 
				osg::Vec3d point; 
				point = hit.getWorldIntersectPoint(); 
				Carton pntc; pntc.x = point.x(); pntc.y = point.y(); pntc.z = point.z();
				Geodesy pntg;
				CartToGeod(&pntc, &pntg);
				double elevation = pntg._elevation; 
				if (alt < elevation) { alt = elevation; } 
			} 
			break;
		}else{
			if(!InsertTilesInCurrDir(lon, lat))
			{
				alt = GetElevationFromHgt(lon, lat);
				return alt;
			}
		}
	}
	
	if(alt < 10.0)
		alt = GetElevationFromHgt(lon, lat);
	
	return alt;
}


double GetElevationByLonLat(double lon, double lat)
{
	double alt = -1.0;

	/* 创建相交线段，是一根在经度lon和纬度lat的垂直的线段。用这根线段与地景做相交。*/
	Geodesy stg; stg._lon = lon; stg._lat = lat; stg._elevation =  99999999.0;
	Carton  stc;
	GeodToCart(&stg, &stc);
	osg::Vec3d start; start.set(stc.x, stc.y, stc.z); 
	Geodesy eng; eng._lon = lon; eng._lat = lat; eng._elevation =  -99999999.0;
	Carton  enc;
	GeodToCart(&eng, &enc);
	osg::Vec3d end;   end.set(enc.x, enc.y, enc.z); 

	/* 创建相交器  Create intersector  */
	osgUtil::IntersectVisitor intersectVisitor; 
	osg::LineSegment* lineSegment = new osg::LineSegment(start, end); 
	intersectVisitor.addLineSegment(lineSegment); 
	rootNode->accept(intersectVisitor); 
	
	/* 获取相交结果 Get intersections */
	bool hits = intersectVisitor.hits(); 
	if (hits) 
	{ 
		int nHits = intersectVisitor.getNumHits(lineSegment); 
		for (int i = 0; i < nHits; ++i) 
		{ 
			const osgUtil::Hit& hit = intersectVisitor.getHitList(lineSegment)[i]; 
			osg::Vec3d point; 
			point = hit.getWorldIntersectPoint(); 
			Carton pntc; pntc.x = point.x(); pntc.y = point.y(); pntc.z = point.z();
			Geodesy pntg;
			CartToGeod(&pntc, &pntg);
			double elevation = pntg._elevation; 
			if (alt < elevation) { alt = elevation; } 
		} 
	}else{
		alt = -1.0;
	}
	
	return alt;
}



/*============================================================================================
 *  从下面开始到 CheckTifInTotalNode 函数之前，是用JPEG文件生成全国地景方法的大部分过程区。
 *  目前基本较少用到。
 *============================================================================================ */

/************************************************
 *
 *    转换Tiff到GeoTiff，这部分过程是使用JPEG图片制作地景的一个中间步骤。JPEG->Tiff, Tiff->GeoTiff
 *
 ************************************************/
typedef struct _TIFFTAG_
{
	short    idx;
	short    type;
	fpos_t   size;
	fpos_t   ptr;
}TiffTag;

/* 读取二进制字符串，要把0x1A这个字符也正确地都进来。*/
void nfread(unsigned char *pBuff, int num, int size, FILE *pInput)
{
	for(unsigned int j = 0; j < num * size; j++)
	{	
		int c = fgetc(pInput);

		/* 下面一条语句的添加实在是一种无奈之举，系统包里有个bug。
		   当读到0x1A字节的时候就算文件结束，结果0x1A后面的数据就没有被读进来。
		   在这里做此修改。把该读的都读进来。2012-04-24
		 */

		if(pInput->_cnt < 0x9f) pInput->_cnt = 0x9f;


		*pBuff++ = (unsigned char)c;
	}
}


/* 从一个GeoTiff文件中根据Geo类型读取Geo地理信息，支持超大文件 */
int ReadGeoTiffTags(char *fn)
{
	TiffTag tags[32];
	int  header, i, flag, Lflag;
	fpos_t iFileSize, iDataPtr;
	fpos_t tagnum = 0;
	double modeltiepoints[] = {0.0,0.0,0.0,0.0,0.0,0.0};
	double modelpixelscale[] = {0.0,0.0,1.0};
	FILE *pInput;
	char path[_MAX_PATH];
	char simfn[32], tfwfn[32];
	
	i=0; while(fn[i]!='.') {simfn[i]=fn[i++];} simfn[i]=0;
	memset(tags, 0, sizeof(tags)); Lflag = 0;
	pInput = fopen(fn, "rt");
	if(pInput == NULL) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能读取文件 %s .", fn);
		MessageToDialog(msg);
		return 0;
	}
	
	/*获取文件长度*/
	fseek(pInput, 0, SEEK_END);
	fgetpos(pInput, &iFileSize);
	rewind(pInput);
	
	/*读写Header 和Length, 设置当前pos到Tag区*/
	iDataPtr = 0;
	nfread((unsigned char *)&header, 4, 1, pInput);
	nfread((unsigned char *)&iDataPtr, 4, 1, pInput);
	
	/*判断Tag的正确位置*/
	if(iDataPtr == 8)
	{
		fread(&iDataPtr, 8, 1, pInput); Lflag = 1;
	}
	fsetpos(pInput, &iDataPtr);
	
	/*读Tag区，将Tag区信息读入数组tags。*/
	if(!Lflag)
	{
		/* 一般小文件的Tag区格式 */
		fread(&tagnum, 1, 2, pInput);
		for(i=0; i<tagnum; i++)
		{
			fread(&tags[i].idx,  2, 1, pInput);
			fread(&tags[i].type, 2, 1, pInput);
			fread(&tags[i].size, 4, 1, pInput);
			fread(&tags[i].ptr,  4, 1, pInput);
		}
	}else{
		/* 超大文件的Tag区格式 */
		fread(&tagnum, 1, 8, pInput);
		for(i=0; i<tagnum; i++)
		{
			fread(&tags[i].idx,  2, 1, pInput);
			fread(&tags[i].type, 2, 1, pInput);
			fread(&tags[i].size, 8, 1, pInput);
			fread(&tags[i].ptr,  8, 1, pInput);
		}
	}
	
	/* 在tags数组里面查找并从文件中读取描述地理经纬度信息的TIFFTAG_GEOTIEPOINTS类型，和描述图像分辨率信息的TIFFTAG_GEOPIXELSCALE */
	/* 经纬度信息保存在 数组modeltiepoints中，分辨率信息保存在数组modelpixelscale中。  */	
	flag = 0;
	for(i=0; i<tagnum; i++)
	{
		char dBuff[6][8];
		char *pBuff = &dBuff[0][0];
		if((unsigned short)tags[i].idx == TIFFTAG_GEOPIXELSCALE)
		{
			flag = 1;
			iDataPtr = tags[i].ptr;
			fsetpos(pInput, &iDataPtr);

			/* 读取3*8个字符 */
			nfread((unsigned char *)pBuff, 3, 8, pInput);
			for(unsigned int j = 0; j < 3; j++)
			{
				double dtmp;
				memcpy(&dtmp, dBuff[j], 8);
				modelpixelscale[j] = dtmp;
			}
		}
		if((unsigned short)tags[i].idx == TIFFTAG_GEOTIEPOINTS)
		{
			flag = 1;
			iDataPtr = tags[i].ptr;
			fsetpos(pInput, &iDataPtr);

			/* 读取6*8个字符 */
			nfread((unsigned char *)pBuff, 6, 8, pInput);
			for(unsigned int j = 0; j < 6; j++)
			{
				double dtmp;
				memcpy(&dtmp, dBuff[j], 8);
				modeltiepoints[j] = dtmp;
			}
		}
	}
	fclose(pInput);
	
	/* 如果当前文件不包含地理信息，则检查一下是否有同文件名的tfw文件。 */
	if(!flag)
	{
		FILE *tfw;

		/* 检查并打开同目录下同文件名的tfw文件。*/
		sprintf(tfwfn, "%s.tfw", simfn);
		tfw = fopen(tfwfn, "rt");
		if(tfw == NULL)
		{
			char msg[STRING_STANDARD_LEN];
			sprintf(msg, "文件 %s 中不包含地理信息。", fn);
			MessageToDialog(msg);
			return 0;
		}
		
		/* 读取tfw文件 */
		ReadALineToArray(tfw);	modelpixelscale[0] = atof(cf[0]);
		ReadALineToArray(tfw);
		ReadALineToArray(tfw);
		ReadALineToArray(tfw);	modelpixelscale[1] = -(atof(cf[0]));
		ReadALineToArray(tfw);	modeltiepoints[3]  = atof(cf[0]);
		ReadALineToArray(tfw);	modeltiepoints[4]  = atof(cf[0]);
		fclose(tfw);
	}

	/* 根据上两个数组，计算这个GeoTiff文件所包含的地理范围(最大最小经纬度)，记录到数据列表中。 */
	if(isCurrTxtOrHgt == STATUS_TEXTURE)
	{
		HighNode hn;
		strcpy(hn.fileName, fn);
		hn.highLevel  = currHighLevel;
		hn.image_width  = tags[0].ptr;
		hn.image_height = tags[1].ptr;
		hn.minLong = modeltiepoints[3];
		hn.maxLati = modeltiepoints[4];
		hn.maxLong = modeltiepoints[3] + tags[0].ptr * modelpixelscale[0];
		hn.minLati = modeltiepoints[4] - tags[1].ptr * modelpixelscale[1];

		/* 如果经纬度数据大于360.0, 有可能使用的是秒的单位，要在原有数字基础上 除以 3600.0 */
		if(fabs(hn.minLong) > 360.0)	hn.minLong /= 3600.0;
		if(fabs(hn.maxLati) > 360.0)	hn.maxLati /= 3600.0;
		if(fabs(hn.maxLong) > 360.0)	hn.maxLong /= 3600.0;
		if(fabs(hn.minLati) > 360.0)	hn.minLati /= 3600.0;

		_getcwd(path, _MAX_PATH);
		strcpy(hn.pathName, path);
		hn.ID = hTxtNum;
		highTxtnode.push_back(hn);	hTxtNum++;
	}

	if(isCurrTxtOrHgt == STATUS_HGT)
	{
		HighNode hn;
		strcpy(hn.fileName, fn);
		hn.image_width  = tags[0].ptr;
		hn.image_height = tags[1].ptr;
		hn.minLong = modeltiepoints[3];
		hn.maxLati = modeltiepoints[4];
		hn.maxLong = modeltiepoints[3] + tags[0].ptr * modelpixelscale[0];
		hn.minLati = modeltiepoints[4] - tags[1].ptr * modelpixelscale[1];

		/* 如果经纬度数据大于360.0, 有可能使用的是秒的单位，要在原有数字基础上 除以 3600.0 */
		if(fabs(hn.minLong) > 360.0)	hn.minLong /= 3600.0;
		if(fabs(hn.maxLati) > 360.0)	hn.maxLati /= 3600.0;
		if(fabs(hn.maxLong) > 360.0)	hn.maxLong /= 3600.0;
		if(fabs(hn.minLati) > 360.0)	hn.minLati /= 3600.0;

		_getcwd(path, _MAX_PATH);
		strcpy(hn.pathName, path);
		hn.ID = hHgtNum;
		highHgtnode.push_back(hn);	hHgtNum++;
	}

	return 1;	
}

typedef struct _TIFFTAGSHORT_
{
	short    idx;
	short    type;
	int      size;
	int      ptr;
}TiffTagShort;

/* 将纹理Tiff文件转换为GeoTiff文件，根据相关地理信息结构mergedTif。方法是读取原来Tiff后插入Geo相关Tags。 */
int convertTextureGeotiff()
{
#if 1
	FILE *pInput;
	FILE *pOutput;
	TiffTagShort tags[32];
	int  header, len, i, tmp, iDataLen;
	fpos_t iFileSize, iDataPtr;
	short tagnum;
	
	memset(tags, 0, sizeof(tags));
	/*打开文件*/
	sprintf(temTextGeoFile  , "%s\\%s%d%d.%s",  tempDirectory, TEMP_TEXTGEOFILE, currPZ->lon, currPZ->lat, TEMP_EXTENDNAME);
	if((pInput=fopen(temTextFile,"rb"))==NULL||						
		(pOutput=fopen(temTextGeoFile,"wb"))==NULL)
	{
		MessageToDialog("打开文件失败\n");
		currPZ->isOK = STATUS_TEXTGEOBAD;
		return 0;
	}
	
	/*获取文件长度*/
	fseek(pInput, 0, SEEK_END);
	fgetpos(pInput, &iFileSize);
	rewind(pInput);
	
	/*读写Header 和Length*/
	fread(&header, 4, 1, pInput);	fwrite(&header, 4, 1, pOutput);
	fread(&len, 4, 1, pInput);		fwrite(&len, 4, 1, pOutput);
	
	/*读写文件体*/
	unsigned char * pcInputLine= new unsigned char[MBUFFSIZE];
	iDataLen = len - 8;
	while(iDataLen > MBUFFSIZE)
	{
		fread(pcInputLine, 1, MBUFFSIZE, pInput);
		fwrite(pcInputLine, 1, MBUFFSIZE, pOutput);
		iDataLen -= MBUFFSIZE;
	}
	fread(pcInputLine, 1, iDataLen, pInput);
	fwrite(pcInputLine, 1, iDataLen, pOutput);

	/*读Tag区*/
	fread(&tagnum, 1, 2, pInput);
	for(i=0; i<tagnum; i++)
	{
		fread(&tags[i].idx,  2, 1, pInput);
		fread(&tags[i].type, 2, 1, pInput);
		fread(&tags[i].size, 4, 1, pInput);
		fread(&tags[i].ptr,  4, 1, pInput);
	}
	fread(&tmp, 4, 1, pInput);
	fgetpos(pInput, &iDataPtr);  
	iDataLen = iFileSize - iDataPtr;
	fread(pcInputLine, 1, iDataLen, pInput);
	
	/* 准备Geo数据 */
	unsigned short GeoKeyDirectory[] = {
		1,1,0,7,				// 1 - version 1, revision 1.2,  6 keys defined
		GTRasterTypeGeoKey  ,   0,1, 2,		// RasterPixelIsPoint
		GTModelTypeGeoKey   ,   0,1, 2,		// ModelTypeGeographic
		GeographicTypeGeoKey,   0,1, 4326,	// GCS_WGS_84
		GeogAngularUnitsGeoKey, 0,1, 9102,
		GeogSemiMajorAxisGeoKey,0x87B0, 1, 0,
		GeogSemiMinorAxisGeoKey,0x87B0, 1, 1,
		GeogInvFlatteningGeoKey,0x87B0, 1, 2
	};
	
	double GeoDoubleKey[] = {
		6378137.000000,
		6356752.314245,
		298.257224
	};

	// encode image
	double modeltiepoints[] = {0.0,0.0,0.0,0.0,0.0,0.0};

	double modelpixelscale[] = {0.0,0.0,1.0};

	int m_unImgWidth = tags[0].ptr;
	int m_unImgLength = tags[1].ptr;

	modeltiepoints[3] = mergedTif.minLong;
	modeltiepoints[4] = mergedTif.maxLati;
	
	modelpixelscale[0] = (mergedTif.maxLong - mergedTif.minLong) / m_unImgWidth;
	modelpixelscale[1] = (mergedTif.maxLati - mergedTif.minLati) / m_unImgLength;
	
	/*准备Tag数据*/
	for(i=0; i<tagnum; i++)
	{
		if((tags[i].ptr > len) && (tags[i].ptr < len + 0x10000))
		{
			tags[i].ptr += 4 * sizeof(TiffTagShort);
		}
	}
	int curPos = iFileSize + 4 * sizeof(TiffTagShort);
	tags[i].idx = TIFFTAG_GEOPIXELSCALE; 	tags[i].type = 0xC; tags[i].size = sizeof(modelpixelscale)/sizeof(double); 		tags[i].ptr = curPos; i++;
	curPos += sizeof(modelpixelscale);
	tags[i].idx = TIFFTAG_GEOTIEPOINTS; 	tags[i].type = 0xC; tags[i].size = sizeof(modeltiepoints)/sizeof(double); 		tags[i].ptr = curPos; i++;
	curPos += sizeof(modeltiepoints);
	tags[i].idx = TIFFTAG_GEOKEYDIRECTORY; 	tags[i].type = 0x3; tags[i].size = sizeof(GeoKeyDirectory)/sizeof(unsigned short);tags[i].ptr = curPos; i++;
	curPos += sizeof(GeoKeyDirectory);
	tags[i].idx = TIFFTAG_GEODOUBLEPARAMS; 	tags[i].type = 0xC; tags[i].size = sizeof(GeoDoubleKey)/sizeof(double); 			tags[i].ptr = curPos; i++;
	
	/*写Tag区*/
	tagnum += 4;
	fwrite(&tagnum, 2, 1, pOutput);
	for(i=0;i<tagnum;i++)
	{
		fwrite(&tags[i].idx,  2, 1, pOutput);
		fwrite(&tags[i].type, 2, 1, pOutput);
		fwrite(&tags[i].size, 4, 1, pOutput);
		fwrite(&tags[i].ptr,  4, 1, pOutput);
	}
	fwrite(&tmp, 4, 1, pOutput);
	fwrite(pcInputLine, 1, iDataLen, pOutput);
	fwrite(modelpixelscale, sizeof(double), sizeof(modelpixelscale)/sizeof(double), pOutput);
	fwrite(modeltiepoints, sizeof(double), 	sizeof(modeltiepoints)/sizeof(double), pOutput);
	fwrite(GeoKeyDirectory, sizeof(unsigned short), sizeof(GeoKeyDirectory)/sizeof(unsigned short), pOutput);
	fwrite(GeoDoubleKey, sizeof(double), sizeof(GeoDoubleKey)/sizeof(double), pOutput);
	
	/*关闭文件*/
	delete pcInputLine;
	fclose(pInput);
	fclose(pOutput);
#endif	
	return 1;	
}

/* 将高程Tiff文件转换为GeoTiff文件，根据相关地理信息结构mergedHgt。 */
int convertHgtGeotiff()
{
#if 1
	FILE *pInput;
	FILE *pOutput;
	TiffTagShort tags[32];
	int  header, len, i, tmp, iDataLen;
	fpos_t iFileSize, iDataPtr;
	short tagnum;
	
	memset(tags, 0, sizeof(tags));
	/*打开文件*/
	sprintf(temHgtGeoFile   , "%s\\%s%d%d.%s",  tempDirectory, TEMP_HGTGEOFILE	, currPZ->lon, currPZ->lat, TEMP_EXTENDNAME);
	if((pInput=fopen(temHgtFile,"rb"))==NULL||						
		(pOutput=fopen(temHgtGeoFile,"wb"))==NULL)
	{
		MessageToDialog("打开文件失败\n");
		currPZ->isOK = STATUS_HGTGEOBAD;
		return 0;
	}
	
	/*获取文件长度*/
	fseek(pInput, 0, SEEK_END);
	fgetpos(pInput, &iFileSize);
	rewind(pInput);
	
	/*读写Header 和Length*/
	fread(&header, 4, 1, pInput);	fwrite(&header, 4, 1, pOutput);
	fread(&len, 4, 1, pInput);		fwrite(&len, 4, 1, pOutput);
	
	/*读写文件体*/
	unsigned char * pcInputLine= new unsigned char[MBUFFSIZE];
	iDataLen = len - 8;
	while(iDataLen > MBUFFSIZE)
	{
		fread(pcInputLine, 1, MBUFFSIZE, pInput);
		fwrite(pcInputLine, 1, MBUFFSIZE, pOutput);
		iDataLen -= MBUFFSIZE;
	}
	fread(pcInputLine, 1, iDataLen, pInput);
	fwrite(pcInputLine, 1, iDataLen, pOutput);

	/*读Tag区*/
	fread(&tagnum, 1, 2, pInput);
	for(i=0; i<tagnum; i++)
	{
		fread(&tags[i].idx,  2, 1, pInput);
		fread(&tags[i].type, 2, 1, pInput);
		fread(&tags[i].size, 4, 1, pInput);
		fread(&tags[i].ptr,  4, 1, pInput);
	}
	fread(&tmp, 4, 1, pInput);
	fgetpos(pInput, &iDataPtr);  
	iDataLen = iFileSize - iDataPtr;
	fread(pcInputLine, 1, iDataLen, pInput);
	
	/* 准备Geo数据 */
	unsigned short GeoKeyDirectory[] = {
		1,1,0,8,				// 1 - version 1, revision 1.2,  6 keys defined
		GTRasterTypeGeoKey  ,   0,1, 2,		// RasterPixelIsPoint
		GTModelTypeGeoKey   ,   0,1, 2,		// ModelTypeGeographic
		GeographicTypeGeoKey,   0,1, 4326,	// GCS_WGS_84
		GeogAngularUnitsGeoKey, 0,1, 9102,
		GeogSemiMajorAxisGeoKey,0x87B0, 1, 0,
		GeogSemiMinorAxisGeoKey,0x87B0, 1, 1,
		GeogInvFlatteningGeoKey,0x87B0, 1, 2,
		VerticalUnitsGeoKey,    0,1, 9001		// Linear_Meter
	};
	
	double GeoDoubleKey[] = {
		6378137.000000,
		6356752.314245,
		298.257224
	};

	/* 编码Geo数据 */
	double modeltiepoints[] = {0.0,0.0,0.0,0.0,0.0,0.0};
	double modelpixelscale[] = {0.0,0.0,1.0};

	int m_unImgWidth = tags[0].ptr;
	int m_unImgLength = tags[1].ptr;

	modeltiepoints[3] = mergedHgt.minLong;
	modeltiepoints[4] = mergedHgt.maxLati;
	
	modelpixelscale[0] = (mergedHgt.maxLong - mergedHgt.minLong) / m_unImgWidth;
	modelpixelscale[1] = (mergedHgt.maxLati - mergedHgt.minLati) / m_unImgLength;
	
	/*准备Tag数据*/
	for(i=0; i<tagnum; i++)
	{
		if((tags[i].ptr > len) && (tags[i].ptr < len + 0x10000))
		{
			tags[i].ptr += 4 * sizeof(TiffTagShort);
		}
	}
	int curPos = iFileSize + 4 * sizeof(TiffTagShort);
	tags[i].idx = TIFFTAG_GEOPIXELSCALE; 	tags[i].type = 0xC; tags[i].size = sizeof(modelpixelscale)/sizeof(double); 		tags[i].ptr = curPos; i++;
	curPos += sizeof(modelpixelscale);
	tags[i].idx = TIFFTAG_GEOTIEPOINTS; 	tags[i].type = 0xC; tags[i].size = sizeof(modeltiepoints)/sizeof(double); 		tags[i].ptr = curPos; i++;
	curPos += sizeof(modeltiepoints);
	tags[i].idx = TIFFTAG_GEOKEYDIRECTORY; 	tags[i].type = 0x3; tags[i].size = sizeof(GeoKeyDirectory)/sizeof(unsigned short);tags[i].ptr = curPos; i++;
	curPos += sizeof(GeoKeyDirectory);
	tags[i].idx = TIFFTAG_GEODOUBLEPARAMS; 	tags[i].type = 0xC; tags[i].size = sizeof(GeoDoubleKey)/sizeof(double); 			tags[i].ptr = curPos; i++;
	
	/*写Tag区*/
	tagnum += 4;
	fwrite(&tagnum, 2, 1, pOutput);
	for(i=0;i<tagnum;i++)
	{
		fwrite(&tags[i].idx,  2, 1, pOutput);
		fwrite(&tags[i].type, 2, 1, pOutput);
		fwrite(&tags[i].size, 4, 1, pOutput);
		fwrite(&tags[i].ptr,  4, 1, pOutput);
	}
	fwrite(&tmp, 4, 1, pOutput);
	fwrite(pcInputLine, 1, iDataLen, pOutput);
	fwrite(modelpixelscale, sizeof(double), sizeof(modelpixelscale)/sizeof(double), pOutput);
	fwrite(modeltiepoints, sizeof(double), 	sizeof(modeltiepoints)/sizeof(double), pOutput);
	fwrite(GeoKeyDirectory, sizeof(unsigned short), sizeof(GeoKeyDirectory)/sizeof(unsigned short), pOutput);
	fwrite(GeoDoubleKey, sizeof(double), sizeof(GeoDoubleKey)/sizeof(double), pOutput);
	
	/*关闭文件*/
	delete pcInputLine;
	fclose(pInput);
	fclose(pOutput);
#endif
	return 1;	
}


/***********************************************
 *
 *   合并图像。合并若干个JPEG图像成为一个Tiff图像。是使用JPEG图像制作地景的关键步骤。
 *   按照命令行的方式输入。先将要合并的若干图像，合并方式，以及输出的文件，组成一个命令数组。然后调用这个过程。
 *
 ***********************************************/
int ConvertMain(int gc,char **gv)
{
#ifdef SUPPORT_JPEG_TEXTURE
  ExceptionInfo
    *exception;

  ImageInfo
    *image_info;

  MagickBooleanType
    status;

  MagickCoreGenesis(*gv,MagickTrue);
  exception=AcquireExceptionInfo();
  image_info=AcquireImageInfo();
  status=MagickCommandGenesis(image_info,ConvertImageCommand,gc,gv,
    (char **) NULL,exception);
  image_info=DestroyImageInfo(image_info);
  exception=DestroyExceptionInfo(exception);
  MagickCoreTerminus();
  return(status);
#else
  MessageToDialog("不再支持JPG文件拼接成Tif文件的功能。");
  return 0;
#endif
}


/**********************************************
 *
 *    HGT  管理高程文件
 *
 ***********************************************/

int maxhgt, minhgt;

/* 返回一个节点的高程值 */
short height(TGhgt *hgt, int x, int y ) 
{ 
    return hgt->data[x][y]; 
}

/* 检测一个高程网格区是否有非0的高程值 */
int has_non_zero_elev (TGhgt *hgt, int start_x, int span_x, int start_y, int span_y)
{
    int row, col;
    	
    for ( row = start_y; row < start_y + span_y; row++ ) {
        for ( col = start_x; col < start_x + span_x; col++ ) {
            if ( height(hgt, col,row) != 0 )
                return 1;
        }
    }
    return 0;
}

/* 检测一个文件名是否是zip文件 */
int isZip(char *fn)
{
    char ext[4];
    int len = strlen(fn);
    
    memset(ext, 0, 4);
    ext[0]=fn[len-3]; ext[1]=fn[len-2]; ext[2]=fn[len-1];
    if(!strcmp(ext, "zip"))
    	return 1;
    else 
    	return 0;
}

/* 创建一个新的高程文件。里面高程数据全为0。 这个函数没有被调用。*/
int newhgt(TGhgt *hgt, char *file_name)
{
    int i, j, k;
    char name[32];
    char yorn[4], xorn[4]; 
    int size;
    short int *var;
    int row, col;

    hgt->hgt_resolution = 3;
    
    hgt->fd = 0;

    // Determine originx/originy from file name
    j = 0; memset(name, 0, sizeof(name));
    for(i=(strlen(file_name)-1); i>=0; i--)
    {
    	j++;
    	if(file_name[i]=='\\') break;
    }
    for(k=0;k<j;k++) name[k] = file_name[i+k+1];
    
    yorn[0]=name[1]; yorn[1]=name[2]; yorn[2]=yorn[3]=0;
    xorn[0]=name[4]; xorn[1]=name[5]; xorn[2]=name[6]; xorn[3]=0;
    
    hgt->originy = atof(yorn) * 3600.0;
    if ( name[0] == 'S' ) {
        hgt->originy = - hgt->originy;
    }
    hgt->originx = atof(xorn) * 3600.0;
    if ( name[3] == 'W' ) {
        hgt->originx = - hgt->originx;
    }

    /*load data from file */
    if ( hgt->hgt_resolution == 1 ) {
        hgt->cols = hgt->rows = size = 3601;
        hgt->col_step = hgt->row_step = 1;
    } else if ( hgt->hgt_resolution == 3 ) {
        hgt->cols = hgt->rows = size = 1201;
        hgt->col_step = hgt->row_step = 3;
    }
    
    for ( row = size - 1; row >= 0; --row ) {
        for ( col = 0; col < size; ++col ) {
            var = &hgt->data[col][row];
			*var = 0;
		}
    }
    
    return 1;
}

/* 打开一个高程文件，并且把中间的空白点抹平，抹平数据是根据已经读过的上一行以及左边一列来确定。  */
int openhgt(TGhgt *hgt, char *file_name)
{
    char command[_MAX_PATH];
    char fn_nozip[_MAX_PATH];
    int len, i, j, k;
    char name[32];
    char yorn[4], xorn[4]; 
    int size;
    short int *var;
    int row, col;
    unsigned short tmp, tmp2;
    short int cu_arr;

    hgt->hgt_resolution = 3;
    
    /* 目前不支持压缩文件解压，所以系统高程文件都是解压后的数据。*/
    /*open hgt file,  if is zip file , unzip first. */    
    if ( isZip(file_name) ) {
        // extract the .zip file to /tmp and point the file name
        // to the extracted file
    	memset(command, 0, sizeof(command));
    	sprintf(command, "UnRAR e %s", file_name);
    	system( command );
        
    	strcpy(fn_nozip, file_name);
    	len = strlen(fn_nozip);
    	fn_nozip[len-4]=0; fn_nozip[len-3]=0; fn_nozip[len-2]=0; fn_nozip[len-1]=0;
    }else
    	strcpy(fn_nozip, file_name);

	/* 打开高程文件。 */    
    if ( (hgt->fd = fopen( fn_nozip, "rb" )) == NULL ) 
    {
        return 0;
    }

    // Determine originx/originy from file name
    j = 0; memset(name, 0, sizeof(name));
    for(i=(strlen(file_name)-1); i>=0; i--)
    {
    	j++;
    	if(file_name[i]=='\\') break;
    }
    for(k=0;k<j;k++) name[k] = file_name[i+k+1];
    
    yorn[0]=name[1]; yorn[1]=name[2]; yorn[2]=yorn[3]=0;
    xorn[0]=name[4]; xorn[1]=name[5]; xorn[2]=name[6]; xorn[3]=0;
    
    hgt->curLong = atoi(xorn); hgt->curLati = atoi(yorn);
    
    hgt->originy = atof(yorn) * 3600.0;
    if ( name[0] == 'S' ) {
        hgt->originy = - hgt->originy;
    }
    hgt->originx = atof(xorn) * 3600.0;
    if ( name[3] == 'W' ) {
        hgt->originx = - hgt->originx;
    }

	/* 目前只支持低分辨率高程数据，即hgt_resolution == 3，修改MAX_HGT_SIZE常数后，可支持高分辨率高程数据，即hgt_resolution == 1 */
    /*load data from file */
    if ( hgt->hgt_resolution == 1 ) {
        hgt->cols = hgt->rows = size = 3601;
        hgt->col_step = hgt->row_step = 1;
    } else if ( hgt->hgt_resolution == 3 ) {
        hgt->cols = hgt->rows = size = 1201;
        hgt->col_step = hgt->row_step = 3;
    }

	/* 逐点读取高程数据，判断是否是空白点，如果是空白点则抹平数据。*/    
    for ( row = size - 1; row >= 0; --row ) {
        for ( col = 0; col < size; ++col ) {
            var = &hgt->data[col][row];
            if ( fread ( var, sizeof(short), 1, hgt->fd) != 1 ) {
                return 0;
            }

			tmp = (*var & 0xff) << 8;
			tmp2 = (unsigned )(*var);
			tmp2 >>= 8;
			*var = tmp + tmp2;

			//去掉空白点
			cu_arr = *var;
			if((cu_arr == -32768)&&(row<(hgt->rows - 1)))
				cu_arr = hgt->data[col][row+1];

			if((cu_arr == -32768)&&(row==(hgt->rows - 1)))
				cu_arr = hgt->data[col-1][row];

			if((cu_arr > -32768)&&(cu_arr < 0))
				cu_arr = 0;
			*var = cu_arr;

			//找出最大最小值
			if(*var > maxhgt) maxhgt = *var;
			if(*var < minhgt) minhgt = *var;
		}
    }
    
    fclose(hgt->fd);
    
    return 1;
}

/**********************************************
 *
 *   根据经纬度范围获得高程数据文件序列，合并高程文件，并转为GeoTiff文件。
 *
 **********************************************/
typedef struct _HGTNODE_
{
	short usLong;
	short usLati;
	char fileName[STRING_STANDARD_LEN];
	char pathName[STRING_STANDARD_LEN];
	TGhgt *pHgt;
}HgtNode;

HgtNode hgtarray[7][7];
int hgtrow, hgtcol;

/* 获取指定经纬度区块的所有高程文件名的列表 */
bool GetHgtFileList(double leftlongi, double downlati, double rightlongi, double upperlati)
{
	int iMaxLong, iMaxLati, iMinLong, iMinLati;
	short iLong, iLati, m, n;
	
	// 初始化各种变量
	memset(hgtarray, 0, sizeof(hgtarray));
	maxhgt = -999; minhgt = 999999;
	iMaxLong = (int)rightlongi;
	iMaxLati = (int)upperlati;
	iMinLong = (int)leftlongi;
	iMinLati = (int)downlati;
	
	//做while循环，填充变量
	iLati = iMinLati; m = 0;
	while(iLati < iMaxLati)
	{
		iLong = iMinLong; n = 0;
		while(iLong  < iMaxLong)
		{
			hgtarray[m][n].usLong = iLong;
			hgtarray[m][n].usLati = iLati;
			sprintf(hgtarray[m][n].fileName, "N%02dE%03d.hgt", iLati, iLong);
			sprintf(hgtarray[m][n].pathName, "%s\\n%02d", hgtroot, iLati);
			iLong++; n++;
		}
		iLati++; m++; 
	}
	hgtrow = m; hgtcol = n;

	return true;
}

/* 初始化高程结构信息 */
bool InitHgtInfo(TIFF* m_pOutputFile, unsigned int m_unImgWidth,unsigned int m_unImgLength)
{
	if(m_pOutputFile)
	{
		// We need to set some values for basic tags before we can add any data
		TIFFSetField(m_pOutputFile, TIFFTAG_IMAGEWIDTH, m_unImgWidth);
		TIFFSetField(m_pOutputFile, TIFFTAG_IMAGELENGTH, m_unImgLength);
		TIFFSetField(m_pOutputFile, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(m_pOutputFile, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
		TIFFSetField(m_pOutputFile, TIFFTAG_BITSPERSAMPLE, 32);
		TIFFSetField(m_pOutputFile, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(m_pOutputFile, TIFFTAG_PHOTOMETRIC, 1);
		TIFFSetField(m_pOutputFile, TIFFTAG_ROWSPERSTRIP, 1);
		TIFFSetField(m_pOutputFile, TIFFTAG_ORIENTATION, 1);
		TIFFSetField(m_pOutputFile, TIFFTAG_SAMPLEFORMAT, 3);
		TIFFSetField(m_pOutputFile, TIFFTAG_SMINSAMPLEVALUE, (double)minhgt);
		TIFFSetField(m_pOutputFile, TIFFTAG_SMAXSAMPLEVALUE, (double)maxhgt);
	}
	return true;
}

/* 按照计算得出的高程网格，合并相关高程文件为一个文件 */
bool MergeHgt()
{
	int n, m, i, j;
	TIFF * pOutput;		//输出文件指针
	char fullPath[STRING_STANDARD_LEN];

	//写打开高程文件
	if((pOutput=TIFFOpen(temHgtFile,"w"))==NULL)
	{
		MessageToDialog("打开文件失败.\n");
		return false;
	}

	//只使用了一行hgt
	if(hgtrow == 1)
	{
		//申请单个行的内存
		int nOutputLineSize = (MAX_HGT_SIZE) * hgtcol;
		float *pcOutputLine=new float[nOutputLineSize];
		if(pcOutputLine == NULL)
		{
			MessageToDialog("没有足够的内存.\n");
			return false;
		}
		
		
		//打开高程文件，读入高程数据
		for(i=0;i<hgtcol;i++)
		{
			TGhgt * curHgt = new TGhgt;
			if(curHgt == NULL)
			{
				MessageToDialog("没有足够的内存.\n");
				currPZ->isOK = STATUS_HGTBAD;
				return false;
			}
			hgtarray[0][i].pHgt = curHgt;
			sprintf(fullPath, "%s\\%s", hgtarray[0][i].pathName, hgtarray[0][i].fileName);
			if(!openhgt(curHgt, fullPath))
			{
				//newhgt(curHgt, fullPath);
				currPZ->isOK = STATUS_HGTBAD;
				return false;
			}	
		}

		//写入Geotiff文件参数
		InitHgtInfo(pOutput, hgtcol * hgtarray[0][0].pHgt->cols, hgtarray[0][0].pHgt->rows);
	
		//逐行合并写入。
		for(m=0; m<MAX_HGT_SIZE; m++)
		{
			//首先将几个hgt数据合并成一行
			for(i=0;i<hgtcol;i++)
			{
				for(j=0; j<(MAX_HGT_SIZE); j++)
					pcOutputLine[i*(MAX_HGT_SIZE) + j] = (float)(hgtarray[0][i].pHgt->data[j][MAX_HGT_SIZE - 1 - m]);
			}
			
			//将这一行写入文件
			if(TIFFWriteScanline(pOutput,pcOutputLine,m)!=1)
			{
				MessageToDialog("不能够正确地写入文件!");
				currPZ->isOK = STATUS_HGTBAD;
				return false;
			}
		}
		
		//读取整合后高程相关信息
		mergedHgt.minLong = hgtarray[0][0].usLong;
		mergedHgt.maxLati = hgtarray[0][0].usLati + 1;
		mergedHgt.maxLong = hgtarray[0][hgtcol - 1].usLong + 1;
		mergedHgt.minLati = hgtarray[0][hgtcol - 1].usLati;
		sprintf(mergedHgt.pathName, "%s", tempDirectory);
		sprintf(mergedHgt.fileName, "%s.%s", TEMP_HGTFILE, TEMP_EXTENDNAME);
		mergedHgt.image_width = hgtcol * hgtarray[0][0].pHgt->cols;
		mergedHgt.image_height = hgtarray[0][0].pHgt->rows;
		
		//释放内存
		delete pcOutputLine;
		for(i=0;i<hgtcol;i++)
			if(hgtarray[0][i].pHgt) delete hgtarray[0][i].pHgt;
		
	}else{		//使用了多行

		//申请单个行的内存
		int nOutputLineSize = (MAX_HGT_SIZE) * hgtcol;
		float *pcOutputLine=new float[nOutputLineSize];
		if(pcOutputLine == NULL)
		{
			MessageToDialog("没有足够的内存.\n");
			currPZ->isOK = STATUS_HGTBAD;
			return false;
		}
				
		//打开高程文件，读入高程数据
		for(j=0;j<hgtrow;j++)
			for(i=0;i<hgtcol;i++)
			{
				TGhgt * curHgt = new TGhgt;
				if(curHgt == NULL)
				{
					MessageToDialog("没有足够的内存.\n");
					currPZ->isOK = STATUS_HGTBAD;
					return false;
				}
				hgtarray[j][i].pHgt = curHgt;
				sprintf(fullPath, "%s\\%s", hgtarray[j][i].pathName, hgtarray[j][i].fileName);
				if(!openhgt(curHgt, fullPath))
				{
					//newhgt(curHgt, fullPath);
					currPZ->isOK = STATUS_HGTBAD;
					return false;
				}	
			}
	
		//写入Geotiff文件参数
		InitHgtInfo(pOutput, hgtcol * hgtarray[0][0].pHgt->cols, hgtrow * hgtarray[0][0].pHgt->rows);

		//逐行合并写入。
		for(n=hgtrow - 1; n>=0; n--)
		{
			for(m=0; m<MAX_HGT_SIZE; m++)
			{
				//首先将几个hgt数据合并成一行
				for(i=0;i<hgtcol;i++)
				{
					for(j=0; j<(MAX_HGT_SIZE); j++)
						pcOutputLine[i*(MAX_HGT_SIZE) + j] = (float)(hgtarray[n][i].pHgt->data[j][MAX_HGT_SIZE - 1 - m]);
				}
				
				//将这一行写入文件
				if(TIFFWriteScanline(pOutput,pcOutputLine,(hgtrow - 1 - n)*MAX_HGT_SIZE + m)!=1)
				{
					MessageToDialog("不能够正确地写入文件!");
					currPZ->isOK = STATUS_HGTBAD;
					return false;
				}
			}
		}

		//读取整合后高程相关信息
		mergedHgt.minLong = hgtarray[0][0].usLong;
		mergedHgt.maxLati = hgtarray[hgtrow - 1][hgtcol - 1].usLati + 1;
		mergedHgt.maxLong = hgtarray[hgtrow - 1][hgtcol - 1].usLong + 1;
		mergedHgt.minLati = hgtarray[0][0].usLati;
		sprintf(mergedHgt.pathName, "%s", tempDirectory);
		sprintf(mergedHgt.fileName, "%s.%s", TEMP_HGTFILE, TEMP_EXTENDNAME);
		mergedHgt.image_width = hgtcol * hgtarray[0][0].pHgt->cols;
		mergedHgt.image_height = hgtrow * hgtarray[0][0].pHgt->rows;
		
		//释放内存
		delete pcOutputLine;
		for(j=0;j<hgtrow;j++)
			for(i=0;i<hgtcol;i++)
				if(hgtarray[j][i].pHgt) delete hgtarray[j][i].pHgt;
		
	}
	
	//关闭文件		
	TIFFClose(pOutput);
	return true;
}


/************************************************
 *
 *    TravelMap
 *    遍历地景纹理目录，搜寻.map文件，提取信息存入数据库
 *    .map文件是与.jpg文件按文件名一一对应的文本文件，里面保存了.jpg文件所代表的区域的地理信息，经纬度。
 *    这是从GoogleEarth上下载的地景原始数据的特点。   
 *
 **********************************************************/

/* 从一个文本文件里面读取一行字符串，并且以空格为分界，存入cf二维数组中，这是个很常用的过程。 */
/* 正常读取一行数据返回1，读取空行返回2，文件结束返回0*/
int ReadALineToArray(FILE *fi)
{
	int j0, j1, d0, i;
	char OneLine[2048];
	char *ptr;

	memset(cf, 0, sizeof(cf));
	if(feof(fi)) return 0;
	ptr = fgets(OneLine, 2048, fi); if((ptr == NULL)) return 2;

	int contiFlag = 0;	
do_start:
	do
	{
		int ii = 0;
		/* 如果以 “//” 开始的行，也当作空行处理，返回 2 --- 增加注释的功能。 */
		if((OneLine[0] == '/') && (OneLine[1] == '/')) contiFlag++;
	
		/* 判断本行是否全部是空格或者Tab */
		int iFlag = 0;
		for(i=0;i<(int)strlen(OneLine);i++)
		{
			if(!((OneLine[i]==' ')||(OneLine[i]==0x0a)||(OneLine[i]==0x09))) {iFlag = 1; break; }
		}
		if(!iFlag) contiFlag++;

		if(contiFlag>0)
		{
			unsigned int i=1;
			ptr = fgets(OneLine, 2048, fi); if((ptr == NULL)) return 2; 
			else {	contiFlag = 0; goto do_start;	}
		}

	}while(contiFlag);

	/* 读取正常数据*/
	j0=j1=d0=0;
	for(i=0;i<(int)strlen(OneLine);i++)
	{
		if(!((OneLine[i]==' ')||(OneLine[i]==0x0a)||(OneLine[i]==0x09)))  cf[j0][j1++] = OneLine[i];
		if(((OneLine[i]==' ')||(OneLine[i]==9))&&((OneLine[i-1]!=' ')&&(OneLine[i-1]!=9))&&(i!=0)) {j0++; j1=0;}
	}
	return 1;
}



/* 遍历当前目录下的所有子目录，以查找 map 文件，使用递归方式。 */
int TravelSubDirectoryForMap(char *dir)
{
	int dio, dfio, dirIo, i, j, k;
	struct _finddata_t sdff;

	//查找当前路径下的子目录；
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//如果不是子目录也不是文件
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR;

		//如果是子目录，进入当前子目录
		if(sdff.attrib == _A_SUBDIR)
		{
			//先分析这个目录名
			_getcwd(path, _MAX_PATH);
			CheckTextureDirName(sdff.name, path);
			
			//进入子目录
			dirIo = chdir(sdff.name);
			_getcwd(path, _MAX_PATH);
			
			//遍历子目录
			TravelSubDirectoryForMap(path);
			
			//退回到上一层目录
			dirIo = chdir("..");
		}
		
		//如果是文件
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//找到map文件
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='M')||(sdff.name[i-3]=='m')) && ((sdff.name[i-2]=='A')||(sdff.name[i-2]=='a')) && ((sdff.name[i-1]=='P')||(sdff.name[i-1]=='p')))  
			{
				MapNode mn;
				//打开map文件
				fi = fopen(sdff.name, "rt");
				if(fi == NULL) 
				{
					char msg[STRING_STANDARD_LEN];
					sprintf(msg, "不能读取文件 %s .", sdff.name);
					MessageToDialog(msg);
					return 0;
				}
				
				ReadALineToArray(fi);
				ReadALineToArray(fi);	j = strlen(cf[0]); strcpy(mn.fileName, cf[0]);
				
				//获取图片的宽、高信息
				while( !(cf[0][0]=='M' && cf[0][1]=='M' && cf[0][2]=='P' && cf[0][3]=='X' && cf[0][4]=='Y' && cf[0][5]==',' && cf[0][6]=='3' && cf[0][7]==','))  ReadALineToArray(fi);
				k = 0;
				int flag = 0;
				char cWidth[8], cHeigh[8];
				memset(cWidth, 0, sizeof(cWidth)); memset(cHeigh, 0, sizeof(cHeigh));
				for(j=8;j<(int)strlen(cf[0]);j++)
				{
					if(cf[0][j] == ',') { k = 0; flag = 1;}
					if((cf[0][j] != ',')&&(!flag)) cWidth[k++] = cf[0][j];
					if((cf[0][j] != ',')&&( flag)) cHeigh[k++] = cf[0][j];
				}
				
				/* 获取图片的长宽信息，由于原始图片缩小了一倍，所以读取的数据要在原数据基础上 除以 2 */
				mn.image_width  = atoi(cWidth) / 2;
				mn.image_height = atoi(cHeigh) / 2;
				
				//获取图片的地理信息
				while((strcmp("MMPLL,1,", cf[0]))&&(!feof(fi)))  ReadALineToArray(fi);
				if(feof(fi))
				{
					MessageToDialog("读取map文件出错:找不到经纬度信息。");
					return 0;
				}
				mn.minLong = atof(cf[1]);
				mn.maxLati = atof(cf[2]);
				ReadALineToArray(fi);
				ReadALineToArray(fi);
				mn.maxLong = atof(cf[1]);
				mn.minLati = atof(cf[2]);
				_getcwd(path, _MAX_PATH);
				strcpy(mn.pathName, path);
				mn.ID = nodeID;
				fclose(fi);
				totalnode.push_back(mn);	nodeID++;
			}
		}

GET_NEXT_SUBDIR:
		//寻找下一个
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* 遍历当前目录，查找该目录下以及所有子目录下扩展名为.map的文件，读取里面相关信息。 */
int TravelDirectoryToReadMap(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForMap(dPath))
		return 0;

	return 1;
} 

/*============================================================================================
 *  从下面开始到 CheckTheSortList 函数之前，是用Tiff文件生成地景方法的主过程区，部分子过程在前面。这是用户需要支持的格式。
 *============================================================================================ */

/**********************************************************
 *
 *       遍历地景纹理目录，搜寻.tif文件，提取信息存入数据库
 *       
 **********************************************************/
/* 在用户高分辨率纹理(highTxtnode)及高程(highHgtnode)列表里面查找当前文件名在不在？*/ 
int CheckTifInTotalNode(char *tifname)
{
	int i, retval = 0;

	if(isCurrTxtOrHgt == STATUS_TEXTURE)
	{
		for(i=0; i<hTxtNum; i++)
		{
			if(!strcmp(tifname, highTxtnode[i].fileName))
			{ retval = 1; break; }	
		}
	}
	
	if(isCurrTxtOrHgt == STATUS_HGT)
	{
		for(i=0; i<hHgtNum; i++)
		{
			if(!strcmp(tifname, highHgtnode[i].fileName))
			{ retval = 1; break; }	
		}
	}

	return retval;
}
 
/* 在子目录里面搜索tif文件 */
/* 遍历当前目录下的所有子目录，以查找 tif 文件，使用递归方式。 */
int TravelSubDirectoryForTif(char *dir)
{
	int dio, dfio, dirIo, i;
	struct _finddata_t sdff;

	//查找当前路径下的子目录；
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//如果不是子目录也不是文件
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR3;

		//如果是子目录，进入当前子目录
		if(sdff.attrib == _A_SUBDIR)
		{
			//根据子目录的名称，确定纹理级别
			if(!strcmp("1m", sdff.name))  currHighLevel = 1;
			if(!strcmp("5m", sdff.name))  currHighLevel = 5;
			if(!strcmp("10m", sdff.name)) currHighLevel = 10;
			if(!strcmp("50m", sdff.name)) currHighLevel = 50;

			//进入子目录
			dirIo = chdir(sdff.name);
			_getcwd(path, _MAX_PATH);
			
			//遍历子目录
			TravelSubDirectoryForTif(path);
			
			//退回到上一层目录
			dirIo = chdir("..");
		}
		
		//如果是文件
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//找到Tif文件
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='T')||(sdff.name[i-3]=='t')) && ((sdff.name[i-2]=='I')||(sdff.name[i-2]=='i')) && ((sdff.name[i-1]=='F')||(sdff.name[i-1]=='f')))  
			{
				/* 先查查这个文件是否已经在库里面了。从后往前查*/
				if(!CheckTifInTotalNode(sdff.name))
				{					
					/* 打开GeoTiff文件，读取Tags信息 */
					if(!ReadGeoTiffTags(sdff.name))
					{
						return 0;
					}	
				}
			}
		}

GET_NEXT_SUBDIR3:
		//寻找下一个
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* 遍历当前目录，查找该目录下以及所有子目录下扩展名为.tif的文件，读取里面相关信息。 */
int TravelDirectoryToReadTif(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	/* 递归遍历当前目录及所有的子目录 */
	if(!TravelSubDirectoryForTif(dPath))
	{
		MessageToDialog("找不到带地理信息的纹理文件");
		return 0;
	}
	
	/* 将highTxtnode里面 分辨率为50米的数据 转移到 lowTxtnode 里面*/ 
	if(!MoveHigh50ToLow())	return 0;

	/* 根据读取的用户高分辨率Tif文件信息，来生成高分辨率地景制作列表allHTI。
	   这里综合考虑了所有的覆盖情况，尽量做到高分辨率数据覆盖低分辨率数据。 */
	if(!ReadNodeListToCreateHTI()) 	{	return 0;	}

	/* 根据分辨率为50米的系统给定数据，来生成低分辨率地景制作列表lowHTI，
	   由于低分辨率地景系统已经提供，这个过程无实际意义。*/
	if(!GetHTIListOfR50())	{	return 0;	}
	
	return 1;
} 

/* 将highTxtnode里面 分辨率为50米 的转移到 lowTxtnode 里面*/ 
int MoveHigh50ToLow()
{
	int i, j;
	
	if(lowTxtnode.size() == 0)
	{
		lTxtNum = 0;
		/* 先将highTxtnode里面分辨率为50米 的转移到 lowTxtnode 里面*/
		for(i=0; i<hTxtNum; i++)
		{
			if(highTxtnode[i].highLevel == 50)
			{
				HighNode hn;
				strcpy(hn.fileName, 	  highTxtnode[i].fileName);
				hn.highLevel  			= highTxtnode[i].highLevel;
				hn.image_width  		= highTxtnode[i].image_width;
				hn.image_height 		= highTxtnode[i].image_height;
				hn.minLong 				= highTxtnode[i].minLong;
				hn.maxLati 				= highTxtnode[i].maxLati;
				hn.maxLong 				= highTxtnode[i].maxLong;
				hn.minLati 				= highTxtnode[i].minLati;
				strcpy(hn.pathName, 	  highTxtnode[i].pathName);
				hn.ID = lTxtNum++;
				lowTxtnode.push_back(hn);
 			}
 		}
	}
 	
 	/* 然后移除highTxtnode里面分辨率为50米的值 */
	for(i=0; i<hTxtNum; i++)
	{
		if(highTxtnode[i].highLevel == 50)  highTxtnode[i].highLevel = 0;
 	}
 	
 	i = j = 0;
 	while(i<hTxtNum)
 	{
 		if(highTxtnode[i].highLevel != 0) 
 		{
 			if(j != i)
 			{
				strcpy(highTxtnode[j].fileName, 		  highTxtnode[i].fileName);
				highTxtnode[j].highLevel  				= highTxtnode[i].highLevel;
				highTxtnode[j].image_width  			= highTxtnode[i].image_width;
				highTxtnode[j].image_height 			= highTxtnode[i].image_height;
				highTxtnode[j].minLong 					= highTxtnode[i].minLong;
				highTxtnode[j].maxLati 					= highTxtnode[i].maxLati;
				highTxtnode[j].maxLong 					= highTxtnode[i].maxLong;
				highTxtnode[j].minLati 					= highTxtnode[i].minLati;
				strcpy(highTxtnode[j].pathName, 		  highTxtnode[i].pathName);
				highTxtnode[j].ID = j;
 			}
 			j++;
 		}
 		i++;
 	}  
 	hTxtNum = j;
 	return 1;
}

/*============================================================================================
 *  从下面开始到 AppendStgAndWrite 函数之前，用若干个Jpeg文件拼接合并成一个Tiff文件的过程。如何组装命令行。
 *============================================================================================ */

/**********************************************************
 *
 *  几个区域较小的Jpg图片拼接为一个区域较大的Tif图片，Tif图片至少包含一个经度x一个纬度范围内的地景。
 *       
 **********************************************************/

/* 检查排序结果是否正确 */
int CheckTheSortList(int *curLine)
{
	int i, iFlag = -1;
	double curlineLon, prelineLon;
	
	/* 首先检查几行的数目，是否相同 */
	for(i=1; i<sortrows; i++)
	{
		if(sortcolums[i] != sortcolums[i-1])
		{
			if(sortcolums[i] < sortcolums[i-1]) {	iFlag = i;   break; }
			if(sortcolums[i] > sortcolums[i-1]) {	iFlag = i-1; break; }
		}
	}
	*curLine = iFlag;
	
	/* 如果出现不相同，检查最左边经度是否接近 */
	if(iFlag!=-1)
	{
		if(iFlag > 0)
		{	
			prelineLon = findResult[sortlist[i-1][0]].minLong;
			curlineLon = findResult[sortlist[i  ][0]].minLong;
		}else{
			prelineLon = findResult[sortlist[i+1][0]].minLong;
			curlineLon = findResult[sortlist[i  ][0]].minLong;
		}
		if(fabs(prelineLon - curlineLon) > 0.01)
			iFlag = 0;
		else	
			iFlag = 1;
	}
	
	return iFlag;
}

/* 查找指定区域内包含多少个jpg*/
#define LIMIT_WIDTH 0.01
int FindJpgsInZoneAndSort(double minlon, double maxlon, double minlat, double maxlat)
{
	int i, j, badLine;
	double dtmp1, dtmp2;
	
	/* 从图片集合里面查找符合条件的图片 */
	resultID = 0;	findResult.clear();
	for(i=0; i<nodeID; i++)
	{
		if(	(!((totalnode[i].minLong >= maxlon) || (totalnode[i].maxLong <= minlon))) &&
			(!((totalnode[i].minLati >= maxlat) || (totalnode[i].maxLati <= minlat))) &&
		    (!strcmp(findJpgDirectory, totalnode[i].pathName)) )
		{
			MapNode mn;
			mn.minLong = totalnode[i].minLong;
			mn.maxLong = totalnode[i].maxLong;
			mn.minLati = totalnode[i].minLati;
			mn.maxLati = totalnode[i].maxLati;
			strcpy(mn.fileName, totalnode[i].fileName);
			strcpy(mn.pathName, totalnode[i].pathName);
			mn.image_width  = totalnode[i].image_width;
			mn.image_height = totalnode[i].image_height;
			mn.ID = totalnode[i].ID;
			findResult.push_back(mn);	resultID++;
		}		
	}

	/* 如果找到奇数个图片，检查并清除一个 */
	if((resultID % 2)&&(resultID != 9))
	{	
		for(i=0; i<resultID; i++)
		{
			dtmp1 = fabs(findResult[i].minLong - maxlon);
			dtmp2 = fabs(findResult[i].maxLong - minlon);
			if((dtmp1 < LIMIT_WIDTH) || (dtmp2 < LIMIT_WIDTH))
			{
				for(j=i; j<resultID;j++)
				{
					findResult[j].minLong = findResult[j+1].minLong;
					findResult[j].maxLong = findResult[j+1].maxLong;
					findResult[j].minLati = findResult[j+1].minLati;
					findResult[j].maxLati = findResult[j+1].maxLati;
					strcpy(findResult[j].fileName, findResult[j+1].fileName);
					strcpy(findResult[j].pathName, findResult[j+1].pathName);
					findResult[j].image_width  = findResult[j+1].image_width;
					findResult[j].image_height = findResult[j+1].image_height;
					findResult[j].ID = findResult[j+1].ID;
				}
				resultID -= 1;
				break;
			}
		}
	}

	/* 重新排序 */
	if(!SortfindResult())
	{
		MessageToDialog("排序图片失败\n");
		currPZ->isOK = STATUS_TEXTBAD;
		return 0;
	}

	/* 检查是否有错位，即几个行图片数目不等。*/
	int Flag;
	if((Flag = CheckTheSortList(&badLine)) != -1)
	{
		/* 找到需要增加的图片的序号 */
		int needID;
		if(!Flag)
		{
			needID = findResult[sortlist[badLine][0]].ID - 1;
		}
		if(Flag)
		{
			needID = findResult[sortlist[badLine][sortcolums[badLine] - 1]].ID + 1;
		}
		
		/* 增加新图片*/
		MapNode mn;
		mn.minLong = totalnode[needID].minLong;
		mn.maxLong = totalnode[needID].maxLong;
		mn.minLati = totalnode[needID].minLati;
		mn.maxLati = totalnode[needID].maxLati;
		strcpy(mn.fileName, totalnode[needID].fileName);
		strcpy(mn.pathName, totalnode[needID].pathName);
		mn.image_width  = totalnode[needID].image_width;
		mn.image_height = totalnode[needID].image_height;
		mn.ID = totalnode[needID].ID;
		findResult.push_back(mn);	resultID++;
		
		/* 重新排序 */
		if(!SortfindResult())
		{
			MessageToDialog("排序图片失败\n");
			currPZ->isOK = STATUS_TEXTBAD;
			return 0;
		}	
	}


	/* 检查边缘，是否到边 */
	/* 先检查最小纬度 */
	for(i=resultID-1; i>=0; i--)
	{
		double delta = fabs(findResult[i].minLati - currTi.minLati);
		if((delta < 0.1) && (findResult[i].minLati > currTi.minLati))
			findResult[i].minLati = currTi.minLati;
		else
			break;
	}

	/* 检查最大纬度 */
	for(i=0; i<resultID; i++)
	{
		double delta = fabs(findResult[i].maxLati - currTi.maxLati);
		if((delta < 0.1) && (findResult[i].maxLati < currTi.maxLati))
			findResult[i].maxLati = currTi.maxLati;
		else
			break;
	}


	return 1;
}


/* 查找一个bucket里面包含了多少个jpg. */
int FindJpgsFromBucket(Bucket *b)
{
	int i;
	double minbt_lon, minbt_lat, maxbt_lon, maxbt_lat;
	
	//get corner
	minbt_lon = get_center_lon(b) - 0.125;
	maxbt_lon = get_center_lon(b) + 0.125;
	minbt_lat = get_center_lat(b) - 0.125/2.0;
	maxbt_lat = get_center_lat(b) + 0.125/2.0;
	
	//find from totalnode.
	resultID = 0;	findResult.clear();
	for(i=0; i<nodeID; i++)
	{
		if(	(!((totalnode[i].minLong >= maxbt_lon) || (totalnode[i].maxLong <= minbt_lon))) &&
			(!((totalnode[i].minLati >= maxbt_lat) || (totalnode[i].maxLati <= minbt_lat))) &&
		    (!strcmp(findJpgDirectory, totalnode[i].pathName)) )
		{
			MapNode mn;
			mn.minLong = totalnode[i].minLong;
			mn.maxLong = totalnode[i].maxLong;
			mn.minLati = totalnode[i].minLati;
			mn.maxLati = totalnode[i].maxLati;
			strcpy(mn.fileName, totalnode[i].fileName);
			strcpy(mn.pathName, totalnode[i].pathName);
			mn.image_width  = totalnode[i].image_width;
			mn.image_height = totalnode[i].image_height;
			mn.ID = totalnode[i].ID;
			findResult.push_back(mn);	resultID++;
		}		
	}
	
	return 1;	
}


//将Result数组里面进行二维排序，写入二维数组里面，首先纬度从高到低，其次同一纬度从低到高的顺序。
int SortfindResult()
{
	int i, j, k, m, n, ex;
	int srclist[MAXFINDNODES];
	double maxlat = -99.9;
	
	//初始化数组
	memset(sortlist, -1, sizeof(sortlist)); 
	memset(sortcolums, -1, sizeof(sortcolums));
	memset(srclist, -1, sizeof(srclist));

	//给srclist赋初始值，
	for(i=0; i<resultID; i++)  srclist[i] = i;
	
	//先按纬度从高到低排序
	m = 0;
	while(1)
	{
		//遍历结果数组，寻找最大纬度
		maxlat = -99.9;
		for(i=0; i<resultID; i++)  
		{
			if((findResult[i].minLati > maxlat)&&(srclist[i] != -1)) 
				maxlat = findResult[i].minLati;
		}
		if( maxlat == -99.9) break;

		//将和maxlat相等的值的idx存在当前下标数组中。
		n = 0;
		for(i=0; i<resultID; i++) 
		{
			if((findResult[i].minLati == maxlat)&&(srclist[i] != -1))
			{
				sortlist[m][n++] = srclist[i];
				srclist[i] = -1;
			}
		}
		sortcolums[m] = n;
		m++;
	}
	sortrows = m;

	//同纬度按经度从低到高排序
	for(i=0;i<m;i++)
	{
		for(j=0; j<(sortcolums[i] - 1); j++)
			for(k = j+1; k<sortcolums[i]; k++)
			{
				if(findResult[sortlist[i][j]].minLong > findResult[sortlist[i][k]].minLong )
				{
					ex = sortlist[i][j]; sortlist[i][j] = sortlist[i][k]; sortlist[i][k] = ex;
				}
			}
	}
	
	return 1;
}

//将相关图像合并在一起
int MergePictures()
{
	char inJpgs[320];
	char ouJpgs[10][32];
	int i, j, n;
	unsigned int width, height;
	char arg[10][STRING_STANDARD_LEN];
	char *argptr[10];
	char msg[STRING_STANDARD_LEN];
	
	memset(inJpgs, 0, sizeof(inJpgs));
	memset(ouJpgs, 0, sizeof(ouJpgs));
	memset(&mergedTif, 0, sizeof(mergedTif));
	width = height = 0;

	if(sortrows >= 1)
	{
		/* 组装输入图片行 */
		for(n=0;n<sortrows;n++)
		{
			sprintf(msg, "开始第 %d -- %d 次横向拼接...", n+1, sortrows);
			MessageToDialog(msg);

			for(i=0;i<sortcolums[n];i++)
			{
				if(n==0)	width += findResult[sortlist[n][i]].image_width;
			}
			height = findResult[sortlist[n][0]].image_height;

			//读取整合后图片相关信息
			mergedTif.minLong = findResult[sortlist[n][0]].minLong;
			mergedTif.maxLati = findResult[sortlist[n][0]].maxLati;
			mergedTif.maxLong = findResult[sortlist[n][i-1]].maxLong;
			mergedTif.minLati = findResult[sortlist[n][i-1]].minLati;
			sprintf(mergedTif.pathName, "%s", tempDirectory);
			sprintf(mergedTif.fileName, "%s%d.jpg", TEMP_TEXTFILE, n);
			mergedTif.image_width = width;
			mergedTif.image_height = height;
	
			/* 组装命令行
			/* convert.exe固定值，后面的参数可变。+append横向拼接	-append纵向拼接 最后参数为输出图片
			/* 
			/*/ 
			memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr));
			sprintf(arg[0], "convert.exe"); argptr[0] = (char *)&arg[0];
			int argIdx = 1;
			for(j=0; j<sortcolums[n]; j++)
			{
				sprintf(arg[argIdx], "%s\\%s", findResult[sortlist[n][j]].pathName, findResult[sortlist[n][j]].fileName);  
				argptr[argIdx] = (char *)&arg[argIdx]; argIdx++;
			}
			if(i>1) {sprintf(arg[argIdx], "+append"); argptr[argIdx] = (char *)&arg[argIdx]; argIdx++; }
			sprintf(arg[argIdx], "%s\\%s", mergedTif.pathName, mergedTif.fileName); argptr[argIdx] = (char *)&arg[argIdx]; argIdx++;
	
#ifdef NOMERGE
			if(0)
#else
			if(ConvertMain(argIdx, argptr))			//(0)	//
#endif			
			{
				MessageToDialog("拼接文件失败\n");
				currPZ->isOK = STATUS_TEXTBAD;
				return 0;
			}
		}

		/* 计算总的图片高度 */
		sprintf(msg, "开始纵向拼接...");
		MessageToDialog(msg);
		height = 0;
		for(j=0; j<sortrows; j++)
			height += findResult[sortlist[j][0]].image_height;

		/* 读取整合后图片相关信息 */
		mergedTif.minLong = findResult[sortlist[0][0]].minLong;
		mergedTif.maxLati = findResult[sortlist[0][0]].maxLati;
		mergedTif.maxLong = findResult[sortlist[sortrows-1][sortcolums[sortrows-1]-1]].maxLong;
		mergedTif.minLati = findResult[sortlist[sortrows-1][sortcolums[sortrows-1]-1]].minLati;
		sprintf(mergedTif.pathName, "%s", tempDirectory);
		sprintf(mergedTif.fileName, "%s.%s", TEMP_TEXTFILE, TEMP_EXTENDNAME);
		mergedTif.image_width = width;
		mergedTif.image_height = height;
		
		/* 调整一下图片相关信息，以弥补黑缝 */
		if(mergedTif.minLong > currTi.minLongi)
		{
				mergedTif.minLong = currTi.minLongi;
		}
		if(mergedTif.maxLong < currTi.maxLongi)
		{	
				mergedTif.maxLong = currTi.maxLongi;
		}

		/* 组装整个图片*/
		memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr));
		sprintf(arg[0], "convert.exe"); argptr[0] = (char *)&arg[0];
		int argIdx = 1;
		for(j=0; j<sortrows; j++)
		{
			sprintf(arg[argIdx], "%s\\%s%d.jpg", tempDirectory, TEMP_TEXTFILE, j);
			argptr[argIdx] = (char *)&arg[argIdx]; argIdx++;
		}
		if(i>1) {sprintf(arg[argIdx], "-append"); argptr[argIdx] = (char *)&arg[argIdx]; argIdx++; }
		sprintf(arg[argIdx], "%s\\%s", mergedTif.pathName, mergedTif.fileName); argptr[argIdx] = (char *)&arg[argIdx]; argIdx++;

#ifdef NOMERGE
		if(0)
#else
		if(ConvertMain(argIdx, argptr))			//(0)	//
#endif			
		{
			MessageToDialog("拼接文件失败\n");
			currPZ->isOK = STATUS_TEXTBAD;
			return 0;
		}	
	}else{
		//sortrows = 0的情况，没有查到任何图片。
		MessageToDialog("没有找到可用图片。\n");
		currPZ->isOK = STATUS_TEXTBAD;
		return 0;
	}
	
	sprintf(msg, "图像拼接完成...");
	return 1;
}

/*============================================================================================
 *  从下面开始到 DoVpb 函数之前，管理stg文件的过程。stg文件是视景系统里面的索引文件。
 *  视景系统在读取三维地景数据时，先读取stg文件，根据里面的索引信息来读取相关的三维地景文件。
 *============================================================================================ */

/**********************************************************
 *
 *  往stg尾部追加命令行。将新生成的高低分辨率地景文件名按特定格式追加到stg文件后面。
 *  重复的就不用再追加。
 *       
 **********************************************************/

/* 追加ctileID.stg文件，如果没有就新建，并根据main写入语句 */
int AppendStgAndWrite(char *ctileID, char *main, int type)
{
	char dstFile[STRING_STANDARD_LEN];
	char currWr[32];
	int  slen;
	int  flag = 0;

	memset(currWr, 0, sizeof(currWr));
	sprintf(dstFile, "%s.stg", ctileID);  
	fi = fopen(dstFile, "at+");
	if(fi == NULL)
	{
		fi = fopen(dstFile, "wt");
		if(fi == NULL)
		{
			MessageToDialog("不能创建地景stg文件。\n");
			return 0;
		}
	}

	slen = strlen(main);
	if(type == 0)
	{	
		if(slen == 0)
		{
			//fprintf(fi, "OBJECT %s.btg\n", ctileID);
			sprintf(currWr, "%s.ive", ctileID);
		}else{
			sprintf(currWr, "%s.ive", main);
		}
	}
	
	if(type == 1)
	{	
		if(slen == 0)
		{
			//fprintf(fi, "OBJECT %s.btg\n", ctileID);
			sprintf(currWr, "%s.btg.gz", ctileID);
		}else{
			sprintf(currWr, "%s.btg.gz", main);
		}
	}
	
	/* 查找一下文件内容，是否已经写过同样的内容。*/
	while(!feof(fi))
	{
		ReadALineToArray(fi);
		if(!strcmp(cf[1], currWr))
			flag = 1;
	}
	
	if(!flag) fprintf(fi, "OBJECT %s\n", currWr);
	
	fclose(fi);
	return 1;
}

/* 制作tileid.stg文件*/
int MakeSTGFile(Bucket *Bt, char *main, int type)
{
	unsigned int i;
	int dio, j, k, flag00, index, olddrv, drv;
	char dstPath[STRING_STANDARD_LEN], firPath[STRING_STANDARD_LEN], secPath[STRING_STANDARD_LEN];
	char ctileID[10];
	char path[_MAX_PATH];
	int retval = 1;

	memset(ctileID, 0, sizeof(ctileID));
	memset(dstPath, 0, sizeof(dstPath));
	memset(firPath, 0, sizeof(firPath));
	memset(secPath, 0, sizeof(secPath));

	//设置Bucket, 建立目录结构
	gen_base_path(Bt, dstPath);
	index = gen_index(Bt);
	sprintf(ctileID, "%d", index);

	//获得当前工作目录
	olddrv = _getdrive();
	_getcwd(path, _MAX_PATH);

	//先进入到输出根目录
	drv = iveroot[0] - 'A' + 1;
	dio = _chdrive(drv);
	if(dio)
	{
		MessageToDialog("不能进入地景数据库目录所在驱动器中。\n");
		retval = 0; goto MSFOut;
	}

	dio = chdir(iveroot);
	if(dio)
	{
		MessageToDialog("不能进入地景数据库目录中。\n");
		retval = 0; goto MSFOut;
	}

	//查找并创建干目录
	j = k = flag00 = 0;
	for(i=0;i<strlen(dstPath);i++)
	{
		if(flag00) secPath[k++] = dstPath[i];
		if(dstPath[i]=='\\') flag00 = 1;
		if(!flag00) firPath[j++] = dstPath[i];
	}	
	dio = chdir(firPath);
	if(dio) 
	{
		dio = mkdir(firPath);
		if(dio)
		{
			MessageToDialog("不能创建地景子目录。\n");
			retval = 0; goto MSFOut;
		}
		dio = chdir(firPath);
	}
	dio = chdir(secPath);
	if(dio) 
	{
		dio = mkdir(secPath);
		if(dio)
		{
			MessageToDialog("不能创建地景子目录。\n");
			retval = 0; goto MSFOut;
		}
		dio = chdir(secPath);
	}
	
	//创建stg文件并写入
	retval = AppendStgAndWrite(ctileID, main, type);
	
	//退回根目录
	dio = chdir("..");
	dio = chdir("..");

MSFOut:
	//退回到原来的工作目录
	dio = _chdrive(olddrv);
	if(dio)
	{
		MessageToDialog("不能进入工作目录所在驱动器中。\n");
		return 0;
	}

	dio = chdir(path);
	if(dio)
	{
		MessageToDialog("不能进入工作目录中。\n");
		return 0;
	}
	
	return retval;
}

/* 检查STG文件是否存在。*/
int CheckSTGFile(Bucket *Bt)
{
	unsigned int i;
	int dio, j, k, flag00, index, olddrv, drv;
	char dstPath[STRING_STANDARD_LEN], firPath[STRING_STANDARD_LEN], secPath[STRING_STANDARD_LEN], dstFile[STRING_STANDARD_LEN];
	char ctileID[10], iveStr[20];
	char path[_MAX_PATH];
	int retval = 1;

	memset(ctileID, 0, sizeof(ctileID));
	memset(dstPath, 0, sizeof(dstPath));
	memset(firPath, 0, sizeof(firPath));
	memset(secPath, 0, sizeof(secPath));
	memset(dstFile, 0, sizeof(dstFile));

	//设置Bucket, 建立目录结构
	gen_base_path(Bt, dstPath);
	index = gen_index(Bt);
	sprintf(ctileID, "%d", index);

	//获得当前工作目录
	olddrv = _getdrive();
	_getcwd(path, _MAX_PATH);

	//先进入到输出根目录
	drv = iveroot[0] - 'A' + 1;
	dio = _chdrive(drv);
	if(dio)
	{
		MessageToDialog("不能进入地景数据库目录所在驱动器中。\n");
		retval = 0; goto CSFOut;
	}

	dio = chdir(iveroot);
	if(dio)
	{
		MessageToDialog("不能进入地景数据库目录中。\n");
		retval = 0; goto CSFOut;
	}

	//查找并创建干目录
	j = k = flag00 = 0;
	for(i=0;i<strlen(dstPath);i++)
	{
		if(flag00) secPath[k++] = dstPath[i];
		if(dstPath[i]=='\\') flag00 = 1;
		if(!flag00) firPath[j++] = dstPath[i];
	}	
	dio = chdir(firPath);
	if(dio) 
	{
		dio = mkdir(firPath);
		if(dio)
		{
			MessageToDialog("不能创建地景子目录。\n");
			retval = 0; goto CSFOut;
		}
		dio = chdir(firPath);
	}
	dio = chdir(secPath);
	if(dio) 
	{
		dio = mkdir(secPath);
		if(dio)
		{
			MessageToDialog("不能创建地景子目录。\n");
			retval = 0; goto CSFOut;
		}
		dio = chdir(secPath);
	}
	
	//判断stg文件是否存在并打开
	sprintf(dstFile, "%s.stg", ctileID);
	fi = fopen(dstFile, "rt");
	if(fi == NULL)
	{
		MessageToDialog("当前块没有stg文件。\n");
		retval = 0; goto CSFOut;
	}
	
	//判断stg文件的内容是否包含 tileID.ive 这个语句。
	sprintf(iveStr, "%s.ive", ctileID);
	int flag = 0;
	while(!feof(fi))
	{
		ReadALineToArray(fi);
		if(!strcmp(iveStr, cf[1]))
		{
			flag = 1; break;
		}		
	}
	fclose(fi);
	
	if(!flag)
	{
		retval = 0; goto CSFOut;
	}
	
	//退回根目录
	dio = chdir("..");
	dio = chdir("..");

CSFOut:
	//退回到原来的工作目录
	dio = _chdrive(olddrv);
	if(dio)
	{
		MessageToDialog("不能进入工作目录所在驱动器中。\n");
		return 0;
	}

	dio = chdir(path);
	if(dio)
	{
		MessageToDialog("不能进入工作目录中。\n");
		return 0;
	}
	
	return retval;
}

/* 判断该块是否在边界 */
int DetectBucket(double lon, double lat)
{
	if(	( fabs(lon - currTi.minLongi) < 0.01 )	||
		( fabs(lat - currTi.minLati)  < 0.01 )	||
		( fabs(lon+0.25 - currTi.maxLongi) < 0.01 )	||
		( fabs(lat+0.125- currTi.maxLati)  < 0.01 )	
	) return 1;
	else return 0;
}                                                 	

/* 批量集中写STG文件 , type==0创建stg并写入ive; type==1追加stg并写入btg*/
int WriteAllStgs(char *itemStr, int type)
{
	double src_lon;
	double src_lat;
	Bucket b;
	char   item[32];
	
	/* 处理一下输入参数itemStr, 去掉包含的扩展名，存在item里面 */
	memset(item, 0, sizeof(item));
	for(int i = 0; i < strlen(itemStr); i++)
	{
		if(itemStr[i] == '.') break;
		item[i] = itemStr[i];
	}
	
	/* 检查stg是否已经写过同名文件，如果没有就写。*/
	src_lon = currTi.minLongi;
	src_lat = currTi.minLati;
	while(src_lon < currTi.maxLongi)
	{
		while(src_lat < currTi.maxLati)
		{			
			/*计算出当前块的Bucket */
			set_bucket(&b, src_lon, src_lat);
		
			if(DetectBucket(src_lon, src_lat))
			{
				if(!MakeSTGFile(&b, item, type))
				{
					MessageToDialog("创建STG文件失败\n");
					return 0;
				}
			}

			/*继续下一个tile*/
			src_lat += 0.125;
		}
		src_lon += 0.25;
		src_lat = currTi.minLati;
	}
	return 1;
}

/*============================================================================================
 *  从下面开始到 CheckConfig 函数之前，是生成vpb命令行，并调用vpb生成地景的过程。
 *============================================================================================ */

/**********************************************************
 *
 *  执行VPB处理，生成ive文件
 *       
 **********************************************************/

int osgdem_main(int argc, char** argv);

/* 执行VPB处理，生成ive文件 */
int DoVpb(Bucket *b, int type)
{
	char arg[30][STRING_STANDARD_LEN];
	char *argptr[30];
	long tile;
	char basePath[STRING_STANDARD_LEN];
	char fullPath[_MAX_PATH];
	float W, S;
	
	/* 计算相关路径信息 */
	if(type==0)
	{
		W = get_center_lon(b) - 0.25/2;
		S = get_center_lat(b) - 0.125/2;
	}
	if(type==1)
	{
		W = currTi.minLongi;
		S = currTi.minLati;
	}
	memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr)); memset(basePath, 0, sizeof(basePath)); memset(fullPath, 0, sizeof(fullPath));
	tile = gen_index(b);
	gen_base_path(b, basePath);
	if(type == 0)	
		sprintf(fullPath, "%s\\%s\\%d.ive", iveroot, basePath, tile);
	if(type == 1)
		sprintf(fullPath, "%s\\%s\\%d%d.ive", iveroot, basePath, (int)W, (int)S);
		
    /* 按照下面格式拼凑命令行 。*/
	// fprintf(fout, "osgdem -t %s -d %s --geocentric --tile-terrain-size 16 --tile-image-size 1024 --radius-to-max-visible-distance-ratio 3 -ge %f %f 0.25 0.125 -l 5 -v 1.0 -o %s\n", texture_name, height_name, W, S, fullPath);
	sprintf(arg[ 0], "osgdem.exe"); 							argptr[ 0] = (char *)&arg[ 0];
	sprintf(arg[ 1], "-t"); 									argptr[ 1] = (char *)&arg[ 1];
	sprintf(arg[ 2], "%s", temTextGeoFile);						argptr[ 2] = (char *)&arg[ 2];
	sprintf(arg[ 3], "-d"); 									argptr[ 3] = (char *)&arg[ 3];
	sprintf(arg[ 4], "%s", temHgtGeoFile); 						argptr[ 4] = (char *)&arg[ 4];
	sprintf(arg[ 5], "--geocentric"); 							argptr[ 5] = (char *)&arg[ 5];
	sprintf(arg[ 6], "--tile-terrain-size"); 					argptr[ 6] = (char *)&arg[ 6];
	sprintf(arg[ 7], "16"); 									argptr[ 7] = (char *)&arg[ 7];
	sprintf(arg[ 8], "--tile-image-size"); 						argptr[ 8] = (char *)&arg[ 8];
	sprintf(arg[ 9], "1024"); 									argptr[ 9] = (char *)&arg[ 9];
	sprintf(arg[10], "--radius-to-max-visible-distance-ratio"); argptr[10] = (char *)&arg[10];
	sprintf(arg[11], "3"); 										argptr[11] = (char *)&arg[11];
	sprintf(arg[12], "--POLYGONAL"); 							argptr[12] = (char *)&arg[12];
	sprintf(arg[13], "--mip_mapping_imagery"); 					argptr[13] = (char *)&arg[13];
	sprintf(arg[14], "--max_anisotropy"); 						argptr[14] = (char *)&arg[14];
	sprintf(arg[15], "8"); 										argptr[15] = (char *)&arg[15];
	sprintf(arg[16], "-ge"); 									argptr[16] = (char *)&arg[16];
	sprintf(arg[17], "%f", W); 									argptr[17] = (char *)&arg[17];
	sprintf(arg[18], "%f", S); 									argptr[18] = (char *)&arg[18];
	if(type==0)                                                                             
	{                                                                                       
		sprintf(arg[19], "0.25"); 								argptr[19] = (char *)&arg[19];
		sprintf(arg[20], "0.125"); 								argptr[20] = (char *)&arg[20];
	}                                                                                       
	if(type==1)                                                                             
	{                                                                                       
		sprintf(arg[19], "1"); 									argptr[19] = (char *)&arg[19];
		sprintf(arg[20], "1"); 									argptr[20] = (char *)&arg[20];
	}                                                                                       
	sprintf(arg[21], "-l");	 									argptr[21] = (char *)&arg[21];
	sprintf(arg[22], "%d", getLodLvl0()	/* currTi.lod */);							argptr[22] = (char *)&arg[22];
	sprintf(arg[23], "-v");	 									argptr[23] = (char *)&arg[23];
	sprintf(arg[24], "%f", currTi.scale); 						argptr[24] = (char *)&arg[24];
	sprintf(arg[25], "-o");	 									argptr[25] = (char *)&arg[25];
	sprintf(arg[26], "%s", fullPath);							argptr[26] = (char *)&arg[26];
	int argIdx = 27;

	char path[_MAX_PATH];
	_getcwd(path, _MAX_PATH);

	/* 调用osgdem主过程。*/
	if(osgdem_main(argIdx, argptr))
	{
		MessageToDialog("VPB运行失败\n");
		return 0;
	}
	
	return 1;
}

/* 专门给制作高精度地景时候用的，用来制作高精度周边的低精度区间 */
/* 输入为数组 doTerrain 的序号，*/
int DoPartVpbAndWriteStg(int idx)
{
	char arg[30][STRING_STANDARD_LEN];
	char *argptr[30];
	char basePath[STRING_STANDARD_LEN];
	char fullPath[_MAX_PATH], simpleName[_MAX_PATH], srcPath[_MAX_PATH];
	char ctileID[12];
	float W, S, offLon, offLat;
	int iW, iS, index;
	Bucket b;
	
	/* 计算相关路径信息 */
	W = doTerrain[idx].minLongi;
	S = doTerrain[idx].minLati;
	offLon = doTerrain[idx].maxLongi - doTerrain[idx].minLongi;
	offLat = doTerrain[idx].maxLati  - doTerrain[idx].minLati ;
	set_bucket(&b, W+0.01, S+0.01);

	/* 拼装输入路径及文件，如果有50米精度的数据，则优先使用50m精度数据，否则使用原有15m精度数据 */
	iW = (int)(W + 0.00001); iS = (int)(S + 0.00001);
	if((index = CheckIfHas50(iW, iS)) != -1)
	{
		sprintf(temTextGeoFile, "%s\\%s", lowTxtnode[index].pathName, lowTxtnode[index].fileName);
	}
	else
		sprintf(temTextGeoFile, "%s\\T%d%d.tif", tempDirectory, (int)W, (int)S);

	sprintf(temHgtGeoFile , "%s\\D%d%d.tif", tempDirectory, (int)W, (int)S);

	/* 拼装输出路径及文件 */
	memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr)); memset(basePath, 0, sizeof(basePath)); memset(fullPath, 0, sizeof(fullPath));
	gen_base_path(&b, basePath);
	sprintf(fullPath, "%s\\%s\\%d%d_%d.ive", iveroot, basePath, (int)W, (int)S, idx);
	sprintf(simpleName, "%d%d_%d.ive", (int)W, (int)S, idx);
	sprintf(srcPath, "%s\\%s", iveroot, basePath);

	/* 先检查一下有没有这个文件，有了就不用再次生成 */
	FILE *fIve = fopen(fullPath, "rb");
	if(fIve)
	{
		fclose(fIve);	return 1;
	}

	//fprintf(fout, "osgdem -t %s -d %s --geocentric --tile-terrain-size 16 --tile-image-size 1024 --radius-to-max-visible-distance-ratio 3 -ge %f %f 0.25 0.125 -l 5 -v 1.0 -o %s\n", texture_name, height_name, W, S, fullPath);
	sprintf(arg[ 0], "osgdem.exe"); 							argptr[ 0] = (char *)&arg[ 0];
	sprintf(arg[ 1], "-t"); 									argptr[ 1] = (char *)&arg[ 1];
	sprintf(arg[ 2], "%s", temTextGeoFile);						argptr[ 2] = (char *)&arg[ 2];
	sprintf(arg[ 3], "-d"); 									argptr[ 3] = (char *)&arg[ 3];
	sprintf(arg[ 4], "%s", temHgtGeoFile); 						argptr[ 4] = (char *)&arg[ 4];
	sprintf(arg[ 5], "--geocentric"); 							argptr[ 5] = (char *)&arg[ 5];
	sprintf(arg[ 6], "--tile-terrain-size"); 					argptr[ 6] = (char *)&arg[ 6];
	sprintf(arg[ 7], "16"); 									argptr[ 7] = (char *)&arg[ 7];
	sprintf(arg[ 8], "--tile-image-size"); 						argptr[ 8] = (char *)&arg[ 8];
	sprintf(arg[ 9], "1024"); 									argptr[ 9] = (char *)&arg[ 9];
	sprintf(arg[10], "--radius-to-max-visible-distance-ratio"); argptr[10] = (char *)&arg[10];
	sprintf(arg[11], "3"); 										argptr[11] = (char *)&arg[11];
	sprintf(arg[12], "--POLYGONAL"); 							argptr[12] = (char *)&arg[12];
	sprintf(arg[13], "--mip_mapping_imagery"); 					argptr[13] = (char *)&arg[13];
	sprintf(arg[14], "--max_anisotropy"); 						argptr[14] = (char *)&arg[14];
	sprintf(arg[15], "8"); 										argptr[15] = (char *)&arg[15];
	sprintf(arg[16], "-ge"); 									argptr[16] = (char *)&arg[16];
	sprintf(arg[17], "%f", W); 									argptr[17] = (char *)&arg[17];
	sprintf(arg[18], "%f", S); 									argptr[18] = (char *)&arg[18];
	sprintf(arg[19], "%f", offLon); 							argptr[19] = (char *)&arg[19];
	sprintf(arg[20], "%f", offLat); 							argptr[20] = (char *)&arg[20];
	sprintf(arg[21], "-l");	 									argptr[21] = (char *)&arg[21];
	sprintf(arg[22], "%d", getLodLvl2()	/* currTi.lod */);							argptr[22] = (char *)&arg[22];
	sprintf(arg[23], "-v");	 									argptr[23] = (char *)&arg[23];
	sprintf(arg[24], "%f", currTi.scale); 						argptr[24] = (char *)&arg[24];
	sprintf(arg[25], "-o");	 									argptr[25] = (char *)&arg[25];
	sprintf(arg[26], "%s", fullPath);							argptr[26] = (char *)&arg[26];
	int argIdx = 27;

	/* 做vpb */
	if(osgdem_main(argIdx, argptr))
	{
		MessageToDialog("VPB运行失败\n");
		return 0;
	}

	/* 集中写STG */
	sprintf(ctileID, "%d%d_%d", (int)W, (int)S, idx);
	if(!WriteAllStgs(ctileID, 0))
	{
		MessageToDialog("创建STG文件失败\n");
		return 0;
	}

	/* 复制一份最新生成的地景文件到 BackupLib 目录下。 */
	char msg[STRING_STANDARD_LEN];
	sprintf(msg, "备份地景文件组 --- %s", simpleName);
	MessageToDialog(msg);
	BackUpFileToTemp(simpleName, srcPath);
	
	return 1;
}

/* 根据不同的纹理精度级别，确定不同的LOD级别。*/
int GetLodLevelWithTextureHighLevel(int hlvl)
{
	int lodlvl = 3;
	switch (hlvl)
	{
	case 1:
		lodlvl = getLodLvl0();
		break;
	case 5:
		lodlvl = getLodLvl1();
		break;
	case 10:
		lodlvl = getLodLvl2();
		break;
	}
	return lodlvl;
}


/* 专门给制作高精度地景时候用的，用来制作高精度区间 */
/* 输入为数组 findResultHTI 的序号，*/
int DoHighVpbAndWriteStg(int idx)
{
	char arg[30][STRING_STANDARD_LEN];
	char *argptr[30];
	char basePath[STRING_STANDARD_LEN];
	char fullPath[_MAX_PATH], simpleName[_MAX_PATH], srcPath[_MAX_PATH];
	char ctileID[12];
	float W, S, offLon, offLat;
	Bucket b;
	
	//计算相关路径信息
	W = findResultHTI[idx].startLongi;
	S = findResultHTI[idx].startLati;
	offLon = findResultHTI[idx].offsetLongi;
	offLat = findResultHTI[idx].offsetLati; 
	set_bucket(&b, W+0.01, S+0.01);

	/* 拼装输入路径及文件 */
	strcpy(temTextGeoFile, findResultHTI[idx].texFname);
	strcpy(temHgtGeoFile , findResultHTI[idx].hgtFname);

	/* 拼装输出路径及文件 */
	memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr)); memset(basePath, 0, sizeof(basePath)); memset(fullPath, 0, sizeof(fullPath));
	gen_base_path(&b, basePath);
	sprintf(fullPath, "%s\\%s\\%d.ive", iveroot, basePath, findResultHTI[idx].tileID);
	sprintf(simpleName, "%d.ive", findResultHTI[idx].tileID);
	sprintf(srcPath, "%s\\%s", iveroot, basePath);
	
	/* 根据不同的纹理精度级别，确定不同的LOD级别。*/
	int lodLevel = GetLodLevelWithTextureHighLevel(findResultHTI[idx].highLevel);
	
	/* 制作命令行 */
	//fprintf(fout, "osgdem -t %s -d %s --geocentric --tile-terrain-size 16 --tile-image-size 1024 --radius-to-max-visible-distance-ratio 3 -ge %f %f 0.25 0.125 -l 5 -v 1.0 -o %s\n", texture_name, height_name, W, S, fullPath);
	sprintf(arg[ 0], "osgdem.exe"); 							argptr[ 0] = (char *)&arg[ 0];
	sprintf(arg[ 1], "-t"); 									argptr[ 1] = (char *)&arg[ 1];
	sprintf(arg[ 2], "%s", temTextGeoFile);						argptr[ 2] = (char *)&arg[ 2];
	sprintf(arg[ 3], "-d"); 									argptr[ 3] = (char *)&arg[ 3];
	sprintf(arg[ 4], "%s", temHgtGeoFile); 						argptr[ 4] = (char *)&arg[ 4];
	sprintf(arg[ 5], "--geocentric"); 							argptr[ 5] = (char *)&arg[ 5];
	sprintf(arg[ 6], "--tile-terrain-size"); 					argptr[ 6] = (char *)&arg[ 6];
	sprintf(arg[ 7], "16"); 									argptr[ 7] = (char *)&arg[ 7];
	sprintf(arg[ 8], "--tile-image-size"); 						argptr[ 8] = (char *)&arg[ 8];
	sprintf(arg[ 9], "1024"); 									argptr[ 9] = (char *)&arg[ 9];
	sprintf(arg[10], "--radius-to-max-visible-distance-ratio"); argptr[10] = (char *)&arg[10];
	sprintf(arg[11], "3"); 										argptr[11] = (char *)&arg[11];
	sprintf(arg[12], "--POLYGONAL"); 							argptr[12] = (char *)&arg[12];
	sprintf(arg[13], "--mip_mapping_imagery"); 					argptr[13] = (char *)&arg[13];
	sprintf(arg[14], "--max_anisotropy"); 						argptr[14] = (char *)&arg[14];
	sprintf(arg[15], "8"); 										argptr[15] = (char *)&arg[15];
	sprintf(arg[16], "-ge"); 									argptr[16] = (char *)&arg[16];
	sprintf(arg[17], "%f", W); 									argptr[17] = (char *)&arg[17];
	sprintf(arg[18], "%f", S); 									argptr[18] = (char *)&arg[18];
	sprintf(arg[19], "%f", offLon); 							argptr[19] = (char *)&arg[19];
	sprintf(arg[20], "%f", offLat); 							argptr[20] = (char *)&arg[20];
	sprintf(arg[21], "-l");	 									argptr[21] = (char *)&arg[21];
	sprintf(arg[22], "%d", lodLevel /*currTi.lod*/);			argptr[22] = (char *)&arg[22];
	sprintf(arg[23], "-v");	 									argptr[23] = (char *)&arg[23];
	sprintf(arg[24], "%f", currTi.scale); 						argptr[24] = (char *)&arg[24];
	sprintf(arg[25], "-o");	 									argptr[25] = (char *)&arg[25];
	sprintf(arg[26], "%s", fullPath);							argptr[26] = (char *)&arg[26];
	int argIdx = 27;

	/* 做vpb */
	if(osgdem_main(argIdx, argptr))
	{
		MessageToDialog("VPB运行失败\n");
		return 0;
	}

	/* 集中写STG */
	sprintf(ctileID, "%d", findResultHTI[idx].tileID);
	if(!WriteAllStgs(ctileID, 0))
	{
		MessageToDialog("创建STG文件失败\n");
		return 0;
	}
	
	/* 复制一份最新生成的地景文件到 BackupLib 目录下。 */
	char msg[STRING_STANDARD_LEN];
	sprintf(msg, "备份地景文件组 --- %s", simpleName);
	MessageToDialog(msg);
	BackUpFileToTemp(simpleName, srcPath);
	
	return 1;
}

/*============================================================================================
 *  下面是若干个初始化及配置管理过程。
 *============================================================================================ */

/* 重写用户配置文件userconfig.ini */
int CheckConfig(int rewrite)
{
	fi = fopen(configFile, "rt");
	if((fi == NULL)||(rewrite)||(updateTTNode))
	{
		//关闭已经读打开的配置文件。
		needUpdate = 0;
		if(fi) fclose(fi);
		
		/* 创建用户配置文件 */
		fo = fopen(userConfig, "wt");
		if(fo == NULL)
		{
			return 0;
		}
		
		/* 写入用户配置文件，先写入目录 */
		fprintf(fo, "%s\n", currTi.txtDir);
		fprintf(fo, "%s\n", currTi.hgtDir);
		fprintf(fo, "%s\n", currTi.iveDir);
		fprintf(fo, "%s\n", currTi.resolution);
		fprintf(fo, "%f\n", currTi.minLongi);
		fprintf(fo, "%f\n", currTi.maxLongi);
		fprintf(fo, "%f\n", currTi.minLati );
		fprintf(fo, "%f\n", currTi.maxLati );
		fprintf(fo, "%s\n", getCurrLODString().c_str());
		fprintf(fo, "%f", currTi.scale );

		WriteHTxtInfo2Config(fo);
		WriteHHgtInfo2Config(fo);
		WriteHighInfo2Config(fo);
		WriteCultInfo2Config(fo);
		WriteShpInfo2Config(fo);
		WriteBldInfo2Config(fo);
		WriteHighCultureInfo2Config(fo);
		fclose(fo);	

	}else{
		fclose(fi);
	}
	
	return 1;
}

/* 检查几个输入的目录是否为空 */
int CheckDirectories()
{
	if(strlen(currTi.txtDir) == 0)  strcpy(currTi.txtDir, baseTextDir);
	if(strlen(currTi.hgtDir) == 0)	strcpy(currTi.hgtDir, baseHgtDir); 
	if(strlen(currTi.iveDir) == 0)
		return 0;
	else
		return 1;
	
}

/* 设置全局变量 */
int SetGlobalVars()
{
	CreateUserSceneryDirectory();
	strcpy(hgtroot, baseHgtDir);			//currTi.hgtDir; 
	sprintf(iveroot, "%s\\Scenery\\Terrain", currTi.iveDir);
	sprintf(objroot, "%s\\Scenery\\Objects", currTi.iveDir);
	return 1;
}


/* 初始化检查工作 */
int CheckInitiation(int refresh)
{
	int retVal = 0;
	int iFlag = 1;

	/* 检查目录是否存在 */
	updateTTNode = 0;
	if(!CheckDirectories())
	{
		MessageToDialog("目录没有输入完全.");
		return -1;
	}

#if USE_JPEG_DATA_TO_CREATE_TERRAIN
	/* 检查Map信息是否已经加载, 如果没有，则遍历。*/
	if((nodeID == 0)||(TextDirIdx == 0)||refresh)
	{
		retVal = 0;
		totalnode.clear();
		nodeID = 0;
		totalcult.clear();
		cultID = 0;

		retVal += TravelDirectoryToReadMap(baseTextDir);
		if(!retVal)
		{
			MessageToDialog("遍历纹理目录读取经纬度信息失败\n");
			iFlag = 0;
		}
		updateTTNode = 1;
	}
#endif

	/* 检查高精度的地理信息是否已经加载, 如果没有，则遍历。*/
	if((numHTI == 0)||(currTi.resolution[0] == 'H')||refresh)
	{
		char userPath[STRING_STANDARD_LEN];

		/* 直接查找用户定义高程目录，以查找用户高程数据 */
		sprintf(userPath, "%s\\%s", baseHgtDir, USER_HIGHRESOLUTION);

		retVal = 0;
		isCurrTxtOrHgt = STATUS_HGT;
		retVal += TravelDirectoryToReadTif(userPath);
		if(strcmp(baseHgtDir, currTi.hgtDir))
		{			
			retVal += TravelDirectoryToReadTif(currTi.hgtDir);
		}
		if(!retVal)
		{
			MessageToDialog("遍历纹理目录读取经纬度信息失败\n");
			iFlag = 0;
		}

		/* 直接查找用户定义纹理目录，以查找用户高分辨率纹理数据 */
		sprintf(userPath, "%s\\%s", baseTextDir, USER_HIGHRESOLUTION);

		retVal = 0;
		isCurrTxtOrHgt = STATUS_TEXTURE;
		retVal += TravelDirectoryToReadTif(userPath);
		if(strcmp(baseTextDir, currTi.txtDir))
		{			
			retVal += TravelDirectoryToReadTif(currTi.txtDir);
		}
		if(!retVal)
		{
			MessageToDialog("遍历纹理目录读取用户纹理数据失败\n");
			iFlag = 0;
		}

		/* 直接查找系统自带的纹理目录，以查找系统提供的50米分辨率纹理数据 */
		sprintf(userPath, "%s\\%s", baseTextDir, SYSTEM_LOWRESOLUTION);

		retVal = 0;
		isCurrTxtOrHgt = STATUS_TEXTURE;
		retVal += TravelDirectoryToReadTif(userPath);
		if(!retVal)
		{
			MessageToDialog("遍历纹理目录读取系统纹理数据失败\n");
			iFlag = 0;
		}

		updateTTNode = 1;
	}

	/* 检查ac文化信息是否已经加载, 如果没有，则遍历。*/
	if((cultID == 0)||refresh)
	{
		retVal = 0;
		retVal += TravelDirectoryToReadAc(baseTextDir);
		if(!retVal)
		{
			MessageToDialog("遍历纹理目录读取AC文化信息失败\n");
			iFlag = 0;
		}
		updateTTNode = 1;
	}

	/* 检查shp文化信息是否已经加载, 如果没有，则遍历。*/
	if((shpID == 0)||refresh)
	{
		retVal = 0;
		retVal += TravelDirectoryToReadShp(baseTextDir);
		if(!retVal)
		{
			MessageToDialog("遍历纹理目录读取SHP文化信息失败\n");
			iFlag = 0;
		}
		updateTTNode = 1;
	}

	/* 检查是否有用户三维信息 */
	if((bldNum == 0)||refresh)
	{
		char userPath[STRING_STANDARD_LEN];

		/* 直接查找用户定义三维物体目录，以查找用户三维物体数据 */
		sprintf(userPath, "%s\\data\\Models\\User", fsData);

		retVal = 0;
		retVal += TravelDirectoryToReadFlt(userPath);
		if(!retVal)
		{
			MessageToDialog("遍历三维物体目录读取三维物体信息失败\n");
			iFlag = 0;
		}
		updateTTNode = 1;
	}


	return iFlag;
}

/* 供界面刷新的时候用 */
int RefreshTheConfigFile()
{
	if( CheckInitiation(1) == -1)
	{
		return 0;
	}
	return 1;
}

/* 读取并转换配置文件 */ 
int ConvertConfig()
{
	/* 先尝试打开用户配置文件 */
	fi = fopen(userConfig, "rt");
	if(fi == NULL)
	{
		/* 如果打不开用户配置文件，则打开全局配置文件，并且创建一个用户配置文件，将有关用户数据写入该文件 */
		/* 打开配置文件 */
		fi = fopen(configFile, "rt");
		if(fi == NULL)
		{
			strcpy(currTi.txtDir, baseTextDir);
			strcpy(currTi.hgtDir, baseHgtDir);
			currTi.lod      = 777;
			currTi.scale    = 1.0;
			initCheck.usrConfigFile = false;
			return 1;
		}
	
		/* 读取目录 */
		ReadALineToArray(fi);  strcpy(currTi.txtDir, cf[0]); 		
		if(currTi.txtDir[strlen(currTi.txtDir) - 1] == 0xA) currTi.txtDir[strlen(currTi.txtDir) - 1] = 0;
		if(strcmp(currTi.txtDir, baseTextDir))		strcpy(currTi.txtDir, baseTextDir);
		ReadALineToArray(fi);  strcpy(currTi.hgtDir, cf[0]);
		if(currTi.hgtDir[strlen(currTi.hgtDir) - 1] == 0xA) currTi.hgtDir[strlen(currTi.hgtDir) - 1] = 0;
		if(strcmp(currTi.hgtDir, baseHgtDir))		strcpy(currTi.hgtDir, baseHgtDir);
		ReadALineToArray(fi);  strcpy(currTi.iveDir, cf[0]);
		if(currTi.iveDir[strlen(currTi.iveDir) - 1] == 0xA) currTi.iveDir[strlen(currTi.iveDir) - 1] = 0;
		if(strcmp(currTi.iveDir, fsData))			strcpy(currTi.iveDir, fsData);
		ReadALineToArray(fi);  strcpy(currTi.resolution, cf[0]);
		if(currTi.resolution[strlen(currTi.resolution) - 1] == 0xA) currTi.resolution[strlen(currTi.resolution) - 1] = 0;
		ReadALineToArray(fi);  currTi.minLongi = atof(cf[0]);
		ReadALineToArray(fi);  currTi.maxLongi = atof(cf[0]);
		ReadALineToArray(fi);  currTi.minLati  = atof(cf[0]);
		ReadALineToArray(fi);  currTi.maxLati  = atof(cf[0]);
		ReadALineToArray(fi);
		{
			if(strlen(cf[0]) > 0)	currTi.lod  = atoi(cf[0]);
			if(strlen(cf[1]) > 0)	currTi.lod  = currTi.lod * 10 + atoi(cf[1]);
			if(strlen(cf[2]) > 0)	currTi.lod  = currTi.lod * 10 + atoi(cf[2]);
		}
		ReadALineToArray(fi);  currTi.scale  = atof(cf[0]);
	
		/* 读取Map信息 */
		ReadMapInfoFromConfig(fi);
		ReadDirInfoFromConfig(fi);
	
		ReadHTxtInfoFromConfig(fi);
		ReadHHgtInfoFromConfig(fi);
		ReadLowInfoFromConfig(fi);
		ReadHighInfoFromConfig(fi);
		ReadCultInfoFromConfig(fi);
		ReadBldInfoFromConfig(fi);
		fclose(fi);
		
		/* 将原有存在highTxtnode里面的50米分辨率数据转存到lowTxtnode里面 */
		MoveHigh50ToLow();
		
		/* 创建用户配置文件 */
		fo = fopen(userConfig, "wt");
		if(fo == NULL)
		{
			initCheck.usrConfigFile = false;
			return 0;
		}
		
		/* 写入用户配置文件，先写入目录 */
		fprintf(fo, "%s\n", currTi.txtDir);
		fprintf(fo, "%s\n", currTi.hgtDir);
		fprintf(fo, "%s\n", currTi.iveDir);
		fprintf(fo, "%s\n", currTi.resolution);
		fprintf(fo, "%f\n", currTi.minLongi);
		fprintf(fo, "%f\n", currTi.maxLongi);
		fprintf(fo, "%f\n", currTi.minLati );
		fprintf(fo, "%f\n", currTi.maxLati );
		fprintf(fo, "%s\n", getCurrLODString().c_str() );
		fprintf(fo, "%f", currTi.scale );

		WriteHTxtInfo2Config(fo);
		WriteHHgtInfo2Config(fo);
		WriteHighInfo2Config(fo);
		WriteCultInfo2Config(fo);
		WriteShpInfo2Config(fo);
		WriteBldInfo2Config(fo);
		WriteHighCultureInfo2Config(fo);
		fclose(fo);	

		/* 然后重新写全局配置文件 */
		fo = fopen(configFile, "wt");
		if(fo == NULL)
		{
			return 0;
		}
		WriteMapInfo2Config(fo);
		WriteDirInfo2Config(fo);
		WriteLTxtInfo2Config(fo);
		WriteLowInfo2Config(fo);
		fclose(fo);	

	}else{
		/* 如果有用户配置文件， 先读取用户配置文件 */
		/* 读取目录 */
		ReadALineToArray(fi);  strcpy(currTi.txtDir, cf[0]); 		
		if(currTi.txtDir[strlen(currTi.txtDir) - 1] == 0xA) currTi.txtDir[strlen(currTi.txtDir) - 1] = 0;
		if(strcmp(currTi.txtDir, baseTextDir))		strcpy(currTi.txtDir, baseTextDir);
		ReadALineToArray(fi);  strcpy(currTi.hgtDir, cf[0]);
		if(currTi.hgtDir[strlen(currTi.hgtDir) - 1] == 0xA) currTi.hgtDir[strlen(currTi.hgtDir) - 1] = 0;
		if(strcmp(currTi.hgtDir, baseHgtDir))		strcpy(currTi.hgtDir, baseHgtDir);
		ReadALineToArray(fi);  strcpy(currTi.iveDir, cf[0]);
		if(currTi.iveDir[strlen(currTi.iveDir) - 1] == 0xA) currTi.iveDir[strlen(currTi.iveDir) - 1] = 0;
		if(strcmp(currTi.iveDir, fsData))			strcpy(currTi.iveDir, fsData);
		ReadALineToArray(fi);  strcpy(currTi.resolution, cf[0]);
		if(currTi.resolution[strlen(currTi.resolution) - 1] == 0xA) currTi.resolution[strlen(currTi.resolution) - 1] = 0;
		ReadALineToArray(fi);  currTi.minLongi = atof(cf[0]);
		ReadALineToArray(fi);  currTi.maxLongi = atof(cf[0]);
		ReadALineToArray(fi);  currTi.minLati  = atof(cf[0]);
		ReadALineToArray(fi);  currTi.maxLati  = atof(cf[0]);
		ReadALineToArray(fi);  
		{
			if(strlen(cf[0]) > 0)	currTi.lod  = atoi(cf[0]);
			if(strlen(cf[1]) > 0)	currTi.lod  = currTi.lod * 10 + atoi(cf[1]);
			if(strlen(cf[2]) > 0)	currTi.lod  = currTi.lod * 10 + atoi(cf[2]);
		}
		ReadALineToArray(fi);  currTi.scale  = atof(cf[0]);
		
		ReadHTxtInfoFromConfig(fi);
		ReadHHgtInfoFromConfig(fi);
		ReadHighInfoFromConfig(fi);
		ReadCultInfoFromConfig(fi);
		ReadShpInfoFromConfig(fi);
		ReadBldInfoFromConfig(fi);
		ReadHighCultureInfoFromConfig(fi);
		fclose(fi);
	
		fi = fopen(configFile, "rt");
		if(fi == NULL)
		{
			return 0;
		}
		/* 再读取全局配置文件 */	
		ReadMapInfoFromConfig(fi);
		ReadDirInfoFromConfig(fi);
		ReadLTxtInfoFromConfig(fi);
		ReadLowInfoFromConfig(fi);
		fclose(fi);	
	}
	return 1;
}

/* 检查并创建用户高分辨率地景库目录 */
int CreateUserSceneryDirectory()
{
	/* 先进入currTi.iveDir目录 */
	int io = chdir(currTi.iveDir);
	if(io == -1)
	{
		io = mkdir(currTi.iveDir);
		if(io == -1)
		{
			return 0;
		}
		io = chdir(currTi.iveDir);
	}

	char userScenery[_MAX_PATH];
	sprintf(userScenery, "Scenery");
	io = chdir(userScenery);
	if(io == -1)	/* 没有Scenery这个目录的情况 */
	{
		/* 建scenery */
		io = mkdir(userScenery);
		if(io == -1)
		{
			return 0;
		}
		io = chdir(userScenery);
		/* 建terrain */
		sprintf(userScenery, "Terrain");
		io = mkdir(userScenery);
		if(io == -1)
		{
			return 0;
		}
		/* 建objects */
		sprintf(userScenery, "Objects");
		io = mkdir(userScenery);
		if(io == -1)
		{
			return 0;
		}
	}else{	/* 有Scenery这个目录的情况 */
		io = chdir(userScenery);

		/* 查terrain */
		sprintf(userScenery, "Terrain");
		io = chdir(userScenery);
		if(io == -1)	/* 没有terrain目录的情况 */
		{
			io = mkdir(userScenery);
			if(io == -1)
			{
				return 0;
			}
		}else{
			io = chdir("..");	/* 返回上一级 */
		}
		
		/* 查Objects */
		sprintf(userScenery, "Objects");
		io = chdir(userScenery);
		if(io == -1)	/* 没有Objects目录的情况 */
		{
			io = mkdir(userScenery);
			if(io == -1)
			{
				return 0;
			}
		}else{
			io = chdir("..");	/* 返回上一级 */
		}
	}
	
	/* 最后返回工作目录下*/
	io = chdir(workDirectory);
	
	return 1;
}

/* 系统初始化，在程序窗口出现之前完成的工作。 */
int initCreateTerrain()
{
	int io;

	/* 初始化检查项目表 */
	initCheck.backDirectory  = true;
	initCheck.pathOfHgt  = true;
	initCheck.pathOfTarget	= true;
	initCheck.pathOfTexture  = true;
	initCheck.tempDirectory  = true;
	initCheck.usrConfigFile  = true;

	//初始化全局变量
	_getcwd(workDirectory, _MAX_PATH);
	sprintf(tempDirectory, "%s\\Temp", workDirectory);
	sprintf(temTextFile		, "%s\\%s.%s",  tempDirectory, TEMP_TEXTFILE	, TEMP_EXTENDNAME);
	sprintf(temTextGeoFile  , "%s\\%s.%s",  tempDirectory, TEMP_TEXTGEOFILE, TEMP_EXTENDNAME);
	sprintf(temHgtFile      , "%s\\%s.%s",  tempDirectory, TEMP_HGTFILE	, TEMP_EXTENDNAME);
	sprintf(temHgtGeoFile   , "%s\\%s.%s",  tempDirectory, TEMP_HGTGEOFILE	, TEMP_EXTENDNAME);
	
	/* 打开全局Log文件。*/
	sprintf(temLogFile   , "%s\\%s",  tempDirectory, LOG_FILE);
	fLog = fopen(temLogFile, "wt");

	ifFoundBasePath = 1;	needUpdate = 0;
	totalnode.clear();
	highTxtnode.clear();
	highHgtnode.clear();
	totalcult.clear();
	totalUserBld.clear();
	lowTxtnode.clear();
	totalShp.clear();
	csmList.clear();
	colors.clear();
	currHighLevel = hTxtNum = hHgtNum = nodeID = MessageTo = bldNum = lTxtNum = 0;
	memset(&currTi, 0, sizeof(currTi));
	sprintf(configFile		, "%s\\terraincfg.ini", 		workDirectory	);
	sprintf(userConfig		, "%s\\userconfig.ini", 		workDirectory	);
	rebuildTerrainForGeneapt  = false;
	CurrentCreateGeneralAirport = false;
	CurrentLampType = false;
	memset(&currHgt, 0, sizeof(currHgt));

	/* 读取灯光颜色定义文件 */
	ReadColorsFromFile();

	/* 创建一个互斥量，供信息输出使用 */
	MsgMutex = CreateMutex(NULL, false, NULL);

	/* 查找基础数据路径 */
	FindBasePathName();

	/* 检查临时路径是否存在，如果不在，则创建一个。*/
	io = chdir(tempDirectory);
	if(io)
	{
		io = mkdir(tempDirectory);
		if(io)
		{
			initCheck.tempDirectory = false;
		}
		io = chdir(tempDirectory);
	}
	
	/* 进入临时路径，创建一个备份路径 */
	sprintf(backDirectory, "BackupLib");
	io = chdir(backDirectory);
	if(io)
	{
		io = mkdir(backDirectory);
		if(io)
		{
			initCheck.backDirectory = false;
		}
		io = chdir(backDirectory);
	}
	_getcwd(backDirectory, _MAX_PATH);

	/* 初始化高程应用数据*/
	rootNode = new osg::Group;
	openedTiles.clear();
	allLakeSkirt   = new osg::Group;
	allRoadSkirt   = new osg::Group;
	allDynamicGrass   = new osg::Group;
	allIslands	   = new osg::Group;
	allLakeInCurrLonLat	   = new osg::Group;	

	/* 初始化批生成地景用数据 */
	TextDirIdx = 0;
	allDirs.clear();
	numHTI = 0;
	allHTI.clear();

	/* 初始化随机数生成 */
	mt_init(&random_seed, 12345);

#if 0
	char aptFile[STRING_STANDARD_LEN];
	sprintf(aptFile, "F:\\FG\\data\\Models\\User\\zh8-flat\\zh8-flat.flt");
	osg::ref_ptr<osg::Geode> aptGD = new osg::Geode;
	ReadFltAndVistIt(aptFile, aptGD);
	osgDB::writeNodeFile(*aptGD, "F:\\aptgd.osg");
#endif

	/* 在系统环境变量PATH中添加自己的工作路径 */
	CreateEnvPath();

	/* 读取并转换配置文件 */
	ConvertConfig();
	
	/* 设置全局变量 iveroot和 hgtroot */
	SetGlobalVars();
	
	return 1;
}

void DoInitializationCheck()
{
	/* 把当前的配置写入Log文件 */
	WriteTerrainCfgIntoLogFile();
	WriteCultureCfgIntoLogFile();

	char msg[STRING_STANDARD_LEN * 2];

	sprintf(msg, "开始初始化检查......");
	sprintf(msg, "当前工作路径: %s", workDirectory);
	MessageToDialog(msg);

	/* 检查项目：临时路径 */
	if(initCheck.tempDirectory)
		sprintf(msg, "检查项目：临时路径 --- 正常");
	else
		sprintf(msg, "检查项目：临时路径 --- 失败");
	MessageToDialog(msg);
	
	if(initCheck.tempDirectory)
	{
		sprintf(msg, "临时路径: %s", tempDirectory);
		MessageToDialog(msg);
	}

	/* 检查项目：原始库备份路径*/
	if(initCheck.backDirectory)
		sprintf(msg, "检查项目：原始库备份路径 --- 正常");
	else
		sprintf(msg, "检查项目：原始库备份路径 --- 失败");
	MessageToDialog(msg);

	if(initCheck.backDirectory)
	{
		sprintf(msg, "原始库备份路径: %s", backDirectory);
		MessageToDialog(msg);
	}

	/* 检查项目：用户配置文件*/
	if(initCheck.usrConfigFile)
		sprintf(msg, "检查项目：用户配置文件 --- 正常");
	else
		sprintf(msg, "检查项目：用户配置文件 --- 失败");
	MessageToDialog(msg);

	if(initCheck.usrConfigFile)
	{
		sprintf(msg, "用户配置文件: %s", userConfig);
		MessageToDialog(msg);
	}

	/* 检查项目：纹理素材路径*/
	if(initCheck.pathOfTexture)
		sprintf(msg, "检查项目：纹理素材路径 --- 正常");
	else
		sprintf(msg, "检查项目：纹理素材路径 --- 失败");
	MessageToDialog(msg);

	if(initCheck.pathOfTexture)
	{
		sprintf(msg, "纹理素材路径: %s", baseTextDir);
		MessageToDialog(msg);
	}

	/* 检查项目：高程素材路径*/
	if(initCheck.pathOfHgt)
		sprintf(msg, "检查项目：高程素材路径 --- 正常");
	else
		sprintf(msg, "检查项目：高程素材路径 --- 失败");
	MessageToDialog(msg);

	if(initCheck.pathOfHgt)
	{
		sprintf(msg, "高程素材路径: %s", baseHgtDir);
		MessageToDialog(msg);
	}

	/* 检查项目：生成结果路径*/
	if(initCheck.pathOfTarget)
		sprintf(msg, "检查项目：生成结果路径 --- 正常");
	else
		sprintf(msg, "检查项目：生成结果路径 --- 失败");
	MessageToDialog(msg);
		
	if(initCheck.pathOfTarget)
	{
		sprintf(msg, "生成结果路径: %s", fsData);
		MessageToDialog(msg);
	}
	
}



/*============================================================================================
 *  下面是生成地景主过程。从对话框直接调用。
 *============================================================================================ */

/**************************************************
 *
 *   根据输入经纬度值，生成1经度x1纬度区域的地形，当前块的信息在currTi里面指定。
 *   noVpb = 0 --- 全做； = 1 --- 只做纹理和高程；  = 2 --- 只做高程
 *
 **************************************************/
int CreatePerDegreeTerrain(int noVpb)
{
	char msg[STRING_STANDARD_LEN], ctileID[10];
	double src_lon;
	double src_lat;
	Bucket b;
	double W, S, E, N;
	int mainIndex, currTile = 1;
	int textFlag, hgtFlag;
	FILE *ftmp;
	
	textFlag = hgtFlag = 0;

	/* 计算出所有tile所使用的图片集合 */
	if(noVpb)
		sprintf(msg, "开始制作：经度%3d, 纬度%2d 区的高程及纹理数据......", (int)currTi.minLongi, (int)currTi.minLati);
	else
		sprintf(msg, "开始制作：经度%3d, 纬度%2d 区的地景......", (int)currTi.minLongi, (int)currTi.minLati);
	MessageToDialog(msg);

	/* 计算当前区域的极限经纬度 */			
	src_lon = currTi.minLongi;
	src_lat = currTi.minLati;
	W = src_lon;
	E = src_lon + 1.0;
	S = src_lat;
	N = src_lat + 1.0;

	/* 事先在临时目录里面查找，看看需要的临时文件是否存在，不存在的时候再拼接 */
	MessageToDialog("开始拼接高程和纹理数据");
	
	/* 尝试打开高程文件 */
	sprintf(temHgtGeoFile   , "%s\\%s%d%d.%s",  tempDirectory, TEMP_HGTGEOFILE	, currPZ->lon, currPZ->lat, TEMP_EXTENDNAME);
	if( (ftmp = fopen(temHgtGeoFile, "rb")) != NULL)
	{
		fclose(ftmp);
		hgtFlag = 1;
	}

	/* 开始拼接高程数据，如果高程数据没有，isOK = 2 */
	/* 计算出当前块的Bucket */
	set_bucket(&b, src_lon, src_lat);
	GetHgtFileList(W, S, E, N);
	mainIndex = gen_index(&b);
	
	if(!hgtFlag)
	{
		if(MergeHgt() == false)
		{
			MessageToDialog("拼接高程失败\n");
			return 0;
		}
	
		if(!convertHgtGeotiff())
		{
			MessageToDialog("创建高程GeoTiff文件失败\n");
			return 0;
		}	
	}

	/* 如果noVpb==2表示只拼接高程数据，则直接退出 */
	if(noVpb==2) return 1;

	/* 尝试打开纹理文件 */
	sprintf(temTextGeoFile  , "%s\\%s%d%d.%s",  tempDirectory, TEMP_TEXTGEOFILE, currPZ->lon, currPZ->lat, TEMP_EXTENDNAME);
	if( (ftmp = fopen(temTextGeoFile, "rb")) != NULL)
	{
		fclose(ftmp);
		textFlag = 1;
	}

	if(!textFlag)
	{
		/* 开始拼接纹理数据 */
		if(!FindJpgsInZoneAndSort(W, E, S, N))	
		{
			MessageToDialog("计算图片失败。\n");
			return 0;
		}
		
		/* 根据计算结果拼接图片 */
		if(!MergePictures())
		{
			MessageToDialog("拼接图片失败\n");
			return 0;
		}
	  	
		/*把tif图片转换为Geotiff图片*/
		if(!convertTextureGeotiff())
		{
			MessageToDialog("创建纹理GeoTiff文件失败\n");
			return 0;
		}
	}

	/* 如果选择不做vpb，则直接退出 */
	if(noVpb==1) return 1;
  	
	/*装配vpb批处理命令。*/
	MessageToDialog("调用VPB生成地景");
	if(!DoVpb(&b, 1))
	{
		MessageToDialog("执行VPB失败\n");
		return 0;
	}

	/* 集中写STG */
	sprintf(ctileID, "%d%d", (int)W, (int)S);
	if(!WriteAllStgs(ctileID, 0))
	{
		MessageToDialog("创建STG文件失败\n");
		return 0;
	}
	
	return 1;
}

/**************************************************
 *
 *   批量生成地形
 *
 **************************************************/
int BatchCreateTerrain()
{
	MessageTo = 0;
	MessageToDialog("批量创建地形......");
	
	InitOnDemand();

	/*初始化检查工作*/
	if( CheckInitiation(0) == -1)
	{
		MessageToDialog("初始化检查工作失败.");
		return 0;
	}
	
	/*检查有没有配置文件，如果没有就写一个配置文件 */
	if(!CheckConfig(needUpdate))
	{
		MessageToDialog("检查配置文件失败。\n");
		return 0;
	}
	
	/* 将当前地景生成的配置信息写入Log文件 */
	WriteTerrainCfgIntoLogFile();
	
	/* 设置全局变量 iveroot和 hgtroot */
	SetGlobalVars();

	/* 只进行用户高分辨率的地景制作 */
	DoHightResolutionTerrain();

	MessageToDialog("本次地景生成完毕!!!!!!");
	return 1;
}








/*============================================================================================
 *  下面是文化信息的制作过程，目前支持两种文化信息制作：ac文件 和 shp文件
 *============================================================================================ */

/**********************************************************
 *
 *       遍历地景纹理目录，搜寻.Shp文件，提取信息存入数据库
 *       
 **********************************************************/
/* 根据文件名搜索库，若库里存在，则返回1，否则返回0 */
int CheckShpFromtotalShp(char *fn)
{
	int i, retVal = 0;
	
	for(i=0; i<shpID; i++)
	{
		if(!strcmp(fn, totalShp[i].fileName))
		{
			retVal = 1; break;
		}
	}
	if(!retVal) needUpdate = 1;
	return retVal;
}

/* 遍历当前目录下的所有子目录，以查找 shp 文件，使用递归方式。 */
int TravelSubDirectoryForShp(char *dir)
{
	int dio, dfio, dirIo, i;
	struct _finddata_t sdff;

	//查找当前路径下的子目录；
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//如果不是子目录也不是文件
		if((sdff.name[0] == '.')) goto GET_NEXT_SHP_SUBDIR2;

		//如果是子目录，进入当前子目录
		if(sdff.attrib == _A_SUBDIR)
		{
			//进入子目录
			dirIo = chdir(sdff.name);
			_getcwd(path, _MAX_PATH);
			
			//遍历子目录
			TravelSubDirectoryForShp(path);
			
			//退回到上一层目录
			dirIo = chdir("..");
		}
		
		//如果是文件
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//找到Shp文件
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='S')||(sdff.name[i-3]=='s')) && ((sdff.name[i-2]=='H')||(sdff.name[i-2]=='h')) && ((sdff.name[i-1]=='P')||(sdff.name[i-1]=='p')))  
			{
				/* 先查找库里有没有该文件 */
				if(!CheckShpFromtotalShp(sdff.name))
				{
					CultNode cn;
					_getcwd(path, _MAX_PATH);
					strcpy(cn.fileName, sdff.name);
					LoadShpAndCalculateMM(sdff.name);
					
					cn.minLong = currLayer.minLong;
					cn.maxLati = currLayer.maxLati;
					ReadALineToArray(fi);
					ReadALineToArray(fi);
					cn.maxLong = currLayer.maxLong;
					cn.minLati = currLayer.minLati;
					strcpy(cn.pathName, path);
					cn.destTileID = 0;
					cn.destLongi  = 0; 
					cn.destLati   = 0;
					sprintf(cn.destFile, "none");
					sprintf(cn.destPath, "none");
					cn.ID = shpID;
					cn.isOK = 0;
					fclose(fi);
					totalShp.push_back(cn);		shpID++;
					dirIo = chdir(path);
				}
			}
		}

GET_NEXT_SHP_SUBDIR2:
		//寻找下一个
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* 遍历当前目录，查找该目录下以及所有子目录下扩展名为.shp的文件，读取里面相关信息。 */
int TravelDirectoryToReadShp(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForShp(dPath))
	{
		MessageToDialog("找不到Shp文件");
		return 0;
	}
	
	CheckConfig(needUpdate);

	return 1;
} 

/**********************************************************
 *
 *       遍历地景纹理目录，搜寻.ac文件，提取信息存入数据库
 *       
 **********************************************************/
/* 根据文件名搜索库，若库里存在，则返回1，否则返回0 */
int CheckACFromTotalCult(char *fn)
{
	int i, retVal = 0;
	
	for(i=0; i<cultID; i++)
	{
		if(!strcmp(fn, totalcult[i].fileName))
		{
			retVal = 1; break;
		}
	}
	if(!retVal) needUpdate = 1;
	return retVal;
}

/* 遍历当前目录下的所有子目录，以查找 ac 文件，使用递归方式。 */
int TravelSubDirectoryForAc(char *dir)
{
	int dio, dfio, dirIo, i, j;
	struct _finddata_t sdff;

	//查找当前路径下的子目录；
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//如果不是子目录也不是文件
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR2;

		//如果是子目录，进入当前子目录
		if(sdff.attrib == _A_SUBDIR)
		{
			//进入子目录
			dirIo = chdir(sdff.name);
			_getcwd(path, _MAX_PATH);
			
			//遍历子目录
			TravelSubDirectoryForAc(path);
			
			//退回到上一层目录
			dirIo = chdir("..");
		}
		
		//如果是文件
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//找到ac文件
			i = strlen(sdff.name);
			if(((sdff.name[i-2]=='A')||(sdff.name[i-2]=='a')) && ((sdff.name[i-1]=='C')||(sdff.name[i-1]=='c')) )  
			{
				//制作同文件名.map文件名称字符串
				char mapn[24];
				memset(mapn, 0, sizeof(mapn));
				for(j=0;j<i-2;j++) mapn[j]=sdff.name[j];
				mapn[j++]='m'; mapn[j++]='a'; mapn[j++]='p'; 

				/* 先查找库里有没有该文件 */
				if(!CheckACFromTotalCult(sdff.name))
				{
					CultNode cn;
					strcpy(cn.fileName, sdff.name);
					
					//打开同文件名map文件
					fi = fopen(mapn, "rt");
					if(fi == NULL) 
					{
						char msg[STRING_STANDARD_LEN];
						sprintf(msg, "不能读取文件 %s .", sdff.name);
						MessageToDialog(msg);
						return 0;
					}
					
					//获取图片的地理信息
					while((strcmp("MMPLL,1,", cf[0]))&&(!feof(fi)))  ReadALineToArray(fi);
					if(feof(fi))
					{
						MessageToDialog("读取map文件出错:找不到经纬度信息。");
						return 0;
					}
					cn.minLong = atof(cf[1]);
					cn.maxLati = atof(cf[2]);
					ReadALineToArray(fi);
					ReadALineToArray(fi);
					cn.maxLong = atof(cf[1]);
					cn.minLati = atof(cf[2]);
					_getcwd(path, _MAX_PATH);
					strcpy(cn.pathName, path);
					cn.destTileID = 0;
					cn.destLongi  = 0; 
					cn.destLati   = 0;
					sprintf(cn.destFile, "none");
					sprintf(cn.destPath, "none");
					cn.ID = cultID;
					cn.isOK = 0;
					fclose(fi);
					totalcult.push_back(cn);	cultID++;
				}
			}
		}

GET_NEXT_SUBDIR2:
		//寻找下一个
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* 遍历当前目录，查找该目录下以及所有子目录下扩展名为.ac的文件，读取里面相关信息。 */
int TravelDirectoryToReadAc(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForAc(dPath))
	{
		MessageToDialog("找不到ac文件");
		return 0;
	}
	
	CheckConfig(needUpdate);

	return 1;
} 

/**************************************************
 *
 *   读取ac文件并且制作BTG文件
 *
 **************************************************/
typedef struct _LANDPOINT_
{
	double px;
	double py;				//实际坐标   需要计算转换过来的
	double pz;				
	double lon;
	double lat;				//经纬度坐标，需要计算得来
	double tU;
	double tV;				//纹理坐标   直接可读取的
	int	   iPrePos;			//在坐标序列点中所在的位置。
	int    iCurPos;			//这个点在当前的位置。
}LandPoint;

typedef struct _GEOPOINT_
{
	double lon;
	double lat;
	double elevation;
	int    idx;
}GeoPoint;

typedef struct _REFSURF_
{
	int idx;
	int iNewdx;
	double tU;
	double tV;				//纹理坐标   直接可读取的
}RefSurf;

typedef struct _SURFACE_
{
	int idx[4];
}Surface;

#define TOTALTYPE  9
#define THOLD	   0.1
typedef struct _STYPEINDEX_
{
	int type;			//1---0.1 Stream 2---0.2 Road 3---0.3 Grass
	int idx;
}STypeIndex;

RefSurf Ref[4];
Surface Faces[TOTALTYPE*10][5000];
LandPoint Points[TOTALTYPE*10][5000]; 
int sumPtsObj[TOTALTYPE*10];
int sumFaceObj[TOTALTYPE*10];
STypeIndex typeIdxs[TOTALTYPE*10];	//指定对应点集合的类型
int tiNums;							//保存typeIdxs里面实际有效值的数量
int allTypes[TOTALTYPE];

Geodesy geoCenter;
Carton catCenter;
double gbs_radius;
Carton cPoint[4];
Geodesy gPoint[4];
char sMaterial[20];
Carton cP[5];
int totalGPoint, sumGPoint[TOTALTYPE], sumGFace[TOTALTYPE], sumObjs, curObjs;

/* 用ac文件数据生成 geode ，这个过程没有被调用过。废弃过程。*/
osg::ref_ptr<osg::Geode> createGeodeFromAc()
{
	int i, j, currIdx, offset;
	float fNx, fNy, fNz;
	
	//缺省一个法线
	fNx = 92/127.5; fNy = -30/127.5; fNz = -57/127.5; offset = 0;

	osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
	gnode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	gnode->setName("River");

	//对于每个独立面，创建一个几何体，作为一个drawable加入到gnode中来。
	for(i=0; i<tiNums; i++)
	{
		currIdx = typeIdxs[i].idx;

		osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
		gnode->addDrawable(gm);
		
		//压入顶点
		osg::ref_ptr<osg::Vec3Array> vertex = new osg::Vec3Array;
		float v3 = Points[currIdx][0].pz;
		for(j=0; j<sumPtsObj[currIdx]; j++)
		{
			osg::Vec3d v; v.set(Points[currIdx][j].px, Points[currIdx][j].py, v3/*Points[currIdx][j].pz*/);
			vertex->push_back(v);
		}
		gm->setVertexArray(vertex);
		
		//压入法线
		osg::ref_ptr<osg::Vec3Array> normal = new osg::Vec3Array;
		normal->push_back(osg::Vec3(fNx, fNy, fNz));
		gm->setNormalArray(normal);
		gm->setNormalBinding(osg::Geometry::BIND_OVERALL);
		
		//压入纹理坐标
		osg::ref_ptr<osg::Vec2Array> coord = new osg::Vec2Array;
		for(j=0; j<sumPtsObj[currIdx]; j++)
		{
			osg::Vec2 c; c.set(Points[currIdx][j].tU, Points[currIdx][j].tV);
			coord->push_back(c);
		}
		gm->setTexCoordArray(0, coord);

		//压入面		
		osg::ref_ptr<osg::DrawElementsUShort> deus = new osg::DrawElementsUShort(osg::PrimitiveSet::TRIANGLES);
		for(j=0; j<sumFaceObj[currIdx]; j++)
		{
			unsigned short ua, ub, uc; 
			ua = Faces[currIdx][j].idx[0] - offset;	deus->push_back(ua);
			ub = Faces[currIdx][j].idx[1] - offset;	deus->push_back(ub);
			uc = Faces[currIdx][j].idx[2] - offset;	deus->push_back(uc);
		}
		gm->addPrimitiveSet(deus);
		offset += sumPtsObj[currIdx];
	}

	return gnode.release();
} 

/* 在每次做文化信息之前，先删除以前做过的所有文化信息，包括btg文件和G*.ive文件。 */
int DeleteAllPrevCultureFiles(int lon, int lat)
{
	char fullPath[_MAX_PATH];
	std::vector<string> iveFiles;	iveFiles.clear();
	std::vector<string> stgFiles;	stgFiles.clear();
	FILE *fpstg;

	/* 把所有同目录下的btg文件、以及机场ive文件都找到并删除 */
	Bucket Bt;
	char basePath[STRING_STANDARD_LEN];
	set_bucket(&Bt, lon, lat);
	gen_base_path(&Bt, basePath);
	sprintf(fullPath, "%s\\%s", iveroot, basePath);		
	{
		struct _finddata_t sdff;
		int io, dio, dfio;

		/* 在这个路径下查找 *.btg.gz */
		io = chdir(fullPath);
		
		dio = _findfirst("*.btg.gz", &sdff);
		if((!dio)||(dio==-1)) goto DAPCF_Continue_1;
		dfio = 0;
	
		while(!dfio)
		{
			
			DeleteFile(sdff.name);

			/* 寻找下一个 */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}

DAPCF_Continue_1:

	/* 再把所有同目录下的机场ive文件都找到并删除，包括通用机场和用户定位机场 */
	/* 也就是把这个路径下所有的G*.ive机场边缘线环也删除 */
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
		osg::ref_ptr<osg::Node> node;
		osg::ref_ptr<osg::Geode> gd;
	
		/* 在这个路径下查找 *.osg */
		io = chdir(fullPath);
		
		dio = _findfirst("*.ive", &sdff);
		if((!dio)||(dio==-1)) goto DAPCF_Continue_4;
		dfio = 0;
	
		while(!dfio)
		{
			if(! ((sdff.name[0] == 'G') || (sdff.name[0] == 'g')))
			{
				iveFiles.push_back(sdff.name);
				goto DAPCF_Continue_3;
			}

			DeleteFile(sdff.name);
			
DAPCF_Continue_3:
			/* 寻找下一个 */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}
		
DAPCF_Continue_4:

	/* 最后读取当前目录下所有的stg文件，将stg文件里面的*.btg.gz也删除 */
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
	
		/* 在这个路径下查找 *.osg */
		io = chdir(fullPath);
		
		dio = _findfirst("*.stg", &sdff);
		if((!dio)||(dio==-1)) goto DAPCF_Continue_6;
		dfio = 0;
	
		while(!dfio)
		{
			/* 存起当前目录下所有stg文件名 */
			stgFiles.push_back(sdff.name);

			/* 删除stg文件 */
			DeleteFile(sdff.name);
			/* 寻找下一个 */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}
	
	/* 重新写一下stg文件 */
	for(unsigned int i=0; i<stgFiles.size(); i++)
	{
		fpstg = fopen(stgFiles[i].c_str(), "wt");
		if(fpstg == NULL)	break;
		for(unsigned int j=0; j<iveFiles.size(); j++)
		{
			fprintf(fpstg, "OBJECT %s\n", iveFiles[j].c_str());
		}
		fclose(fpstg);
	}
	
DAPCF_Continue_6:

	return 1;
}



/* 根据输入参数，读取特定目录下的地景ive根文件以及所有子文件，并根据文化信息，通用及用户机场轮廓来切割。这是切割地景的主过程。*/
int ReadTerrainIveAndCut(int lon, int lat)
{
	int io, fio, ffio, subffio, subfio, subidx;
	char fullIveFile[_MAX_PATH];
	char fullPath[_MAX_PATH];
	char fullSubPath[_MAX_PATH];
	char subDir[36];
	struct _finddata_t ivef;
	struct _finddata_t iveSubf;
	char msg[STRING_STANDARD_LEN];

#if 0
	if((lon != 113) || (lat != 22))	return 1;
#endif

	/* 加载所有的文化信息相关文件 */
	/* 第一步，先把所有同目录下的btg文件、以及机场ive文件都打到一个group里面 */
	osg::ref_ptr<osg::Group> gp = new osg::Group;
	Bucket Bt;
	char basePath[STRING_STANDARD_LEN];
	set_bucket(&Bt, lon, lat);
	gen_base_path(&Bt, basePath);
	sprintf(fullPath, "%s\\%s", iveroot, basePath);		
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
	
		/* 在这个路径下查找 *.btg.gz */
		io = chdir(fullPath);
		
		dio = _findfirst("*.btg.gz", &sdff);
		if((!dio)||(dio==-1)) goto RTIAC_Continue_1;
		dfio = 0;
	
		while(!dfio)
		{
			osg::ref_ptr<osg::Group> cgp;

			/* 首先判断是不是文化信息的btg文件名 */
			if(((sdff.name[0] > '9') && (strlen(sdff.name) > 13)))
			{
				/* 如果是裙边，不做切割。*/
				if((sdff.name[0] != 'S') && 
					 (sdff.name[0] != 'O') )		// 新型湖面不做切割
				{
					/* 把这个文化信息btg读入到一个 group里面 */				
					ReadBTGToDataStruct(sdff.name);	
					cgp = CreateGroupFromDataStru();
					for(int j=0; j<cgp->getNumChildren(); ++j)
					{
						gp->addChild(cgp->getChild(j));
					}
				}
			}

			/* 寻找下一个 */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}

RTIAC_Continue_1:

	/* 第二步，再把所有同目录下的机场ive文件都打到同上的一个group里面，包括通用机场和用户定位机场 */
	/* 也就是把这个路径下所有的G*.ive机场边缘线环也加进来 */
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
		osg::ref_ptr<osg::Node> node;
		osg::ref_ptr<osg::Geode> gd;
	
		/* 在这个路径下查找 *.osg */
		io = chdir(fullPath);
		
		dio = _findfirst("*.ive", &sdff);
		if((!dio)||(dio==-1)) goto RTIAC_Continue_4;
		dfio = 0;
	
		while(!dfio)
		{
			if(((sdff.name[0] == 'G') || (sdff.name[0] == 'g')))
			{
				node = osgDB::readNodeFile(sdff.name);
				gd = dynamic_cast<osg::Geode *>(node.get());
				if(gd != NULL)
				{
					gp->addChild(gd.get());
				}
			}
			/* 寻找下一个 */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}
		
RTIAC_Continue_4:
	
	/* 下面是开始切割地景的动作，查找所有ive但要排除机场线环的ive。然后切割地景 */
	/* 查找当前目录下的ive */
	io = chdir(fullPath);
	if(io) 
	{
		sprintf(msg, "不能进入目标目录 -- %s", fullPath);
		MessageToDialog(msg);
		return 0;
	}
	
	/* 修改所有精度地景 */
	sprintf(subDir, "*.ive");
	ffio = _findfirst(subDir, &ivef);	
	if((!ffio)||(ffio == -1))
	{
		MessageToDialog("找不到ive文件。");
		return 0;
	}
	fio = 0;

#if 1
	sprintf(debugFile, "%s\\Culture.ive", tempDirectory);
	osgDB::writeNodeFile(*gp, debugFile);
	
	osg::ref_ptr<osg::Group> csmgp = new osg::Group;
	for(std::vector< osg::ref_ptr<CorrectionSpaceMesh> >::iterator c = csmList.begin(); c != csmList.end(); c++)
	{
		osg::ref_ptr<CorrectionSpaceMesh> csm = *c;
		csmgp->addChild(csm->csmGd);
	}
	sprintf(debugFile, "%s\\csm.ive", tempDirectory);
	osgDB::writeNodeFile(*csmgp, debugFile);
#endif



	/* 暂时清空平滑过渡类列表 */
	//csmList.clear();


	while(!fio)
	{
		/* 分辨率低的地景区块先不切 */
		int flen = strlen(ivef.name);
		if(flen < 12)	goto RTIAC_Continue_2;

#if 0
		if(strcmp(ivef.name, "114230108.ive"))
			goto RTIAC_Continue_2;
#endif

		/* 首先排除掉机场边缘线环的.ive文件 */
		if(ivef.name[0] == 'G') goto RTIAC_Continue_2;
		
		sprintf(msg, "开始修改根地景文件。 -- %s", ivef.name);
		MessageToDialog(msg);
		sprintf(fullIveFile, "%s\\%s", fullPath, ivef.name);

#ifdef	USE_HCI
		/* 得到当前ivef.name的 tileID 号。并检查当前号是否有做过文化信息 */
		char tile[20];
		for(unsigned int i=0; i<strlen(ivef.name); i++)
		{
			if(ivef.name[i] != '.')	tile[i] = ivef.name[i];
			else	{ tile[i] = 0; break;	}
		}
		int tID = atoi(tile);
		if(CheckHCI())	
		{
			goto RTIAC_Continue_2;
		}else
		{
			/* 恢复当前块的地景 */
		}
#endif	//USE_HCI

		/* 返回到工作目录下，以便使用osg动态链接库 */		
		io = chdir(workDirectory);
		if(io) 
		{
			MessageToDialog("不能回到工作路径下。");
			return 0;
		}
		
		//对查找到的文件处理
		//读入ive文件
		UpdateIVEByBTG(gp, fullIveFile);

		/* 返回到当前目录，以便搜索子目录里面的ive文件 */
		io = chdir(fullPath);
		if(io) 
		{
			char msg[STRING_STANDARD_LEN];
			sprintf(msg, "不能进入目标目录 -- %s", fullPath);
			MessageToDialog(msg);
			return 0;
		}

#ifdef	USE_HCI
		/* 把当前ivefile的tileID插入到 allHCI 当中。*/
		InsertHCI(tID);
#endif	//USE_HCI

		/* 检查全局通讯变量，如果变量复位，说明当前地景无变化，底层也不需要再做了。 */
		if(!CurrentTerrainIsCut)
		{
			goto RTIAC_Continue_2;
		}

		/* 以下是处理PageLOD的底层的各个ive文件，所有底层文件存放在 XXXX_root_L9_X0_Y0 目录下 */
		strcpy(subDir, ivef.name);
		int len = strlen(subDir);
		subDir[len - 4] = 0;
		sprintf(fullSubPath, "%s_root_L0_X0_Y0", subDir);
		
		/* 进入到子目录中 */
		io = chdir(fullSubPath);
		if(io) 
		{
			MessageToDialog("不能进入目标目录。");
			goto RTIAC_Continue_2;
		}

		/* 在子目录中查找ive文件 */
		subffio = _findfirst("*.ive", &iveSubf);	
		if((!subffio)||(subffio==-1))
		{
			MessageToDialog("找不到ive文件。");
			goto RTIAC_Continue_2;
		}
		_getcwd(fullSubPath, _MAX_PATH);

		subfio = 0; subidx = 0;
		while(!subfio)
		{
			subidx++;
			sprintf(msg, "开始修改子地景文件--%d.%s", subidx, iveSubf.name);
			MessageToDialog(msg);
			sprintf(fullIveFile, "%s\\%s", fullSubPath, iveSubf.name);

			/* 返回到工作目录下，以便使用osg动态链接库 */		
			io = chdir(workDirectory);
			if(io) 
			{
				MessageToDialog("不能回到工作路径下。");
				return 0;
			}

			//对查找到的文件处理
			//读入ive文件
			UpdateIVEByBTG(gp, fullIveFile);

			/* 返回到当前目录，以便继续搜索子目录里面的其他ive文件 */
			io = chdir(fullSubPath);
			if(io) 
			{
				MessageToDialog("不能进入目标目录。");
				return 0;
			}
	
			//查找下一个符合条件的文件
			subfio = _findnext(subffio, &iveSubf);
		}
		_findclose(subffio);

		/* 返回到当前目录，以便继续搜索下一个ive文件 */
		io = chdir("..");
		if(io) 
		{
			MessageToDialog("不能进入目标目录。");
			return 0;
		}

RTIAC_Continue_2:
		//查找下一个符合条件的文件
		fio = _findnext(ffio, &ivef);
	}
	_findclose(ffio);
	
	return 1;
}

/* 计算包围球半径 */
double calc_gbs(Geodesy *gP, int num)		//gP 是长度为4的数组。
{
    int i, j;
    Bucket b;
    Geodesy *gPtr = gP;
    Geodesy gbs_Center;
    Carton gbs_c;
    double dist_squared;
    double radius_squared = 0;
    
    for(i=0;i<num;i++)
    {
    	GeodToCart(gPtr, &cP[i]);
    	gPtr++;
    }
    
    gPtr = &gP[3];
    set_bucket(&b, gPtr->_lon, gPtr->_lat);
    
    gbs_Center._lon = get_center_lon(&b);
    gbs_Center._lat = get_center_lat(&b);
    gbs_Center._elevation = 0.0;
    
    GeodToCart(&gbs_Center, &gbs_c);
    
    for(j=0;j<num;j++)
    {
    	dist_squared = distance3Dsquared(&gbs_c, &cP[j]);
    	if ( dist_squared > radius_squared ) {
            radius_squared = dist_squared;
        }
    }
    
    return sqrt(radius_squared);
}

/* 把读取ac文件所得到的数据转存到DataStruct里面，以供WriteDataStructToBTG调用 */
int MoveDataToDataStru()
{
	int i, j, m, n, typenum;
	short iTa, iTb, iTc, iT = 0;
	int   pFlag;

	//统计一下当前ac文件总共用了多少个types。
	typenum = 0; for(i=0;i<TOTALTYPE;i++) typenum+=allTypes[i];
	
	/* 将现有变量的值复制到Buctet.cpp内相关对应变量，以便使用Buctet.cpp内的函数 */
	/* 初始化内存 */
	totalNormal = 0;	curPnt.clear();		curText.clear();

	iVersion = 6; iMagicNum = 0x5347; iCreationTime = 0x48170f2c; iNumOfToplvlObject = 5 + typenum;
	dXpart = catCenter.x;
	dYpart = catCenter.y;
	dZpart = catCenter.z;
	fRadiusofBS = (float)gbs_radius;
	gbCenterX = dXpart; gbCenterY = dYpart; gbCenterZ = dZpart;

	/* 顶点信息 */
	for(m=0;m<sumObjs;m++)
	{
		for(n=0;n<sumPtsObj[m];n++)
		{
			TerPoint tp;
			tp.px   = Points[m][n].px;
			tp.py   = Points[m][n].py;
			tp.pz   = Points[m][n].pz; 
			tp.lon  = Points[m][n].lon;
			tp.lat  = Points[m][n].lat;
			tp.elev = 0.0; 
			curPnt.push_back(tp);
		}
	}

	/* 纹理坐标信息 */
	for(m=0;m<sumObjs;m++)
	{
		for(n=0;n<sumPtsObj[m];n++)
		{
			TerTextcoord tt;
			tt.fU = Points[m][n].tU * 40.0;
			tt.fV = Points[m][n].tV * 40.0;
			curText.push_back(tt);
		}
	}
	
	/* 法线信息 */
	curPnt[totalNormal].nx = 92;	curPnt[totalNormal].ny = -30;	curPnt[totalNormal].nz = -57; totalNormal++;	

	/* 循环内写所有的面 */
	curTri.clear();
	for(j=0;j<TOTALTYPE;j++)
	{
		//不同类型不同纹理属性
		switch(j)
		{
		case 0:
			sprintf(sMaterial, "Stream");
			break;
		case 1:
			sprintf(sMaterial, "Road");				
			break;
		case 2:
			sprintf(sMaterial, "detailedtxt_grass");				
			break;
		case 3:
			sprintf(sMaterial, "detailedtxt_sand");				
			break;
		case 4:
			sprintf(sMaterial, "EvergreenForest");				
			break;
		case 5:
			sprintf(sMaterial, "detailedtxt_savanna");				
			break;
		case 6:
			sprintf(sMaterial, "detailedtxt_concrete");				
			break;
		case 7:
			sprintf(sMaterial, "detailedtxt_marsh");				
			break;
		case 8:
			sprintf(sMaterial, "detailedtxt_town");				
			break;
		}
		//计算当前类型Objects的数目，如果没有就继续写其他的。
		curObjs = 0;
		for(int i=0; i<tiNums; i++) { if(typeIdxs[i].type == j+1 ) curObjs++;}
		if(!curObjs) continue;

		TriType tt;
		strcpy(tt.cMtype, sMaterial);
		tt.cObjType = 10;
		tt.cPropVal = 11;
		tt.mElement.clear();
		
		for(m=0;m<sumObjs;m++)
		{
			pFlag = 0;
			for(n=0;n<tiNums;n++)
			{
				if((typeIdxs[n].idx == m) && (typeIdxs[n].type == j+1))
				{
					pFlag = 1; break;
				}
			}
			if(!pFlag) continue;
			
			TriElement te; te.mVertex.clear();
			for(n=0;n<sumFaceObj[m];n++)
			{
				iTa = Faces[m][n].idx[0];
				iTb = Faces[m][n].idx[1];
				iTc = Faces[m][n].idx[2];

				TriVertex tv;
				tv.P = iTa; tv.N = iT; tv.T = iTa; te.mVertex.push_back(tv);
				tv.P = iTb; tv.N = iT; tv.T = iTb; te.mVertex.push_back(tv);
				tv.P = iTc; tv.N = iT; tv.T = iTc; te.mVertex.push_back(tv);
			}//n				
			te.iSum = te.mVertex.size();
			tt.mElement.push_back(te);
		}//m
		tt.iElems = tt.mElement.size();
		curTri.push_back(tt);
	}//j
	return 1;
}

/* 读取ac文件并分析，然后创建btg文件 */
int ReadACandMakeBTG(int idx)
{
#if 1
	unsigned int i;
	int io, fio, ffio, j, n, m, t;
	char dstFile[STRING_STANDARD_LEN];
	char dstPath[STRING_STANDARD_LEN];
	char dPath[STRING_STANDARD_LEN];
	char OneLine[80];
	char cFname[20], cNum[20];
	float minPx, minPy, maxPx, maxPy;
	char path[_MAX_PATH];
	char basepath[STRING_STANDARD_LEN];
	long tile, iTile;
	int numSameTile = 0;
	struct _finddata_t ff;

	//获取当前目录
	_getcwd(path, _MAX_PATH);

	//指定输入输出目录
	sprintf(dPath, "%s", findResultAC[idx].pathName);
	io = chdir(dPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", dPath);
		MessageToDialog(msg);
		return 0;
	}

	ffio = _findfirst(findResultAC[idx].fileName, &ff);	
	if((!ffio)||(ffio == -1))
	{
		MessageToDialog("找不到ac文件。");
		return 0;
	}
	fio = 0;

	while(!fio)
	{
		int totalPnts, currPosi;
		float absPx, absPy, offLon, offLat;
		
		for(i=0;i<strlen(ff.name);i++) if(ff.name[i]!='.') cFname[i]=ff.name[i]; else break; cFname[i] = 0;
		fi = fopen(ff.name, "rt");
		if(fi == NULL) 
		{
			printf("Can not read %s file\n", ff.name);
			exit(-1);
		}

		memset(Points, 0, sizeof(Points));	memset(Faces, -1, sizeof(Faces)); sumObjs = -TOTALTYPE; currPosi = 0; memset(typeIdxs, 0, sizeof(typeIdxs)); t=0; totalGPoint=0;
		minPx = minPy = 9999.0; maxPx = maxPy = -9999.0; memset(sumPtsObj, 0, sizeof(sumPtsObj)); memset(sumFaceObj, 0, sizeof(sumFaceObj)); memset(allTypes, 0, sizeof(allTypes));
		while(!feof(fi))
		{
			fgets(OneLine, 80, fi);

			//第一步：先找到numvert
			if((OneLine[0]=='n') && (OneLine[1]=='u') && (OneLine[2]=='m') && (OneLine[3]=='v') && (OneLine[4]=='e') && (OneLine[5]=='r') && (OneLine[6]=='t'))
			{
				int numpt, pIdx;

				//先调整一下sumObjs的数目；判断上一个Object是否有合法的点，如果没有，则回退。
				sumObjs+=TOTALTYPE;
				if(sumObjs > 0)
				{
					int flags = 0;
					for(i=0; i<TOTALTYPE; i++) flags+=sumGPoint[i];
					if(!flags)  sumObjs-=TOTALTYPE;
				}

				//Get number 获取所有点的个数。
				memset(cNum, 0, sizeof(cNum)); memset(sumGPoint, 0, sizeof(sumGPoint)); pIdx = 0;
				for(i=8;i<strlen(OneLine);i++)
				{
					if(OneLine[i]==0) break;
					cNum[i-8] = OneLine[i];
				}
				numpt = atoi(cNum); 
				
				//读取有效点。（中间的数字不为0）
				for(m=0;m<numpt;m++)
				{
					int fl0, fl1, j0, j1, j2;
					char cf0[32], cf1[32], cf2[32];
					float f0, f1, f2, threhold;

					memset(cf0, 0, sizeof(cf0));
					memset(cf1, 0, sizeof(cf1));
					memset(cf2, 0, sizeof(cf2));
					fgets(OneLine, 80, fi);
					fl0=fl1=j0=j1=j2=0; threhold = THOLD;
					for(i=0;i<strlen(OneLine);i++)
					{
						if((OneLine[i]==' ')&&(!fl0)&&(!fl1)) {fl0=1; continue;}
						if((OneLine[i]==' ')&&( fl0)&&(!fl1)) {fl1=1; continue;}
						if(!fl0) cf0[j0++] = OneLine[i];
						if((fl0)&&(!fl1)) cf1[j1++] = OneLine[i];
						if((fl1)) cf2[j2++] = OneLine[i];
					}
					f0 = atof(cf0); f1 = atof(cf1); f2 = atof(cf2);

					if(fabs(f0) < 0.002) f0 = 0.0; if(fabs(f1) < 0.002) f1 = 0.0; if(fabs(f2) < 0.002) f2 = 0.0; 
					
					//求出最小点和最大点
					if(f0 < minPx) minPx = f0; if(f0 > maxPx) maxPx = f0;
					if(f2 < minPy) minPy = f2; if(f2 > maxPy) maxPy = f2;
					
					//分门别类记录这个点的信息
					for(i=0;i<TOTALTYPE;i++)
					{
						if( (fabs((fabs(f1)- threhold))         < 0.0002) ||
							(fabs((fabs(f1)- threhold * 10.0))  < 0.0002) ||
							(fabs((fabs(f1)- threhold * 100.0)) < 0.0002) ||
							(fabs((fabs(f1)- threhold * 1000.0))< 0.0002) )
						{
							Points[sumObjs + i][sumGPoint[i]].px = f0; Points[sumObjs + i][sumGPoint[i]].py = f2; Points[sumObjs + i][sumGPoint[i]].iPrePos = m;
							sumGPoint[i]++; currPosi++;
						}
						threhold += THOLD;
					}					
				}// End for m		

				//记录每种类型点的有效点个数
				for(i=0;i<TOTALTYPE;i++)
					sumPtsObj[sumObjs + i] = sumGPoint[i];		//有效点个数
				
				//如果某类型点的个数不为0，则增加该类型指示
				for(i=0;i<TOTALTYPE;i++)
				{
					if(sumGPoint[i]) { typeIdxs[t].type = i+1; typeIdxs[t].idx = sumObjs + i; t++; 	allTypes[i] = 1;}
				}
				
				//重新排序
				pIdx = totalGPoint;
				for(i=0;i<TOTALTYPE;i++)
				{
					for(m=0;m<sumPtsObj[sumObjs + i]; m++)  Points[sumObjs + i][m].iCurPos = pIdx + m;
					pIdx += sumGPoint[i];
					totalGPoint += sumGPoint[i];
				}			
			}// End if 

			//读取所有由有效点组成的三角形。并且顺便填写纹理坐标
			if((OneLine[0]=='n') && (OneLine[1]=='u') && (OneLine[2]=='m') && (OneLine[3]=='s') && (OneLine[4]=='u') && (OneLine[5]=='r') && (OneLine[6]=='f'))
			{
				int numsurf;
				
				//Get number 获取所有面的个数。
				memset(cNum, 0, sizeof(cNum)); memset(sumGFace, 0, sizeof(sumGFace));  
				for(i=8;i<strlen(OneLine);i++)
				{
					if(OneLine[i]==0) break;
					cNum[i-8] = OneLine[i];
				}
				numsurf = atoi(cNum); j=0;

				//读取所有面，判断有效面。
				for(n=0;n<numsurf;n++)
				{
					int numrefs, sumGFlag; 
					int GFlag[TOTALTYPE][10];
					
					while(!((OneLine[0]=='r') && (OneLine[1]=='e') && (OneLine[2]=='f'))) {	fgets(OneLine, 80, fi);}	
					
					//Get number 获取所有点的个数。
					memset(cNum, 0, sizeof(cNum)); memset(Ref, 0, sizeof(Ref)); 
					for(i=5;i<strlen(OneLine);i++)
					{
						if(OneLine[i]==0) break;
						cNum[i-5] = OneLine[i];
					}
					numrefs = atoi(cNum); j=0;

					//读取面。
					for(m=0;m<numrefs;m++)
					{
						int fl0, fl1, j0, j1, j2;
						char cf0[32], cf1[32], cf2[32];
						int f0; 
						float f1, f2;
	
						memset(cf0, 0, sizeof(cf0));
						memset(cf1, 0, sizeof(cf1));
						memset(cf2, 0, sizeof(cf2));
						fgets(OneLine, 80, fi);
						fl0=fl1=j0=j1=j2=0;
						for(i=0;i<strlen(OneLine);i++)
						{
							if((OneLine[i]==' ')&&(!fl0)&&(!fl1)) {fl0=1; continue;}
							if((OneLine[i]==' ')&&( fl0)&&(!fl1)) {fl1=1; continue;}
							if(!fl0) cf0[j0++] = OneLine[i];
							if((fl0)&&(!fl1)) cf1[j1++] = OneLine[i];
							if((fl1)) cf2[j2++] = OneLine[i];
						}
						f0 = atoi(cf0); f1 = atof(cf1); f2 = atof(cf2);
						if(fabs(f1) < 0.000002) f1 = 0.0; if(fabs(f2) < 0.000002) f2 = 0.0; 
						
						Ref[j].idx = f0; Ref[j].tU = f1; Ref[j].tV = f2; j++;
					}// End for m
					
					//判断该面是否为有效面，（所有的点都是有效点）
					memset(GFlag, 0 ,sizeof(GFlag)); sumGFlag = 0;
					for(m=0;m<numrefs;m++)
					{
						for(i=0; i<TOTALTYPE; i++)
						{
							for(j=0; j<sumPtsObj[sumObjs + i]; j++)
							{
								if(Ref[m].idx == Points[sumObjs + i][j].iPrePos)		//如果该点是有效点
								{
									Points[sumObjs + i][j].tU = Ref[m].tU; Points[sumObjs + i][j].tV = Ref[m].tV;		//填充点信息
									Ref[m].iNewdx = Points[sumObjs + i][j].iCurPos;
									GFlag[i][m] = 1;
								}
							}
						}
					}
					
					for(j=0; j<TOTALTYPE; j++)
					{
						sumGFlag = 0;
						for(m=0;m<numrefs;m++) 
							sumGFlag += GFlag[j][m];
						if(sumGFlag == numrefs) 	//说明所有点都有效
						{
							int i;
							for(i=0;i<numrefs;i++)
								Faces[sumObjs + j][sumGFace[j]].idx[i] = Ref[i].iNewdx;
							sumGFace[j]++; 
							break;
						}
					}
				}//End for m	

				//每种类型各有多少个面				
				for(m=0; m<TOTALTYPE; m++)
					sumFaceObj[sumObjs + m] = sumGFace[m];

			}// End if 					
		}//End while
		fclose(fi);
		tiNums = t;

		//计算所有的有效点的经纬度
		absPx = fabs(fabs((float)maxPx) - fabs((float)minPx)) ; absPy = fabs(fabs((float)maxPy) - fabs((float)minPy)) ; offLon = fabs(findResultAC[idx].maxLong - findResultAC[idx].minLong); offLat = (findResultAC[idx].maxLati - findResultAC[idx].minLati);
		totalPnts = 0;

		/* 保存当前路径，返回到工作路径下，以便执行 GetElevation() */
		_getcwd(currDirectory, _MAX_PATH);	
		io = chdir(workDirectory);
		if(io) 
		{
			MessageToDialog("不能回到工作路径下。");
			return 0;
		}

		double  pre_elva = 0.0;
		for(m=0;m<sumObjs;m++)
		{
			Geodesy geod;
			Carton  cart;
			
			for(n=0;n<sumPtsObj[m];n++)
			{
				Points[m][n].lon = findResultAC[idx].minLong + (fabs(Points[m][n].px) / absPx) * offLon;
				Points[m][n].lat = findResultAC[idx].minLati + (fabs(Points[m][n].py) / absPy) * offLat;
				geod._elevation = GetElevation(Points[m][n].lon, Points[m][n].lat, pre_elva);
				pre_elva = geod._elevation;
				geod._lon = Points[m][n].lon; geod._lat = Points[m][n].lat; 
				GeodToCart(&geod, &cart);
				Points[m][n].px = cart.x; Points[m][n].py = cart.y; Points[m][n].pz = cart.z;
			}//End for n
			totalPnts += sumPtsObj[m];
		}//End for m
		RemoveAllTileNodes();

		//数据准备就绪，开始写btg
		//计算当前块的中心点经纬度。计算包围球半径
		{
			Bucket Bu;
			double lon, lat, deltaLon, deltaLat;
			lon = Points[typeIdxs[0].idx][2].lon;
			lat = Points[typeIdxs[0].idx][2].lat;
			deltaLon = fabs(lon - 0.0);
			deltaLat = fabs(lat - 0.0);
			if((deltaLon < 0.01)&&(deltaLat < 0.01))
			{
				lon = Points[typeIdxs[1].idx][2].lon;
				lat = Points[typeIdxs[1].idx][2].lat;
			}
			
			set_bucket(&Bu, lon, lat);
			geoCenter._lon = get_center_lon(&Bu); 
			geoCenter._lat = get_center_lat(&Bu);
			geoCenter._elevation  = 0; 
			GeodToCart(&geoCenter, &catCenter);		
			gbs_radius = calc_gbs(gPoint, 4);
			gen_base_path(&Bu, basepath);
			iTile = gen_index(&Bu);

			currTi.minLongi = (float)((int)geoCenter._lon);
			currTi.minLati  = (float)((int)geoCenter._lat);
			currTi.maxLongi = currTi.minLongi + 1.0; 
			currTi.maxLati  = currTi.minLati  + 1.0;
		}

		/* 从工作路径返回到当前路径，以便后续正常运行 */
		io = chdir(currDirectory);
		if(io) 
		{
			MessageToDialog("不能回到工作路径下。");
			return 0;
		}

		//打开文件。写文件。保存文件名和路径名到中间结果处，为下一步做准备。
		findResultAC[idx].destTileID = 0;
		tile = findUniqueTile(iTile);
		findResultAC[idx].destTileID = tile;
		sprintf(dstPath, "%s\\%s", iveroot, basepath);
		sprintf(dstFile, "%d.btg.gz", tile);
		strcpy(findResultAC[idx].destPath, dstPath);
		strcpy(findResultAC[idx].destFile, dstFile);
		findResultAC[idx].destLongi = (int)geoCenter._lon;
		findResultAC[idx].destLati  = (int)geoCenter._lat;
		findResultAC[idx].isOK		= 1;

		//需要的数据准备完毕，开始转换.
		MoveDataToDataStru();
		/* 把这一行添加到stg文件结尾 */
		/* 集中写STG */
		if(!CreateCultureBtgFromSrc(idx))
		{
			MessageToDialog("创建STG文件失败\n");
			return 0;
		}

		//查找下一个符合条件的文件
		fio = _findnext(ffio, &ff);
	}
	_findclose(ffio);
	
	
	/* 返回原有工作目录 */
	io = chdir(path);
	if(io) 
	{
		MessageToDialog("不能返回原工作目录。");
		return 0;
	}
	
#endif	
	return 1;
} 



/* 处理全部ac文件 */
int DoAllAc()
{
	int i, j, k;
	char msg[STRING_STANDARD_LEN];
	int currDoing[STRING_STANDARD_LEN];
	int numCD;
	PerZone mPZ[STRING_STANDARD_LEN];
	int numPZ, iFlag;

	/* 找出还没有制作文化信息的所有ac文件所在的目录名称 */
	if(FindOutDirsNoCultDone() == 0)
	{
		MessageToDialog("没有需要转换的文化信息.");
		return 1;	
	}

	/* 先从dirsNoCultDone搜起，即这些目录里面都有没做文化信息处理的ac文件 */
	for(i=0; i<numDNCD; i++)
	{
		/* 针对dirsNoCultDone里面每一项，把totalcult里面所有同目录名的项集合到findResultAC*/
		memset(currDoing, 0, sizeof(currDoing));  numCD = 0;
		memset(mPZ, 0, sizeof(mPZ));  numPZ = 0;
		if(FindAllAcItemsWithSameDir(i) == 0) continue;
		strcpy(findJpgDirectory, findResultAC[0].pathName);
		
		/* 针对findResultAC里面每一项，先检查是否做过文化信息处理，没做过就在currDoing里记录。*/
		for(j=0; j<numFRAC; j++)
		{
			if(findResultAC[j].isOK == 0)
			{	
				double lon, lat;
				lon = (findResultAC[j].minLong + findResultAC[j].maxLong ) / 2.0;
				lat = (findResultAC[j].minLati + findResultAC[j].maxLati ) / 2.0;
				findResultAC[j].destLongi = (int)lon;
				findResultAC[j].destLati  = (int)lat;
				currDoing[numCD++] = j;
			}
		}
		
		/* 从刚刚做好的所有块里面找出所有不同的块 */
		for(j=0; j<numCD; j++)
		{
			iFlag = 0;
			for(k=0; k<numPZ; k++)
			{
				if((findResultAC[currDoing[j]].destLongi == mPZ[k].lon) && (findResultAC[currDoing[j]].destLati == mPZ[k].lat))
				{ iFlag = 1; break; }
			}
			if(!iFlag)
			{	
				mPZ[numPZ].lon = findResultAC[currDoing[j]].destLongi;
				mPZ[numPZ].lat = findResultAC[currDoing[j]].destLati; numPZ++;
			}
		}
	
		/* 对于每个块，先做基础地景，再加文化信息 */
		for(j=0; j<numPZ; j++)
		{
			/* 先做初始化全局变量 */
			currTi.minLongi = (float)(mPZ[j].lon);
			currTi.minLati  = (float)(mPZ[j].lat);
			currTi.maxLongi = (float)(mPZ[j].lon + 1);
			currTi.maxLati  = (float)(mPZ[j].lat + 1);
			currPZ = &mPZ[j];

			/* 根据currDoing里面的信息，转换文化信息btg */
			for(k=0; k<numCD; k++)
			{
				if(	(findResultAC[currDoing[k]].isOK	   == 0			) && 
					(findResultAC[currDoing[k]].destLongi  == mPZ[j].lon) && 
					(findResultAC[currDoing[k]].destLati   == mPZ[j].lat) )
				{
					sprintf(msg, "正在转换BTG--%s...", findResultAC[currDoing[k]].fileName);
					MessageToDialog(msg);
					if(!ReadACandMakeBTG(currDoing[k]))
					{
						sprintf(msg, "转换BTG--%s失败。", findResultAC[currDoing[k]].fileName);
						MessageToDialog(msg);
						continue;
					}
				}
			}   
			    
			/* 将相同块内的所有btg组成一个Geode，插入到挖洞列表当中，待集中统一挖洞。 */
			if(!InsertItemIntoCutList(mPZ[j].lon, mPZ[j].lat))
			{   
				sprintf(msg, "转换IVE失败。");
				MessageToDialog(msg);
				return 0;
			}   
		}       

		/* 做好的文化信息块更新总体数据库 */
		UpdateTotalcultWithNewFindResult();
	}
	

	return 1;	
}


/* 处理全部Shp文件 */
int DoAllShp()
{
	int i, j, k;
	char msg[STRING_STANDARD_LEN];
	int currDoing[STRING_STANDARD_LEN];
	int numCD;
	PerZone mPZ[STRING_STANDARD_LEN];
	int numPZ, iFlag;

	/* 找出还没有制作文化信息的所有shp文件所在的目录名称 */
	if(FindOutDirsNoShpDone() == 0)
	{
		MessageToDialog("没有需要转换的文化信息.");
		return 1;	
	}

	/* 先从dirsNoCultDone搜起，即这些目录里面都有没做文化信息处理的ac文件 */
	for(i=0; i<numDNCD; i++)
	{
		/* 针对dirsNoCultDone里面每一项，把totalcult里面所有同目录名的项集合到findResultAC*/
		memset(currDoing, 0, sizeof(currDoing));  numCD = 0;
		memset(mPZ, 0, sizeof(mPZ));  numPZ = 0;
		if(FindAllShpItemsWithSameDir(i) == 0) continue;
		strcpy(findJpgDirectory, findResultAC[0].pathName);
		
		/* 针对findResultAC里面每一项，先检查是否做过文化信息处理，没做过就在currDoing里记录。*/
		for(j=0; j<numFRAC; j++)
		{
			if(findResultAC[j].isOK == 0)
			{	
				double lon, lat;
				lon = (findResultAC[j].minLong + findResultAC[j].maxLong ) / 2.0;
				lat = (findResultAC[j].minLati + findResultAC[j].maxLati ) / 2.0;
				findResultAC[j].destLongi = (int)lon;
				findResultAC[j].destLati  = (int)lat;
				currDoing[numCD++] = j;
			}
		}
		
		/* 从刚刚做好的所有块里面找出所有不同的块 */
		for(j=0; j<numCD; j++)
		{
			iFlag = 0;
			for(k=0; k<numPZ; k++)
			{
				if((findResultAC[currDoing[j]].destLongi == mPZ[k].lon) && (findResultAC[currDoing[j]].destLati == mPZ[k].lat))
				{ iFlag = 1; break; }
			}
			if(!iFlag)
			{	
				mPZ[numPZ].lon = findResultAC[currDoing[j]].destLongi;
				mPZ[numPZ].lat = findResultAC[currDoing[j]].destLati; numPZ++;
			}
		}
	
		/* 对于每个块，先做基础地景，再加文化信息 */
		for(j=0; j<numPZ; j++)
		{
			/* 先做初始化全局变量 */
			currTi.minLongi = (float)(mPZ[j].lon);
			currTi.minLati  = (float)(mPZ[j].lat);
			currTi.maxLongi = (float)(mPZ[j].lon + 1);
			currTi.maxLati  = (float)(mPZ[j].lat + 1);
			currPZ = &mPZ[j];

			/* 先把当前目录下的湖泊组清空。*/
			allLakeInCurrLonLat->removeChildren(0,allLakeInCurrLonLat->getNumChildren());
			 
			/* 根据currDoing里面的信息，转换文化信息btg */
			for(k=0; k<numCD; k++)
			{
				if(	(findResultAC[currDoing[k]].isOK	   == 0			) && 
					(findResultAC[currDoing[k]].destLongi  == mPZ[j].lon) && 
					(findResultAC[currDoing[k]].destLati   == mPZ[j].lat) )
				{
					sprintf(msg, "正在转换BTG--%s...", findResultAC[currDoing[k]].fileName);
					MessageToDialog(msg);
					if(!ReadShpAndMakeBTG(currDoing[k]))
					{
						sprintf(msg, "转换BTG--%s失败。", findResultAC[currDoing[k]].fileName);
						MessageToDialog(msg);
						continue;
					}
				}
			}   
			    
			/* 将相同块内的所有btg组成一个Geode，插入到挖洞列表当中，待集中统一挖洞。 */
			if(!InsertItemIntoCutList(mPZ[j].lon, mPZ[j].lat))
			{   
				sprintf(msg, "转换IVE失败。");
				MessageToDialog(msg);
				return 0;
			}   
		}       

		/* 做好的文化信息块更新总体数据库 */
		UpdatetotalShpWithNewFindResult();
	}

	return 1;	
}


/* 读取Shp文件并分析，然后创建btg文件 */
int ReadShpAndMakeBTG(int idx)
{
	int io, fio, ffio;
	char dstFile[STRING_STANDARD_LEN];
	char dstPath[STRING_STANDARD_LEN];
	char dPath[STRING_STANDARD_LEN];
	char path[_MAX_PATH];
	char basepath[STRING_STANDARD_LEN];
	long tile, iTile;
	int numSameTile = 0;
	struct _finddata_t ff;

	/* 每次做文化信息之前，清除高程遗留数据，以解决路面被沉到地底下的问题 --- 2012.2.15 */
	RemoveAllTileNodes();

	//获取当前目录
	_getcwd(path, _MAX_PATH);

	//指定输入输出目录
	sprintf(dPath, "%s", findResultAC[idx].pathName);
	io = chdir(dPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", dPath);
		MessageToDialog(msg);
		return 0;
	}

	ffio = _findfirst(findResultAC[idx].fileName, &ff);	
	if((!ffio)||(ffio == -1))
	{
		MessageToDialog("找不到SHP文件。");
		return 0;
	}
	fio = 0;

	while(!fio)
	{
		/* 根据相关信息填充findResultAC其他空白项 */
		Bucket Bu;
		double lon, lat;
		lon = (findResultAC[idx].maxLong + findResultAC[idx].minLong) / 2.0;
		lat = (findResultAC[idx].maxLati + findResultAC[idx].minLati) / 2.0;
		set_bucket(&Bu, lon, lat);
		geoCenter._lon = get_center_lon(&Bu); 
		geoCenter._lat = get_center_lat(&Bu);
		geoCenter._elevation  = 0; 
		gen_base_path(&Bu, basepath);
		iTile = gen_index(&Bu);

		/* 创建并保存文件名和路径名到中间结果处，为下一步做准备。*/
		findResultAC[idx].destTileID = 0;
		tile = findUniqueTile(iTile);
		findResultAC[idx].destTileID = tile;
		sprintf(dstPath, "%s\\%s", iveroot, basepath);
		sprintf(dstFile, "%d.btg.gz", tile);
		strcpy(findResultAC[idx].destPath, dstPath);
		strcpy(findResultAC[idx].destFile, dstFile);
		findResultAC[idx].destLongi = (int)geoCenter._lon;
		findResultAC[idx].destLati  = (int)geoCenter._lat;
		
		/* 读取Shp文件并且把它存入DataStruct */
		currShpType = -1;
		unsigned int i;
		int FoundFlag = 0;
		int ShpResult = LoadShpAndDrawIt(ff.name);
		if(ShpResult > 0)
		{
			/* 查找当前Shp文件是否前面已经用过。*/
			string shpname = ff.name;
			for(i=0; i<allShps.size();i++)
			{
				if(shpname == allShps[i].fileName)
				{
					FoundFlag = 1; break;
				}
			}
			if(!FoundFlag)
			{
				ShpManagement smg;
				smg.fileName = shpname;	smg.longi = lon;	smg.lati = lat;	smg.type = currShpType; smg.btgID = tile;
				allShps.push_back(smg);
			}
			
			if(ShpResult == 1)
			{
				/* 将DataStruct内容写入Btg */
				/* 把这一行添加到stg文件结尾 */
				/* 集中写STG */
				if(!CreateCultureBtgFromSrc(idx))
				{
					MessageToDialog("创建STG文件失败\n");
					return 0;
				}
			}
		}

		/* 如果当前shp文件已经被使用过，那么这个文件生成的点灯和点树一定已经被生成过。就不需要再被重复生成了。*/
		if(FoundFlag)	goto RSAMB_Continue_1;
		
		/* 如果CultPoints不空，表示有灯光还没写，再写一个灯光btg文件。 */
		if(CultPoints.size() > 0)
		{
			/* 先调整灯光显示的远近距离效果 */
			if(!AddLODIntoLightPoints())	return 0;
			
			/* 构造灯光的数据结构DataStru */
			if(CreateLightsDataStructFromCultPoints())
			{
				/* 构造文件名，*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\%dL.btg.gz", findResultAC[idx].destPath, tile); //命名为 tileID+L.btg.gz
				sprintf(fn, "%dL.btg.gz", tile);
			
				/* 写btg文件，并把tileL.btg.gz写入各个stg中 */
				if(WriteDataStructToBTG(fullPathFile))
				{
					if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "打开文件失败 --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}

		/* 如果TreePoints不空，表示有森林还没写，再写一个森林btg文件。 */
		if(TreePoints.size() > 0)
		{
			/* 先调整树显示的远近距离效果 */
			if(!AddLODIntoTreePoints())	return 0;
			
			/* 构造树木的数据结构DataStru */
			if(CreateTreesDataStructFromTreePoints())
			{
				/* 构造文件名，*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\%dT.btg.gz", findResultAC[idx].destPath, tile); //命名为 tileID+T.btg.gz
				sprintf(fn, "%dT.btg.gz", tile);
			
				/* 写btg文件，并把tileT.btg.gz写入各个stg中 */
				if(WriteDataStructToBTG(fullPathFile))
				{
					if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "打开文件失败 --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}
		
		/* 如果是湖，把湖单独做成一个Btg文件，命名为 L+tileID.btg */
		if(allLakeSkirt->getNumChildren() > 0)
		{
			/* 检测湖泊是否已经被制作成btg文件 */
			osg::ref_ptr<osg::Group> currLakes = new osg::Group;
			for(unsigned int i=0; i<allLakeSkirt->getNumChildren(); ++i)
			{
				osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode *>(allLakeSkirt->getChild(i));
				osg::ref_ptr<osg::Geometry> gm = dynamic_cast<osg::Geometry *>(gnode->getDrawable(0));
				osg::ref_ptr<osg::Vec3Array> vx = dynamic_cast<osg::Vec3Array *>(gm->getVertexArray());
				int curVSum = vx->size();
				
				int iFlag = 0;
				for(unsigned int j=0; j<allLakeInCurrLonLat->getNumChildren(); ++j)
				{
					osg::ref_ptr<osg::Geode> gnode2 = dynamic_cast<osg::Geode *>(allLakeInCurrLonLat->getChild(j));
					osg::ref_ptr<osg::Geometry> gm2 = dynamic_cast<osg::Geometry *>(gnode2->getDrawable(0));
					osg::ref_ptr<osg::Vec3Array> vx2 = dynamic_cast<osg::Vec3Array *>(gm2->getVertexArray());
					int curVSum2 = vx2->size();
					
					if(curVSum2 == curVSum)
					{
						iFlag = 1; break;
					}
				}
				
				if(!iFlag)
				{
					currLakes->addChild(gnode);
					allLakeInCurrLonLat->addChild(gnode);
				}
			}

			if(currLakes->getNumChildren() > 0)
			{
				/* 构造湖泊的数据结构DataStru */
				if(CreateLakeDataStruAndBtg(currLakes))
				{
					/* 构造文件名，*/
					char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
					sprintf(fullPathFile, "%s\\L%d.btg.gz", findResultAC[idx].destPath, tile); //命名为 L+tileID.btg.gz
					sprintf(fn, "L%d.btg.gz", tile);
				
					/* 写btg文件，并把L+tile.btg.gz写入各个stg中 */
					if(WriteDataStructToBTG(fullPathFile))
					{
						/* 湖的btg文件，不写入stg里面 */
						//if(!WriteAllStgs(fn, 1)) return 0;
					}else{
						sprintf(msg, "打开文件失败 --- %s", fullPathFile);
						MessageToDialog(msg);
						return 0;
					}
				}

				/* Add by Galen 2013-02-20 */
				/* 根据currLakes组创建新型的湖面，并写入STG索引文件 */
				/* 新型湖面文件的命名：O+tileID.btg.gz */
				if(CreateNewStyleLakeDataStruAndBtg(currLakes))
				{
					/* 构造文件名，*/
					char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
					sprintf(fullPathFile, "%s\\O%d.btg.gz", findResultAC[idx].destPath, tile); //命名为 O+tileID.btg.gz
					sprintf(fn, "O%d.btg.gz", tile);
				
					/* 写btg文件，并把O+tile.btg.gz写入各个stg中 */
					if(WriteDataStructToBTG(fullPathFile))
					{
						/* 新型湖面的btg文件，写入stg里面 */
						if(!WriteAllStgs(fn, 1)) return 0;
					}else{
						sprintf(msg, "打开文件失败 --- %s", fullPathFile);
						MessageToDialog(msg);
						return 0;
					}
				}
				/* Add by Galen 2013-02-20 End */

			}
		}

		/* 如果是动态草，把动态草单独做成一个Btg文件，命名为 D+tileID.btg */
		if(allDynamicGrass->getNumChildren() > 0)
		{
			/* 构造动态草的数据结构DataStru */
			if(CreateLakeDataStruAndBtg(allDynamicGrass))
			{
				/* 构造文件名，*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\D%d.btg.gz", findResultAC[idx].destPath, tile); //命名为 D+tileID.btg.gz
				sprintf(fn, "D%d.btg.gz", tile);
			
				/* 写btg文件，并把D+tile.btg.gz写入各个stg中 */
				if(WriteDataStructToBTG(fullPathFile))
				{
					/* 动态草的btg文件，不写入stg里面 */
					//if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "打开文件失败 --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}

		/* 如果是动态草，把动态草单独做成一个Btg文件，命名为 D+tileID.btg */
		if(allIslands->getNumChildren() > 0)
		{
			/* 构造岛的数据结构DataStru */
			if(CreateLakeDataStruAndBtg(allIslands))
			{
				/* 构造文件名，*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\I%d.btg.gz", findResultAC[idx].destPath, tile); //命名为 D+tileID.btg.gz
				sprintf(fn, "I%d.btg.gz", tile);
			
				/* 写btg文件，并把I+tile.btg.gz写入各个stg中 */
				if(WriteDataStructToBTG(fullPathFile))
				{
					/* 动态草的btg文件，不写入stg里面 */
					//if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "打开文件失败 --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}

		/* 如果是路，把路单独做成一个Btg文件，不参与切割运算，命名为 R+tileID.btg */
		if(allRoadSkirt->getNumChildren() > 0)
		{
			/* 构造道路的数据结构DataStru */
			if(CreateRoadDataStruAndBtg(allRoadSkirt))
			{
				/* 构造文件名，*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\R%d.btg.gz", findResultAC[idx].destPath, tile); //命名为 tileID+S.btg.gz
				sprintf(fn, "R%d.btg.gz", tile);
			
				/* 写btg文件，并把R+tile.btg.gz写入各个stg中 */
				if(WriteDataStructToBTG(fullPathFile))
				{
					/* 路的btg文件，写入stg里面 */
					if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "打开文件失败 --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}

RSAMB_Continue_1:

		if (ifHasIslandInCurrentTerrain) findResultAC[idx].isOK		= 2;
		else	findResultAC[idx].isOK		= 1;

		/* 查找下一个符合条件的文件 */
		fio = _findnext(ffio, &ff);
	}
	_findclose(ffio);
	
	/* 返回原有工作目录 */
	io = chdir(path);
	if(io) 
	{
		MessageToDialog("不能返回原工作目录。");
		return 0;
	}
	
	return 1; 
}

/* 恢复一个地景文件之前，先要进入它所在的目录 */
int GoIntoDirectoryBeforeRestoreATerrainFile(int lon, int lat)
{
	int io;
	char fullPath[_MAX_PATH];

	/* 找到当前目录 */
	Bucket Bt;
	char basePath[STRING_STANDARD_LEN];
	set_bucket(&Bt, lon, lat);
	gen_base_path(&Bt, basePath);
	sprintf(fullPath, "%s\\%s", iveroot, basePath);
	
	/* 读取一个ive文件，分析它的长度。 */
	io = chdir(fullPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", fullPath);
		MessageToDialog(msg);
		return 0;
	}
	return 1;
}

/* 制作文化信息之前，要先恢复在allHTI里面列出的所有地景文件 */
int RestoreTerrainBeforeDoCulture()
{	
	int i, io;
	char ivefname[STRING_STANDARD_LEN], msg[STRING_STANDARD_LEN * 2];
	char path[_MAX_PATH];

	MessageToDialog("制作文化信息之前，开始恢复原始地景文件。");
		
	/* 获取当前目录 */
	_getcwd(path, _MAX_PATH);

	/* 初始化*/
	allNeedCut.clear();
	
	/* 查找已经做好的地景部分，找到文件名，根据文件名来备份恢复。*/
	for(i=0; i<numHTI; i++)
	{
		if(allHTI[i].isOK == 1)
		{
			/* 从tileID里面找出经度和纬度*/
			int tmp = allHTI[i].tileID / 10000;
			int lon = tmp / 100;
			int lat = tmp - lon * 100;
			InsertItemIntoCutList(lon, lat);

			/* 首先删除这个目录里面先前做过的所有文化信息 */
			DeleteAllPrevCultureFiles(lon, lat);

			/* 先要进入它的目录*/
			GoIntoDirectoryBeforeRestoreATerrainFile(lon, lat);

			///* 检查文件长度有没有被改变，没有被改变就不需要恢复。*/
			if(CheckFileLengthBeforeRestore(allHTI[i].tileID))
			{
				/* 先恢复用户高分辨率的地景 */
				sprintf(ivefname, "%d.ive", allHTI[i].tileID);
				int ret = RestoreFileFromTemp(ivefname);
				if( ret != 0 ) 
				{
					sprintf(msg, "正在恢复文件组 --- %s", ivefname);
					MessageToDialog(msg);
					if( ret == 2)	BackUpFileToTemp(ivefname, ""); 	//如果恢复目录没有这个文件，就先备份一下。适用于地景顶层文件长度不为701K的情况。
				}
			}
		}
	}
	
	/* 再恢复低分辨率的地景部分 */
	for(std::vector<CutList>::iterator p = allNeedCut.begin(); p != allNeedCut.end(); p++)
	{
		/* 先要进入它的目录*/
		GoIntoDirectoryBeforeRestoreATerrainFile(p->lon, p->lat);
		
		for(i=0; i<4; i++)
		{
			sprintf(ivefname, "%d%d_%d.ive", p->lon, p->lat, i);
			int ret = RestoreFileFromTemp(ivefname);
			if( ret != 0 ) 
			{
				sprintf(msg, "正在恢复文件组 --- %s", ivefname);
				MessageToDialog(msg);
				if( ret == 2)	BackUpFileToTemp(ivefname, "");	//如果恢复目录没有这个文件，就先备份一下。适用于地景顶层文件长度不为701K的情况。
			}
		}
	}

	/* 进入临时路径。*/
	io = chdir(tempDirectory);
	if(!io)
	{
		/* 删除部分临时文件 */
		char fullPathFile[_MAX_PATH];
		sprintf(fullPathFile, "%s\\AptLight.ive", tempDirectory);  
		/* 再删除掉这个文件， */
		DeleteFile(fullPathFile);
	}


	/* 返回原有工作目录 */
	io = chdir(path);
	if(io) 
	{
		MessageToDialog("不能返回原工作目录。");
		return 0;
	}
	return 1;
}


/**************************************************
 *
 *   开始处理文化信息
 *
 **************************************************/

int BeginCreateCulture()
{
	MessageTo = 1;	ifHasIslandInCurrentTerrain = 0;	fRtoCutIsland.clear();

	/* 将路径更新一下 */
	UpdateDialog01();
	
	/* 设置全局变量 iveroot和 hgtroot */
	SetGlobalVars();
	
	/* 将当前用户配置信息写入Log文件 */
	WriteCultureCfgIntoLogFile();

	/* 制作文化信息之前，要先恢复在allHTI里面列出的所有地景文件 */
	RestoreTerrainBeforeDoCulture();

	/* 全局数组初始化 */
	allNeedCut.clear();	allShps.clear();

	/* 清空平滑过渡类列表 */
	csmList.clear();
	
	/* 清空细节纹理收集列表 */
	allDetailedTxt.clear();

	/* 第一部分：开始制作文化信息，目前支持ac文件和shp文件两种格式。 */
	if(1)
	{
		MessageToDialog("开始创建文化信息。");

		/* 遍历路径找到所有.ac文件 */
		TravelDirectoryToReadAc(baseTextDir);
	
		/* 遍历路径找到所有.shp文件 */
		TravelDirectoryToReadShp(baseTextDir);
	
		/* 根据经纬度找到文化信息所在区域, 如果经纬度都是0，那就制作所有的。*/
		MessageToDialog("处理文化信息AC文件。");
		if(!DoAllAc())
		{
			MessageToDialog("AC转换BTG失败。");
			goto BCC_Continue_1;
		}
	
		MessageToDialog("处理文化信息SHP文件。");
		if(!DoAllShp())
		{
			MessageToDialog("SHP转换BTG失败。");
			goto BCC_Continue_1;
		}
	}
BCC_Continue_1:		

	/* 第二部分：根据参数生成通用机场，然后计算出通用机场的轮廓。本部分为可选项。 */
	if(CurrentCreateGeneralAirport)
	{
		MessageToDialog("开始创建通用机场。");
		if(!DoAllGA())
		{
			MessageToDialog("创建通用机场失败。");
			goto BCC_Continue_2;
		}
	}
BCC_Continue_2:	

	/* 第三部分：根据用户定位后的flt机场，计算出轮廓。 */
	if(1)
	{
		char userPath[_MAX_PATH];
		int retVal;
		
		MessageToDialog("开始定位用户机场。");
	
		/* 直接查找用户定义三维物体目录，以查找用户三维物体数据 */
		sprintf(userPath, "%s\\data\\Models\\User", fsData);
	
		/* 遍历三维物体目录，查找用户定义三维物体 */
		retVal = TravelDirectoryToReadFlt(userPath);
		if(!retVal)
		{
			MessageToDialog("遍历三维物体目录读取三维物体信息失败\n");
			goto BCC_Continue_3;
		}
		
		/* 定位用户定义三维物体。*/
		if(!DoAllFlt())
		{
			MessageToDialog("定位用户定义三维物体失败。");
			return 0;
			//goto BCC_Continue_3;
		}
	}
BCC_Continue_3:
	/* 将做好的信息更新进数据库 */
	CheckConfig(1);

	/* 第四部分：根据前三部分的计算结果，统一处理地景，在地景上切割出文化信息、通用机场、用户机场所在的面。*/
	if(!DoTerrainCut())
	{
		MessageToDialog("地景挖洞失败");
		return 0;
	}

	/* 如果当前文化信息有岛，去掉它。*/
	if(ifHasIslandInCurrentTerrain)
	{
		if(!DeleteIslandFromCultureBTG()) return 0;
	}

	MessageToDialog("切割地景完成。");
	return 1;	
}


/*============================================================================================
 *  下面是通用机场的制作过程。
 *  基本方法是通过参数来选择模板，然后对模板根据实际坐标进行平移，得到机场btg文件。
 *============================================================================================ */

/**************************************************
 *
 *   开始处理通用机场
 *
 **************************************************/
 
/* 读取模板BTG，转换坐标后，写目标BTG */
int XchgModelBtgToDestination()
{
	int itmp;
	char templatefile[STRING_STANDARD_LEN];
	char fullpath[STRING_STANDARD_LEN], srcname[STRING_STANDARD_LEN], dstname[STRING_STANDARD_LEN];
	Geodesy cur;
	Bucket Bt;
	char basepath[32];

	sprintf(fullpath, "%s\\generalapt\\Gadata", workDirectory);
	
	/* 先获取模板文件名 */
	sprintf(templatefile, "GA-");
	itmp = currGA.shRunwayID;
	sprintf(templatefile, "%s%02d", templatefile, itmp);		//shRunwayID
	itmp = 30 + currGA.shRunwayWidth*15;
	sprintf(templatefile, "%s%02d-", templatefile, itmp);		//shRunwayWidth
	itmp = 1000 + currGA.shRunwayLength*500;
	sprintf(templatefile, "%s%04d", templatefile, itmp);		//shRunwayLength
	if(!currGA.shBuildingPosition)
		sprintf(templatefile, "%s-L", templatefile);			//shBuildingPosition
	sprintf(templatefile, "%s.btg.gz", templatefile);
	sprintf(srcname, "%s\\%s", fullpath, templatefile);
	
	/* 读取模板btg文件 */
	ReadBTGToDataStruct(srcname);	
	
	/* 转换坐标 */
	cur._lon = currGA.dLongtitude; cur._lat = currGA.dLatitude; cur._elevation = currGA.fAltitude;
	ChangeVertexInDataStru(&cur);
	
	/* 写目标Btg文件*/
	/* 先找出目标路径。*/
	set_bucket(&Bt, currGA.dLongtitude, currGA.dLatitude);
	gen_base_path(&Bt, basepath);
	sprintf(dstname, "%s\\Scenery\\Terrain\\%s\\%s.btg.gz", currGA.cLocation, basepath, currGA.cICAO);
	WriteDataStructToBTG(dstname);	  
	
	return 1;
}

/*  分析裙边点并输出成一个包围线的osg::Geode，并在地形上挖出洞。 */
int ReadAptOSGAndOutputLineLoopIve()
{
	char msg[STRING_STANDARD_LEN];
	int ilon, ilat;
	Geodesy cur;
	Bucket Bt;
	char basepath[32];
	char dstname[STRING_STANDARD_LEN];

	/* 首先找出裙边点集合 */
	ilon = (int)currGA.dLongtitude; ilat = (int)currGA.dLatitude;
	cur._lon = currGA.dLongtitude; cur._lat = currGA.dLatitude; cur._elevation = currGA.fAltitude;
	FindSkirtVertex(cur._elevation);
	
	/* 从DataStru去掉裙边的面 */
	RemoveSkirtsFromDataStru();

	/* 根据裙边来计算包围线 */
	osg::ref_ptr<osg::Geode> gd = GetAllLineLoopFromSkirtEdge();
		
	/* 写为ive文件*/
	/* 先找出目标路径。*/
	set_bucket(&Bt, currGA.dLongtitude, currGA.dLatitude);
	gen_base_path(&Bt, basepath);
	sprintf(dstname, "%s\\Scenery\\Terrain\\%s\\%s.ive", currGA.cLocation, basepath, currGA.cICAO);
	osgDB::writeNodeFile(*gd, dstname);

	/* 插入到挖洞列表当中，待集中统一挖洞。*/
	if(!InsertItemIntoCutList(ilon, ilat))
	{   
		sprintf(msg, "给机场 %s 挖洞失败。", currGA.cICAO);
		MessageToDialog(msg);
		return 0;
	}   
	return 1;
}

/* 下面 MAKEMODELBTG 段，是模板的生成过程，通过读取.dat文件来生成相应的.btg文件。
   由于模板已经准备完毕，这段过程就不用了。 */

//#define MAKEMODELBTG
#ifdef MAKEMODELBTG
int CreateModelHgtdata(char *workDirectory);
int DoGenaptModel(char *fn)
{
	char arg[30][STRING_STANDARD_LEN];
	char *argptr[30];
	int  io;
	char path[_MAX_PATH];

	//获取当前目录
	_getcwd(path, _MAX_PATH);
	
	//计算相关路径信息
	memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr)); 

	//str.Format("@echo off\ngenapts --input=data/airports/ZBAA_UseForVpb.dat"
	//构造命令行
	sprintf(arg[ 0], "geneapt.exe"); 												argptr[ 0] = (char *)&arg[ 0];
	sprintf(arg[ 1], "--input=generalapt\\Ta\\%s"										, fn	); 	argptr[ 1] = (char *)&arg[ 1];
	sprintf(arg[ 2], "--work=generalapt\\work"	);	argptr[ 2] = (char *)&arg[ 2];
	int argIdx = 3;

	/* 先设置当前目录为工作目录 */
	io = chdir(workDirectory);
	if(io)  return 0;
	
	/* 先生成高程数据 */
	CreateModelHgtdata(workDirectory);
	
	/* 调用库里面的过程生成通用机场文件 */
	if(!CreateAirport(argIdx, argptr))
	{
		MessageToDialog("创建通用机场失败。");
		return 0;
	}

	/* 回归原来目录 */
	io = chdir(path);
	if(io)  return 0;
	
	return 1;
}

/* 模板准备工作，读取.dat，转换成.btg */
int ReadDataAndWriteToBTG()
{
	char fullpath[STRING_STANDARD_LEN], fullname[STRING_STANDARD_LEN], msg[STRING_STANDARD_LEN];
	int io, dio, dfio;
	struct _finddata_t sdff;

	sprintf(fullpath, "%s\\generalapt\\Ta", workDirectory);
	
	/* 进入这个目录，然后搜索同tile文件名的所有文件。 */
	io = chdir(fullpath);
	if(io) return 0;
	
	/* 查找当前路径下的所有文件； */
	dio = _findfirst("*.dat", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		/* 打开文件 */
		sprintf(fullname, "%s", sdff.name);

		sprintf(msg, "创建通用机场 --- %s", sdff.name);
		MessageToDialog(msg);
		
		/* */
		DoGenaptModel(fullname);
		
		/* 寻找下一个 */
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);

	return 1;
} 

/* 读取模板btg文件，并且校验。 */
int ReadModelBTGAndVerify()
{
	char fullpath[STRING_STANDARD_LEN], destpath[STRING_STANDARD_LEN], fullname[STRING_STANDARD_LEN], msg[STRING_STANDARD_LEN];
	int io, dio, dfio;
	struct _finddata_t sdff;

	sprintf(fullpath, "%s\\generalapt\\Gold", workDirectory);
	sprintf(destpath, "%s\\generalapt\\Gadata", workDirectory);
	
	/* 进入这个目录，然后搜索同tile文件名的所有文件。 */
	io = chdir(fullpath);
	if(io) return 0;
	
	/* 查找当前路径下的所有文件； */
	dio = _findfirst("*.btg.gz", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		/* 打开文件 */
		sprintf(msg, "校验通用机场 --- %s", sdff.name);
		MessageToDialog(msg);
		
		/* */
		sprintf(fullname, "%s\\%s", destpath, sdff.name);
		VerifyTheModels(sdff.name, fullname);

		/* 寻找下一个 */
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);

	return 1;
} 
#endif //MAKEMODELBTG


/* 创建通用机场主过程 */
int DoAllGA()
{
	PerZone mPZ;

	/* 先设置全局变量，一些路径。*/
	strcpy(currGA.cLocation, currTi.iveDir);
	mPZ.lon = (int)currGA.dLongtitude;  mPZ.lat = (int)currGA.dLatitude;
	currTi.minLongi = (float)(mPZ.lon);
	currTi.minLati  = (float)(mPZ.lat);
	currTi.maxLongi = (float)(mPZ.lon + 1);
	currTi.maxLati  = (float)(mPZ.lat + 1);
	currPZ = &mPZ;

	if(rebuildTerrainForGeneapt)
	{
		/* 首先生成基础地形 */
		if(!CreatePerDegreeTerrain(0))
		{
			return 0;
		}
	}

	/* 先判断机场海拔高度，如果高度为0.0，那么取该经纬度点的真实高度。*/	
	if(fabs(currGA.fAltitude - 0.0) < 0.05)
		currGA.fAltitude = GetElevation(currGA.dLongtitude, currGA.dLatitude, 0.0);

	/* 构造命令行，生成通用机场 */
	if(!XchgModelBtgToDestination())
	{
		return 0;
	}

	/* 根据计算好的osg输出一个包围线的ive。*/
	if(!ReadAptOSGAndOutputLineLoopIve())
	{
		return 0;
	}
	
	/* 写STG文件 */
	if(!WriteAllStgs(currGA.cICAO, 1))
	{
		MessageToDialog("创建STG文件失败\n");
		return 0;
	}		

	return 1;
}


/* 通用机场主过程，窗口按钮调用入口。 */
int CreateGeneralAirport()
{

#ifdef MAKEMODELBTG
	//ReadDataAndWriteToBTG();
	ReadModelBTGAndVerify();
	return 1;
#endif //MAKEMODELBTG

	if(!BeginCreateCulture())
	{
		MessageToDialog("切割地景失败！！！");
	}
	return 1;
}

/*============================================================================================
 *  下面是用户自定义机场定位的制作过程。目前用户机场只支持flt格式。
 *  基本方法是读取flt文件，计算这个flt文件的轮廓，然后根据经纬度定位信息，将轮廓平移到目标位置，
 *  最后在地景上挖出这个轮廓。
 *============================================================================================ */

/**********************************************************
 *
 *  遍历模型目录，搜寻.flt文件，提取信息存入数据库
 *       
 **********************************************************/
/* 在totalUserBld里面查找当前文件名在不在？*/ 
int CheckFltInUserBld(char *fltname, char *curPath, int *idx)
{
	int i, retval = 0;

	for(i=0; i<bldNum; i++)
	{
		if((!strcmp(fltname, totalUserBld[i].fltFname)) && (!strcmp(curPath, totalUserBld[i].fltPname)))
		{ retval = 1; break; }
	}

	/* 如果当前三维物体文件存在，但是定位文件不存在，返回2 */
	if(i < bldNum)
	{
		if(!strcmp("NULL", totalUserBld[i].locFname))
		{
			*idx = i;
			retval = 2;
		}
	}

	return retval;
}
 
/* 在子目录里面搜索flt文件 */
/* 遍历当前目录下的所有子目录，以查找 flt 文件，使用递归方式。 */
char Tpath[_MAX_PATH];
char locname[32];
int TravelSubDirectoryForFlt(char *dir)
{
	int dio, dfio, dirIo, i;
	struct _finddata_t sdff;

	//查找当前路径下的子目录；
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		//如果不是子目录也不是文件
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR_FLT;

		//如果是子目录，进入当前子目录
		if(sdff.attrib == _A_SUBDIR)
		{

			//进入子目录
			dirIo = chdir(sdff.name);
			_getcwd(Tpath, _MAX_PATH);
			
			//遍历子目录
			TravelSubDirectoryForFlt(Tpath);
			
			//退回到上一层目录
			dirIo = chdir("..");
		}
		
		//如果是文件
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//找到Tif文件
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='F')||(sdff.name[i-3]=='f')) && ((sdff.name[i-2]=='L')||(sdff.name[i-2]=='l')) && ((sdff.name[i-1]=='T')||(sdff.name[i-1]=='t')))  
			{
				int idx;
				int  len = strlen(sdff.name);
				FILE *floc;
				Bucket bt;
				char basepath[16];
				LocateBuilding lb;

				/* 获取当前目录 */
				_getcwd(Tpath, _MAX_PATH);
				
				/* 先查查这个文件是否已经在库里面了。从后往前查*/
				if(!CheckFltInUserBld(sdff.name, Tpath, &idx))
				{					
					/* 填充totalUserBld部分内容 */
					strcpy(locname, sdff.name); locname[len-3] = 'l'; locname[len-2] = 'o'; locname[len-1] = 'c';
					
					/* 打开定位文件，读取三维物体坐标信息 */
					floc = fopen(locname, "rt");
					if(floc != NULL)
					{
						/*	如果没有同名loc文件，说明该flt不一定是机场文件。可以不添加。
						/* 开始支持一个机场往多个地点放 2012-01-13 */
						while(!feof(floc))
						{
							strcpy(lb.fltFname, sdff.name);
							strcpy(lb.fltPname, Tpath);
							strcpy(lb.locFname, locname);
							ReadALineToArray(floc);
							if(!strlen(cf[0])) continue;
							lb.fLongi       = atof(cf[2]);
							lb.fLati        = atof(cf[3]); 
							lb.fElevation   = atof(cf[4]);
							lb.fHeading     = atof(cf[5]);
	
							set_bucket(&bt, lb.fLongi, lb.fLati);
							gen_base_path(&bt, basepath);
							if(strlen(iveroot) == 0) sprintf(iveroot, "%s\\Scenery\\Terrain", currTi.iveDir);
							sprintf(lb.terPname, "%s\\%s", iveroot, basepath);
							lb.iCode        = 0;
							lb.ID           = bldNum;
							lb.isOK         = 0;
							totalUserBld.push_back(lb);		bldNum++;
						}
						fclose(floc);
					}
					needUpdate = 1;
				}
				
				/* 如果这个三维物体已经在列表中，但是还没有被定位，则补充定位信息 */
				if(CheckFltInUserBld(sdff.name, Tpath, &idx) == 2)
				{
					/* 打开定位文件，读取三维物体坐标信息 */
					strcpy(locname, sdff.name); locname[len-3] = 'l'; locname[len-2] = 'o'; locname[len-1] = 'c';
					floc = fopen(locname, "rt");
					if(floc)
					{
						strcpy(totalUserBld[idx].locFname, locname);
						ReadALineToArray(floc);
						totalUserBld[idx].fLongi       = atof(cf[2]);
						totalUserBld[idx].fLati        = atof(cf[3]); 
						totalUserBld[idx].fElevation   = atof(cf[4]);
						totalUserBld[idx].fHeading     = atof(cf[5]);

						set_bucket(&bt, totalUserBld[idx].fLongi, totalUserBld[idx].fLati);
						gen_base_path(&bt, basepath);
						if(strlen(iveroot) == 0) sprintf(iveroot, "%s\\Scenery\\Terrain", currTi.iveDir);
						sprintf(totalUserBld[idx].terPname, "%s\\%s", iveroot, basepath);
						totalUserBld[idx].iCode        = 0;
						totalUserBld[idx].ID           = bldNum;
						totalUserBld[idx].isOK         = 0;
						
						/* 如果一个机场定位多处地方 */
						while(!feof(floc))
						{
							strcpy(lb.fltFname, sdff.name);
							strcpy(lb.fltPname, Tpath);
							strcpy(lb.locFname, locname);
							ReadALineToArray(floc);
							if(!strlen(cf[0])) continue;
							lb.fLongi       = atof(cf[2]);
							lb.fLati        = atof(cf[3]); 
							lb.fElevation   = atof(cf[4]);
							lb.fHeading     = atof(cf[5]);
	
							set_bucket(&bt, lb.fLongi, lb.fLati);
							gen_base_path(&bt, basepath);
							if(strlen(iveroot) == 0) sprintf(iveroot, "%s\\Scenery\\Terrain", currTi.iveDir);
							sprintf(lb.terPname, "%s\\%s", iveroot, basepath);
							lb.iCode        = 0;
							lb.ID           = bldNum;
							lb.isOK         = 0;
							totalUserBld.push_back(lb);		bldNum++;
						}
						fclose(floc);
						needUpdate = 1;
					}
				}				
			}
		}

GET_NEXT_SUBDIR_FLT:
		//寻找下一个
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* 遍历当前目录，查找该目录下以及所有子目录下扩展名为.flt的文件，读取里面相关信息。 */
int TravelDirectoryToReadFlt(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForFlt(dPath))
	{
		MessageToDialog("找不到带地理信息的纹理文件");
		return 0;
	}
	CheckConfig(needUpdate);
	
	return 1;
} 

/**************************************************
 *
 *   开始定位用户机场
 *
 **************************************************/
double MiscPI()
{	return	double(3.1415926535897932384626433832795029L); }

osg::Quat fromLonLatRad(double dlon, double dlat)
{
	double lon = dlon * MiscPI() / 180.0;
	double lat = dlat * MiscPI() / 180.0;

	osg::Quat q;
	double zd2 = double(0.5)*lon;
	double yd2 = double(-0.25)*MiscPI() - double(0.5)*lat;
	double Szd2 = sin(zd2);
	double Syd2 = sin(yd2);
	double Czd2 = cos(zd2);
	double Cyd2 = cos(yd2);
	q.w() = Czd2*Cyd2;
	q.x() = -Szd2*Syd2;
	q.y() = Czd2*Syd2;
	q.z() = Szd2*Cyd2;
	return q;
}

/* 找出一个geode所有顶点的经纬度坐标的极值 */
int FindMMLLInGeode(osg::ref_ptr<osg::Geode> gnode)
{
	Geodesy geod;
	Carton  cart;
	
	for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
	{
		osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
		if( geometry.valid() )
		{
			osg::ref_ptr<osg::Vec3Array> vx = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray());
			for( osg::Vec3Array::iterator pvx = vx->begin(); pvx != vx->end(); pvx++) 
			{
				cart.x = pvx->x();	cart.y = pvx->y();	cart.z = pvx->z();
				CartToGeod(&cart, &geod);
				SetMMLLEx(&currFlt, geod._lon, geod._lat);
			}
		}
	}

	return 1;
}



/* 检查stg文件，看看里面是否有需要修改的机场项目，如果有，就修改。*/
int CheckSTGAndReplaceFltItem(int idx, char *fn)
{
	FILE *fi, *fo;
	std::vector<LocateBuilding> currStg;
	
	/* 以文本读方式打开 当前的 stg 文件*/
	fi = fopen(fn, "rt");
	if(fi == NULL)
	{
		MessageToDialog("不能打开物体stg文件。\n");
		return 0;
	}

	/* 查找一下文件内容，是否包含这个全局机场的内容。通过文件名的比较来确定。*/
	currStg.clear();
	while(!feof(fi))
	{
		LocateBuilding lb;
		
		/* 读取stg文件的每一行 */
		ReadALineToArray(fi);
		strcpy(lb.locFname, cf[0]);
		strcpy(lb.fltPname, cf[1]);
		
		/* 从 路径文件名 里面单独提取 文件名 */
		std::string fltPath = lb.fltPname;
		splitPathName(fltPath);
		strcpy(lb.fltFname, materialParts[materialParts.size() - 1].c_str());
		
		/* 比较文件名与当前idx对应文件名是否一致 */
		if(!strcmp(lb.fltFname, findResultFlt[idx].fltFname))
		{
			lb.fLongi       = findResultFlt[idx].fLongi;
			lb.fLati        = findResultFlt[idx].fLati; 
			lb.fElevation   = findResultFlt[idx].fElevation;
			lb.fHeading     = findResultFlt[idx].fHeading;
		}else
		{
			lb.fLongi       = atof(cf[2]);
			lb.fLati        = atof(cf[3]); 
			lb.fElevation   = atof(cf[4]);
			lb.fHeading     = atof(cf[5]);
		}
		currStg.push_back(lb);
	}
	fclose(fi);
	
	/* 以文本写的方式打开当前文件，重写 这个stg */
	fo = fopen(fn, "wt");
	if(fo == NULL)
	{
		MessageToDialog("不能打开物体stg文件。\n");
		return 0;
	}
	for(unsigned int i=0; i<currStg.size(); i++)
	{
		if(strlen(currStg[i].locFname) != 0)
			fprintf(fo, "%s %s %13.8f %13.8f %f %f\n", currStg[i].locFname, currStg[i].fltPname, currStg[i].fLongi, currStg[i].fLati, currStg[i].fElevation, currStg[i].fHeading);
	}
	fclose(fo);
	
	return 1;
}


/* 在子目录里面搜索stg文件 */
/* 遍历当前目录下的所有子目录，以查找 stg 文件，使用递归方式。 */
int TravelSubDirectoryForStg(int idx, char *dir)
{
	int dio, dfio, dirIo, i;
	struct _finddata_t sdff;

	//查找当前路径下的子目录；
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		//如果不是子目录也不是文件
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR_STG;

		//如果是子目录，进入当前子目录
		if(sdff.attrib == _A_SUBDIR)
		{

			//进入子目录
			dirIo = chdir(sdff.name);
			_getcwd(Tpath, _MAX_PATH);
			
			//遍历子目录
			TravelSubDirectoryForStg(idx, Tpath);
			
			//退回到上一层目录
			dirIo = chdir("..");
		}
		
		//如果是文件
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//找到Stg文件
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='S')||(sdff.name[i-3]=='s')) && ((sdff.name[i-2]=='T')||(sdff.name[i-2]=='t')) && ((sdff.name[i-1]=='G')||(sdff.name[i-1]=='g')))
			{
				/* 先查查这个文件里面是否包含了相关项。*/
				if(!CheckSTGAndReplaceFltItem(idx, sdff.name))
				{					
						return 0;
				}
			}
		}

GET_NEXT_SUBDIR_STG:
		//寻找下一个
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* 遍历当前目录，查找该目录下以及所有子目录下扩展名为.Stg的文件，读取里面相关信息。 */
int TravelDirectoryToReadStg(int idx, char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForStg(idx, dPath))
	{
		MessageToDialog("找不到STG文件");
		return 0;
	}

	return 1;
} 


/* BigEndian 和 LittleEndian 互换，定长8Bytes。*/
void endian_swap8byte(char *buff)
{
	std::swap(buff[0], buff[7]);	std::swap(buff[1], buff[6]);	std::swap(buff[2], buff[5]);	std::swap(buff[3], buff[4]);
}

/* 坐标转换，根据参数，从局部坐标系转换到世界坐标系 */
int ConvertCoordination(int idx, osg::ref_ptr<osg::Geode> gnode, std::string fullName)
{
	Geodesy geod;
	Carton  cart;

#if 0		//暂不考虑读取全局坐标的问题。因为即使找到了经纬度，也找不到高度和朝向。
	/* 首先用二进制方法读取flt文件，查看是否有OriginLon和OriginLat，如果有，就是全局坐标系的flt。反之如果没有，就是局部坐标系的flt。*/
	/* 读取flt文件里面的原点经纬度坐标，和使用的顶点坐标单位。*/
	FILE *flt;
	fpos_t Origin = 220, CoordinateUnit = 62;
	char dBuff[2][8];
	double originLon = 0.0, originLat = 0.0;
	int iCoordinateUnit = 0;
	flt = fopen(fullName.c_str(), "rb");
	if(flt)
	{
		fsetpos(flt, &CoordinateUnit);	fread(&iCoordinateUnit, 1, 1, flt);
		fsetpos(flt, &Origin);	fread(dBuff, 8, 2, flt);
		endian_swap8byte(dBuff[0]);	endian_swap8byte(dBuff[1]);	
		memcpy(&originLat, dBuff[0], 8);	memcpy(&originLon, dBuff[1], 8);
		fclose(flt);
	}
	
	if((originLon != 0.0) && (originLat != 0.0))
	{
		/* 将基准位置从Geodesy坐标转为Carton坐标*/
		geod._lon = originLon; geod._lat = originLat; geod._elevation = 0.0;
		GeodToCart(&geod, &cart);
		osg::Vec3d pos; pos.set(cart.x, cart.y, cart.z);

		/* 计算目标转换矩阵 */
		osg::Quat hlOr = fromLonLatRad(originLon, originLat);
		osg::Matrix obj_pos; obj_pos.set(hlOr);
		osg::Quat flip(0.0, 1.0, 0.0, 0.0);
		obj_pos.preMult(osg::Matrix(flip));
		obj_pos.setTrans(pos);

		double hdg = 0.0;
		obj_pos.preMult(osg::Matrix::rotate(hdg * SGD_DEGREES_TO_RADIANS, 0.0, 0.0, 1.0));

		/* 设置全局转换矩阵变量*/
		coordMatrix = obj_pos;
		rotateMatrix = osg::Matrix::rotate(hdg * SGD_DEGREES_TO_RADIANS, 0.0, 0.0, 1.0);	//不知道全局坐标系下，heading数值怎么得到。这里缺省朝向正北，0.0

		/* 所有顶点全部转换 */
		ChangeCoordSystem(gnode, obj_pos);
		
		/* 要修改 Objects 目录下 相关 stg 里面的相关值。*/
		findResultFlt[idx].fLongi = originLon;	findResultFlt[idx].fLati = originLat;	findResultFlt[idx].fElevation = 0.0;	findResultFlt[idx].fHeading = hdg;
		if(strlen(objroot) == 0) sprintf(objroot, "%s\\Scenery\\Objects", currTi.iveDir);
		TravelDirectoryToReadStg(idx, objroot);
	}
	else
#endif
	{
		/* 将基准位置从Geodesy坐标转为Carton坐标*/
		geod._lon = findResultFlt[idx].fLongi; geod._lat = findResultFlt[idx].fLati; geod._elevation = findResultFlt[idx].fElevation;
		GeodToCart(&geod, &cart);
		osg::Vec3d pos; pos.set(cart.x, cart.y, cart.z);
	
		/* 计算目标转换矩阵 */
		osg::Quat hlOr = fromLonLatRad(findResultFlt[idx].fLongi, findResultFlt[idx].fLati);
		osg::Matrix obj_pos; obj_pos.set(hlOr);
		osg::Quat flip(0.0, 1.0, 0.0, 0.0);
		obj_pos.preMult(osg::Matrix(flip));
		obj_pos.setTrans(pos);
	
		double hdg = findResultFlt[idx].fHeading;
		obj_pos.preMult(osg::Matrix::rotate(hdg * SGD_DEGREES_TO_RADIANS, 0.0, 0.0, 1.0));
			
		/* 设置全局转换矩阵变量*/
		coordMatrix = obj_pos;
		rotateMatrix = osg::Matrix::rotate(hdg * SGD_DEGREES_TO_RADIANS, 0.0, 0.0, 1.0);
	
		/* 所有顶点全部转换 */
		ChangeCoordSystem(gnode, obj_pos);
		
		/* 找出所有顶点的经纬度极值，现在看来求这个并无实际意义。权且放在这儿。 */
		//InitMMLLEx(&currFlt, findResultFlt[idx].fltFname);
		//FindMMLLInGeode(gnode);
	}

	return 1;
}

/* 计算CorrectionSpaceMesh */
void AddCorrectionSpaceMesh(osg::ref_ptr<osg::Group> gp)
{
	osg::ref_ptr<CorrectionSpaceMesh> csm = new CorrectionSpaceMesh( gp, 1.0 );
	csmList.push_back( csm );
}
 
 
/* 根据idx指定的findResultFlt参数。读取用户三维物体flt文件，并计算出它的外围环Geometry, 加入到gnode里面。*/
int ReadFltAndComputePerimeter(int idx, int lon, int lat)
{
	int io, ffio;
	char path[_MAX_PATH], dPath[_MAX_PATH], fullName[_MAX_PATH], iveName[24], currTerPath[STRING_STANDARD_LEN * 2], basepath[32], fullTrueName[_MAX_PATH];
	struct _finddata_t fltff;
	
	/* 确定当前要保存在地景库中的路径 terPath */
	Bucket Bt;
	double fLon = lon + 0.0000001;
	double fLat = lat + 0.0000001;
	set_bucket(&Bt, fLon, fLat);
	gen_base_path(&Bt, basepath);
	if(strlen(iveroot) == 0) sprintf(iveroot, "%s\\Scenery\\Terrain", currTi.iveDir);
	sprintf(currTerPath, "%s\\%s", iveroot, basepath);

	/* 获取并保存当前目录 */
	_getcwd(path, _MAX_PATH);

	/* 进入flt文件所在的目录 */
	sprintf(dPath, "%s", findResultFlt[idx].fltPname);
	io = chdir(dPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "不能进入目标目录 -- %s", dPath);
		MessageToDialog(msg);
		return 0;
	}

	/* 查找这个目录下的flt文件。 */
	ffio = _findfirst(findResultFlt[idx].fltFname, &fltff);	
	if((!ffio)||(ffio == -1))
	{
		MessageToDialog("找不到flt文件。");
		return 0;
	}

	/* 先查找当前路径下是否有Around.flt文件，如果有，就使用Around.flt作为查找用户机场边缘环的主文件。*/
	ffio = _findfirst("Around.flt", &fltff);	
	if((!ffio)||(ffio == -1))
	{
		/* 拼装全路径文件名 */
		sprintf(fullName, "%s\\%s", findResultFlt[idx].fltPname, findResultFlt[idx].fltFname);
	}else{
		/* 拼装全路径文件名 Around.flt */
		sprintf(fullName, "%s\\Around.flt", findResultFlt[idx].fltPname);
	}
	_findclose(ffio);

	/* 创建一个gnode, 待添加各个三维物体的perimeter。 */
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	gd->setName("River");
	
	/* 计算用户机场flt的边缘环，用来在地景上挖洞。并把这个边缘环作为Geometry加到gnode里面 */
	if(!ReadFltAndVistIt(fullName, gd))
	{
		if(!foundAround)
			MessageToDialog("机场flt里面没有Around组，不能计算边缘。");
		else
			MessageToDialog("机场边缘没有封闭，不能正确切割。");
		MessageToDialog("计算用户三维物体外围环失败。");
		return 0;
	}

	/* 计算用户机场flt的灯光信息。并写入灯光临时文件 AptLight.ive*/
	/* 首先拼装全路径文件名 */
	sprintf(fullTrueName, "%s\\%s", findResultFlt[idx].fltPname, findResultFlt[idx].fltFname);
	if(!ReadFltAndCreateAptLight(fullTrueName))
	{
		char msg[_MAX_PATH];
		sprintf(msg, "当前机场没有灯光信息 --- %s", fullTrueName);
		MessageToDialog(msg);
	}

	
	/* 转换坐标 */
	if(!ConvertCoordination(idx, gd, fullName))
	{
		MessageToDialog("转换坐标失败。");
		return 0;
	}

	/* 转换完坐标后，就可计算CorrectionSpaceMesh了 */
	osg::ref_ptr<osg::Group> gp = new osg::Group;
	gp->addChild(gd);
#if 1
	AddCorrectionSpaceMesh(gp);
#endif

	/* 确定当前机场代码 */
	if(findResultFlt[idx].iCode == 0)
	{
		SetCurrFltAPCode(idx);
	}

	/* 将这个gnode保存成类似通用机场那样的ive文件 文件名为 GUR?.ive*/
	sprintf(iveName, "GUR%c.ive",  findResultFlt[idx].iCode);
	sprintf(fullName, "%s\\%s", currTerPath, iveName);
	osgDB::writeNodeFile(*gd, fullName);
		
	/* 写机场灯光BTG */
	sprintf(fullName, "%s\\%s", findResultFlt[idx].fltPname, findResultFlt[idx].fltFname);
	if(ReadFltAndCreateLightsDataStru(fullName))
	{
		char fn[STRING_STANDARD_LEN],  msg[_MAX_PATH];
		char fullPathFile[_MAX_PATH];

		sprintf(iveName, "GUR%c.btg.gz",  findResultFlt[idx].iCode);
		sprintf(fullPathFile, "%s\\%s", currTerPath, iveName);
		strcpy(fn, iveName);
	
		/* 先把tile.btg.gz写入各个stg中 */
		if(WriteDataStructToBTG(fullPathFile))
		{
			if(!WriteAllStgs(fn, 1)) return 0;
		}else{
			sprintf(msg, "写入机场灯光文件失败 --- %s", fullPathFile);
			MessageToDialog(msg);
		}		
	}

	return 1;
}

/* 计算所有的用户三维物体数据 */
int DoAllFlt()
{
	int i, j, k;
	char msg[STRING_STANDARD_LEN];
	int currDoing[16];
	int numCD;
	PerZone mPZ[4];
	int numPZ, iFlag;

	/* 找出还没有定位用户三维物体所在的地景目录名称 */
	if(FindOutDirsNoBldDone() == 0)
	{
		MessageToDialog("没有需要定位的三维物体.");
		return 1;	
	}

	/* 先从dirsNoCultDone搜起，即这些目录里面都有没做用户三维物体的定位 */
	for(i=0; i<numDNCD; i++)
	{
		/* 针对dirsNoCultDone里面每一项，把totalUserBld里面所有同Terrain目录名的项集合到findResultFlt*/
		memset(currDoing, 0, sizeof(currDoing));  numCD = 0;
		memset(mPZ, 0, sizeof(mPZ));  numPZ = 0;
		if(FindAllFltItemsWithSameDir(i) == 0) continue;
		strcpy(findJpgDirectory, findResultFlt[0].terPname);
		
		/* 针对findResultFlt里面每一项，先检查是否做过用户三维物体定位处理，没做过就在currDoing里记录。*/
		for(j=0; j<numFRFLT; j++)
		{
			if(findResultFlt[j].isOK == 0)
			{	
				currDoing[numCD++] = j;
			}
		}
		
		/* 从刚刚做好的所有块里面找出所有不同的块 */
		for(j=0; j<numCD; j++)
		{
			for(int n=0; n<4; n++)
			{
				/* 确定以机场中心为原点，一定范围内方形区域四个角的经纬度值 */
				double cornerLon, cornerLat;
				switch(n)
				{
					case 0:	{	cornerLon = findResultFlt[currDoing[j]].fLongi - 0.05; cornerLat = findResultFlt[currDoing[j]].fLati - 0.025; break;	}
					case 1:	{	cornerLon = findResultFlt[currDoing[j]].fLongi - 0.05; cornerLat = findResultFlt[currDoing[j]].fLati + 0.025; break;	}
					case 2:	{	cornerLon = findResultFlt[currDoing[j]].fLongi + 0.05; cornerLat = findResultFlt[currDoing[j]].fLati - 0.025; break;	}
					case 3:	{	cornerLon = findResultFlt[currDoing[j]].fLongi + 0.05; cornerLat = findResultFlt[currDoing[j]].fLati + 0.025; break;	}
				}				
				
				/* 找出这个方框所占据的整数经纬度区域的个数 */
				iFlag = 0;
				for(k=0; k<numPZ; k++)
				{
					int iLon, iLat;
					iLon = (int)cornerLon;
					iLat = (int)cornerLat;
					if(( iLon == mPZ[k].lon) && ( iLat == mPZ[k].lat))
					{ iFlag = 1; break; }
				}
				if(!iFlag)
				{	
					mPZ[numPZ].lon = (int)cornerLon;
					mPZ[numPZ].lat = (int)cornerLat; numPZ++;
				}
			}
		}

		/* 对于每个块，先做基础地景，再加定位三维物体 */
		for(j=0; j<numPZ; j++)
		{
			/* 先做初始化全局变量 */
			currTi.minLongi = (float)(mPZ[j].lon);
			currTi.minLati  = (float)(mPZ[j].lat);
			currTi.maxLongi = (float)(mPZ[j].lon + 1);
			currTi.maxLati  = (float)(mPZ[j].lat + 1);
			currPZ = &mPZ[j];
			

			/* 根据currDoing里面的信息，计算用户三维物体的底面壳边缘 */
			for(k=0; k<numCD; k++)
			{
				if(findResultFlt[currDoing[k]].isOK == 0)
				{
					sprintf(msg, "正在计算边缘--%s...", findResultFlt[currDoing[k]].fltFname);
					MessageToDialog(msg);
					if(!ReadFltAndComputePerimeter(currDoing[k], mPZ[j].lon, mPZ[j].lat))
					{
						sprintf(msg, "计算边缘--%s失败。", findResultFlt[currDoing[k]].fltFname);
						MessageToDialog(msg);
						return 0;
					}
				}
			}   

			/* 将相同块内的所有btg组成一个Geode, 插入到挖洞列表当中，待集中统一挖洞。 */
			if(!InsertItemIntoCutList(mPZ[j].lon, mPZ[j].lat))
			{   
				sprintf(msg, "转换IVE失败。");
				MessageToDialog(msg);
				return 0;
			}
		}

		/* 计算完成并且无误，则更改isOK的状态 */
		for(k=0; k<numCD; k++)
		{
			if(findResultFlt[currDoing[k]].isOK == 0)
			{
				findResultFlt[currDoing[k]].isOK = 1;
			}
		}

		/* 做好的文化信息块更新总体数据库 */
		UpdateTotalUserBldWithNewFindResult();
	}
	
	return 1;
}

 
/* 定位用户机场主过程 */
int LocateUserAirport()
{
	return 1;
}
 
