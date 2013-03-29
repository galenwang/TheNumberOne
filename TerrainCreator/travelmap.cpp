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

/* ���岢��¼��ǰ״̬ */
#define STATUS_UNWORK				0
#define STATUS_SUCCESS				1
#define STATUS_HGTBAD				2
#define STATUS_TEXTBAD				3
#define STATUS_HGTGEOBAD			4
#define STATUS_TEXTGEOBAD			5

/* ��ǰ�������Ǹ߳� */
#define STATUS_TEXTURE				1
#define STATUS_HGT					2
#define MBUFFSIZE 					0x100000

/* �����ַ����ռ䶨�壬��ֹ�����������*/
#define STRING_STANDARD_LEN	128

using std::string;

/************************************************
 *
 *    �����ⲿ���̼�����
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
 *    ����˵���ö�
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
 *    ȫ�ֱ���
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
int currHighLevel;			//��ǰ�������𣬸���Ŀ¼������ȷ��
int updateTTNode;
char cf[64][STRING_STANDARD_LEN];
FILE *fi, *fo, *fLog;
gzFile gzfo;
MapNode mergedTif;						//���Ϻ��ͼƬ�����Ϣ
MapNode mergedHgt;						//���Ϻ�ĸ߳������Ϣ

/* ȫ��·�� */
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

/* ��ȡ�̵߳�ȫ�ֱ���*/
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

/* ������ive��������������������ݽṹ */

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

/* �����ؾ���߶ȣ�ʹ�õؾ������ӻ���ƽ������ */
std::vector< osg::ref_ptr<CorrectionSpaceMesh> > csmList;
	
/* ��ȡcolors.txt������������ */
std::vector<ColorName> colors;

/* flt��������任ȫ�ֱ��� */
osg::Matrix coordMatrix;
osg::Matrix rotateMatrix;

/* flt����������ռ����ļ�ֵ */
MaxMinLonLat currFlt;

/* ϸ�������ռ� */
std::vector<std::string> allDetailedTxt;
	
/* Shp�ļ����� */
std::vector<ShpManagement> allShps;
	
/* ��ʼ�������Ŀ */
InitializationCheck initCheck;	

/* �������е��Ļ���Ϣ������Ŀ */
std::vector<HighCultureInfo> allHCI;

/* ��·��������ȹ�� */
osg::ref_ptr<osg::Group> allLakeSkirt;
osg::ref_ptr<osg::Group> allRoadSkirt;
osg::ref_ptr<osg::Group> allDynamicGrass;
osg::ref_ptr<osg::Group> allIslands;
osg::ref_ptr<osg::Group> allLakeInCurrLonLat;


/**************************************************
 *
 *   �� currTi.lod �ֽ�Ϊ 3 �� LOD �ַ���
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
 *   ����tileID�ţ��Ӹ���ؾ������б����ҵ������Ŀ����ֵallHCI
 *
 **************************************************/
void InsertHCI(int tileID)
{
	/* ����tileID����allHCI�в��ң����Ƿ��Ѿ���ӹ���*/
	for(unsigned int i = 0; i<allHCI.size(); i++)
	{
		if(allHCI[i].tileID == tileID)	return;
	}
	
	/* ��allHTI�в���ͬtileID����Ŀ��*/
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

/* ����tileID�ţ���鵱ǰ�Ļ���Ϣ�Ƿ�������*/
int CheckHCI(int tileID)
{
	/* ����tileID����allHCI�в��ң����Ƿ��Ѿ���ӹ���*/
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

/* д���ã����û��߷ֱ����Ļ���Ϣ�������д�������ļ��� */
int WriteHighCultureInfo2Config(FILE *fo)
{
	unsigned int i;
	
	fprintf(fo, "\n%d",  allHCI.size());
	//��д��Map��
	for(i=0; i<allHCI.size(); i++)
	{
		fprintf(fo, "\n%s %11.6f %11.6f %11.6f %11.6f %9d %2d", allHCI[i].texFname, allHCI[i].startLongi, allHCI[i].startLati, allHCI[i].offsetLongi, allHCI[i].offsetLati, allHCI[i].tileID, allHCI[i].isOK);
	}
	return 1;
}

/* �����ã��������ļ��ж�ȡ�û��߷ֱ����Ļ���Ϣ��������� */
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
 *   ���ؾ������ļ�д��Log
 *
 **************************************************/
void WriteTerrainCfgIntoLogFile()
{
	/* д�û��������ݽṹ����Log */
	fprintf(fLog, "\n��ǰ�ؾ�������Ϣ��\n");
	fprintf(fLog, "��ǰ�ؾ�����LOD����: %d\n", currTi.lod );
	WriteHTxtInfo2Config(fLog);
	WriteHHgtInfo2Config(fLog);
	WriteHighInfo2Config(fLog);
	fprintf(fLog, "\n");
}

void WriteCultureCfgIntoLogFile()
{
	/* д�û��������ݽṹ����Log */
	fprintf(fLog, "\n");
	fprintf(fLog, "\n��ǰ�Ļ���Ϣ������Ϣ��\n");
	WriteShpInfo2Config(fLog);
	WriteBldInfo2Config(fLog);
	WriteHighCultureInfo2Config(fLog);
	fprintf(fLog, "\n");
}

/************************************************
 *
 *    ��ȡ�ļ�����
 *
 ************************************************/

fpos_t GetFileSize(char *fn)
{
	FILE *ft;
	fpos_t iFileSize;
	
		
	ft = fopen(fn, "rb");
	if(ft == NULL)	return 0;

	/*��ȡ�ļ�����*/
	fseek(ft, 0, SEEK_END);
	fgetpos(ft, &iFileSize);
	fclose(ft);
	
	return iFileSize;
}

/************************************************
 *
 *    ��ָ����γ��Ŀ¼�µ�����btg�ļ����Ϊtxt�ļ����Թ�����
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
	
		/* �����·���²��� *.btg.gz */
		io = chdir(fullPath);
		
		dio = _findfirst("*.btg.gz", &sdff);
		if((!dio)||(dio==-1)) goto CBTT_Continue_1;
		dfio = 0;
	
		while(!dfio)
		{
			/* ������Ļ���Ϣbtg���뵽һ�� group���� */				
			AnalysisBTG(sdff.name);	

			/* Ѱ����һ�� */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}
CBTT_Continue_1:
	return;
}


/************************************************
 *
 *    �Ӹ߳��ļ��ж�ȡ��ǰ��ĸ߳�����
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
	
	/* ��õ�ǰλ�ø߳��ļ�·�������ļ��� */
	iCurrLong = (int) lon;
	iCurrLat  = (int) lat;

	sprintf(hgtFilename, "N%02dE%03d.hgt", iCurrLat, iCurrLong);
	sprintf(hgtPathname, "%s\\n%02d", hgtroot, iCurrLat);
	
	/* �򿪲���ȡ�߳��ļ� */
	TGhgt * curHgt = &currHgt;
	if((curHgt->curLong != iCurrLong) || (curHgt->curLati != iCurrLat))
	{
		sprintf(fullPath, "%s\\%s", hgtPathname, hgtFilename);
		if(!openhgt(curHgt, fullPath))
		{
			return 1.0;
		}
	}
	
	/* ��ȡ�߳����� */
	clon = lon * 3600.0; clat = lat * 3600.0;
	cElev = altitude_from_grid(curHgt, clon, clat);
	
	return cElev;
}
	


/************************************************
 *
 *    ��һ���ļ��ж�ȡcolors����������Ϣ
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

/* ����һ��color����ɫֵ����color���в�������λ�� */
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
 *    ��һ��totalShp���棬���ݲ�ͬ��������γ�ȷֳ�һ������findResultAC
 *    ���ڴ����Ļ���Ϣ��Խ��γ�ȵ������
 *
 ************************************************/
int SplitByLonLatShp(int idx)
{
	int minLon, maxLon, minLat, maxLat, i, j;

	/* ���ҳ�findResultAC��ռ�ݵ������С��γ������ֵ */
	minLon = (int)totalShp[idx].minLong;
	maxLon = (int)totalShp[idx].maxLong;
	minLat = (int)totalShp[idx].minLati;
	maxLat = (int)totalShp[idx].maxLati;
	
	/* ��˫��ѭ����������ռ�����������findAllItemsAC */
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
			
			/* ���һ�£���ǰ���Ƿ�ֵ���� ��*/
			if(fabs(cn.minLong - cn.maxLong) < 0.001) continue;
			if(fabs(cn.minLati - cn.maxLati) < 0.001) continue;
			
			/* ���ֵ��������������Ҫ��������Ŀ ��*/
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
 *    ��һ��totalcult���棬���ݲ�ͬ��������γ�ȷֳ�һ������findResultAC
 *    ���ڴ����Ļ���Ϣ��Խ��γ�ȵ������
 *
 ************************************************/
int SplitByLonLatAC(int idx)
{
	int minLon, maxLon, minLat, maxLat, i, j;

	/* ���ҳ�findResultAC��ռ�ݵ������С��γ������ֵ */
	minLon = (int)totalcult[idx].minLong;
	maxLon = (int)totalcult[idx].maxLong;
	minLat = (int)totalcult[idx].minLati;
	maxLat = (int)totalcult[idx].maxLati;
	
	/* ��˫��ѭ����������ռ�����������findAllItemsAC */
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
			
			/* ���һ�£���ǰ���Ƿ�ֵ���� ��*/
			if(fabs(cn.minLong - cn.maxLong) < 0.001) continue;
			if(fabs(cn.minLati - cn.maxLati) < 0.001) continue;
			
			/* ���ֵ��������������Ҫ��������Ŀ ��*/
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
 *    ���û������룬����findResultFltǰ���Ѿ�����Ĵ�����������õ�ǰ���룬���ܺ�ǰ����ظ��ˡ�
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
 *    ����tile��findResultAC������ң��Ƿ�����ͬ��tile���ڣ�ֱ��ȷ��һ����һ�޶��ġ�
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
 *    ��ȡ��������Path�������������Լ��Ļ�������
 *
 ************************************************/
static void CreateEnvPath()
{
	std::vector<string> paths;
	char onepath[STRING_STANDARD_LEN];
	int i, j = 0;
	
	/* ��ȡ��ǰ��������PATH������paths�б����� */
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
    
    /* paths�б��뵱ǰ����Ŀ¼��һ�Ƚϣ��������ͬ�ľ��˳������û����ͬ�ľͼ����� */
    string currwork = workDirectory;
    for(i=0; i<paths.size(); i++)
    {
    	if(paths[i] == currwork)
    	{
    		return;
    	}
    }
    
    /* ����ǰ����Ŀ¼���뵽�б��У������µ�path��д�뻷���������� */
	char newpath[2048];
	sprintf(newpath, "set path=%s;%s", envp, workDirectory);
	::putenv(newpath);
}

/************************************************
 *
 *    ���ݻ��߻ָ������Ļ���Ϣ�ĵؾ����ļ�
 *
 ************************************************/

/* �����ļ�����src���Ƶ�dst���� */
int CopyFile(char *src, char *dst)
{
	FILE *pInput;
	FILE *pOutput;
	fpos_t iFileSize, iDataLen;

	/* ���ļ� */	
	if((pInput=fopen(src,"rb"))==NULL||						
		(pOutput=fopen(dst,"wb"))==NULL)
	{
		return 0;
	}

	/*��ȡ�ļ�����*/
	fseek(pInput, 0, SEEK_END);
	fgetpos(pInput, &iFileSize);
	rewind(pInput);

	/*��д�ļ���*/
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
	
	/* �ر��ļ� */
	delete pcInputLine;
	fclose(pInput);
	fclose(pOutput);
	return 1;
}

/* �ѵ�ǰ�ļ��ӵ�ǰĿ¼���Ƶ�TempĿ¼�£��������ļ��������͸�����Ŀ¼�µ�����ļ� */
int BackUpFileToTemp(char *fn, char *fnpath)
{
	char fullSubPath[_MAX_PATH], currDir[_MAX_PATH], backSubDir[_MAX_PATH];
	char dstFile[_MAX_PATH];
	char subDir[_MAX_PATH], currRoot[_MAX_PATH];
	struct _finddata_t iveSubf;

	/* ���ȱ��ݵ�ǰ��·��, Ȼ����뵱ǰ������·�� */
	if(strlen(fnpath) != 0)
	{
		_getcwd(currRoot, _MAX_PATH);
		int ioCur = chdir(fnpath);
		if(ioCur) 
		{
			char msg[STRING_STANDARD_LEN];
			sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", fnpath);
			MessageToDialog(msg);
			return 0;
		}
	}
	
	/* ���ƶ��� */
	sprintf(dstFile, "%s\\%s", backDirectory, fn);
	
	/* ���һ�±����ļ��ڲ��ڣ�����Ѿ����ݹ��ˣ��Ͳ����ٱ����� --- �ݲ���Ҫ */
	if(!CopyFile(fn, dstFile))	return 0;
	
	/* �����Ǵ���PageLOD�ĵײ�ĸ���ive�ļ������еײ��ļ������ XXXX_root_L0_X0_Y0 Ŀ¼�� */
	strcpy(subDir, fn);
	int len = strlen(subDir);
	subDir[len - 4] = 0;
	sprintf(fullSubPath, "%s_root_L0_X0_Y0", subDir);

	/* �ڱ���·���´�����Ŀ¼ */
	_getcwd(currDir, _MAX_PATH);
	int io = chdir(backDirectory);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", backDirectory);
		MessageToDialog(msg);
		return 0;
	}
	
	io = chdir(fullSubPath);
	if(io) 
	{
		int dio = mkdir(fullSubPath);
		if(dio)
		{
			MessageToDialog("���ܴ����ؾ���Ŀ¼��\n");
			return 0;
		}
	}
	sprintf(backSubDir, "%s\\%s", backDirectory, fullSubPath);

	/* ����Դ����Ŀ¼�����뵽��Ŀ¼�� */
	io = chdir(currDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", currDir);
		MessageToDialog(msg);
		return 0;
	}

	io = chdir(fullSubPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", fullSubPath);
		MessageToDialog(msg);
		return 0;
	}

	/* ����Ŀ¼�в���ive�ļ� */
	int subffio = _findfirst("*.ive", &iveSubf);	
	if((!subffio)||(subffio == -1))
	{
		MessageToDialog("�Ҳ���ive�ļ���");
		return 0;
	}

	int subfio = 0, subidx = 0;
	while(!subfio)
	{
		subidx++;
		sprintf(dstFile, "%s\\%s", backSubDir, iveSubf.name);
		
		/* ���ø��ƺ��� */
		if(!CopyFile(iveSubf.name, dstFile))	return 0;

		/* ������һ�������������ļ� */
		subfio = _findnext(subffio, &iveSubf);
	}
	_findclose(subffio);

	/* �˳�ǰ����Դ����Ŀ¼ */
	io = chdir(currDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", currDir);
		MessageToDialog(msg);
		return 0;
	}
	
	/* ���й������֮�󣬷��ظ�·�� */
	if(strlen(fnpath) != 0)
	{
		int ioCur = chdir(currRoot);
		if(ioCur) 
		{
			char msg[STRING_STANDARD_LEN];
			sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", currRoot);
			MessageToDialog(msg);
			return 0;
		}
	}
	
	return 1;
}

/* ��TempĿ¼�ָ��ļ�����ǰĿ¼�£��������ļ����ָ���Ŀ¼�µ�����ļ� */
int RestoreFileFromTemp(char *fn)
{
	char fullSubPath[_MAX_PATH], currDir[_MAX_PATH], backSubDir[_MAX_PATH];
	char dstFile[_MAX_PATH];
	char subDir[24];
	struct _finddata_t iveSubf;

	sprintf(dstFile, "%s\\%s", backDirectory, fn);

	/* �ȱȽ϶���������ļ�����,����������,˵�������ļ�û�仯,�Ǿ�û��Ҫ�ٸ�����.*/
	fpos_t srcsize, tgtsize;
	if((srcsize = GetFileSize(fn)) == 0) return 0;
	if((tgtsize = GetFileSize(dstFile)) == 0) return 2;
	if(srcsize == tgtsize) return 0;

	/* ���ƶ��� */
	if(!CopyFile(dstFile, fn))	return 2;
	
	/* �����Ǵ���PageLOD�ĵײ�ĸ���ive�ļ������еײ��ļ������ XXXX_root_L9_X0_Y0 Ŀ¼�� */
	strcpy(subDir, fn);
	int len = strlen(subDir);
	subDir[len - 4] = 0;
	sprintf(fullSubPath, "%s_root_L0_X0_Y0", subDir);
	sprintf(backSubDir, "%s\\%s", backDirectory, fullSubPath);

	/* �ڱ���·���´�����Ŀ¼ */
	_getcwd(currDir, _MAX_PATH);
	int io = chdir(backSubDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", backSubDir);
		MessageToDialog(msg);
		return 0;
	}
	
	/* ����Դ����Ŀ¼�����뵽��Ŀ¼�� */
	io = chdir(currDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", currDir);
		MessageToDialog(msg);
		return 0;
	}

	io = chdir(fullSubPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", fullSubPath);
		MessageToDialog(msg);
		return 0;
	}

	/* ����Ŀ¼�в���ive�ļ� */
	int subffio = _findfirst("*.ive", &iveSubf);	
	if((!subffio)||(subffio == -1))
	{
		MessageToDialog("�Ҳ���ive�ļ���");
		return 0;
	}

	int subfio = 0, subidx = 0;
	while(!subfio)
	{
		subidx++;
		sprintf(dstFile, "%s\\%s", backSubDir, iveSubf.name);
		
		/* ���ø��ƺ��� */
		if(!CopyFile(dstFile, iveSubf.name))	return 0;

		/* ������һ�������������ļ� */
		subfio = _findnext(subffio, &iveSubf);
	}
	_findclose(subffio);

	/* �˳�ǰ����Դ����Ŀ¼ */
	io = chdir(currDir);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", currDir);
		MessageToDialog(msg);
		return 0;
	}
	return 1;
}

/* �����ļ��������������ݻ��ǻָ����ļ�������FILELENMIN��FILELENMAX֮�䣬�򱸷ݣ�����ָ��� */
#define FILELENMAX 717300
#define FILELENMIN 717200

/* �ڻָ��ؾ��ļ�֮ǰ���ȼ��һ�±��ݹ���ive���ļ����ȣ�
   ��������ļ������ڣ�˵����û�б����ݡ�
   ������ļ�����û�仯�����б����ļ���˵���ļ�û�б��Ķ�������Ҫ�����ݻ��߻ָ�����
   �ļ�����û�仯�򷵻� 0�� �����ں��б仯�򷵻� 1��
  */
int CheckFileLengthBeforeRestore(int tID)
{
	char backIveFile[_MAX_PATH];
	char fullIveFile[_MAX_PATH];
	fpos_t iFileSizeSrc, iFileSizeBkp;
	int iFlag = 1;

	/* ��ȡ��ǰ��һ��ive�ļ��� */
	sprintf(fullIveFile, "%d.ive", tID);
	iFileSizeSrc = GetFileSize(fullIveFile);

	/* ��ȡ���ݹ���һ��ive�ļ��� */
	sprintf(backIveFile, "%s\\%d.ive",backDirectory, tID);
	iFileSizeBkp = GetFileSize(backIveFile);
	if(iFileSizeBkp == 0)	return iFlag;

	/* �Ƚ������ļ��ĳ��� */
	if(iFileSizeSrc == iFileSizeBkp) iFlag = 0;

	return iFlag;
}


/************************************************
 *
 *    ����һ��tile.btg.gz��Ȼ��ֻ�������cultures����д��һ���µ�Btg���棬�ļ��� A+tile.btg.gz
 *
 ************************************************/
int CreateCultureBtgFromSrc(int idx)
{
	char fn[STRING_STANDARD_LEN], afn[STRING_STANDARD_LEN], msg[_MAX_PATH];

	char fullPathFile[_MAX_PATH];
	sprintf(fullPathFile, "%s\\%s", findResultAC[idx].destPath, findResultAC[idx].destFile);
	strcpy(fn, findResultAC[idx].destFile);

	/* �Ȱ�tile.btg.gzд�����stg�� */
	if(WriteDataStructToBTG(fullPathFile))
	{
		if(!WriteAllStgs(fn, 1)) return 0;
	}else{
		sprintf(msg, "���ļ�ʧ�� --- %s", fullPathFile);
		MessageToDialog(msg);
		return 0;
	}

	/* ����tile.btg.gz��ȥ�����Ļ���Ϣ���֣����A+tile.btg.gz*/
	sprintf(fullPathFile, "%s\\A%s", findResultAC[idx].destPath, findResultAC[idx].destFile);
	sprintf(afn, "A%s", fn);
	if(WriteDataStructToBTG_CulturePlusIsland(fullPathFile))
	{
		/* ��A+tile.btg.gzҲд�����stg�� */
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
			
			/* ����tile.btg.gz��ȥ�����Ļ���Ϣ���֣����A+tile.btg.gz*/
			sprintf(fullPathFile, "%s\\A%s", fRtoCutIsland[i].destPath, fRtoCutIsland[i].destFile);  

			/* �ȶ�������ļ����ڴ�. */
			if(!ReadBTGToDataStruct(fullPathFile)) continue;

			/* ��ɾ��������ļ��� */
			DeleteFile(fullPathFile);

			/* Ȼ���ٸ����ڴ��е�������������*/
			if(!WriteDataStructToBTG_OnlyCulture(fullPathFile)) continue;
		}
	}
	fRtoCutIsland.clear();

	return 1;
}

/************************************************
 *
 *    �����Ļ���Ϣ�������
 *
 ************************************************/
typedef struct _CUTLIST_
{
	int lon;
	int lat;
}CutList;
std::vector<CutList> allNeedCut;			//ͳһ�ڶ����б�ÿһ���ʾ��һ����γ���������ڶ���

/* �ѵ�ǰ�ľ�γ����ϲ��뵽CutList�У������������ͬ��ֵ�������ظ����� */
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

/* ����CutList�����ݣ��Եؾ�������ڶ����� */
int DoTerrainCut()
{
	char msg[STRING_STANDARD_LEN];

	for(std::vector<CutList>::iterator p = allNeedCut.begin(); p != allNeedCut.end(); p++)
	{
		/* �Ȱѵ�ǰ��γ��Ŀ¼�µ�btg���Ϊtxt */
		ConvertBTGToTXT(p->lon, p->lat);
		
		/* ����ǰ�飨1������x1��γ���������Ŀռ䣩�ڵ����е��Ļ���ϢBtg�ļ���Ȼ�����ͨ�û������û�����������ive�ļ�����ͬ���һ����Group���������Group�Եؾ������ڶ��� */
		if(!ReadTerrainIveAndCut(p->lon, p->lat))
		{   
			sprintf(msg, "ת��IVEʧ�ܡ����򣺾��� %d γ�� %d", p->lon, p->lat);
			MessageToDialog(msg);
		}
	}
	allNeedCut.clear();
	return 1;
}

/************************************************
 *
 *    �����û��Զ�����ά�����б�
 *
 ************************************************/
/* ��totalUserBld�����棬�ҳ�����dirsNoCultDone����һ��terPname�����м�¼���������õĺ�û���õ� */
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

/* �ҳ�����δ��λ�û���ά����ĵؾ�Ŀ¼���� */
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

/* �Ѿ���λ��ɵ���ά������µ�������ά�����б��� */
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

/* д���ã���������ά�����б�д�������ļ� */
int WriteBldInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  bldNum);
	//��д���
	for(i=0; i<bldNum; i++)
	{
		fprintf(fo, "\n%s %s %s %s %13.8f %13.8f %13.8f %13.8f %d %5d %1d", totalUserBld[i].fltPname, totalUserBld[i].fltFname, totalUserBld[i].locFname, totalUserBld[i].terPname, totalUserBld[i].fLongi, totalUserBld[i].fLati, totalUserBld[i].fElevation, totalUserBld[i].fHeading, totalUserBld[i].iCode, totalUserBld[i].ID, totalUserBld[i].isOK);
	}

	return 1;
}

/* �����ã��������ļ��ж�ȡ������ά�����б� */
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
 *    ���Ӱ������ɵĹ���
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

/* ��ʼ�����������������*/
int InitOnDemand()
{
	inputParam.maxLongi = 0.0;
	inputParam.maxLati  = 0.0; 
	inputParam.minLongi = 0.0;
	inputParam.minLati  = 0.0; 
	return 1;
}

/* ���ð������ɵ����� */
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

/* �жϵ�ǰ���������Ƿ��ڰ������ɷ�Χ�� */
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
 *  �Ӵ��������ҵ��ض���Ŀ¼����
 *  Ŀǰ�Ӿ�ϵͳ������Ҫ���ض�Ŀ¼���£�
 *  BaseTexture: ��������ؾ�����Դ���ݣ�����Դ���ݷ�Ϊ���֣�һ����JPEG+Map���ݣ�һ����GeoTiff���ݣ���User\50mĿ¼�¡�
 *               �û��߷ֱ��ʵؾ����������User��Ŀ¼�¡��û�Ӧ�����ݸ��Ƶ� User\1m  User\5m  User\10m Ŀ¼��
 *  BaseHgt:     ��������߳�Դ���ݡ��û��ĸ߷ֱ��ʸ߳�����Ӧ���Ƶ� User ���档
 *
 ************************************************/
char PathFound[_MAX_PATH];
int  PathIsFound;

/* ����Ŀ¼�������pathname���֡�*/
int TravelSubDirectoryForPathname(char *dir,  char *pathname, int level)
{
	int dio, dfio, dirIo;
	struct _finddata_t sdff;

	//����Ѿ��ҵ��ˣ��Ͳ���Ҫ������
	if(PathIsFound) return 1;

	//���ҵ�ǰ·���µ���Ŀ¼��
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//���������Ŀ¼Ҳ�����ļ�
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR_PATHNAME;

		//�������Ŀ¼�����뵱ǰ��Ŀ¼
		if(sdff.attrib == _A_SUBDIR)
		{
			//����ҵ��ˡ�
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
				//������Ŀ¼
				if(level >= 4)
				{

				}
				else
				{
					dirIo = chdir(sdff.name);
					if(!dirIo)
					{
						_getcwd(path, _MAX_PATH);
				
						//������Ŀ¼
						level++;
						TravelSubDirectoryForPathname(path, pathname, level);

						//�˻ص���һ��Ŀ¼
						dirIo = chdir("..");
						level--;
					}
				}
			}
		}

GET_NEXT_SUBDIR_PATHNAME:
		//Ѱ����һ��
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* �Ӵ����ϲ����ض�Ŀ¼ */
int FindSpecDirFromDisk(char *dirname)
{
	int i, drv, dio;
	char path[_MAX_PATH];
	
	memset(PathFound, 0, sizeof(PathFound)); PathIsFound = 0;
	//�Ƚ��뵽�����Ŀ¼
	drv = 'C' - 'A' + 1;
	for(i=0; i<8; i++)
	{
		//����Ѿ��ҵ��ˣ��Ͳ���Ҫ������
		if(PathIsFound) break;
		
		dio = _chdrive(drv);
		if(!dio)
		{
			sprintf(path, "%c:\\", drv + 'A' - 1);
			dio = chdir(path);
			
			//������Ŀ¼
			sprintf(path, "%c:", drv + 'A' - 1);
			TravelSubDirectoryForPathname(path, dirname, 1);
		}
		drv++;
	}
	
	return 1;
}

/* �ڹ���Ŀ¼�¶�ȡĿ¼�����ļ���������������ʹӴ����в��ң�������ҵ�����д���Ŀ¼�����ļ�
   ������Ҳ���������0��
 */
int FindBasePathName()
{
	FILE *fpath;
	int io, isValid = 1;

	/* ����Ҫȷ���ڹ���Ŀ¼��*/
	io = chdir(workDirectory);
	if(io)  return 0;
	
	/* ��ȡPathName�ļ� */
	fpath = fopen("secret.pth", "rt");
	if(fpath != NULL )
	{
		fgets(baseTextDir, _MAX_PATH, fpath);
		fgets(baseHgtDir,  _MAX_PATH, fpath);
		fgets(fsData,  _MAX_PATH, fpath);
		fclose(fpath);
		
		/* �������·���Ƿ���Ч */
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
	
	/* ���û�ҵ�����ļ� */
	int iFlag = 1;
	memset(baseTextDir, 0, sizeof(baseTextDir)); memset(baseHgtDir, 0, sizeof(baseHgtDir));	memset(fsData, 0, sizeof(fsData));
	FindSpecDirFromDisk("BaseTexture");
	if(PathIsFound) {	strcpy(baseTextDir, PathFound);}	else {	iFlag *= 0;	initCheck.pathOfTexture = false;	}
	FindSpecDirFromDisk("BaseHgt");
	if(PathIsFound) {	strcpy(baseHgtDir,  PathFound);}	else {	iFlag *= 0;	initCheck.pathOfHgt = false;	}
	FindSpecDirFromDisk("FG");
	if(PathIsFound) {	strcpy(fsData,		PathFound);}	else {	iFlag *= 0;	initCheck.pathOfTarget = false;	}
	ifFoundBasePath = iFlag;
	
	/* ���û���ҵ���������Ŀ¼������һ�������˳���*/
	if(!iFlag)	return 0;
	
	/* ����Ҫȷ���ص�����Ŀ¼��*/
	io = chdir(workDirectory);
	if(io)  return 0;

	/* �������Ŀ¼���ҵ��ˣ�д�ļ�����*/
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
 *    �ҳ�����δ�����Ļ���Ϣ��shp�ļ���Ŀ¼����
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

/* ��totalShp�����棬�ҳ�����dirsNoCultDone����һ��pathName�����м�¼���������õĺ�û���õ� */
int FindAllShpItemsWithSameDir(int idx)
{
	int i;
	
	if(idx >= numDNCD) return 0;
	findResultAC.clear();	numFRAC = 0;
	for(i=0; i<shpID; i++)
	{
		if(!strcmp(dirsNoCultDone[idx].c_str(), totalShp[i].pathName))
		{
			/* ���һ��shp�Ļ���Ϣ��Խ�����γ�ȣ���ÿ����γ��ϸ��һ�¡�*/
			SplitByLonLatShp(i);
		}
	}
	return numFRAC;
}

/* �Ѿ�������ɵ�shp�Ļ���Ϣ���µ�����shp�б��� */
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
		
		/* ���isOK == 2��˵��Ҫ���һ�����죬���һ��fRtoCutIsland */
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
 *    �ҳ�����δ�����Ļ���Ϣ��ac�ļ���Ŀ¼����
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

/* ��totalcult�����棬�ҳ�����dirsNoCultDone����һ��pathName�����м�¼���������õĺ�û���õ� */
int FindAllAcItemsWithSameDir(int idx)
{
	int i;
	
	if(idx >= numDNCD) return 0;
	findResultAC.clear();	numFRAC = 0;
	for(i=0; i<cultID; i++)
	{
		if(!strcmp(dirsNoCultDone[idx].c_str(), totalcult[i].pathName))
		{
			/* ���һ��AC�Ļ���Ϣ��Խ�����γ�ȣ���ÿ����γ��ϸ��һ�¡�*/
			SplitByLonLatAC(i);
		}
	}
	return numFRAC;
}

/* �Ѿ�������ɵ�ac�Ļ����µ�����ac�б��� */
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

		/* ���isOK == 2��˵��Ҫ���һ�����죬���һ��fRtoCutIsland */
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
 *    ����ȫ�ֱ�����������Ϣ���������ĶԻ���
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
	/* ͬʱҪдLog�ļ���*/
	/* ��Ϣǰ�����ʱ��ͷ */
	struct tm *local;  
	time_t t;  
	t = time(NULL);  
	local = localtime(&t); 
	
	/* ƴװ��������Ϣ */
	char fullmsg[128];
	sprintf(fullmsg, "%02d:%02d:%02d: %s", local->tm_hour, local->tm_min, local->tm_sec, msg);
	
	/* д���¼�ļ���*/
	fprintf(fLog, "%s\n", fullmsg);
#endif
}


/************************************************
 *
 *    ����ÿһ��Ŀ¼����
 *
 ************************************************/
int CheckTextureDirName(char *dirname, char *currPath)
{
	int i, slen, ilat, ilon, jlat, jlon, flag;
	char clat[2][5], clon[2][5];
	int  lat[2], lon[2], minlat, minlon;
	
	/* ���ַ��������ȡ�ĸ����� */
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
	
	/* �жϵ�ǰ·���Ƿ�Ϸ� */
	if((strlen(clon[0])==0)||(strlen(clat[0])==0)||flag==-1)
		return 1;
	
	/* ת��������, ȡСֵ */
	lat[0] = atoi(clat[0]);  lat[1] = atoi(clat[1]); lon[0] = atoi(clon[0]); lon[1] = atoi(clon[1]);
	if(lat[0]<lat[1]) minlat = lat[0]; else minlat = lat[1]; if(lon[0]<lon[1]) minlon = lon[0]; else minlon = lon[1]; 
	
	/* ��д���ݽṹ */
	TextDirInfo tdi;
	sprintf(tdi.dirName, "%s\\%s", currPath, dirname);
	tdi.zone[0].isOK = STATUS_UNWORK;	tdi.zone[0].lat = minlat;	tdi.zone[0].lon = minlon;
	tdi.zone[1].isOK = STATUS_UNWORK;	tdi.zone[1].lat = minlat;	tdi.zone[1].lon = minlon+1;
	tdi.zone[2].isOK = STATUS_UNWORK;	tdi.zone[2].lat = minlat+1;	tdi.zone[2].lon = minlon;
	tdi.zone[3].isOK = STATUS_UNWORK;	tdi.zone[3].lat = minlat+1;	tdi.zone[3].lon = minlon+1;
	allDirs.push_back(tdi);		TextDirIdx++;
	
	return 1;
}

/*  ��鵱ǰ���Ƿ��Ѿ����꣬���귵��0�� û������currTi, ����1 */
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

/* д���ã���Ŀ¼��Ϣд�������ļ��� */
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

/* �����ã��������ļ��ж�ȡ����Ŀ¼��Ϣ������ͼƬƴ�ӷ�ʽ�ĵؾ����� */
int ReadDirInfoFromConfig(FILE *fi)
{
	int i, j, k, blen;
	char currPath[STRING_STANDARD_LEN];
	
	blen = strlen(baseTextDir);	
	ReadALineToArray(fi);	TextDirIdx = atoi(cf[0]);	allDirs.clear();
	for(i=0; i<TextDirIdx; i++)
	{
		ReadALineToArray(fi);
		
		/*���Ŀ¼��ǰ�벿����baseTextDir��ͬ������baseTextDir���� */
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
 *    ��дMap��Ϣ���������ļ�����
 *
 ************************************************/
/* д���ã���Map��Ϣд�������ļ��� */
int WriteMapInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  nodeID);
	//��д��Map��
	for(i=0; i<nodeID; i++)
	{
		fprintf(fo, "\n%11.6f %11.6f %11.6f %11.6f %5d %5d %s %s", totalnode[i].minLong, totalnode[i].maxLong, totalnode[i].minLati, totalnode[i].maxLati, totalnode[i].image_width, totalnode[i].image_height, totalnode[i].fileName, totalnode[i].pathName);
	}
	return 1;
}

/* �����ã��������ļ��ж�ȡ����Map��Ϣ������ͼƬƴ�ӷ�ʽ�ĵؾ����� */
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
		
		/*���Ŀ¼��ǰ�벿����baseTextDir��ͬ������baseTextDir���� */
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
 *    ��дHighTxt��Ϣ���������ļ�����
 *
 ************************************************/
/* д���ã����û��߷ֱ��ʵؾ�������Ϣд�������ļ��� */
int WriteHTxtInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  hTxtNum);
	//��д��HighTxt��
	for(i=0; i<hTxtNum; i++)
	{
		fprintf(fo, "\n%11.6f %11.6f %11.6f %11.6f %5d %5d %1d %s %s", highTxtnode[i].minLong, highTxtnode[i].maxLong, highTxtnode[i].minLati, highTxtnode[i].maxLati, highTxtnode[i].image_width, highTxtnode[i].image_height, highTxtnode[i].highLevel, highTxtnode[i].fileName, highTxtnode[i].pathName);
	}
	return 1;
}

/* �����ã��������ļ��ж�ȡ�����û��߷ֱ��ʵؾ�������Ϣ������Tif��ʽ�ĵؾ����� */
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

		/*���Ŀ¼��ǰ�벿����baseTextDir��ͬ������baseTextDir���� */
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
 *    ��дLowTxt��Ϣ���������ļ�����
 *
 ************************************************/
/* д���ã���ϵͳȱʡ�ͷֱ��ʵؾ�������Ϣд�������ļ��� */
int WriteLTxtInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  lTxtNum);
	//��д��LowTxt��
	for(i=0; i<lTxtNum; i++)
	{
		fprintf(fo, "\n%11.6f %11.6f %11.6f %11.6f %5d %5d %1d %s %s", lowTxtnode[i].minLong, lowTxtnode[i].maxLong, lowTxtnode[i].minLati, lowTxtnode[i].maxLati, lowTxtnode[i].image_width, lowTxtnode[i].image_height, lowTxtnode[i].highLevel, lowTxtnode[i].fileName, lowTxtnode[i].pathName);
	}
	return 1;
}

/* �����ã��������ļ��ж�ȡ����ϵͳȱʡ�ֱ��ʵؾ�������Ϣ������Tif��ʽ�ĵؾ����� */
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

		/*���Ŀ¼��ǰ�벿����baseTextDir��ͬ������baseTextDir���� */
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
 *    ��дHighHgt��Ϣ���������ļ�����
 *
 ************************************************/
/* д���ã����û��߷ֱ��ʵؾ��߳���Ϣд�������ļ��� */
int WriteHHgtInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  hHgtNum);
	//��д��HighTxt��
	for(i=0; i<hHgtNum; i++)
	{
		fprintf(fo, "\n%11.6f %11.6f %11.6f %11.6f %5d %5d %1d %s %s", highHgtnode[i].minLong, highHgtnode[i].maxLong, highHgtnode[i].minLati, highHgtnode[i].maxLati, highHgtnode[i].image_width, highHgtnode[i].image_height, highHgtnode[i].highLevel, highHgtnode[i].fileName, highHgtnode[i].pathName);
	}
	return 1;
}

/* �����ã��������ļ��ж�ȡ�����û��߷ֱ��ʵؾ��߳���Ϣ������Tif��ʽ�ĵؾ����� */
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

		/*���Ŀ¼��ǰ�벿����baseTextDir��ͬ������baseTextDir���� */
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
 *    ��дShp�Ļ���ʼ��Ϣ���������ļ�����
 *
 ************************************************/
/* д���ã����û�shp��ʽ�Ļ���Ϣд�������ļ��� */
int WriteShpInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  shpID);
	//��д��Map��
	for(i=0; i<shpID; i++)
	{
		fprintf(fo, "\n%s %s %11.6f %11.6f %11.6f %11.6f %s %s %3d %2d %d %d %d", totalShp[i].fileName, totalShp[i].pathName, totalShp[i].minLong, totalShp[i].maxLong, totalShp[i].minLati, totalShp[i].maxLati, totalShp[i].destFile, totalShp[i].destPath, totalShp[i].destLongi, totalShp[i].destLati, totalShp[i].destTileID, totalShp[i].ID, totalShp[i].isOK );
	}
	return 1;
}

/* �����ã��������ļ��ж�ȡ�����û�shp��ʽ��Ϣ�� */
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
 *    ��дAc�Ļ���ʼ��Ϣ���������ļ�����
 *
 ************************************************/
/* д���ã����û�ac��ʽ�Ļ���Ϣд�������ļ��� */
int WriteCultInfo2Config(FILE *fo)
{
	int i;
	
	fprintf(fo, "\n%d",  cultID);
	//��д��Map��
	for(i=0; i<cultID; i++)
	{
		fprintf(fo, "\n%s %s %11.6f %11.6f %11.6f %11.6f %s %s %3d %2d %d %d %d", totalcult[i].fileName, totalcult[i].pathName, totalcult[i].minLong, totalcult[i].maxLong, totalcult[i].minLati, totalcult[i].maxLati, totalcult[i].destFile, totalcult[i].destPath, totalcult[i].destLongi, totalcult[i].destLati, totalcult[i].destTileID, totalcult[i].ID, totalcult[i].isOK );
	}
	return 1;
}

/* �����ã��������ļ��ж�ȡ�����û�ac��ʽ��Ϣ�� */
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
 *    ����߾��ȵؾ��������б�
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
	double py;				//ʵ��ֱ������ϵ����   ��Ҫ����ת��������
	double pz;				
	double lon;
	double lat;				//��γ�����꣬��Ҫ�������
	double elev;
	double tU;
	double tV;				//��������   ֱ�ӿɶ�ȡ��
}TerrainPoint;

std::vector<TerrainOnDemand> doTerrain;
int numDT;


/**************************************************
 *
 *   �жϵ�ǰ���tileID�ڵؾ������Ƿ��еؾ����ڣ����ڷ���1�� ���򷵻�0
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

	//�������·����Ϣ
	W = hti->startLongi;
	S = hti->startLati;
	offLon = hti->offsetLongi;
	offLat = hti->offsetLati; 
	set_bucket(&b, W+0.01, S+0.01);

	/* ƴװ���·�����ļ� */
	gen_base_path(&b, basePath);
	sprintf(fullPath, "%s\\%s\\%d.ive", iveroot, basePath, hti->tileID);
	
	/* �����Ƿ��ܹ����ļ� */
	bin = fopen(fullPath, "rb");
	if(bin != NULL)
	{
		iFlag = 1;	fclose(bin);
	}
	return iFlag;
}


/* ��allHTI��������ռ���˶���ƽ����  */
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

/* ��allHTI�����ҳ���ͬһ��Zone(ƽ����)�е�������Ŀ�������findResultHTI���� */
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


/* ���ĳһ����ĸ߾��ȵؾ��Ƿ��Ѿ�ȫ��������� */
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
 *  ����������߾��ȵؾ������㷨��֧�ָ��Ӹ��ӵ���� ---- 2012.01.20
 *****************************************************************************************************/
/* ��findResultHTI������������о�γ����Ϣ���ֱ���д�С��������� */
double LongiSortList[32], LatiSortList[32];
int numSLon, numSLat;
#define DELTA_THREAHOLD 0.0001
int SortTheFRHTI(TerrainOnDemand *curr)
{
	int i, j;
	
	/* �Ȱѱ߽���Ϣ�ӽ�ȥ */
	memset(LongiSortList, 0, sizeof(LongiSortList)); memset(LatiSortList, 0, sizeof(LatiSortList));
	LongiSortList[0] = curr->minLongi; LongiSortList[1] = curr->maxLongi;  numSLon = 2;
	LatiSortList[0]  = curr->minLati;  LatiSortList[1]  = curr->maxLati;   numSLat = 2;
	
	/* �ٰ�findResultHTI����������С��γ��ֵ���ӽ�ȥ��������ظ��ģ�����Բ��� */
	for(i=0; i<numFRHTI; i++)
	{
		double minFrLon, minFrLat, maxFrLon, maxFrLat;
		int iFlag = 0;
		minFrLon = findResultHTI[i].startLongi; maxFrLon = findResultHTI[i].startLongi + findResultHTI[i].offsetLongi;
		minFrLat = findResultHTI[i].startLati;  maxFrLat = findResultHTI[i].startLati  + findResultHTI[i].offsetLati;

		/* ����߾��ȿ�ȵ�ǰ���ȿ���������ڵ�ǰ���ȿ��� */
		if(minFrLon < curr->minLongi) minFrLon = curr->minLongi;
		if(minFrLat < curr->minLati ) minFrLat = curr->minLati ;
		if(maxFrLon > curr->maxLongi) maxFrLon = curr->maxLongi;
		if(maxFrLat > curr->maxLati ) maxFrLat = curr->maxLati ;
		
		/* ��� minFrLon ���ظ��󣬲��������� */
		iFlag = 0;
		for(j=0; j<numSLon; j++)
		{
			double delta = fabs(LongiSortList[j] - minFrLon);
			if(delta < DELTA_THREAHOLD)	{ iFlag = 1; break; }
		}
		if(!iFlag) { LongiSortList[numSLon++] = minFrLon; }

		/* ��� maxFrlon ���ظ��󣬲��������� */
		iFlag = 0;
		for(j=0; j<numSLon; j++)
		{
			double delta = fabs(LongiSortList[j] - maxFrLon);
			if(delta < DELTA_THREAHOLD)	{ iFlag = 1; break; }
		}
		if(!iFlag) { LongiSortList[numSLon++] = maxFrLon; }

		/* ��� minFrLat ���ظ��󣬲��������� */
		iFlag = 0;
		for(j=0; j<numSLat; j++)
		{
			double delta = fabs(LatiSortList[j] - minFrLat);
			if(delta < DELTA_THREAHOLD)	{ iFlag = 1; break; }
		}
		if(!iFlag) { LatiSortList[numSLat++] = minFrLat; }

		/* ��� maxFrLat ���ظ��󣬲��������� */
		iFlag = 0;
		for(j=0; j<numSLat; j++)
		{
			double delta = fabs(LatiSortList[j] - maxFrLat);
			if(delta < DELTA_THREAHOLD)	{ iFlag = 1; break; }
		}
		if(!iFlag) { LatiSortList[numSLat++] = maxFrLat; }
	}
	
	/* ������������������ LongiSortList,LatiSortList ���дӵ͵��ߵ����� */
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

/* ����һ�����ݽṹ����¼��ά����ÿһ������㣨��㣩����Ϣ */
typedef struct _PERGRIDINFO_
{
	int 	RowIndex;
	int		ColIndex;			//�ڶ�ά�����е�λ�ã�
	double	MinLongi;
	double	MaxLongi;
	double	MinLati;
	double	MaxLati;			//����������������
	double	MidLongi;
	double	MidLati;			//��������ĵ㾭γ��
	int		Info;				//�����״̬	0��δ�����߾��ȸ��ǣ� 1���ѱ����߾��ȸ���
}PerGridInfo;
PerGridInfo AllGrids[32][32];

/* ������õ�LongiSortList��LatiSortList��ɶ�ά���񣬲�����¼ÿ��������״̬ */
int MakeGridFindInfo()
{
	int i, j, k, iFlag;
	memset(AllGrids, 0, sizeof(AllGrids));
	
	/* �ȶ�������и�ֵ ��Ϊ���ȣ���Ϊγ�� */
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
			
			/* �����е��жϱ�ռ��*/
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

/* ���������������Ϣ�����ռ��;��ȵؾ�������������doTerrain */
typedef struct _RECORDGRIDINFO_
{
	int 	Row;
	int		Col;			//�ڶ�ά�����е�λ�ã�
}RecordGridInfo;

/* ���ݶ�ά�����ڸ�����ķֲ������������Info>0�Ǳ��߷ֱ��ʵؾ����ǣ�Info=0û�б����ǣ���������󷽿�ԭ�������ͷֱ��ʵؾ���*/
int CalcDTnByGrids(HighTerrainInfo *curr)
{
	int i, j, k;
	RecordGridInfo RecGrid[32][32];
	int numGRow, numGCol;
	
	while(1)
	{
		int StartRow, StartCol, iFlag, NextRow, NextCol, Lines, bkFlag ;
		
		/* ���Ȳ��ҵ�ǰ�����Ͻǵ� Info = 0 �Ŀ飬��¼��ʼ��StartRow�Ϳ�ʼ��StartCol */
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
		
		/* ����Ҳ���Info == 0 �Ŀ飬���ʾ�Ѿ����꣬�˳�ѭ�� */
		if(!iFlag) break;
		
		/* ���ݿ�ʼ��StartRow�Ϳ�ʼ��StartCol�����²��ң��ҵ���һ�� Info != 0�Ŀ飬���߾�ͷ����¼�� RecGrid��*/
		memset(RecGrid, -1, sizeof(RecGrid)); RecGrid[0][0].Row = StartRow; RecGrid[0][0].Col = StartCol; 	numGRow = numGCol = 1;
		NextRow = StartRow + 1; Lines = 1;
		AllGrids[StartRow][StartCol].Info = 2;
		while((AllGrids[NextRow][StartCol].Info == 0) && (NextRow < (numSLat-1)))
		{
			AllGrids[NextRow][StartCol].Info = 2;
			RecGrid[numGRow][0].Row = NextRow; RecGrid[numGRow][0].Col = StartCol; NextRow++; Lines++; numGRow++;
		}
		
		/* �ٴ� StartCol+1 ��ʼ���Ҳ��ң�һ��һ�еز飬�е�������Ԫ��Info��Ϊ0����ȷ���� */
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
		
		/* ��DoTerrain���鸳ֵ */
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

/* ��������� */
/* ���1����γ�������ڣ�û�б��߾��ȵؾ����ǵĲ��� */
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

/* ����һ���߷ֱ�������鱻��һ�����߷ֱ��ʵ������ȫ���� */
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
 *  ���㷨������� ---- 2012.01.21
 *****************************************************************************************************/

//�����ռ�����doTerrain�����;��ȵؾ���
int DoLowRTerrain()
{
	int i;
	char msg[STRING_STANDARD_LEN];
	
	for(i=0; i<numDT; i++)
	{
		sprintf(msg, "������ %d ��;������顣", i+1);
		MessageToDialog(msg);

		if(!DoPartVpbAndWriteStg(i))
		{
			MessageToDialog("ִ��VPBʧ��\n");
			return 0;			
		}
	}
	return 1;
}

/* д�ص���¼�ܼ����� */
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

/* ���ݾ��Ⱥ�γ�ȵ�����ֵ������findJpgDirectory */
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


/* �����߷ֱ��ʵؾ������� */
int DoHightResolutionTerrain()
{
	int n, j;
	PerZone mPZ;
	char msg[STRING_STANDARD_LEN];
	
	/* ���ȼ�����и߷ֱ��������ܹ�ռ���˶��پ���xγ�ȣ�ƽ���ȣ�*/
	CheckZoneInAllHti();
	
	/* ��ÿ����ռ�õ�ƽ������ѭ��*/
	for(n=0; n<numZWHTI; n++)
	{
		int lowEnable = 1;		//����γ�ȿ��Ƿ������ͷֱ��ʵؾ����ݣ� 2012-11-1

		/* ���ռ���ƽ���������еĸ߾������飬����ŵ� findResultHTI���档 */
		FindAllHTIInSameZone(n);
		mPZ.lon = ZoneWithHTI[n].lon; mPZ.lat = ZoneWithHTI[n].lat; mPZ.isOK = 0;
		
		/* ���������ƽ�����ڵĻ����ͷֱ��ʵؾ� *//* ���ɻ������� */
		sprintf(msg, "��ʼ��������%d γ��%d ����Ļ����ؾ�.", mPZ.lon, mPZ.lat);
		MessageToDialog(msg);
		currTi.minLongi = (float)(mPZ.lon);
		currTi.minLati  = (float)(mPZ.lat);
		currTi.maxLongi = (float)(mPZ.lon + 1);
		currTi.maxLati  = (float)(mPZ.lat + 1);
		currPZ = &mPZ;

		/* ����Ҫƴ���;��ȵؾ��ĸ߳����ݺ��������� */
		int iW = (int)currTi.minLongi, iS = (int)currTi.minLati;
		if(CheckIfHas50(iW, iS) != -1)
		{
			if(!CreatePerDegreeTerrain(2))
			{
				sprintf(msg, "��������%d γ��%d ����Ļ����ؾ�����.", mPZ.lon, mPZ.lat);
				MessageToDialog(msg);
				lowEnable = 0;
			}
		}else{
			if(!CreatePerDegreeTerrain(1))
			{
				sprintf(msg, "��������%d γ��%d ����Ļ����ؾ�����.", mPZ.lon, mPZ.lat);
				MessageToDialog(msg);
				lowEnable = 0;
			}
		}

		/* �ѱ�����û�����߾��ȵؾ���tileȫ������ */
		for(j=0;j<numFRHTI;j++)
		{
			/* �Թ��Ѿ������ */
			if(findResultHTI[j].isOK) continue;
			
			/* ��ʼ����һ��tile�ĸ߾��ȵؾ� */
			sprintf(msg, "��ʼ����Tile %d ����ĸ߾��ȵؾ�.", findResultHTI[j].tileID);
			MessageToDialog(msg);

			/* װ��vpb���������*/
			MessageToDialog("����VPB���ɵؾ�");

			/* װ��vpb���������������vpb��������ǰ��߾��ȵؾ�*/
			if(!DoHighVpbAndWriteStg(j))
			{
				MessageToDialog("ִ��VPBʧ��\n");
				continue;		//return 0;
			}
			
			/* ��¼������Ѿ����� */
			findResultHTI[j].isOK = 1;
		}

		/* ������Ҫ���ɵĵ;�������ϵ�� doTerrain[]*/
		if(!CheckNoHighInZoneEx()) return 0;
		
		/* ���numDTΪ0��˵���õؿ�ȫ�������˸߾��ȵؾ�������Ҫ�����;��ȵؾ��� */
		if((numDT != 0)&&(lowEnable))
		{
			/* �ֱ����������;������� */
			if(!DoLowRTerrain())
			{   
				sprintf(msg, "�����;�������--���� %d γ�� %d ʧ�ܡ�", mPZ.lon, mPZ.lat);
				MessageToDialog(msg);
				return 0;
			}
		}

		/* �������Ժ󣬻�Ҫ�����õĽ��д�ص�allHTI */
		WriteBackToallHTI();

		/* ��allHTIд�������ļ� */
		CheckConfig(1);
		
		sprintf(msg, "���� %d γ�� %d �����ڸ߾��ȵؾ�������ɡ�", mPZ.lon, mPZ.lat);
		MessageToDialog(msg);
	}
	return 1;
}

/************************************************
 *
 *    ���߾�������ϵ��allHTI
 *
 ************************************************/

/* ���õ�ǰHTI, �����ǰ����ˣ������� */
void SetCurrHTI(HighTerrainInfo *curr)
{
	int i, iFlag, iNum;
	
	/* ���ݾ�γ��������Ƿ�֮ǰ����Ŀ���ظ���*/
	iFlag = 0;
	for(i=0;i<numHTI;i++)
	{
		if(	(fabs(allHTI[i].startLongi  - curr->startLongi ) < 0.001) &&
			(fabs(allHTI[i].startLati   - curr->startLati  ) < 0.001) && 
			(fabs(allHTI[i].offsetLongi - curr->offsetLongi) < 0.001) && 
			(fabs(allHTI[i].offsetLati  - curr->offsetLati ) < 0.001)  )
		{ iFlag = 1; break; }
	}
	
	/* ���û���ظ���������Ŀ��ӵ��б����һ�*/
	if(!iFlag)
	{
		HighTerrainInfo hti;

		/* ��֤tileID���ظ� */
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
		hti.tileID      = (curr->tileID * 100) + iNum; 	/* ���赱ǰ����ؿ�һ�����ظ���ID�� */
		strcpy(hti.texFname, curr->texFname);
		strcpy(hti.hgtFname, curr->hgtFname);
		hti.ID          = numHTI;
		hti.isOK		= IfTileIDExistInLib(&hti);
		allHTI.push_back(hti);	numHTI++;
	}
}

/* ���ñȽ�HTI */
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

/* �жϵ�ǰ�ĸ߷ֱ��ʿ飬�Ƿ��������һ��ֱ��ʵ��ӿ� */
int DetectHTContainHigherHT(HighTerrainInfo *curr, int higherLevel)
{
	int i;

	/* �жϸ���һ���ֱ��ʵĿ� �� ��ǰ�ֱ��ʵĿ� �Ƿ����ཻ������о� SetHigherHTI */
	float curW = curr->startLongi, curE = curr->startLongi + curr->offsetLongi, curS = curr->startLati, curN = curr->startLati + curr->offsetLati;
	findResultHTI.clear(); 	numFRHTI=0;
	
	for(i=0; i<numHTI; i++)
	{
		if(allHTI[i].highLevel  <= higherLevel)
		{
			float tmpW = allHTI[i].startLongi, tmpE = allHTI[i].startLongi + allHTI[i].offsetLongi, tmpS = allHTI[i].startLati, tmpN = allHTI[i].startLati + allHTI[i].offsetLati;
			
			/* �߾��ȵؾ������뵱ǰ���ȵؾ������н��棬���ߵ�ǰ���ȵؾ�������ȫ���Ǹ߾��ȵؾ����� */
			if(	(((tmpW > curW) && (tmpW < curE))	|| ((tmpE > curW) && (tmpE < curE))	|| ((curW > tmpW) && (curW < tmpE))	|| ((curE > tmpW) && (curE < tmpE))) &&
				  (((tmpS > curS) && (tmpS < curN))	|| ((tmpN > curS) && (tmpN < curN)) || ((curS > tmpS) && (curS < tmpN))	|| ((curN > tmpS) && (curN < tmpN))) )
			{
				SetHigherHTIEx(&allHTI[i]);
			}
			
			/* �߾��ȵؾ�������ȫ�����˵�ǰ���ȵؾ����ݣ������ξͲ�����ֱ���˳� */
			else if( (tmpW <= curW) && (tmpE >= curE) && (tmpS <= curS) && (tmpN >= curN) )
			{
				return 1;
			}
			
			/* �߾��ȵؾ������뵱ǰ���ȵؾ�����û���غ� */
			else
			{
			}
			
		}
	}
	
	/* ������า�ǵ�������������ηָǰ��curr */
	if(numFRHTI==0)
	{
		/* û�в��ҵ�����������˳� */
		return 0;
	}else{
		/* �в��ҵ�������� ���㸲�Ƿ��� */
		if(!CalcuHighOverlayByHigherEx(curr))	return 0;			
	}
	return 1;
}



/* ���㲻ͬ�����������ݵ��� */
/*  
   ������ӵ�����൱���ӡ�
   �������������ࣺ��ȫ���ǵģ��а븲�ǵģ�
   �����ȼ������֣��ж��㸲�ǵģ������㸲�ǵģ�
   �������ڸ��������֣��е������ǵģ��ж������ǵģ�
   �������㷨ԭ����֧�����е������
*/
int CalculateTextOverlay(HighTerrainInfo *curr)
{
	/* �����ǰ�Ѿ�����߾��ȵ��ˣ����ü��㣬ֱ�Ӹ�ֵ�˳� */
	if(curr->highLevel == 1)
	{
		SetCurrHTI(curr);
		return 1;
	}
	
	/* ���Ҹ��ڵ�ǰ���ȵģ�������ͬһ��γ�������ڵ�����HTI ������ 5�׾��� ���� 1�׾��� */
	if(curr->highLevel == 5)
	{
		if(!DetectHTContainHigherHT(curr, 1))
		{
			SetCurrHTI(curr);
			return 1;
		}
	}

	/* �� 10 �׾��� ���� 1�� �� 5�׾���*/
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


/* �����������Ⱥ�γ����Ϣ�����allHTI�������tileIDs */
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
			/* ���㵱ǰ��ռ���� */
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

			/* ������̫С�Ͳ�ֵ��һ���� */
			/* 2012-11-1 ���̫СҲ����������������з죬��仰Ҫע�͵���*/
			//////////////////////////	if((currHTI.offsetLongi < 0.01) || (currHTI.offsetLati < 0.01)) goto CHFL_NextBlock;

			/* ���û��߳�����������ң��Ƿ���ƥ��ĸ߳��ļ������û��ƥ����û��߳��ļ�����ʹ��ϵͳ�߳��ļ� */
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

			/* ���㲻ͬ�����������ݵ��� */
			CalculateTextOverlay(&currHTI);
//////////CHFL_NextBlock:
			lat++;
		}//while(lat)
		lon++;
	}//while(lon)

	return 1;
}

/* ����totalnode���������ļ���չ��Ϊtif��������Ϣ�����allHTI���*/
int ReadNodeListToCreateHTI()
{
	
	/*  ����Ṇ̀߳����˳� */
	if(isCurrTxtOrHgt == STATUS_HGT) return 1;
	
	/* ���highTxtnode���棬��Щ�Ǹ߾����������� */
	/* ������1m���ȵĵؾ��������� */
	for(int i=0;i<hTxtNum;i++)
	{
		if(highTxtnode[i].highLevel != 1) continue;

		/* ����һ�������е�tile��������� */
		CalculateHTIFromLonLat(i);
		
	}//for(i)

	/* ������5m���ȵĵؾ��������� */
	for(int i=0;i<hTxtNum;i++)
	{
		if(highTxtnode[i].highLevel != 5) continue;

		/* ����һ�������е�tile��������� */
		CalculateHTIFromLonLat(i);
		
	}//for(i)

	/* ������10m���ȵĵؾ��������� */
	for(int i=0;i<hTxtNum;i++)
	{
		if(highTxtnode[i].highLevel != 10) continue;

		/* ����һ�������е�tile��������� */
		CalculateHTIFromLonLat(i);
		
	}//for(i)

	return 1;	
}

/* д���ã����û��߷ֱ��ʵؾ��������д�������ļ��� */
int WriteHighInfo2Config(FILE *fo)
{
	unsigned int i;
	
	fprintf(fo, "\n%d",  numHTI);
	//��д��Map��
	for(i=0; i<numHTI; i++)
	{
		fprintf(fo, "\n%s %s %11.6f %11.6f %11.6f %11.6f %d %8d %2d", allHTI[i].texFname, allHTI[i].hgtFname, allHTI[i].startLongi, allHTI[i].startLati, allHTI[i].offsetLongi, allHTI[i].offsetLati, allHTI[i].highLevel, allHTI[i].tileID, allHTI[i].isOK);
	}
	return 1;
}

/* �����ã��������ļ��ж�ȡ�û��߷ֱ��ʵؾ���������� */
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
 *    ���ݷֱ���Ϊ50m��Geotiff�����������;��ȵؾ�
 *
 ************************************************/
/* ���õ�ǰLowHTI, �����ǰ����ˣ������� */
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

/* ��У��50m���ȵ��������ݣ������ļ�����У�� */
int CheckTextureNameTif50()
{
	int n, i, slen, ilat, ilon, jlat, jlon, flag;
	char clat[2][5], clon[2][5], dirname[24];
	int  lat[2], lon[2], minlat, minlon, maxlon, maxlat;

	for(n=0; n<lTxtNum; n++)
	{
		if(lowTxtnode[n].highLevel != 50) continue;
			
		/* �Ȼ�ȡ�ļ�����������չ���ġ� */
		memset(dirname, 0, sizeof(dirname));
		for(i=0; i<strlen(lowTxtnode[n].fileName); i++) 
		{
			if(lowTxtnode[n].fileName[i] == '.') break;
			dirname[i] = lowTxtnode[n].fileName[i];
		}
		
		/* ���ַ��������ȡ�ĸ����� */
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
		
		/* �жϵ�ǰ·���Ƿ�Ϸ� */
		if((strlen(clon[0])==0)||(strlen(clat[0])==0)||flag==-1)
			return 1;
		
		/* ת��������, ȡСֵ */
		lat[0] = atoi(clat[0]);  lat[1] = atoi(clat[1]); lon[0] = atoi(clon[0]); lon[1] = atoi(clon[1]);
		if(lat[0]<lat[1]) minlat = lat[0]; else minlat = lat[1]; if(lon[0]<lon[1]) minlon = lon[0]; else minlon = lon[1]; 
		if(lat[0]<lat[1]) maxlat = lat[1]; else maxlat = lat[0]; if(lon[0]<lon[1]) maxlon = lon[1]; else maxlon = lon[0]; 
		
		/* ��д���ݽṹ */
		lowTxtnode[n].minLong = minlon;
		lowTxtnode[n].minLati = minlat;
		lowTxtnode[n].maxLong = maxlon;
		lowTxtnode[n].maxLati = maxlat;
	}

	return 1;
}

/* ����Ƿ����50m���ȵĵͷֱ������� */
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

/* ���Ѽ����� */
int GetHTIListOfR50()
{
	double minLon, minLat;
	char fnPath[_MAX_PATH], tileStr[32];
	HighTerrainInfo currHTI;
	int  iMinLon, iMinLat;
	
	/* ��������У��һ�£�������� */
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
					
					/* ������ǰHTI*/
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
					
					/* д��ͷֱ��ʵؾ���������б���*/
					SetCurrHTIToLowHTI(&currHTI);
					
					/* ��һ��γ�� */					
					minLat += 1.0;
				}
				/* ��һ������ */
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
	//��д��Map��
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

/* �����ͷֱ��ʵؾ������� */
int DoLowResolutionTerrain()
{
	int n;
	PerZone mPZ;
	char msg[STRING_STANDARD_LEN];
	
	/* ������LowHTI��ѭ��*/
	for(n=0; n<numLowHTI; n++)
	{
		/* ����Ƿ������������ľ����� */
		if(lowHTI[n].isOK) continue;
		
		/* ���ռ���ƽ���������еĸ߾������飬����ŵ� findResultHTI���档 */
		mPZ.lon = lowHTI[n].startLongi; mPZ.lat = lowHTI[n].startLati; mPZ.isOK = 0;

		/* ���������ƽ�����ڵĻ����ͷֱ��ʵؾ� *//* ���ɻ������� */
		sprintf(msg, "��ʼ��������%d γ��%d ����Ļ����ؾ�.", mPZ.lon, mPZ.lat);
		MessageToDialog(msg);
		currTi.minLongi = (float)(mPZ.lon);
		currTi.minLati  = (float)(mPZ.lat);
		currTi.maxLongi = (float)(mPZ.lon + 1);
		currTi.maxLati  = (float)(mPZ.lat + 1);
		currPZ = &mPZ;

		/* ����Ҫƴ���;��ȵؾ��ĸ߳����� */
		if(!CreatePerDegreeTerrain(2))
		{
			sprintf(msg, "��������%d γ��%d ����ĸ̳߳���.", mPZ.lon, mPZ.lat);
			MessageToDialog(msg);
			lowHTI[n].isOK = 2;
			continue;
		}

		/* ƴװһ��findResultHTI */
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

		/* װ��vpb���������*/
		MessageToDialog("����VPB���ɵؾ�");
		if(!DoHighVpbAndWriteStg(0))
		{
			MessageToDialog("ִ��VPBʧ��\n");
			lowHTI[n].isOK = 4;
			continue;
		}

		/* ��¼������Ѿ����� */
		lowHTI[n].isOK = 1;

		/* ��allHTIд�������ļ� */
		CheckConfig(1);
	}
	return 1;
}

/*============================================================================================
 *  ��������˳���������һ��߶ȵļ�����̡�
 *============================================================================================ */

/************************************************
 *
 *    ����openedTiles���顣����Ϊ����ؾ���ĳһ��ĸ߳����á������ڱ�����һ����γ���ڵ�������ά�ؾ�����
 *	  ����̵߳ķ����ǣ�������һƬ�ؾ�����ά���ݣ�Ȼ���ڼ���߳�ȷ��ĳһ��(�ɾ��ȡ�γ�Ⱦ���)����һ����ֱ����ؾ�����ײ��⣬
 *    ��ⷵ�صĽ��������һ��ĸ̡߳�
 *
 ************************************************/
 
 /* ����һ���ؾ��ļ��� �̼߳�����rootNode ���� */
int InsertATile(int tile, osg::ref_ptr<osg::Node> node)
{
	OpenTiles ots;
	
	ots.tile = tile;
	ots.node = node;
	openedTiles.push_back(ots);
	
	rootNode->addChild(node);
	return 1;
}


/* ����ǰ�ؾ�Ŀ¼�����еĸ��ؾ��ļ����뵽 �̼߳�����rootNode ���档*/
int InsertTilesInCurrDir(double lon, double lat)
{
	char fullpath[_MAX_PATH], head[24], lowPart[5], iveHead[48];
	int io, dio, dfio, dirIo, j, iio, iveio, iFlag = 0;
	struct _finddata_t sdff;
	struct _finddata_t iveff;

	/* ���Ȼ�ȡ��ǰĿ¼ */
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
	
	/* �������Ŀ¼��Ȼ��������ǰĿ¼����չ��Ϊive�������ļ��� */
	io = chdir(fullpath);
	if(io) return 0;
	
	/* ���ҵ�ǰ·���µ������ļ��� */
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;
	
	/* ��õ�ǰtileID */
	int currTileID = 0;
	if(openedTiles.size()>0) currTileID = openedTiles[0].tile;

	while(!dfio)
	{

		//���������Ŀ¼Ҳ�����ļ�
		if((sdff.name[0] == '.')) goto ITICD_GET_NEXT_SUBDIR;

		//�������Ŀ¼�����뵱ǰ��Ŀ¼
		if(sdff.attrib == _A_SUBDIR)
		{
			_getcwd(fullpath, _MAX_PATH);

			//������Ŀ¼
			dirIo = chdir(sdff.name);
			if(dirIo == -1) goto ITICD_GET_NEXT_SUBDIR;
			
			for(j=0; j<strlen(sdff.name); j++)
			{ 
				if(sdff.name[j] != '_')	head[j] = sdff.name[j];
				else	break;
			}
			head[j++] = 0; lowPart[0] = sdff.name[j]; lowPart[1] = 0;
			tileID = atoi(head);

			/* �����жϵ;��ȵ��εĹ���	*/
			if(strlen(head) < 8)
			{
				sprintf(head, "%s_%s", head, lowPart);
			}
			
			/* �Ƚϵ�ǰĿ¼ͷ�������Ƿ�͵�ǰtileIDһ�������һ����˵�����Ŀ¼��Ч��ֱ�Ӷ�ȡ��һ��Ŀ¼ */
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

			/* ����鲻���κ��ļ� */
			if(isSucess == 0)	{	_findclose(iio);	break;	}
			
			/* �ڽ�����Ŀ¼�µ����еؾ��ļ����뵽��֮ǰ���ȳ�ʼ���������� */
			RemoveAllTileNodes();
			
			iveio = 0;
			while(!iveio)
			{
				/* ���ļ� */
				osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(iveff.name);
				if(node == NULL) goto ITICD_Continue;

				/* �����´򿪵�tile */
				InsertATile(tileID, node);
ITICD_Continue:
				/* Ѱ����һ�� */
				iveio = _findnext(iio, &iveff);
			}
			_findclose(iio);
			
			double elev = GetElevationByLonLat(lon, lat);
			if(elev != -1.0)
			{
				iFlag = 1;
			}
ITICD_OUT_TO_UPPER_DIR:			
			//�˻ص���һ��Ŀ¼
			dirIo = chdir("..");
			
			if(iFlag) break;
		}
ITICD_GET_NEXT_SUBDIR:

		if(iFlag) break;

		/* Ѱ����һ�� */
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);

	if(iFlag == 0)	return 0; else	return 1;
}

/* �� openedTiles �̼߳�������������Ƴ����е����ݡ� */
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

/* ����һ���㣬��������ӦrootNode�ĵؾ��� */
int LoadRootNodeByVertex(osg::Vec3 v)
{
	Carton ct;	ct.x = v[0]; ct.y = v[1]; ct.z = v[2];
	Geodesy gd;
	CartToGeod(&ct, &gd);
	
	/* ֱ�Ӳ��� ����γ��.ive�ļ���*/
	if(!InsertTilesInCurrDir(gd._lon, gd._lat))
		return 0;

	return 1;
}


/**************************************************
 *
 *   ��ȡ�߶ȡ���ȡ�ؾ���ĳһ�㣨�ɾ�γ�Ⱦ������ĸ߶ȣ�������openedTiles�ĵط��н��ܡ�
 *
 **************************************************/
double GetElevation(double lon, double lat, double pre)
{
	double alt = 1.0; 

#if 0
	alt = GetElevationFromHgt(lon, lat);
	return alt;
#endif

	/* ֱ�Ӳ��� ����γ��.ive�ļ�	������Ҳ������ؾ����ͷ���pre����Ϊȱʡ�߶ȡ�*/
	if(rootNode->getNumChildren() == 0)
	{
		if(!InsertTilesInCurrDir(lon, lat))
		{
			alt = GetElevationFromHgt(lon, lat);
			return alt;
		}
	}
	
	/* �����ཻ�߶Σ���һ���ھ���lon��γ��lat�Ĵ�ֱ���߶Ρ�������߶���ؾ����ཻ��*/
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
		/* �����ཻ��  Create intersector  */
		osgUtil::IntersectVisitor intersectVisitor; 
		osg::LineSegment* lineSegment = new osg::LineSegment(start, end); 
		intersectVisitor.addLineSegment(lineSegment); 
		rootNode->accept(intersectVisitor); 
		
		/* ��ȡ�ཻ��� Get intersections */
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

	/* �����ཻ�߶Σ���һ���ھ���lon��γ��lat�Ĵ�ֱ���߶Ρ�������߶���ؾ����ཻ��*/
	Geodesy stg; stg._lon = lon; stg._lat = lat; stg._elevation =  99999999.0;
	Carton  stc;
	GeodToCart(&stg, &stc);
	osg::Vec3d start; start.set(stc.x, stc.y, stc.z); 
	Geodesy eng; eng._lon = lon; eng._lat = lat; eng._elevation =  -99999999.0;
	Carton  enc;
	GeodToCart(&eng, &enc);
	osg::Vec3d end;   end.set(enc.x, enc.y, enc.z); 

	/* �����ཻ��  Create intersector  */
	osgUtil::IntersectVisitor intersectVisitor; 
	osg::LineSegment* lineSegment = new osg::LineSegment(start, end); 
	intersectVisitor.addLineSegment(lineSegment); 
	rootNode->accept(intersectVisitor); 
	
	/* ��ȡ�ཻ��� Get intersections */
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
 *  �����濪ʼ�� CheckTifInTotalNode ����֮ǰ������JPEG�ļ�����ȫ���ؾ������Ĵ󲿷ֹ�������
 *  Ŀǰ���������õ���
 *============================================================================================ */

/************************************************
 *
 *    ת��Tiff��GeoTiff���ⲿ�ֹ�����ʹ��JPEGͼƬ�����ؾ���һ���м䲽�衣JPEG->Tiff, Tiff->GeoTiff
 *
 ************************************************/
typedef struct _TIFFTAG_
{
	short    idx;
	short    type;
	fpos_t   size;
	fpos_t   ptr;
}TiffTag;

/* ��ȡ�������ַ�����Ҫ��0x1A����ַ�Ҳ��ȷ�ض�������*/
void nfread(unsigned char *pBuff, int num, int size, FILE *pInput)
{
	for(unsigned int j = 0; j < num * size; j++)
	{	
		int c = fgetc(pInput);

		/* ����һ���������ʵ����һ������֮�٣�ϵͳ�����и�bug��
		   ������0x1A�ֽڵ�ʱ������ļ����������0x1A��������ݾ�û�б���������
		   �����������޸ġ��Ѹö��Ķ���������2012-04-24
		 */

		if(pInput->_cnt < 0x9f) pInput->_cnt = 0x9f;


		*pBuff++ = (unsigned char)c;
	}
}


/* ��һ��GeoTiff�ļ��и���Geo���Ͷ�ȡGeo������Ϣ��֧�ֳ����ļ� */
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
		sprintf(msg, "���ܶ�ȡ�ļ� %s .", fn);
		MessageToDialog(msg);
		return 0;
	}
	
	/*��ȡ�ļ�����*/
	fseek(pInput, 0, SEEK_END);
	fgetpos(pInput, &iFileSize);
	rewind(pInput);
	
	/*��дHeader ��Length, ���õ�ǰpos��Tag��*/
	iDataPtr = 0;
	nfread((unsigned char *)&header, 4, 1, pInput);
	nfread((unsigned char *)&iDataPtr, 4, 1, pInput);
	
	/*�ж�Tag����ȷλ��*/
	if(iDataPtr == 8)
	{
		fread(&iDataPtr, 8, 1, pInput); Lflag = 1;
	}
	fsetpos(pInput, &iDataPtr);
	
	/*��Tag������Tag����Ϣ��������tags��*/
	if(!Lflag)
	{
		/* һ��С�ļ���Tag����ʽ */
		fread(&tagnum, 1, 2, pInput);
		for(i=0; i<tagnum; i++)
		{
			fread(&tags[i].idx,  2, 1, pInput);
			fread(&tags[i].type, 2, 1, pInput);
			fread(&tags[i].size, 4, 1, pInput);
			fread(&tags[i].ptr,  4, 1, pInput);
		}
	}else{
		/* �����ļ���Tag����ʽ */
		fread(&tagnum, 1, 8, pInput);
		for(i=0; i<tagnum; i++)
		{
			fread(&tags[i].idx,  2, 1, pInput);
			fread(&tags[i].type, 2, 1, pInput);
			fread(&tags[i].size, 8, 1, pInput);
			fread(&tags[i].ptr,  8, 1, pInput);
		}
	}
	
	/* ��tags����������Ҳ����ļ��ж�ȡ��������γ����Ϣ��TIFFTAG_GEOTIEPOINTS���ͣ�������ͼ��ֱ�����Ϣ��TIFFTAG_GEOPIXELSCALE */
	/* ��γ����Ϣ������ ����modeltiepoints�У��ֱ�����Ϣ����������modelpixelscale�С�  */	
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

			/* ��ȡ3*8���ַ� */
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

			/* ��ȡ6*8���ַ� */
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
	
	/* �����ǰ�ļ�������������Ϣ������һ���Ƿ���ͬ�ļ�����tfw�ļ��� */
	if(!flag)
	{
		FILE *tfw;

		/* ��鲢��ͬĿ¼��ͬ�ļ�����tfw�ļ���*/
		sprintf(tfwfn, "%s.tfw", simfn);
		tfw = fopen(tfwfn, "rt");
		if(tfw == NULL)
		{
			char msg[STRING_STANDARD_LEN];
			sprintf(msg, "�ļ� %s �в�����������Ϣ��", fn);
			MessageToDialog(msg);
			return 0;
		}
		
		/* ��ȡtfw�ļ� */
		ReadALineToArray(tfw);	modelpixelscale[0] = atof(cf[0]);
		ReadALineToArray(tfw);
		ReadALineToArray(tfw);
		ReadALineToArray(tfw);	modelpixelscale[1] = -(atof(cf[0]));
		ReadALineToArray(tfw);	modeltiepoints[3]  = atof(cf[0]);
		ReadALineToArray(tfw);	modeltiepoints[4]  = atof(cf[0]);
		fclose(tfw);
	}

	/* �������������飬�������GeoTiff�ļ��������ĵ���Χ(�����С��γ��)����¼�������б��С� */
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

		/* �����γ�����ݴ���360.0, �п���ʹ�õ�����ĵ�λ��Ҫ��ԭ�����ֻ����� ���� 3600.0 */
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

		/* �����γ�����ݴ���360.0, �п���ʹ�õ�����ĵ�λ��Ҫ��ԭ�����ֻ����� ���� 3600.0 */
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

/* ������Tiff�ļ�ת��ΪGeoTiff�ļ���������ص�����Ϣ�ṹmergedTif�������Ƕ�ȡԭ��Tiff�����Geo���Tags�� */
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
	/*���ļ�*/
	sprintf(temTextGeoFile  , "%s\\%s%d%d.%s",  tempDirectory, TEMP_TEXTGEOFILE, currPZ->lon, currPZ->lat, TEMP_EXTENDNAME);
	if((pInput=fopen(temTextFile,"rb"))==NULL||						
		(pOutput=fopen(temTextGeoFile,"wb"))==NULL)
	{
		MessageToDialog("���ļ�ʧ��\n");
		currPZ->isOK = STATUS_TEXTGEOBAD;
		return 0;
	}
	
	/*��ȡ�ļ�����*/
	fseek(pInput, 0, SEEK_END);
	fgetpos(pInput, &iFileSize);
	rewind(pInput);
	
	/*��дHeader ��Length*/
	fread(&header, 4, 1, pInput);	fwrite(&header, 4, 1, pOutput);
	fread(&len, 4, 1, pInput);		fwrite(&len, 4, 1, pOutput);
	
	/*��д�ļ���*/
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

	/*��Tag��*/
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
	
	/* ׼��Geo���� */
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
	
	/*׼��Tag����*/
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
	
	/*дTag��*/
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
	
	/*�ر��ļ�*/
	delete pcInputLine;
	fclose(pInput);
	fclose(pOutput);
#endif	
	return 1;	
}

/* ���߳�Tiff�ļ�ת��ΪGeoTiff�ļ���������ص�����Ϣ�ṹmergedHgt�� */
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
	/*���ļ�*/
	sprintf(temHgtGeoFile   , "%s\\%s%d%d.%s",  tempDirectory, TEMP_HGTGEOFILE	, currPZ->lon, currPZ->lat, TEMP_EXTENDNAME);
	if((pInput=fopen(temHgtFile,"rb"))==NULL||						
		(pOutput=fopen(temHgtGeoFile,"wb"))==NULL)
	{
		MessageToDialog("���ļ�ʧ��\n");
		currPZ->isOK = STATUS_HGTGEOBAD;
		return 0;
	}
	
	/*��ȡ�ļ�����*/
	fseek(pInput, 0, SEEK_END);
	fgetpos(pInput, &iFileSize);
	rewind(pInput);
	
	/*��дHeader ��Length*/
	fread(&header, 4, 1, pInput);	fwrite(&header, 4, 1, pOutput);
	fread(&len, 4, 1, pInput);		fwrite(&len, 4, 1, pOutput);
	
	/*��д�ļ���*/
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

	/*��Tag��*/
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
	
	/* ׼��Geo���� */
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

	/* ����Geo���� */
	double modeltiepoints[] = {0.0,0.0,0.0,0.0,0.0,0.0};
	double modelpixelscale[] = {0.0,0.0,1.0};

	int m_unImgWidth = tags[0].ptr;
	int m_unImgLength = tags[1].ptr;

	modeltiepoints[3] = mergedHgt.minLong;
	modeltiepoints[4] = mergedHgt.maxLati;
	
	modelpixelscale[0] = (mergedHgt.maxLong - mergedHgt.minLong) / m_unImgWidth;
	modelpixelscale[1] = (mergedHgt.maxLati - mergedHgt.minLati) / m_unImgLength;
	
	/*׼��Tag����*/
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
	
	/*дTag��*/
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
	
	/*�ر��ļ�*/
	delete pcInputLine;
	fclose(pInput);
	fclose(pOutput);
#endif
	return 1;	
}


/***********************************************
 *
 *   �ϲ�ͼ�񡣺ϲ����ɸ�JPEGͼ���Ϊһ��Tiffͼ����ʹ��JPEGͼ�������ؾ��Ĺؼ����衣
 *   ���������еķ�ʽ���롣�Ƚ�Ҫ�ϲ�������ͼ�񣬺ϲ���ʽ���Լ�������ļ������һ���������顣Ȼ�����������̡�
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
  MessageToDialog("����֧��JPG�ļ�ƴ�ӳ�Tif�ļ��Ĺ��ܡ�");
  return 0;
#endif
}


/**********************************************
 *
 *    HGT  ����߳��ļ�
 *
 ***********************************************/

int maxhgt, minhgt;

/* ����һ���ڵ�ĸ߳�ֵ */
short height(TGhgt *hgt, int x, int y ) 
{ 
    return hgt->data[x][y]; 
}

/* ���һ���߳��������Ƿ��з�0�ĸ߳�ֵ */
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

/* ���һ���ļ����Ƿ���zip�ļ� */
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

/* ����һ���µĸ߳��ļ�������߳�����ȫΪ0�� �������û�б����á�*/
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

/* ��һ���߳��ļ������Ұ��м�Ŀհ׵�Ĩƽ��Ĩƽ�����Ǹ����Ѿ���������һ���Լ����һ����ȷ����  */
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
    
    /* Ŀǰ��֧��ѹ���ļ���ѹ������ϵͳ�߳��ļ����ǽ�ѹ������ݡ�*/
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

	/* �򿪸߳��ļ��� */    
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

	/* Ŀǰֻ֧�ֵͷֱ��ʸ߳����ݣ���hgt_resolution == 3���޸�MAX_HGT_SIZE�����󣬿�֧�ָ߷ֱ��ʸ߳����ݣ���hgt_resolution == 1 */
    /*load data from file */
    if ( hgt->hgt_resolution == 1 ) {
        hgt->cols = hgt->rows = size = 3601;
        hgt->col_step = hgt->row_step = 1;
    } else if ( hgt->hgt_resolution == 3 ) {
        hgt->cols = hgt->rows = size = 1201;
        hgt->col_step = hgt->row_step = 3;
    }

	/* ����ȡ�߳����ݣ��ж��Ƿ��ǿհ׵㣬����ǿհ׵���Ĩƽ���ݡ�*/    
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

			//ȥ���հ׵�
			cu_arr = *var;
			if((cu_arr == -32768)&&(row<(hgt->rows - 1)))
				cu_arr = hgt->data[col][row+1];

			if((cu_arr == -32768)&&(row==(hgt->rows - 1)))
				cu_arr = hgt->data[col-1][row];

			if((cu_arr > -32768)&&(cu_arr < 0))
				cu_arr = 0;
			*var = cu_arr;

			//�ҳ������Сֵ
			if(*var > maxhgt) maxhgt = *var;
			if(*var < minhgt) minhgt = *var;
		}
    }
    
    fclose(hgt->fd);
    
    return 1;
}

/**********************************************
 *
 *   ���ݾ�γ�ȷ�Χ��ø߳������ļ����У��ϲ��߳��ļ�����תΪGeoTiff�ļ���
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

/* ��ȡָ����γ����������и߳��ļ������б� */
bool GetHgtFileList(double leftlongi, double downlati, double rightlongi, double upperlati)
{
	int iMaxLong, iMaxLati, iMinLong, iMinLati;
	short iLong, iLati, m, n;
	
	// ��ʼ�����ֱ���
	memset(hgtarray, 0, sizeof(hgtarray));
	maxhgt = -999; minhgt = 999999;
	iMaxLong = (int)rightlongi;
	iMaxLati = (int)upperlati;
	iMinLong = (int)leftlongi;
	iMinLati = (int)downlati;
	
	//��whileѭ����������
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

/* ��ʼ���߳̽ṹ��Ϣ */
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

/* ���ռ���ó��ĸ߳����񣬺ϲ���ظ߳��ļ�Ϊһ���ļ� */
bool MergeHgt()
{
	int n, m, i, j;
	TIFF * pOutput;		//����ļ�ָ��
	char fullPath[STRING_STANDARD_LEN];

	//д�򿪸߳��ļ�
	if((pOutput=TIFFOpen(temHgtFile,"w"))==NULL)
	{
		MessageToDialog("���ļ�ʧ��.\n");
		return false;
	}

	//ֻʹ����һ��hgt
	if(hgtrow == 1)
	{
		//���뵥���е��ڴ�
		int nOutputLineSize = (MAX_HGT_SIZE) * hgtcol;
		float *pcOutputLine=new float[nOutputLineSize];
		if(pcOutputLine == NULL)
		{
			MessageToDialog("û���㹻���ڴ�.\n");
			return false;
		}
		
		
		//�򿪸߳��ļ�������߳�����
		for(i=0;i<hgtcol;i++)
		{
			TGhgt * curHgt = new TGhgt;
			if(curHgt == NULL)
			{
				MessageToDialog("û���㹻���ڴ�.\n");
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

		//д��Geotiff�ļ�����
		InitHgtInfo(pOutput, hgtcol * hgtarray[0][0].pHgt->cols, hgtarray[0][0].pHgt->rows);
	
		//���кϲ�д�롣
		for(m=0; m<MAX_HGT_SIZE; m++)
		{
			//���Ƚ�����hgt���ݺϲ���һ��
			for(i=0;i<hgtcol;i++)
			{
				for(j=0; j<(MAX_HGT_SIZE); j++)
					pcOutputLine[i*(MAX_HGT_SIZE) + j] = (float)(hgtarray[0][i].pHgt->data[j][MAX_HGT_SIZE - 1 - m]);
			}
			
			//����һ��д���ļ�
			if(TIFFWriteScanline(pOutput,pcOutputLine,m)!=1)
			{
				MessageToDialog("���ܹ���ȷ��д���ļ�!");
				currPZ->isOK = STATUS_HGTBAD;
				return false;
			}
		}
		
		//��ȡ���Ϻ�߳������Ϣ
		mergedHgt.minLong = hgtarray[0][0].usLong;
		mergedHgt.maxLati = hgtarray[0][0].usLati + 1;
		mergedHgt.maxLong = hgtarray[0][hgtcol - 1].usLong + 1;
		mergedHgt.minLati = hgtarray[0][hgtcol - 1].usLati;
		sprintf(mergedHgt.pathName, "%s", tempDirectory);
		sprintf(mergedHgt.fileName, "%s.%s", TEMP_HGTFILE, TEMP_EXTENDNAME);
		mergedHgt.image_width = hgtcol * hgtarray[0][0].pHgt->cols;
		mergedHgt.image_height = hgtarray[0][0].pHgt->rows;
		
		//�ͷ��ڴ�
		delete pcOutputLine;
		for(i=0;i<hgtcol;i++)
			if(hgtarray[0][i].pHgt) delete hgtarray[0][i].pHgt;
		
	}else{		//ʹ���˶���

		//���뵥���е��ڴ�
		int nOutputLineSize = (MAX_HGT_SIZE) * hgtcol;
		float *pcOutputLine=new float[nOutputLineSize];
		if(pcOutputLine == NULL)
		{
			MessageToDialog("û���㹻���ڴ�.\n");
			currPZ->isOK = STATUS_HGTBAD;
			return false;
		}
				
		//�򿪸߳��ļ�������߳�����
		for(j=0;j<hgtrow;j++)
			for(i=0;i<hgtcol;i++)
			{
				TGhgt * curHgt = new TGhgt;
				if(curHgt == NULL)
				{
					MessageToDialog("û���㹻���ڴ�.\n");
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
	
		//д��Geotiff�ļ�����
		InitHgtInfo(pOutput, hgtcol * hgtarray[0][0].pHgt->cols, hgtrow * hgtarray[0][0].pHgt->rows);

		//���кϲ�д�롣
		for(n=hgtrow - 1; n>=0; n--)
		{
			for(m=0; m<MAX_HGT_SIZE; m++)
			{
				//���Ƚ�����hgt���ݺϲ���һ��
				for(i=0;i<hgtcol;i++)
				{
					for(j=0; j<(MAX_HGT_SIZE); j++)
						pcOutputLine[i*(MAX_HGT_SIZE) + j] = (float)(hgtarray[n][i].pHgt->data[j][MAX_HGT_SIZE - 1 - m]);
				}
				
				//����һ��д���ļ�
				if(TIFFWriteScanline(pOutput,pcOutputLine,(hgtrow - 1 - n)*MAX_HGT_SIZE + m)!=1)
				{
					MessageToDialog("���ܹ���ȷ��д���ļ�!");
					currPZ->isOK = STATUS_HGTBAD;
					return false;
				}
			}
		}

		//��ȡ���Ϻ�߳������Ϣ
		mergedHgt.minLong = hgtarray[0][0].usLong;
		mergedHgt.maxLati = hgtarray[hgtrow - 1][hgtcol - 1].usLati + 1;
		mergedHgt.maxLong = hgtarray[hgtrow - 1][hgtcol - 1].usLong + 1;
		mergedHgt.minLati = hgtarray[0][0].usLati;
		sprintf(mergedHgt.pathName, "%s", tempDirectory);
		sprintf(mergedHgt.fileName, "%s.%s", TEMP_HGTFILE, TEMP_EXTENDNAME);
		mergedHgt.image_width = hgtcol * hgtarray[0][0].pHgt->cols;
		mergedHgt.image_height = hgtrow * hgtarray[0][0].pHgt->rows;
		
		//�ͷ��ڴ�
		delete pcOutputLine;
		for(j=0;j<hgtrow;j++)
			for(i=0;i<hgtcol;i++)
				if(hgtarray[j][i].pHgt) delete hgtarray[j][i].pHgt;
		
	}
	
	//�ر��ļ�		
	TIFFClose(pOutput);
	return true;
}


/************************************************
 *
 *    TravelMap
 *    �����ؾ�����Ŀ¼����Ѱ.map�ļ�����ȡ��Ϣ�������ݿ�
 *    .map�ļ�����.jpg�ļ����ļ���һһ��Ӧ���ı��ļ������汣����.jpg�ļ������������ĵ�����Ϣ����γ�ȡ�
 *    ���Ǵ�GoogleEarth�����صĵؾ�ԭʼ���ݵ��ص㡣   
 *
 **********************************************************/

/* ��һ���ı��ļ������ȡһ���ַ����������Կո�Ϊ�ֽ磬����cf��ά�����У����Ǹ��ܳ��õĹ��̡� */
/* ������ȡһ�����ݷ���1����ȡ���з���2���ļ���������0*/
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
		/* ����� ��//�� ��ʼ���У�Ҳ�������д������� 2 --- ����ע�͵Ĺ��ܡ� */
		if((OneLine[0] == '/') && (OneLine[1] == '/')) contiFlag++;
	
		/* �жϱ����Ƿ�ȫ���ǿո����Tab */
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

	/* ��ȡ��������*/
	j0=j1=d0=0;
	for(i=0;i<(int)strlen(OneLine);i++)
	{
		if(!((OneLine[i]==' ')||(OneLine[i]==0x0a)||(OneLine[i]==0x09)))  cf[j0][j1++] = OneLine[i];
		if(((OneLine[i]==' ')||(OneLine[i]==9))&&((OneLine[i-1]!=' ')&&(OneLine[i-1]!=9))&&(i!=0)) {j0++; j1=0;}
	}
	return 1;
}



/* ������ǰĿ¼�µ�������Ŀ¼���Բ��� map �ļ���ʹ�õݹ鷽ʽ�� */
int TravelSubDirectoryForMap(char *dir)
{
	int dio, dfio, dirIo, i, j, k;
	struct _finddata_t sdff;

	//���ҵ�ǰ·���µ���Ŀ¼��
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//���������Ŀ¼Ҳ�����ļ�
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR;

		//�������Ŀ¼�����뵱ǰ��Ŀ¼
		if(sdff.attrib == _A_SUBDIR)
		{
			//�ȷ������Ŀ¼��
			_getcwd(path, _MAX_PATH);
			CheckTextureDirName(sdff.name, path);
			
			//������Ŀ¼
			dirIo = chdir(sdff.name);
			_getcwd(path, _MAX_PATH);
			
			//������Ŀ¼
			TravelSubDirectoryForMap(path);
			
			//�˻ص���һ��Ŀ¼
			dirIo = chdir("..");
		}
		
		//������ļ�
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//�ҵ�map�ļ�
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='M')||(sdff.name[i-3]=='m')) && ((sdff.name[i-2]=='A')||(sdff.name[i-2]=='a')) && ((sdff.name[i-1]=='P')||(sdff.name[i-1]=='p')))  
			{
				MapNode mn;
				//��map�ļ�
				fi = fopen(sdff.name, "rt");
				if(fi == NULL) 
				{
					char msg[STRING_STANDARD_LEN];
					sprintf(msg, "���ܶ�ȡ�ļ� %s .", sdff.name);
					MessageToDialog(msg);
					return 0;
				}
				
				ReadALineToArray(fi);
				ReadALineToArray(fi);	j = strlen(cf[0]); strcpy(mn.fileName, cf[0]);
				
				//��ȡͼƬ�Ŀ�����Ϣ
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
				
				/* ��ȡͼƬ�ĳ�����Ϣ������ԭʼͼƬ��С��һ�������Զ�ȡ������Ҫ��ԭ���ݻ����� ���� 2 */
				mn.image_width  = atoi(cWidth) / 2;
				mn.image_height = atoi(cHeigh) / 2;
				
				//��ȡͼƬ�ĵ�����Ϣ
				while((strcmp("MMPLL,1,", cf[0]))&&(!feof(fi)))  ReadALineToArray(fi);
				if(feof(fi))
				{
					MessageToDialog("��ȡmap�ļ�����:�Ҳ�����γ����Ϣ��");
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
		//Ѱ����һ��
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* ������ǰĿ¼�����Ҹ�Ŀ¼���Լ�������Ŀ¼����չ��Ϊ.map���ļ�����ȡ���������Ϣ�� */
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
 *  �����濪ʼ�� CheckTheSortList ����֮ǰ������Tiff�ļ����ɵؾ����������������������ӹ�����ǰ�档�����û���Ҫ֧�ֵĸ�ʽ��
 *============================================================================================ */

/**********************************************************
 *
 *       �����ؾ�����Ŀ¼����Ѱ.tif�ļ�����ȡ��Ϣ�������ݿ�
 *       
 **********************************************************/
/* ���û��߷ֱ�������(highTxtnode)���߳�(highHgtnode)�б�������ҵ�ǰ�ļ����ڲ��ڣ�*/ 
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
 
/* ����Ŀ¼��������tif�ļ� */
/* ������ǰĿ¼�µ�������Ŀ¼���Բ��� tif �ļ���ʹ�õݹ鷽ʽ�� */
int TravelSubDirectoryForTif(char *dir)
{
	int dio, dfio, dirIo, i;
	struct _finddata_t sdff;

	//���ҵ�ǰ·���µ���Ŀ¼��
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//���������Ŀ¼Ҳ�����ļ�
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR3;

		//�������Ŀ¼�����뵱ǰ��Ŀ¼
		if(sdff.attrib == _A_SUBDIR)
		{
			//������Ŀ¼�����ƣ�ȷ��������
			if(!strcmp("1m", sdff.name))  currHighLevel = 1;
			if(!strcmp("5m", sdff.name))  currHighLevel = 5;
			if(!strcmp("10m", sdff.name)) currHighLevel = 10;
			if(!strcmp("50m", sdff.name)) currHighLevel = 50;

			//������Ŀ¼
			dirIo = chdir(sdff.name);
			_getcwd(path, _MAX_PATH);
			
			//������Ŀ¼
			TravelSubDirectoryForTif(path);
			
			//�˻ص���һ��Ŀ¼
			dirIo = chdir("..");
		}
		
		//������ļ�
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//�ҵ�Tif�ļ�
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='T')||(sdff.name[i-3]=='t')) && ((sdff.name[i-2]=='I')||(sdff.name[i-2]=='i')) && ((sdff.name[i-1]=='F')||(sdff.name[i-1]=='f')))  
			{
				/* �Ȳ������ļ��Ƿ��Ѿ��ڿ������ˡ��Ӻ���ǰ��*/
				if(!CheckTifInTotalNode(sdff.name))
				{					
					/* ��GeoTiff�ļ�����ȡTags��Ϣ */
					if(!ReadGeoTiffTags(sdff.name))
					{
						return 0;
					}	
				}
			}
		}

GET_NEXT_SUBDIR3:
		//Ѱ����һ��
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* ������ǰĿ¼�����Ҹ�Ŀ¼���Լ�������Ŀ¼����չ��Ϊ.tif���ļ�����ȡ���������Ϣ�� */
int TravelDirectoryToReadTif(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	/* �ݹ������ǰĿ¼�����е���Ŀ¼ */
	if(!TravelSubDirectoryForTif(dPath))
	{
		MessageToDialog("�Ҳ�����������Ϣ�������ļ�");
		return 0;
	}
	
	/* ��highTxtnode���� �ֱ���Ϊ50�׵����� ת�Ƶ� lowTxtnode ����*/ 
	if(!MoveHigh50ToLow())	return 0;

	/* ���ݶ�ȡ���û��߷ֱ���Tif�ļ���Ϣ�������ɸ߷ֱ��ʵؾ������б�allHTI��
	   �����ۺϿ��������еĸ�����������������߷ֱ������ݸ��ǵͷֱ������ݡ� */
	if(!ReadNodeListToCreateHTI()) 	{	return 0;	}

	/* ���ݷֱ���Ϊ50�׵�ϵͳ�������ݣ������ɵͷֱ��ʵؾ������б�lowHTI��
	   ���ڵͷֱ��ʵؾ�ϵͳ�Ѿ��ṩ�����������ʵ�����塣*/
	if(!GetHTIListOfR50())	{	return 0;	}
	
	return 1;
} 

/* ��highTxtnode���� �ֱ���Ϊ50�� ��ת�Ƶ� lowTxtnode ����*/ 
int MoveHigh50ToLow()
{
	int i, j;
	
	if(lowTxtnode.size() == 0)
	{
		lTxtNum = 0;
		/* �Ƚ�highTxtnode����ֱ���Ϊ50�� ��ת�Ƶ� lowTxtnode ����*/
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
 	
 	/* Ȼ���Ƴ�highTxtnode����ֱ���Ϊ50�׵�ֵ */
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
 *  �����濪ʼ�� AppendStgAndWrite ����֮ǰ�������ɸ�Jpeg�ļ�ƴ�Ӻϲ���һ��Tiff�ļ��Ĺ��̡������װ�����С�
 *============================================================================================ */

/**********************************************************
 *
 *  ���������С��JpgͼƬƴ��Ϊһ������ϴ��TifͼƬ��TifͼƬ���ٰ���һ������xһ��γ�ȷ�Χ�ڵĵؾ���
 *       
 **********************************************************/

/* ����������Ƿ���ȷ */
int CheckTheSortList(int *curLine)
{
	int i, iFlag = -1;
	double curlineLon, prelineLon;
	
	/* ���ȼ�鼸�е���Ŀ���Ƿ���ͬ */
	for(i=1; i<sortrows; i++)
	{
		if(sortcolums[i] != sortcolums[i-1])
		{
			if(sortcolums[i] < sortcolums[i-1]) {	iFlag = i;   break; }
			if(sortcolums[i] > sortcolums[i-1]) {	iFlag = i-1; break; }
		}
	}
	*curLine = iFlag;
	
	/* ������ֲ���ͬ���������߾����Ƿ�ӽ� */
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

/* ����ָ�������ڰ������ٸ�jpg*/
#define LIMIT_WIDTH 0.01
int FindJpgsInZoneAndSort(double minlon, double maxlon, double minlat, double maxlat)
{
	int i, j, badLine;
	double dtmp1, dtmp2;
	
	/* ��ͼƬ����������ҷ���������ͼƬ */
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

	/* ����ҵ�������ͼƬ����鲢���һ�� */
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

	/* �������� */
	if(!SortfindResult())
	{
		MessageToDialog("����ͼƬʧ��\n");
		currPZ->isOK = STATUS_TEXTBAD;
		return 0;
	}

	/* ����Ƿ��д�λ����������ͼƬ��Ŀ���ȡ�*/
	int Flag;
	if((Flag = CheckTheSortList(&badLine)) != -1)
	{
		/* �ҵ���Ҫ���ӵ�ͼƬ����� */
		int needID;
		if(!Flag)
		{
			needID = findResult[sortlist[badLine][0]].ID - 1;
		}
		if(Flag)
		{
			needID = findResult[sortlist[badLine][sortcolums[badLine] - 1]].ID + 1;
		}
		
		/* ������ͼƬ*/
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
		
		/* �������� */
		if(!SortfindResult())
		{
			MessageToDialog("����ͼƬʧ��\n");
			currPZ->isOK = STATUS_TEXTBAD;
			return 0;
		}	
	}


	/* ����Ե���Ƿ񵽱� */
	/* �ȼ����Сγ�� */
	for(i=resultID-1; i>=0; i--)
	{
		double delta = fabs(findResult[i].minLati - currTi.minLati);
		if((delta < 0.1) && (findResult[i].minLati > currTi.minLati))
			findResult[i].minLati = currTi.minLati;
		else
			break;
	}

	/* ������γ�� */
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


/* ����һ��bucket��������˶��ٸ�jpg. */
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


//��Result����������ж�ά����д���ά�������棬����γ�ȴӸߵ��ͣ����ͬһγ�ȴӵ͵��ߵ�˳��
int SortfindResult()
{
	int i, j, k, m, n, ex;
	int srclist[MAXFINDNODES];
	double maxlat = -99.9;
	
	//��ʼ������
	memset(sortlist, -1, sizeof(sortlist)); 
	memset(sortcolums, -1, sizeof(sortcolums));
	memset(srclist, -1, sizeof(srclist));

	//��srclist����ʼֵ��
	for(i=0; i<resultID; i++)  srclist[i] = i;
	
	//�Ȱ�γ�ȴӸߵ�������
	m = 0;
	while(1)
	{
		//����������飬Ѱ�����γ��
		maxlat = -99.9;
		for(i=0; i<resultID; i++)  
		{
			if((findResult[i].minLati > maxlat)&&(srclist[i] != -1)) 
				maxlat = findResult[i].minLati;
		}
		if( maxlat == -99.9) break;

		//����maxlat��ȵ�ֵ��idx���ڵ�ǰ�±������С�
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

	//ͬγ�Ȱ����ȴӵ͵�������
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

//�����ͼ��ϲ���һ��
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
		/* ��װ����ͼƬ�� */
		for(n=0;n<sortrows;n++)
		{
			sprintf(msg, "��ʼ�� %d -- %d �κ���ƴ��...", n+1, sortrows);
			MessageToDialog(msg);

			for(i=0;i<sortcolums[n];i++)
			{
				if(n==0)	width += findResult[sortlist[n][i]].image_width;
			}
			height = findResult[sortlist[n][0]].image_height;

			//��ȡ���Ϻ�ͼƬ�����Ϣ
			mergedTif.minLong = findResult[sortlist[n][0]].minLong;
			mergedTif.maxLati = findResult[sortlist[n][0]].maxLati;
			mergedTif.maxLong = findResult[sortlist[n][i-1]].maxLong;
			mergedTif.minLati = findResult[sortlist[n][i-1]].minLati;
			sprintf(mergedTif.pathName, "%s", tempDirectory);
			sprintf(mergedTif.fileName, "%s%d.jpg", TEMP_TEXTFILE, n);
			mergedTif.image_width = width;
			mergedTif.image_height = height;
	
			/* ��װ������
			/* convert.exe�̶�ֵ������Ĳ����ɱ䡣+append����ƴ��	-append����ƴ�� ������Ϊ���ͼƬ
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
				MessageToDialog("ƴ���ļ�ʧ��\n");
				currPZ->isOK = STATUS_TEXTBAD;
				return 0;
			}
		}

		/* �����ܵ�ͼƬ�߶� */
		sprintf(msg, "��ʼ����ƴ��...");
		MessageToDialog(msg);
		height = 0;
		for(j=0; j<sortrows; j++)
			height += findResult[sortlist[j][0]].image_height;

		/* ��ȡ���Ϻ�ͼƬ�����Ϣ */
		mergedTif.minLong = findResult[sortlist[0][0]].minLong;
		mergedTif.maxLati = findResult[sortlist[0][0]].maxLati;
		mergedTif.maxLong = findResult[sortlist[sortrows-1][sortcolums[sortrows-1]-1]].maxLong;
		mergedTif.minLati = findResult[sortlist[sortrows-1][sortcolums[sortrows-1]-1]].minLati;
		sprintf(mergedTif.pathName, "%s", tempDirectory);
		sprintf(mergedTif.fileName, "%s.%s", TEMP_TEXTFILE, TEMP_EXTENDNAME);
		mergedTif.image_width = width;
		mergedTif.image_height = height;
		
		/* ����һ��ͼƬ�����Ϣ�����ֲ��ڷ� */
		if(mergedTif.minLong > currTi.minLongi)
		{
				mergedTif.minLong = currTi.minLongi;
		}
		if(mergedTif.maxLong < currTi.maxLongi)
		{	
				mergedTif.maxLong = currTi.maxLongi;
		}

		/* ��װ����ͼƬ*/
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
			MessageToDialog("ƴ���ļ�ʧ��\n");
			currPZ->isOK = STATUS_TEXTBAD;
			return 0;
		}	
	}else{
		//sortrows = 0�������û�в鵽�κ�ͼƬ��
		MessageToDialog("û���ҵ�����ͼƬ��\n");
		currPZ->isOK = STATUS_TEXTBAD;
		return 0;
	}
	
	sprintf(msg, "ͼ��ƴ�����...");
	return 1;
}

/*============================================================================================
 *  �����濪ʼ�� DoVpb ����֮ǰ������stg�ļ��Ĺ��̡�stg�ļ����Ӿ�ϵͳ����������ļ���
 *  �Ӿ�ϵͳ�ڶ�ȡ��ά�ؾ�����ʱ���ȶ�ȡstg�ļ������������������Ϣ����ȡ��ص���ά�ؾ��ļ���
 *============================================================================================ */

/**********************************************************
 *
 *  ��stgβ��׷�������С��������ɵĸߵͷֱ��ʵؾ��ļ������ض���ʽ׷�ӵ�stg�ļ����档
 *  �ظ��ľͲ�����׷�ӡ�
 *       
 **********************************************************/

/* ׷��ctileID.stg�ļ������û�о��½���������mainд����� */
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
			MessageToDialog("���ܴ����ؾ�stg�ļ���\n");
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
	
	/* ����һ���ļ����ݣ��Ƿ��Ѿ�д��ͬ�������ݡ�*/
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

/* ����tileid.stg�ļ�*/
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

	//����Bucket, ����Ŀ¼�ṹ
	gen_base_path(Bt, dstPath);
	index = gen_index(Bt);
	sprintf(ctileID, "%d", index);

	//��õ�ǰ����Ŀ¼
	olddrv = _getdrive();
	_getcwd(path, _MAX_PATH);

	//�Ƚ��뵽�����Ŀ¼
	drv = iveroot[0] - 'A' + 1;
	dio = _chdrive(drv);
	if(dio)
	{
		MessageToDialog("���ܽ���ؾ����ݿ�Ŀ¼�����������С�\n");
		retval = 0; goto MSFOut;
	}

	dio = chdir(iveroot);
	if(dio)
	{
		MessageToDialog("���ܽ���ؾ����ݿ�Ŀ¼�С�\n");
		retval = 0; goto MSFOut;
	}

	//���Ҳ�������Ŀ¼
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
			MessageToDialog("���ܴ����ؾ���Ŀ¼��\n");
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
			MessageToDialog("���ܴ����ؾ���Ŀ¼��\n");
			retval = 0; goto MSFOut;
		}
		dio = chdir(secPath);
	}
	
	//����stg�ļ���д��
	retval = AppendStgAndWrite(ctileID, main, type);
	
	//�˻ظ�Ŀ¼
	dio = chdir("..");
	dio = chdir("..");

MSFOut:
	//�˻ص�ԭ���Ĺ���Ŀ¼
	dio = _chdrive(olddrv);
	if(dio)
	{
		MessageToDialog("���ܽ��빤��Ŀ¼�����������С�\n");
		return 0;
	}

	dio = chdir(path);
	if(dio)
	{
		MessageToDialog("���ܽ��빤��Ŀ¼�С�\n");
		return 0;
	}
	
	return retval;
}

/* ���STG�ļ��Ƿ���ڡ�*/
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

	//����Bucket, ����Ŀ¼�ṹ
	gen_base_path(Bt, dstPath);
	index = gen_index(Bt);
	sprintf(ctileID, "%d", index);

	//��õ�ǰ����Ŀ¼
	olddrv = _getdrive();
	_getcwd(path, _MAX_PATH);

	//�Ƚ��뵽�����Ŀ¼
	drv = iveroot[0] - 'A' + 1;
	dio = _chdrive(drv);
	if(dio)
	{
		MessageToDialog("���ܽ���ؾ����ݿ�Ŀ¼�����������С�\n");
		retval = 0; goto CSFOut;
	}

	dio = chdir(iveroot);
	if(dio)
	{
		MessageToDialog("���ܽ���ؾ����ݿ�Ŀ¼�С�\n");
		retval = 0; goto CSFOut;
	}

	//���Ҳ�������Ŀ¼
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
			MessageToDialog("���ܴ����ؾ���Ŀ¼��\n");
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
			MessageToDialog("���ܴ����ؾ���Ŀ¼��\n");
			retval = 0; goto CSFOut;
		}
		dio = chdir(secPath);
	}
	
	//�ж�stg�ļ��Ƿ���ڲ���
	sprintf(dstFile, "%s.stg", ctileID);
	fi = fopen(dstFile, "rt");
	if(fi == NULL)
	{
		MessageToDialog("��ǰ��û��stg�ļ���\n");
		retval = 0; goto CSFOut;
	}
	
	//�ж�stg�ļ��������Ƿ���� tileID.ive �����䡣
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
	
	//�˻ظ�Ŀ¼
	dio = chdir("..");
	dio = chdir("..");

CSFOut:
	//�˻ص�ԭ���Ĺ���Ŀ¼
	dio = _chdrive(olddrv);
	if(dio)
	{
		MessageToDialog("���ܽ��빤��Ŀ¼�����������С�\n");
		return 0;
	}

	dio = chdir(path);
	if(dio)
	{
		MessageToDialog("���ܽ��빤��Ŀ¼�С�\n");
		return 0;
	}
	
	return retval;
}

/* �жϸÿ��Ƿ��ڱ߽� */
int DetectBucket(double lon, double lat)
{
	if(	( fabs(lon - currTi.minLongi) < 0.01 )	||
		( fabs(lat - currTi.minLati)  < 0.01 )	||
		( fabs(lon+0.25 - currTi.maxLongi) < 0.01 )	||
		( fabs(lat+0.125- currTi.maxLati)  < 0.01 )	
	) return 1;
	else return 0;
}                                                 	

/* ��������дSTG�ļ� , type==0����stg��д��ive; type==1׷��stg��д��btg*/
int WriteAllStgs(char *itemStr, int type)
{
	double src_lon;
	double src_lat;
	Bucket b;
	char   item[32];
	
	/* ����һ���������itemStr, ȥ����������չ��������item���� */
	memset(item, 0, sizeof(item));
	for(int i = 0; i < strlen(itemStr); i++)
	{
		if(itemStr[i] == '.') break;
		item[i] = itemStr[i];
	}
	
	/* ���stg�Ƿ��Ѿ�д��ͬ���ļ������û�о�д��*/
	src_lon = currTi.minLongi;
	src_lat = currTi.minLati;
	while(src_lon < currTi.maxLongi)
	{
		while(src_lat < currTi.maxLati)
		{			
			/*�������ǰ���Bucket */
			set_bucket(&b, src_lon, src_lat);
		
			if(DetectBucket(src_lon, src_lat))
			{
				if(!MakeSTGFile(&b, item, type))
				{
					MessageToDialog("����STG�ļ�ʧ��\n");
					return 0;
				}
			}

			/*������һ��tile*/
			src_lat += 0.125;
		}
		src_lon += 0.25;
		src_lat = currTi.minLati;
	}
	return 1;
}

/*============================================================================================
 *  �����濪ʼ�� CheckConfig ����֮ǰ��������vpb�����У�������vpb���ɵؾ��Ĺ��̡�
 *============================================================================================ */

/**********************************************************
 *
 *  ִ��VPB��������ive�ļ�
 *       
 **********************************************************/

int osgdem_main(int argc, char** argv);

/* ִ��VPB��������ive�ļ� */
int DoVpb(Bucket *b, int type)
{
	char arg[30][STRING_STANDARD_LEN];
	char *argptr[30];
	long tile;
	char basePath[STRING_STANDARD_LEN];
	char fullPath[_MAX_PATH];
	float W, S;
	
	/* �������·����Ϣ */
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
		
    /* ���������ʽƴ�������� ��*/
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

	/* ����osgdem�����̡�*/
	if(osgdem_main(argIdx, argptr))
	{
		MessageToDialog("VPB����ʧ��\n");
		return 0;
	}
	
	return 1;
}

/* ר�Ÿ������߾��ȵؾ�ʱ���õģ����������߾����ܱߵĵ;������� */
/* ����Ϊ���� doTerrain ����ţ�*/
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
	
	/* �������·����Ϣ */
	W = doTerrain[idx].minLongi;
	S = doTerrain[idx].minLati;
	offLon = doTerrain[idx].maxLongi - doTerrain[idx].minLongi;
	offLat = doTerrain[idx].maxLati  - doTerrain[idx].minLati ;
	set_bucket(&b, W+0.01, S+0.01);

	/* ƴװ����·�����ļ��������50�׾��ȵ����ݣ�������ʹ��50m�������ݣ�����ʹ��ԭ��15m�������� */
	iW = (int)(W + 0.00001); iS = (int)(S + 0.00001);
	if((index = CheckIfHas50(iW, iS)) != -1)
	{
		sprintf(temTextGeoFile, "%s\\%s", lowTxtnode[index].pathName, lowTxtnode[index].fileName);
	}
	else
		sprintf(temTextGeoFile, "%s\\T%d%d.tif", tempDirectory, (int)W, (int)S);

	sprintf(temHgtGeoFile , "%s\\D%d%d.tif", tempDirectory, (int)W, (int)S);

	/* ƴװ���·�����ļ� */
	memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr)); memset(basePath, 0, sizeof(basePath)); memset(fullPath, 0, sizeof(fullPath));
	gen_base_path(&b, basePath);
	sprintf(fullPath, "%s\\%s\\%d%d_%d.ive", iveroot, basePath, (int)W, (int)S, idx);
	sprintf(simpleName, "%d%d_%d.ive", (int)W, (int)S, idx);
	sprintf(srcPath, "%s\\%s", iveroot, basePath);

	/* �ȼ��һ����û������ļ������˾Ͳ����ٴ����� */
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

	/* ��vpb */
	if(osgdem_main(argIdx, argptr))
	{
		MessageToDialog("VPB����ʧ��\n");
		return 0;
	}

	/* ����дSTG */
	sprintf(ctileID, "%d%d_%d", (int)W, (int)S, idx);
	if(!WriteAllStgs(ctileID, 0))
	{
		MessageToDialog("����STG�ļ�ʧ��\n");
		return 0;
	}

	/* ����һ���������ɵĵؾ��ļ��� BackupLib Ŀ¼�¡� */
	char msg[STRING_STANDARD_LEN];
	sprintf(msg, "���ݵؾ��ļ��� --- %s", simpleName);
	MessageToDialog(msg);
	BackUpFileToTemp(simpleName, srcPath);
	
	return 1;
}

/* ���ݲ�ͬ�������ȼ���ȷ����ͬ��LOD����*/
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


/* ר�Ÿ������߾��ȵؾ�ʱ���õģ����������߾������� */
/* ����Ϊ���� findResultHTI ����ţ�*/
int DoHighVpbAndWriteStg(int idx)
{
	char arg[30][STRING_STANDARD_LEN];
	char *argptr[30];
	char basePath[STRING_STANDARD_LEN];
	char fullPath[_MAX_PATH], simpleName[_MAX_PATH], srcPath[_MAX_PATH];
	char ctileID[12];
	float W, S, offLon, offLat;
	Bucket b;
	
	//�������·����Ϣ
	W = findResultHTI[idx].startLongi;
	S = findResultHTI[idx].startLati;
	offLon = findResultHTI[idx].offsetLongi;
	offLat = findResultHTI[idx].offsetLati; 
	set_bucket(&b, W+0.01, S+0.01);

	/* ƴװ����·�����ļ� */
	strcpy(temTextGeoFile, findResultHTI[idx].texFname);
	strcpy(temHgtGeoFile , findResultHTI[idx].hgtFname);

	/* ƴװ���·�����ļ� */
	memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr)); memset(basePath, 0, sizeof(basePath)); memset(fullPath, 0, sizeof(fullPath));
	gen_base_path(&b, basePath);
	sprintf(fullPath, "%s\\%s\\%d.ive", iveroot, basePath, findResultHTI[idx].tileID);
	sprintf(simpleName, "%d.ive", findResultHTI[idx].tileID);
	sprintf(srcPath, "%s\\%s", iveroot, basePath);
	
	/* ���ݲ�ͬ�������ȼ���ȷ����ͬ��LOD����*/
	int lodLevel = GetLodLevelWithTextureHighLevel(findResultHTI[idx].highLevel);
	
	/* ���������� */
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

	/* ��vpb */
	if(osgdem_main(argIdx, argptr))
	{
		MessageToDialog("VPB����ʧ��\n");
		return 0;
	}

	/* ����дSTG */
	sprintf(ctileID, "%d", findResultHTI[idx].tileID);
	if(!WriteAllStgs(ctileID, 0))
	{
		MessageToDialog("����STG�ļ�ʧ��\n");
		return 0;
	}
	
	/* ����һ���������ɵĵؾ��ļ��� BackupLib Ŀ¼�¡� */
	char msg[STRING_STANDARD_LEN];
	sprintf(msg, "���ݵؾ��ļ��� --- %s", simpleName);
	MessageToDialog(msg);
	BackUpFileToTemp(simpleName, srcPath);
	
	return 1;
}

/*============================================================================================
 *  ���������ɸ���ʼ�������ù�����̡�
 *============================================================================================ */

/* ��д�û������ļ�userconfig.ini */
int CheckConfig(int rewrite)
{
	fi = fopen(configFile, "rt");
	if((fi == NULL)||(rewrite)||(updateTTNode))
	{
		//�ر��Ѿ����򿪵������ļ���
		needUpdate = 0;
		if(fi) fclose(fi);
		
		/* �����û������ļ� */
		fo = fopen(userConfig, "wt");
		if(fo == NULL)
		{
			return 0;
		}
		
		/* д���û������ļ�����д��Ŀ¼ */
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

/* ��鼸�������Ŀ¼�Ƿ�Ϊ�� */
int CheckDirectories()
{
	if(strlen(currTi.txtDir) == 0)  strcpy(currTi.txtDir, baseTextDir);
	if(strlen(currTi.hgtDir) == 0)	strcpy(currTi.hgtDir, baseHgtDir); 
	if(strlen(currTi.iveDir) == 0)
		return 0;
	else
		return 1;
	
}

/* ����ȫ�ֱ��� */
int SetGlobalVars()
{
	CreateUserSceneryDirectory();
	strcpy(hgtroot, baseHgtDir);			//currTi.hgtDir; 
	sprintf(iveroot, "%s\\Scenery\\Terrain", currTi.iveDir);
	sprintf(objroot, "%s\\Scenery\\Objects", currTi.iveDir);
	return 1;
}


/* ��ʼ����鹤�� */
int CheckInitiation(int refresh)
{
	int retVal = 0;
	int iFlag = 1;

	/* ���Ŀ¼�Ƿ���� */
	updateTTNode = 0;
	if(!CheckDirectories())
	{
		MessageToDialog("Ŀ¼û��������ȫ.");
		return -1;
	}

#if USE_JPEG_DATA_TO_CREATE_TERRAIN
	/* ���Map��Ϣ�Ƿ��Ѿ�����, ���û�У��������*/
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
			MessageToDialog("��������Ŀ¼��ȡ��γ����Ϣʧ��\n");
			iFlag = 0;
		}
		updateTTNode = 1;
	}
#endif

	/* ���߾��ȵĵ�����Ϣ�Ƿ��Ѿ�����, ���û�У��������*/
	if((numHTI == 0)||(currTi.resolution[0] == 'H')||refresh)
	{
		char userPath[STRING_STANDARD_LEN];

		/* ֱ�Ӳ����û�����߳�Ŀ¼���Բ����û��߳����� */
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
			MessageToDialog("��������Ŀ¼��ȡ��γ����Ϣʧ��\n");
			iFlag = 0;
		}

		/* ֱ�Ӳ����û���������Ŀ¼���Բ����û��߷ֱ����������� */
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
			MessageToDialog("��������Ŀ¼��ȡ�û���������ʧ��\n");
			iFlag = 0;
		}

		/* ֱ�Ӳ���ϵͳ�Դ�������Ŀ¼���Բ���ϵͳ�ṩ��50�׷ֱ����������� */
		sprintf(userPath, "%s\\%s", baseTextDir, SYSTEM_LOWRESOLUTION);

		retVal = 0;
		isCurrTxtOrHgt = STATUS_TEXTURE;
		retVal += TravelDirectoryToReadTif(userPath);
		if(!retVal)
		{
			MessageToDialog("��������Ŀ¼��ȡϵͳ��������ʧ��\n");
			iFlag = 0;
		}

		updateTTNode = 1;
	}

	/* ���ac�Ļ���Ϣ�Ƿ��Ѿ�����, ���û�У��������*/
	if((cultID == 0)||refresh)
	{
		retVal = 0;
		retVal += TravelDirectoryToReadAc(baseTextDir);
		if(!retVal)
		{
			MessageToDialog("��������Ŀ¼��ȡAC�Ļ���Ϣʧ��\n");
			iFlag = 0;
		}
		updateTTNode = 1;
	}

	/* ���shp�Ļ���Ϣ�Ƿ��Ѿ�����, ���û�У��������*/
	if((shpID == 0)||refresh)
	{
		retVal = 0;
		retVal += TravelDirectoryToReadShp(baseTextDir);
		if(!retVal)
		{
			MessageToDialog("��������Ŀ¼��ȡSHP�Ļ���Ϣʧ��\n");
			iFlag = 0;
		}
		updateTTNode = 1;
	}

	/* ����Ƿ����û���ά��Ϣ */
	if((bldNum == 0)||refresh)
	{
		char userPath[STRING_STANDARD_LEN];

		/* ֱ�Ӳ����û�������ά����Ŀ¼���Բ����û���ά�������� */
		sprintf(userPath, "%s\\data\\Models\\User", fsData);

		retVal = 0;
		retVal += TravelDirectoryToReadFlt(userPath);
		if(!retVal)
		{
			MessageToDialog("������ά����Ŀ¼��ȡ��ά������Ϣʧ��\n");
			iFlag = 0;
		}
		updateTTNode = 1;
	}


	return iFlag;
}

/* ������ˢ�µ�ʱ���� */
int RefreshTheConfigFile()
{
	if( CheckInitiation(1) == -1)
	{
		return 0;
	}
	return 1;
}

/* ��ȡ��ת�������ļ� */ 
int ConvertConfig()
{
	/* �ȳ��Դ��û������ļ� */
	fi = fopen(userConfig, "rt");
	if(fi == NULL)
	{
		/* ����򲻿��û������ļ������ȫ�������ļ������Ҵ���һ���û������ļ������й��û�����д����ļ� */
		/* �������ļ� */
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
	
		/* ��ȡĿ¼ */
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
	
		/* ��ȡMap��Ϣ */
		ReadMapInfoFromConfig(fi);
		ReadDirInfoFromConfig(fi);
	
		ReadHTxtInfoFromConfig(fi);
		ReadHHgtInfoFromConfig(fi);
		ReadLowInfoFromConfig(fi);
		ReadHighInfoFromConfig(fi);
		ReadCultInfoFromConfig(fi);
		ReadBldInfoFromConfig(fi);
		fclose(fi);
		
		/* ��ԭ�д���highTxtnode�����50�׷ֱ�������ת�浽lowTxtnode���� */
		MoveHigh50ToLow();
		
		/* �����û������ļ� */
		fo = fopen(userConfig, "wt");
		if(fo == NULL)
		{
			initCheck.usrConfigFile = false;
			return 0;
		}
		
		/* д���û������ļ�����д��Ŀ¼ */
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

		/* Ȼ������дȫ�������ļ� */
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
		/* ������û������ļ��� �ȶ�ȡ�û������ļ� */
		/* ��ȡĿ¼ */
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
		/* �ٶ�ȡȫ�������ļ� */	
		ReadMapInfoFromConfig(fi);
		ReadDirInfoFromConfig(fi);
		ReadLTxtInfoFromConfig(fi);
		ReadLowInfoFromConfig(fi);
		fclose(fi);	
	}
	return 1;
}

/* ��鲢�����û��߷ֱ��ʵؾ���Ŀ¼ */
int CreateUserSceneryDirectory()
{
	/* �Ƚ���currTi.iveDirĿ¼ */
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
	if(io == -1)	/* û��Scenery���Ŀ¼����� */
	{
		/* ��scenery */
		io = mkdir(userScenery);
		if(io == -1)
		{
			return 0;
		}
		io = chdir(userScenery);
		/* ��terrain */
		sprintf(userScenery, "Terrain");
		io = mkdir(userScenery);
		if(io == -1)
		{
			return 0;
		}
		/* ��objects */
		sprintf(userScenery, "Objects");
		io = mkdir(userScenery);
		if(io == -1)
		{
			return 0;
		}
	}else{	/* ��Scenery���Ŀ¼����� */
		io = chdir(userScenery);

		/* ��terrain */
		sprintf(userScenery, "Terrain");
		io = chdir(userScenery);
		if(io == -1)	/* û��terrainĿ¼����� */
		{
			io = mkdir(userScenery);
			if(io == -1)
			{
				return 0;
			}
		}else{
			io = chdir("..");	/* ������һ�� */
		}
		
		/* ��Objects */
		sprintf(userScenery, "Objects");
		io = chdir(userScenery);
		if(io == -1)	/* û��ObjectsĿ¼����� */
		{
			io = mkdir(userScenery);
			if(io == -1)
			{
				return 0;
			}
		}else{
			io = chdir("..");	/* ������һ�� */
		}
	}
	
	/* ��󷵻ع���Ŀ¼��*/
	io = chdir(workDirectory);
	
	return 1;
}

/* ϵͳ��ʼ�����ڳ��򴰿ڳ���֮ǰ��ɵĹ����� */
int initCreateTerrain()
{
	int io;

	/* ��ʼ�������Ŀ�� */
	initCheck.backDirectory  = true;
	initCheck.pathOfHgt  = true;
	initCheck.pathOfTarget	= true;
	initCheck.pathOfTexture  = true;
	initCheck.tempDirectory  = true;
	initCheck.usrConfigFile  = true;

	//��ʼ��ȫ�ֱ���
	_getcwd(workDirectory, _MAX_PATH);
	sprintf(tempDirectory, "%s\\Temp", workDirectory);
	sprintf(temTextFile		, "%s\\%s.%s",  tempDirectory, TEMP_TEXTFILE	, TEMP_EXTENDNAME);
	sprintf(temTextGeoFile  , "%s\\%s.%s",  tempDirectory, TEMP_TEXTGEOFILE, TEMP_EXTENDNAME);
	sprintf(temHgtFile      , "%s\\%s.%s",  tempDirectory, TEMP_HGTFILE	, TEMP_EXTENDNAME);
	sprintf(temHgtGeoFile   , "%s\\%s.%s",  tempDirectory, TEMP_HGTGEOFILE	, TEMP_EXTENDNAME);
	
	/* ��ȫ��Log�ļ���*/
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

	/* ��ȡ�ƹ���ɫ�����ļ� */
	ReadColorsFromFile();

	/* ����һ��������������Ϣ���ʹ�� */
	MsgMutex = CreateMutex(NULL, false, NULL);

	/* ���һ�������·�� */
	FindBasePathName();

	/* �����ʱ·���Ƿ���ڣ�������ڣ��򴴽�һ����*/
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
	
	/* ������ʱ·��������һ������·�� */
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

	/* ��ʼ���߳�Ӧ������*/
	rootNode = new osg::Group;
	openedTiles.clear();
	allLakeSkirt   = new osg::Group;
	allRoadSkirt   = new osg::Group;
	allDynamicGrass   = new osg::Group;
	allIslands	   = new osg::Group;
	allLakeInCurrLonLat	   = new osg::Group;	

	/* ��ʼ�������ɵؾ������� */
	TextDirIdx = 0;
	allDirs.clear();
	numHTI = 0;
	allHTI.clear();

	/* ��ʼ����������� */
	mt_init(&random_seed, 12345);

#if 0
	char aptFile[STRING_STANDARD_LEN];
	sprintf(aptFile, "F:\\FG\\data\\Models\\User\\zh8-flat\\zh8-flat.flt");
	osg::ref_ptr<osg::Geode> aptGD = new osg::Geode;
	ReadFltAndVistIt(aptFile, aptGD);
	osgDB::writeNodeFile(*aptGD, "F:\\aptgd.osg");
#endif

	/* ��ϵͳ��������PATH������Լ��Ĺ���·�� */
	CreateEnvPath();

	/* ��ȡ��ת�������ļ� */
	ConvertConfig();
	
	/* ����ȫ�ֱ��� iveroot�� hgtroot */
	SetGlobalVars();
	
	return 1;
}

void DoInitializationCheck()
{
	/* �ѵ�ǰ������д��Log�ļ� */
	WriteTerrainCfgIntoLogFile();
	WriteCultureCfgIntoLogFile();

	char msg[STRING_STANDARD_LEN * 2];

	sprintf(msg, "��ʼ��ʼ�����......");
	sprintf(msg, "��ǰ����·��: %s", workDirectory);
	MessageToDialog(msg);

	/* �����Ŀ����ʱ·�� */
	if(initCheck.tempDirectory)
		sprintf(msg, "�����Ŀ����ʱ·�� --- ����");
	else
		sprintf(msg, "�����Ŀ����ʱ·�� --- ʧ��");
	MessageToDialog(msg);
	
	if(initCheck.tempDirectory)
	{
		sprintf(msg, "��ʱ·��: %s", tempDirectory);
		MessageToDialog(msg);
	}

	/* �����Ŀ��ԭʼ�ⱸ��·��*/
	if(initCheck.backDirectory)
		sprintf(msg, "�����Ŀ��ԭʼ�ⱸ��·�� --- ����");
	else
		sprintf(msg, "�����Ŀ��ԭʼ�ⱸ��·�� --- ʧ��");
	MessageToDialog(msg);

	if(initCheck.backDirectory)
	{
		sprintf(msg, "ԭʼ�ⱸ��·��: %s", backDirectory);
		MessageToDialog(msg);
	}

	/* �����Ŀ���û������ļ�*/
	if(initCheck.usrConfigFile)
		sprintf(msg, "�����Ŀ���û������ļ� --- ����");
	else
		sprintf(msg, "�����Ŀ���û������ļ� --- ʧ��");
	MessageToDialog(msg);

	if(initCheck.usrConfigFile)
	{
		sprintf(msg, "�û������ļ�: %s", userConfig);
		MessageToDialog(msg);
	}

	/* �����Ŀ�������ز�·��*/
	if(initCheck.pathOfTexture)
		sprintf(msg, "�����Ŀ�������ز�·�� --- ����");
	else
		sprintf(msg, "�����Ŀ�������ز�·�� --- ʧ��");
	MessageToDialog(msg);

	if(initCheck.pathOfTexture)
	{
		sprintf(msg, "�����ز�·��: %s", baseTextDir);
		MessageToDialog(msg);
	}

	/* �����Ŀ���߳��ز�·��*/
	if(initCheck.pathOfHgt)
		sprintf(msg, "�����Ŀ���߳��ز�·�� --- ����");
	else
		sprintf(msg, "�����Ŀ���߳��ز�·�� --- ʧ��");
	MessageToDialog(msg);

	if(initCheck.pathOfHgt)
	{
		sprintf(msg, "�߳��ز�·��: %s", baseHgtDir);
		MessageToDialog(msg);
	}

	/* �����Ŀ�����ɽ��·��*/
	if(initCheck.pathOfTarget)
		sprintf(msg, "�����Ŀ�����ɽ��·�� --- ����");
	else
		sprintf(msg, "�����Ŀ�����ɽ��·�� --- ʧ��");
	MessageToDialog(msg);
		
	if(initCheck.pathOfTarget)
	{
		sprintf(msg, "���ɽ��·��: %s", fsData);
		MessageToDialog(msg);
	}
	
}



/*============================================================================================
 *  ���������ɵؾ������̡��ӶԻ���ֱ�ӵ��á�
 *============================================================================================ */

/**************************************************
 *
 *   �������뾭γ��ֵ������1����x1γ������ĵ��Σ���ǰ�����Ϣ��currTi����ָ����
 *   noVpb = 0 --- ȫ���� = 1 --- ֻ������͸̣߳�  = 2 --- ֻ���߳�
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

	/* ���������tile��ʹ�õ�ͼƬ���� */
	if(noVpb)
		sprintf(msg, "��ʼ����������%3d, γ��%2d ���ĸ̼߳���������......", (int)currTi.minLongi, (int)currTi.minLati);
	else
		sprintf(msg, "��ʼ����������%3d, γ��%2d ���ĵؾ�......", (int)currTi.minLongi, (int)currTi.minLati);
	MessageToDialog(msg);

	/* ���㵱ǰ����ļ��޾�γ�� */			
	src_lon = currTi.minLongi;
	src_lat = currTi.minLati;
	W = src_lon;
	E = src_lon + 1.0;
	S = src_lat;
	N = src_lat + 1.0;

	/* ��������ʱĿ¼������ң�������Ҫ����ʱ�ļ��Ƿ���ڣ������ڵ�ʱ����ƴ�� */
	MessageToDialog("��ʼƴ�Ӹ̺߳���������");
	
	/* ���Դ򿪸߳��ļ� */
	sprintf(temHgtGeoFile   , "%s\\%s%d%d.%s",  tempDirectory, TEMP_HGTGEOFILE	, currPZ->lon, currPZ->lat, TEMP_EXTENDNAME);
	if( (ftmp = fopen(temHgtGeoFile, "rb")) != NULL)
	{
		fclose(ftmp);
		hgtFlag = 1;
	}

	/* ��ʼƴ�Ӹ߳����ݣ�����߳�����û�У�isOK = 2 */
	/* �������ǰ���Bucket */
	set_bucket(&b, src_lon, src_lat);
	GetHgtFileList(W, S, E, N);
	mainIndex = gen_index(&b);
	
	if(!hgtFlag)
	{
		if(MergeHgt() == false)
		{
			MessageToDialog("ƴ�Ӹ߳�ʧ��\n");
			return 0;
		}
	
		if(!convertHgtGeotiff())
		{
			MessageToDialog("�����߳�GeoTiff�ļ�ʧ��\n");
			return 0;
		}	
	}

	/* ���noVpb==2��ʾֻƴ�Ӹ߳����ݣ���ֱ���˳� */
	if(noVpb==2) return 1;

	/* ���Դ������ļ� */
	sprintf(temTextGeoFile  , "%s\\%s%d%d.%s",  tempDirectory, TEMP_TEXTGEOFILE, currPZ->lon, currPZ->lat, TEMP_EXTENDNAME);
	if( (ftmp = fopen(temTextGeoFile, "rb")) != NULL)
	{
		fclose(ftmp);
		textFlag = 1;
	}

	if(!textFlag)
	{
		/* ��ʼƴ���������� */
		if(!FindJpgsInZoneAndSort(W, E, S, N))	
		{
			MessageToDialog("����ͼƬʧ�ܡ�\n");
			return 0;
		}
		
		/* ���ݼ�����ƴ��ͼƬ */
		if(!MergePictures())
		{
			MessageToDialog("ƴ��ͼƬʧ��\n");
			return 0;
		}
	  	
		/*��tifͼƬת��ΪGeotiffͼƬ*/
		if(!convertTextureGeotiff())
		{
			MessageToDialog("��������GeoTiff�ļ�ʧ��\n");
			return 0;
		}
	}

	/* ���ѡ����vpb����ֱ���˳� */
	if(noVpb==1) return 1;
  	
	/*װ��vpb���������*/
	MessageToDialog("����VPB���ɵؾ�");
	if(!DoVpb(&b, 1))
	{
		MessageToDialog("ִ��VPBʧ��\n");
		return 0;
	}

	/* ����дSTG */
	sprintf(ctileID, "%d%d", (int)W, (int)S);
	if(!WriteAllStgs(ctileID, 0))
	{
		MessageToDialog("����STG�ļ�ʧ��\n");
		return 0;
	}
	
	return 1;
}

/**************************************************
 *
 *   �������ɵ���
 *
 **************************************************/
int BatchCreateTerrain()
{
	MessageTo = 0;
	MessageToDialog("������������......");
	
	InitOnDemand();

	/*��ʼ����鹤��*/
	if( CheckInitiation(0) == -1)
	{
		MessageToDialog("��ʼ����鹤��ʧ��.");
		return 0;
	}
	
	/*�����û�������ļ������û�о�дһ�������ļ� */
	if(!CheckConfig(needUpdate))
	{
		MessageToDialog("��������ļ�ʧ�ܡ�\n");
		return 0;
	}
	
	/* ����ǰ�ؾ����ɵ�������Ϣд��Log�ļ� */
	WriteTerrainCfgIntoLogFile();
	
	/* ����ȫ�ֱ��� iveroot�� hgtroot */
	SetGlobalVars();

	/* ֻ�����û��߷ֱ��ʵĵؾ����� */
	DoHightResolutionTerrain();

	MessageToDialog("���εؾ��������!!!!!!");
	return 1;
}








/*============================================================================================
 *  �������Ļ���Ϣ���������̣�Ŀǰ֧�������Ļ���Ϣ������ac�ļ� �� shp�ļ�
 *============================================================================================ */

/**********************************************************
 *
 *       �����ؾ�����Ŀ¼����Ѱ.Shp�ļ�����ȡ��Ϣ�������ݿ�
 *       
 **********************************************************/
/* �����ļ��������⣬��������ڣ��򷵻�1�����򷵻�0 */
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

/* ������ǰĿ¼�µ�������Ŀ¼���Բ��� shp �ļ���ʹ�õݹ鷽ʽ�� */
int TravelSubDirectoryForShp(char *dir)
{
	int dio, dfio, dirIo, i;
	struct _finddata_t sdff;

	//���ҵ�ǰ·���µ���Ŀ¼��
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//���������Ŀ¼Ҳ�����ļ�
		if((sdff.name[0] == '.')) goto GET_NEXT_SHP_SUBDIR2;

		//�������Ŀ¼�����뵱ǰ��Ŀ¼
		if(sdff.attrib == _A_SUBDIR)
		{
			//������Ŀ¼
			dirIo = chdir(sdff.name);
			_getcwd(path, _MAX_PATH);
			
			//������Ŀ¼
			TravelSubDirectoryForShp(path);
			
			//�˻ص���һ��Ŀ¼
			dirIo = chdir("..");
		}
		
		//������ļ�
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//�ҵ�Shp�ļ�
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='S')||(sdff.name[i-3]=='s')) && ((sdff.name[i-2]=='H')||(sdff.name[i-2]=='h')) && ((sdff.name[i-1]=='P')||(sdff.name[i-1]=='p')))  
			{
				/* �Ȳ��ҿ�����û�и��ļ� */
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
		//Ѱ����һ��
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* ������ǰĿ¼�����Ҹ�Ŀ¼���Լ�������Ŀ¼����չ��Ϊ.shp���ļ�����ȡ���������Ϣ�� */
int TravelDirectoryToReadShp(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForShp(dPath))
	{
		MessageToDialog("�Ҳ���Shp�ļ�");
		return 0;
	}
	
	CheckConfig(needUpdate);

	return 1;
} 

/**********************************************************
 *
 *       �����ؾ�����Ŀ¼����Ѱ.ac�ļ�����ȡ��Ϣ�������ݿ�
 *       
 **********************************************************/
/* �����ļ��������⣬��������ڣ��򷵻�1�����򷵻�0 */
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

/* ������ǰĿ¼�µ�������Ŀ¼���Բ��� ac �ļ���ʹ�õݹ鷽ʽ�� */
int TravelSubDirectoryForAc(char *dir)
{
	int dio, dfio, dirIo, i, j;
	struct _finddata_t sdff;

	//���ҵ�ǰ·���µ���Ŀ¼��
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		char path[_MAX_PATH];

		//���������Ŀ¼Ҳ�����ļ�
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR2;

		//�������Ŀ¼�����뵱ǰ��Ŀ¼
		if(sdff.attrib == _A_SUBDIR)
		{
			//������Ŀ¼
			dirIo = chdir(sdff.name);
			_getcwd(path, _MAX_PATH);
			
			//������Ŀ¼
			TravelSubDirectoryForAc(path);
			
			//�˻ص���һ��Ŀ¼
			dirIo = chdir("..");
		}
		
		//������ļ�
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//�ҵ�ac�ļ�
			i = strlen(sdff.name);
			if(((sdff.name[i-2]=='A')||(sdff.name[i-2]=='a')) && ((sdff.name[i-1]=='C')||(sdff.name[i-1]=='c')) )  
			{
				//����ͬ�ļ���.map�ļ������ַ���
				char mapn[24];
				memset(mapn, 0, sizeof(mapn));
				for(j=0;j<i-2;j++) mapn[j]=sdff.name[j];
				mapn[j++]='m'; mapn[j++]='a'; mapn[j++]='p'; 

				/* �Ȳ��ҿ�����û�и��ļ� */
				if(!CheckACFromTotalCult(sdff.name))
				{
					CultNode cn;
					strcpy(cn.fileName, sdff.name);
					
					//��ͬ�ļ���map�ļ�
					fi = fopen(mapn, "rt");
					if(fi == NULL) 
					{
						char msg[STRING_STANDARD_LEN];
						sprintf(msg, "���ܶ�ȡ�ļ� %s .", sdff.name);
						MessageToDialog(msg);
						return 0;
					}
					
					//��ȡͼƬ�ĵ�����Ϣ
					while((strcmp("MMPLL,1,", cf[0]))&&(!feof(fi)))  ReadALineToArray(fi);
					if(feof(fi))
					{
						MessageToDialog("��ȡmap�ļ�����:�Ҳ�����γ����Ϣ��");
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
		//Ѱ����һ��
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* ������ǰĿ¼�����Ҹ�Ŀ¼���Լ�������Ŀ¼����չ��Ϊ.ac���ļ�����ȡ���������Ϣ�� */
int TravelDirectoryToReadAc(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForAc(dPath))
	{
		MessageToDialog("�Ҳ���ac�ļ�");
		return 0;
	}
	
	CheckConfig(needUpdate);

	return 1;
} 

/**************************************************
 *
 *   ��ȡac�ļ���������BTG�ļ�
 *
 **************************************************/
typedef struct _LANDPOINT_
{
	double px;
	double py;				//ʵ������   ��Ҫ����ת��������
	double pz;				
	double lon;
	double lat;				//��γ�����꣬��Ҫ�������
	double tU;
	double tV;				//��������   ֱ�ӿɶ�ȡ��
	int	   iPrePos;			//���������е������ڵ�λ�á�
	int    iCurPos;			//������ڵ�ǰ��λ�á�
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
	double tV;				//��������   ֱ�ӿɶ�ȡ��
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
STypeIndex typeIdxs[TOTALTYPE*10];	//ָ����Ӧ�㼯�ϵ�����
int tiNums;							//����typeIdxs����ʵ����Чֵ������
int allTypes[TOTALTYPE];

Geodesy geoCenter;
Carton catCenter;
double gbs_radius;
Carton cPoint[4];
Geodesy gPoint[4];
char sMaterial[20];
Carton cP[5];
int totalGPoint, sumGPoint[TOTALTYPE], sumGFace[TOTALTYPE], sumObjs, curObjs;

/* ��ac�ļ��������� geode ���������û�б����ù����������̡�*/
osg::ref_ptr<osg::Geode> createGeodeFromAc()
{
	int i, j, currIdx, offset;
	float fNx, fNy, fNz;
	
	//ȱʡһ������
	fNx = 92/127.5; fNy = -30/127.5; fNz = -57/127.5; offset = 0;

	osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
	gnode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	gnode->setName("River");

	//����ÿ�������棬����һ�������壬��Ϊһ��drawable���뵽gnode������
	for(i=0; i<tiNums; i++)
	{
		currIdx = typeIdxs[i].idx;

		osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
		gnode->addDrawable(gm);
		
		//ѹ�붥��
		osg::ref_ptr<osg::Vec3Array> vertex = new osg::Vec3Array;
		float v3 = Points[currIdx][0].pz;
		for(j=0; j<sumPtsObj[currIdx]; j++)
		{
			osg::Vec3d v; v.set(Points[currIdx][j].px, Points[currIdx][j].py, v3/*Points[currIdx][j].pz*/);
			vertex->push_back(v);
		}
		gm->setVertexArray(vertex);
		
		//ѹ�뷨��
		osg::ref_ptr<osg::Vec3Array> normal = new osg::Vec3Array;
		normal->push_back(osg::Vec3(fNx, fNy, fNz));
		gm->setNormalArray(normal);
		gm->setNormalBinding(osg::Geometry::BIND_OVERALL);
		
		//ѹ����������
		osg::ref_ptr<osg::Vec2Array> coord = new osg::Vec2Array;
		for(j=0; j<sumPtsObj[currIdx]; j++)
		{
			osg::Vec2 c; c.set(Points[currIdx][j].tU, Points[currIdx][j].tV);
			coord->push_back(c);
		}
		gm->setTexCoordArray(0, coord);

		//ѹ����		
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

/* ��ÿ�����Ļ���Ϣ֮ǰ����ɾ����ǰ�����������Ļ���Ϣ������btg�ļ���G*.ive�ļ��� */
int DeleteAllPrevCultureFiles(int lon, int lat)
{
	char fullPath[_MAX_PATH];
	std::vector<string> iveFiles;	iveFiles.clear();
	std::vector<string> stgFiles;	stgFiles.clear();
	FILE *fpstg;

	/* ������ͬĿ¼�µ�btg�ļ����Լ�����ive�ļ����ҵ���ɾ�� */
	Bucket Bt;
	char basePath[STRING_STANDARD_LEN];
	set_bucket(&Bt, lon, lat);
	gen_base_path(&Bt, basePath);
	sprintf(fullPath, "%s\\%s", iveroot, basePath);		
	{
		struct _finddata_t sdff;
		int io, dio, dfio;

		/* �����·���²��� *.btg.gz */
		io = chdir(fullPath);
		
		dio = _findfirst("*.btg.gz", &sdff);
		if((!dio)||(dio==-1)) goto DAPCF_Continue_1;
		dfio = 0;
	
		while(!dfio)
		{
			
			DeleteFile(sdff.name);

			/* Ѱ����һ�� */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}

DAPCF_Continue_1:

	/* �ٰ�����ͬĿ¼�µĻ���ive�ļ����ҵ���ɾ��������ͨ�û������û���λ���� */
	/* Ҳ���ǰ����·�������е�G*.ive������Ե�߻�Ҳɾ�� */
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
		osg::ref_ptr<osg::Node> node;
		osg::ref_ptr<osg::Geode> gd;
	
		/* �����·���²��� *.osg */
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
			/* Ѱ����һ�� */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}
		
DAPCF_Continue_4:

	/* ����ȡ��ǰĿ¼�����е�stg�ļ�����stg�ļ������*.btg.gzҲɾ�� */
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
	
		/* �����·���²��� *.osg */
		io = chdir(fullPath);
		
		dio = _findfirst("*.stg", &sdff);
		if((!dio)||(dio==-1)) goto DAPCF_Continue_6;
		dfio = 0;
	
		while(!dfio)
		{
			/* ����ǰĿ¼������stg�ļ��� */
			stgFiles.push_back(sdff.name);

			/* ɾ��stg�ļ� */
			DeleteFile(sdff.name);
			/* Ѱ����һ�� */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}
	
	/* ����дһ��stg�ļ� */
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



/* ���������������ȡ�ض�Ŀ¼�µĵؾ�ive���ļ��Լ��������ļ����������Ļ���Ϣ��ͨ�ü��û������������и�����и�ؾ��������̡�*/
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

	/* �������е��Ļ���Ϣ����ļ� */
	/* ��һ�����Ȱ�����ͬĿ¼�µ�btg�ļ����Լ�����ive�ļ�����һ��group���� */
	osg::ref_ptr<osg::Group> gp = new osg::Group;
	Bucket Bt;
	char basePath[STRING_STANDARD_LEN];
	set_bucket(&Bt, lon, lat);
	gen_base_path(&Bt, basePath);
	sprintf(fullPath, "%s\\%s", iveroot, basePath);		
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
	
		/* �����·���²��� *.btg.gz */
		io = chdir(fullPath);
		
		dio = _findfirst("*.btg.gz", &sdff);
		if((!dio)||(dio==-1)) goto RTIAC_Continue_1;
		dfio = 0;
	
		while(!dfio)
		{
			osg::ref_ptr<osg::Group> cgp;

			/* �����ж��ǲ����Ļ���Ϣ��btg�ļ��� */
			if(((sdff.name[0] > '9') && (strlen(sdff.name) > 13)))
			{
				/* �����ȹ�ߣ������и*/
				if((sdff.name[0] != 'S') && 
					 (sdff.name[0] != 'O') )		// ���ͺ��治���и�
				{
					/* ������Ļ���Ϣbtg���뵽һ�� group���� */				
					ReadBTGToDataStruct(sdff.name);	
					cgp = CreateGroupFromDataStru();
					for(int j=0; j<cgp->getNumChildren(); ++j)
					{
						gp->addChild(cgp->getChild(j));
					}
				}
			}

			/* Ѱ����һ�� */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}

RTIAC_Continue_1:

	/* �ڶ������ٰ�����ͬĿ¼�µĻ���ive�ļ�����ͬ�ϵ�һ��group���棬����ͨ�û������û���λ���� */
	/* Ҳ���ǰ����·�������е�G*.ive������Ե�߻�Ҳ�ӽ��� */
	{
		struct _finddata_t sdff;
		int io, dio, dfio;
		osg::ref_ptr<osg::Node> node;
		osg::ref_ptr<osg::Geode> gd;
	
		/* �����·���²��� *.osg */
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
			/* Ѱ����һ�� */
			dfio = _findnext(dio, &sdff);
		}
		_findclose(dio);		
	}
		
RTIAC_Continue_4:
	
	/* �����ǿ�ʼ�и�ؾ��Ķ�������������ive��Ҫ�ų������߻���ive��Ȼ���и�ؾ� */
	/* ���ҵ�ǰĿ¼�µ�ive */
	io = chdir(fullPath);
	if(io) 
	{
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", fullPath);
		MessageToDialog(msg);
		return 0;
	}
	
	/* �޸����о��ȵؾ� */
	sprintf(subDir, "*.ive");
	ffio = _findfirst(subDir, &ivef);	
	if((!ffio)||(ffio == -1))
	{
		MessageToDialog("�Ҳ���ive�ļ���");
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



	/* ��ʱ���ƽ���������б� */
	//csmList.clear();


	while(!fio)
	{
		/* �ֱ��ʵ͵ĵؾ������Ȳ��� */
		int flen = strlen(ivef.name);
		if(flen < 12)	goto RTIAC_Continue_2;

#if 0
		if(strcmp(ivef.name, "114230108.ive"))
			goto RTIAC_Continue_2;
#endif

		/* �����ų���������Ե�߻���.ive�ļ� */
		if(ivef.name[0] == 'G') goto RTIAC_Continue_2;
		
		sprintf(msg, "��ʼ�޸ĸ��ؾ��ļ��� -- %s", ivef.name);
		MessageToDialog(msg);
		sprintf(fullIveFile, "%s\\%s", fullPath, ivef.name);

#ifdef	USE_HCI
		/* �õ���ǰivef.name�� tileID �š�����鵱ǰ���Ƿ��������Ļ���Ϣ */
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
			/* �ָ���ǰ��ĵؾ� */
		}
#endif	//USE_HCI

		/* ���ص�����Ŀ¼�£��Ա�ʹ��osg��̬���ӿ� */		
		io = chdir(workDirectory);
		if(io) 
		{
			MessageToDialog("���ܻص�����·���¡�");
			return 0;
		}
		
		//�Բ��ҵ����ļ�����
		//����ive�ļ�
		UpdateIVEByBTG(gp, fullIveFile);

		/* ���ص���ǰĿ¼���Ա�������Ŀ¼�����ive�ļ� */
		io = chdir(fullPath);
		if(io) 
		{
			char msg[STRING_STANDARD_LEN];
			sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", fullPath);
			MessageToDialog(msg);
			return 0;
		}

#ifdef	USE_HCI
		/* �ѵ�ǰivefile��tileID���뵽 allHCI ���С�*/
		InsertHCI(tID);
#endif	//USE_HCI

		/* ���ȫ��ͨѶ���������������λ��˵����ǰ�ؾ��ޱ仯���ײ�Ҳ����Ҫ�����ˡ� */
		if(!CurrentTerrainIsCut)
		{
			goto RTIAC_Continue_2;
		}

		/* �����Ǵ���PageLOD�ĵײ�ĸ���ive�ļ������еײ��ļ������ XXXX_root_L9_X0_Y0 Ŀ¼�� */
		strcpy(subDir, ivef.name);
		int len = strlen(subDir);
		subDir[len - 4] = 0;
		sprintf(fullSubPath, "%s_root_L0_X0_Y0", subDir);
		
		/* ���뵽��Ŀ¼�� */
		io = chdir(fullSubPath);
		if(io) 
		{
			MessageToDialog("���ܽ���Ŀ��Ŀ¼��");
			goto RTIAC_Continue_2;
		}

		/* ����Ŀ¼�в���ive�ļ� */
		subffio = _findfirst("*.ive", &iveSubf);	
		if((!subffio)||(subffio==-1))
		{
			MessageToDialog("�Ҳ���ive�ļ���");
			goto RTIAC_Continue_2;
		}
		_getcwd(fullSubPath, _MAX_PATH);

		subfio = 0; subidx = 0;
		while(!subfio)
		{
			subidx++;
			sprintf(msg, "��ʼ�޸��ӵؾ��ļ�--%d.%s", subidx, iveSubf.name);
			MessageToDialog(msg);
			sprintf(fullIveFile, "%s\\%s", fullSubPath, iveSubf.name);

			/* ���ص�����Ŀ¼�£��Ա�ʹ��osg��̬���ӿ� */		
			io = chdir(workDirectory);
			if(io) 
			{
				MessageToDialog("���ܻص�����·���¡�");
				return 0;
			}

			//�Բ��ҵ����ļ�����
			//����ive�ļ�
			UpdateIVEByBTG(gp, fullIveFile);

			/* ���ص���ǰĿ¼���Ա����������Ŀ¼���������ive�ļ� */
			io = chdir(fullSubPath);
			if(io) 
			{
				MessageToDialog("���ܽ���Ŀ��Ŀ¼��");
				return 0;
			}
	
			//������һ�������������ļ�
			subfio = _findnext(subffio, &iveSubf);
		}
		_findclose(subffio);

		/* ���ص���ǰĿ¼���Ա����������һ��ive�ļ� */
		io = chdir("..");
		if(io) 
		{
			MessageToDialog("���ܽ���Ŀ��Ŀ¼��");
			return 0;
		}

RTIAC_Continue_2:
		//������һ�������������ļ�
		fio = _findnext(ffio, &ivef);
	}
	_findclose(ffio);
	
	return 1;
}

/* �����Χ��뾶 */
double calc_gbs(Geodesy *gP, int num)		//gP �ǳ���Ϊ4�����顣
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

/* �Ѷ�ȡac�ļ����õ�������ת�浽DataStruct���棬�Թ�WriteDataStructToBTG���� */
int MoveDataToDataStru()
{
	int i, j, m, n, typenum;
	short iTa, iTb, iTc, iT = 0;
	int   pFlag;

	//ͳ��һ�µ�ǰac�ļ��ܹ����˶��ٸ�types��
	typenum = 0; for(i=0;i<TOTALTYPE;i++) typenum+=allTypes[i];
	
	/* �����б�����ֵ���Ƶ�Buctet.cpp����ض�Ӧ�������Ա�ʹ��Buctet.cpp�ڵĺ��� */
	/* ��ʼ���ڴ� */
	totalNormal = 0;	curPnt.clear();		curText.clear();

	iVersion = 6; iMagicNum = 0x5347; iCreationTime = 0x48170f2c; iNumOfToplvlObject = 5 + typenum;
	dXpart = catCenter.x;
	dYpart = catCenter.y;
	dZpart = catCenter.z;
	fRadiusofBS = (float)gbs_radius;
	gbCenterX = dXpart; gbCenterY = dYpart; gbCenterZ = dZpart;

	/* ������Ϣ */
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

	/* ����������Ϣ */
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
	
	/* ������Ϣ */
	curPnt[totalNormal].nx = 92;	curPnt[totalNormal].ny = -30;	curPnt[totalNormal].nz = -57; totalNormal++;	

	/* ѭ����д���е��� */
	curTri.clear();
	for(j=0;j<TOTALTYPE;j++)
	{
		//��ͬ���Ͳ�ͬ��������
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
		//���㵱ǰ����Objects����Ŀ�����û�оͼ���д�����ġ�
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

/* ��ȡac�ļ���������Ȼ�󴴽�btg�ļ� */
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

	//��ȡ��ǰĿ¼
	_getcwd(path, _MAX_PATH);

	//ָ���������Ŀ¼
	sprintf(dPath, "%s", findResultAC[idx].pathName);
	io = chdir(dPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", dPath);
		MessageToDialog(msg);
		return 0;
	}

	ffio = _findfirst(findResultAC[idx].fileName, &ff);	
	if((!ffio)||(ffio == -1))
	{
		MessageToDialog("�Ҳ���ac�ļ���");
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

			//��һ�������ҵ�numvert
			if((OneLine[0]=='n') && (OneLine[1]=='u') && (OneLine[2]=='m') && (OneLine[3]=='v') && (OneLine[4]=='e') && (OneLine[5]=='r') && (OneLine[6]=='t'))
			{
				int numpt, pIdx;

				//�ȵ���һ��sumObjs����Ŀ���ж���һ��Object�Ƿ��кϷ��ĵ㣬���û�У�����ˡ�
				sumObjs+=TOTALTYPE;
				if(sumObjs > 0)
				{
					int flags = 0;
					for(i=0; i<TOTALTYPE; i++) flags+=sumGPoint[i];
					if(!flags)  sumObjs-=TOTALTYPE;
				}

				//Get number ��ȡ���е�ĸ�����
				memset(cNum, 0, sizeof(cNum)); memset(sumGPoint, 0, sizeof(sumGPoint)); pIdx = 0;
				for(i=8;i<strlen(OneLine);i++)
				{
					if(OneLine[i]==0) break;
					cNum[i-8] = OneLine[i];
				}
				numpt = atoi(cNum); 
				
				//��ȡ��Ч�㡣���м�����ֲ�Ϊ0��
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
					
					//�����С�������
					if(f0 < minPx) minPx = f0; if(f0 > maxPx) maxPx = f0;
					if(f2 < minPy) minPy = f2; if(f2 > maxPy) maxPy = f2;
					
					//���ű����¼��������Ϣ
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

				//��¼ÿ�����͵����Ч�����
				for(i=0;i<TOTALTYPE;i++)
					sumPtsObj[sumObjs + i] = sumGPoint[i];		//��Ч�����
				
				//���ĳ���͵�ĸ�����Ϊ0�������Ӹ�����ָʾ
				for(i=0;i<TOTALTYPE;i++)
				{
					if(sumGPoint[i]) { typeIdxs[t].type = i+1; typeIdxs[t].idx = sumObjs + i; t++; 	allTypes[i] = 1;}
				}
				
				//��������
				pIdx = totalGPoint;
				for(i=0;i<TOTALTYPE;i++)
				{
					for(m=0;m<sumPtsObj[sumObjs + i]; m++)  Points[sumObjs + i][m].iCurPos = pIdx + m;
					pIdx += sumGPoint[i];
					totalGPoint += sumGPoint[i];
				}			
			}// End if 

			//��ȡ��������Ч����ɵ������Ρ�����˳����д��������
			if((OneLine[0]=='n') && (OneLine[1]=='u') && (OneLine[2]=='m') && (OneLine[3]=='s') && (OneLine[4]=='u') && (OneLine[5]=='r') && (OneLine[6]=='f'))
			{
				int numsurf;
				
				//Get number ��ȡ������ĸ�����
				memset(cNum, 0, sizeof(cNum)); memset(sumGFace, 0, sizeof(sumGFace));  
				for(i=8;i<strlen(OneLine);i++)
				{
					if(OneLine[i]==0) break;
					cNum[i-8] = OneLine[i];
				}
				numsurf = atoi(cNum); j=0;

				//��ȡ�����棬�ж���Ч�档
				for(n=0;n<numsurf;n++)
				{
					int numrefs, sumGFlag; 
					int GFlag[TOTALTYPE][10];
					
					while(!((OneLine[0]=='r') && (OneLine[1]=='e') && (OneLine[2]=='f'))) {	fgets(OneLine, 80, fi);}	
					
					//Get number ��ȡ���е�ĸ�����
					memset(cNum, 0, sizeof(cNum)); memset(Ref, 0, sizeof(Ref)); 
					for(i=5;i<strlen(OneLine);i++)
					{
						if(OneLine[i]==0) break;
						cNum[i-5] = OneLine[i];
					}
					numrefs = atoi(cNum); j=0;

					//��ȡ�档
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
					
					//�жϸ����Ƿ�Ϊ��Ч�棬�����еĵ㶼����Ч�㣩
					memset(GFlag, 0 ,sizeof(GFlag)); sumGFlag = 0;
					for(m=0;m<numrefs;m++)
					{
						for(i=0; i<TOTALTYPE; i++)
						{
							for(j=0; j<sumPtsObj[sumObjs + i]; j++)
							{
								if(Ref[m].idx == Points[sumObjs + i][j].iPrePos)		//����õ�����Ч��
								{
									Points[sumObjs + i][j].tU = Ref[m].tU; Points[sumObjs + i][j].tV = Ref[m].tV;		//������Ϣ
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
						if(sumGFlag == numrefs) 	//˵�����е㶼��Ч
						{
							int i;
							for(i=0;i<numrefs;i++)
								Faces[sumObjs + j][sumGFace[j]].idx[i] = Ref[i].iNewdx;
							sumGFace[j]++; 
							break;
						}
					}
				}//End for m	

				//ÿ�����͸��ж��ٸ���				
				for(m=0; m<TOTALTYPE; m++)
					sumFaceObj[sumObjs + m] = sumGFace[m];

			}// End if 					
		}//End while
		fclose(fi);
		tiNums = t;

		//�������е���Ч��ľ�γ��
		absPx = fabs(fabs((float)maxPx) - fabs((float)minPx)) ; absPy = fabs(fabs((float)maxPy) - fabs((float)minPy)) ; offLon = fabs(findResultAC[idx].maxLong - findResultAC[idx].minLong); offLat = (findResultAC[idx].maxLati - findResultAC[idx].minLati);
		totalPnts = 0;

		/* ���浱ǰ·�������ص�����·���£��Ա�ִ�� GetElevation() */
		_getcwd(currDirectory, _MAX_PATH);	
		io = chdir(workDirectory);
		if(io) 
		{
			MessageToDialog("���ܻص�����·���¡�");
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

		//����׼����������ʼдbtg
		//���㵱ǰ������ĵ㾭γ�ȡ������Χ��뾶
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

		/* �ӹ���·�����ص���ǰ·�����Ա������������ */
		io = chdir(currDirectory);
		if(io) 
		{
			MessageToDialog("���ܻص�����·���¡�");
			return 0;
		}

		//���ļ���д�ļ��������ļ�����·�������м�������Ϊ��һ����׼����
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

		//��Ҫ������׼����ϣ���ʼת��.
		MoveDataToDataStru();
		/* ����һ����ӵ�stg�ļ���β */
		/* ����дSTG */
		if(!CreateCultureBtgFromSrc(idx))
		{
			MessageToDialog("����STG�ļ�ʧ��\n");
			return 0;
		}

		//������һ�������������ļ�
		fio = _findnext(ffio, &ff);
	}
	_findclose(ffio);
	
	
	/* ����ԭ�й���Ŀ¼ */
	io = chdir(path);
	if(io) 
	{
		MessageToDialog("���ܷ���ԭ����Ŀ¼��");
		return 0;
	}
	
#endif	
	return 1;
} 



/* ����ȫ��ac�ļ� */
int DoAllAc()
{
	int i, j, k;
	char msg[STRING_STANDARD_LEN];
	int currDoing[STRING_STANDARD_LEN];
	int numCD;
	PerZone mPZ[STRING_STANDARD_LEN];
	int numPZ, iFlag;

	/* �ҳ���û�������Ļ���Ϣ������ac�ļ����ڵ�Ŀ¼���� */
	if(FindOutDirsNoCultDone() == 0)
	{
		MessageToDialog("û����Ҫת�����Ļ���Ϣ.");
		return 1;	
	}

	/* �ȴ�dirsNoCultDone���𣬼���ЩĿ¼���涼��û���Ļ���Ϣ�����ac�ļ� */
	for(i=0; i<numDNCD; i++)
	{
		/* ���dirsNoCultDone����ÿһ���totalcult��������ͬĿ¼������ϵ�findResultAC*/
		memset(currDoing, 0, sizeof(currDoing));  numCD = 0;
		memset(mPZ, 0, sizeof(mPZ));  numPZ = 0;
		if(FindAllAcItemsWithSameDir(i) == 0) continue;
		strcpy(findJpgDirectory, findResultAC[0].pathName);
		
		/* ���findResultAC����ÿһ��ȼ���Ƿ������Ļ���Ϣ����û��������currDoing���¼��*/
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
		
		/* �Ӹո����õ����п������ҳ����в�ͬ�Ŀ� */
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
	
		/* ����ÿ���飬���������ؾ����ټ��Ļ���Ϣ */
		for(j=0; j<numPZ; j++)
		{
			/* ������ʼ��ȫ�ֱ��� */
			currTi.minLongi = (float)(mPZ[j].lon);
			currTi.minLati  = (float)(mPZ[j].lat);
			currTi.maxLongi = (float)(mPZ[j].lon + 1);
			currTi.maxLati  = (float)(mPZ[j].lat + 1);
			currPZ = &mPZ[j];

			/* ����currDoing�������Ϣ��ת���Ļ���Ϣbtg */
			for(k=0; k<numCD; k++)
			{
				if(	(findResultAC[currDoing[k]].isOK	   == 0			) && 
					(findResultAC[currDoing[k]].destLongi  == mPZ[j].lon) && 
					(findResultAC[currDoing[k]].destLati   == mPZ[j].lat) )
				{
					sprintf(msg, "����ת��BTG--%s...", findResultAC[currDoing[k]].fileName);
					MessageToDialog(msg);
					if(!ReadACandMakeBTG(currDoing[k]))
					{
						sprintf(msg, "ת��BTG--%sʧ�ܡ�", findResultAC[currDoing[k]].fileName);
						MessageToDialog(msg);
						continue;
					}
				}
			}   
			    
			/* ����ͬ���ڵ�����btg���һ��Geode�����뵽�ڶ��б��У�������ͳһ�ڶ��� */
			if(!InsertItemIntoCutList(mPZ[j].lon, mPZ[j].lat))
			{   
				sprintf(msg, "ת��IVEʧ�ܡ�");
				MessageToDialog(msg);
				return 0;
			}   
		}       

		/* ���õ��Ļ���Ϣ������������ݿ� */
		UpdateTotalcultWithNewFindResult();
	}
	

	return 1;	
}


/* ����ȫ��Shp�ļ� */
int DoAllShp()
{
	int i, j, k;
	char msg[STRING_STANDARD_LEN];
	int currDoing[STRING_STANDARD_LEN];
	int numCD;
	PerZone mPZ[STRING_STANDARD_LEN];
	int numPZ, iFlag;

	/* �ҳ���û�������Ļ���Ϣ������shp�ļ����ڵ�Ŀ¼���� */
	if(FindOutDirsNoShpDone() == 0)
	{
		MessageToDialog("û����Ҫת�����Ļ���Ϣ.");
		return 1;	
	}

	/* �ȴ�dirsNoCultDone���𣬼���ЩĿ¼���涼��û���Ļ���Ϣ�����ac�ļ� */
	for(i=0; i<numDNCD; i++)
	{
		/* ���dirsNoCultDone����ÿһ���totalcult��������ͬĿ¼������ϵ�findResultAC*/
		memset(currDoing, 0, sizeof(currDoing));  numCD = 0;
		memset(mPZ, 0, sizeof(mPZ));  numPZ = 0;
		if(FindAllShpItemsWithSameDir(i) == 0) continue;
		strcpy(findJpgDirectory, findResultAC[0].pathName);
		
		/* ���findResultAC����ÿһ��ȼ���Ƿ������Ļ���Ϣ����û��������currDoing���¼��*/
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
		
		/* �Ӹո����õ����п������ҳ����в�ͬ�Ŀ� */
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
	
		/* ����ÿ���飬���������ؾ����ټ��Ļ���Ϣ */
		for(j=0; j<numPZ; j++)
		{
			/* ������ʼ��ȫ�ֱ��� */
			currTi.minLongi = (float)(mPZ[j].lon);
			currTi.minLati  = (float)(mPZ[j].lat);
			currTi.maxLongi = (float)(mPZ[j].lon + 1);
			currTi.maxLati  = (float)(mPZ[j].lat + 1);
			currPZ = &mPZ[j];

			/* �Ȱѵ�ǰĿ¼�µĺ�������ա�*/
			allLakeInCurrLonLat->removeChildren(0,allLakeInCurrLonLat->getNumChildren());
			 
			/* ����currDoing�������Ϣ��ת���Ļ���Ϣbtg */
			for(k=0; k<numCD; k++)
			{
				if(	(findResultAC[currDoing[k]].isOK	   == 0			) && 
					(findResultAC[currDoing[k]].destLongi  == mPZ[j].lon) && 
					(findResultAC[currDoing[k]].destLati   == mPZ[j].lat) )
				{
					sprintf(msg, "����ת��BTG--%s...", findResultAC[currDoing[k]].fileName);
					MessageToDialog(msg);
					if(!ReadShpAndMakeBTG(currDoing[k]))
					{
						sprintf(msg, "ת��BTG--%sʧ�ܡ�", findResultAC[currDoing[k]].fileName);
						MessageToDialog(msg);
						continue;
					}
				}
			}   
			    
			/* ����ͬ���ڵ�����btg���һ��Geode�����뵽�ڶ��б��У�������ͳһ�ڶ��� */
			if(!InsertItemIntoCutList(mPZ[j].lon, mPZ[j].lat))
			{   
				sprintf(msg, "ת��IVEʧ�ܡ�");
				MessageToDialog(msg);
				return 0;
			}   
		}       

		/* ���õ��Ļ���Ϣ������������ݿ� */
		UpdatetotalShpWithNewFindResult();
	}

	return 1;	
}


/* ��ȡShp�ļ���������Ȼ�󴴽�btg�ļ� */
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

	/* ÿ�����Ļ���Ϣ֮ǰ������߳��������ݣ��Խ��·�汻�����ص��µ����� --- 2012.2.15 */
	RemoveAllTileNodes();

	//��ȡ��ǰĿ¼
	_getcwd(path, _MAX_PATH);

	//ָ���������Ŀ¼
	sprintf(dPath, "%s", findResultAC[idx].pathName);
	io = chdir(dPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", dPath);
		MessageToDialog(msg);
		return 0;
	}

	ffio = _findfirst(findResultAC[idx].fileName, &ff);	
	if((!ffio)||(ffio == -1))
	{
		MessageToDialog("�Ҳ���SHP�ļ���");
		return 0;
	}
	fio = 0;

	while(!fio)
	{
		/* ���������Ϣ���findResultAC�����հ��� */
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

		/* �����������ļ�����·�������м�������Ϊ��һ����׼����*/
		findResultAC[idx].destTileID = 0;
		tile = findUniqueTile(iTile);
		findResultAC[idx].destTileID = tile;
		sprintf(dstPath, "%s\\%s", iveroot, basepath);
		sprintf(dstFile, "%d.btg.gz", tile);
		strcpy(findResultAC[idx].destPath, dstPath);
		strcpy(findResultAC[idx].destFile, dstFile);
		findResultAC[idx].destLongi = (int)geoCenter._lon;
		findResultAC[idx].destLati  = (int)geoCenter._lat;
		
		/* ��ȡShp�ļ����Ұ�������DataStruct */
		currShpType = -1;
		unsigned int i;
		int FoundFlag = 0;
		int ShpResult = LoadShpAndDrawIt(ff.name);
		if(ShpResult > 0)
		{
			/* ���ҵ�ǰShp�ļ��Ƿ�ǰ���Ѿ��ù���*/
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
				/* ��DataStruct����д��Btg */
				/* ����һ����ӵ�stg�ļ���β */
				/* ����дSTG */
				if(!CreateCultureBtgFromSrc(idx))
				{
					MessageToDialog("����STG�ļ�ʧ��\n");
					return 0;
				}
			}
		}

		/* �����ǰshp�ļ��Ѿ���ʹ�ù�����ô����ļ����ɵĵ�ƺ͵���һ���Ѿ������ɹ����Ͳ���Ҫ�ٱ��ظ������ˡ�*/
		if(FoundFlag)	goto RSAMB_Continue_1;
		
		/* ���CultPoints���գ���ʾ�еƹ⻹ûд����дһ���ƹ�btg�ļ��� */
		if(CultPoints.size() > 0)
		{
			/* �ȵ����ƹ���ʾ��Զ������Ч�� */
			if(!AddLODIntoLightPoints())	return 0;
			
			/* ����ƹ�����ݽṹDataStru */
			if(CreateLightsDataStructFromCultPoints())
			{
				/* �����ļ�����*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\%dL.btg.gz", findResultAC[idx].destPath, tile); //����Ϊ tileID+L.btg.gz
				sprintf(fn, "%dL.btg.gz", tile);
			
				/* дbtg�ļ�������tileL.btg.gzд�����stg�� */
				if(WriteDataStructToBTG(fullPathFile))
				{
					if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "���ļ�ʧ�� --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}

		/* ���TreePoints���գ���ʾ��ɭ�ֻ�ûд����дһ��ɭ��btg�ļ��� */
		if(TreePoints.size() > 0)
		{
			/* �ȵ�������ʾ��Զ������Ч�� */
			if(!AddLODIntoTreePoints())	return 0;
			
			/* ������ľ�����ݽṹDataStru */
			if(CreateTreesDataStructFromTreePoints())
			{
				/* �����ļ�����*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\%dT.btg.gz", findResultAC[idx].destPath, tile); //����Ϊ tileID+T.btg.gz
				sprintf(fn, "%dT.btg.gz", tile);
			
				/* дbtg�ļ�������tileT.btg.gzд�����stg�� */
				if(WriteDataStructToBTG(fullPathFile))
				{
					if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "���ļ�ʧ�� --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}
		
		/* ����Ǻ����Ѻ���������һ��Btg�ļ�������Ϊ L+tileID.btg */
		if(allLakeSkirt->getNumChildren() > 0)
		{
			/* �������Ƿ��Ѿ���������btg�ļ� */
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
				/* ������������ݽṹDataStru */
				if(CreateLakeDataStruAndBtg(currLakes))
				{
					/* �����ļ�����*/
					char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
					sprintf(fullPathFile, "%s\\L%d.btg.gz", findResultAC[idx].destPath, tile); //����Ϊ L+tileID.btg.gz
					sprintf(fn, "L%d.btg.gz", tile);
				
					/* дbtg�ļ�������L+tile.btg.gzд�����stg�� */
					if(WriteDataStructToBTG(fullPathFile))
					{
						/* ����btg�ļ�����д��stg���� */
						//if(!WriteAllStgs(fn, 1)) return 0;
					}else{
						sprintf(msg, "���ļ�ʧ�� --- %s", fullPathFile);
						MessageToDialog(msg);
						return 0;
					}
				}

				/* Add by Galen 2013-02-20 */
				/* ����currLakes�鴴�����͵ĺ��棬��д��STG�����ļ� */
				/* ���ͺ����ļ���������O+tileID.btg.gz */
				if(CreateNewStyleLakeDataStruAndBtg(currLakes))
				{
					/* �����ļ�����*/
					char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
					sprintf(fullPathFile, "%s\\O%d.btg.gz", findResultAC[idx].destPath, tile); //����Ϊ O+tileID.btg.gz
					sprintf(fn, "O%d.btg.gz", tile);
				
					/* дbtg�ļ�������O+tile.btg.gzд�����stg�� */
					if(WriteDataStructToBTG(fullPathFile))
					{
						/* ���ͺ����btg�ļ���д��stg���� */
						if(!WriteAllStgs(fn, 1)) return 0;
					}else{
						sprintf(msg, "���ļ�ʧ�� --- %s", fullPathFile);
						MessageToDialog(msg);
						return 0;
					}
				}
				/* Add by Galen 2013-02-20 End */

			}
		}

		/* ����Ƕ�̬�ݣ��Ѷ�̬�ݵ�������һ��Btg�ļ�������Ϊ D+tileID.btg */
		if(allDynamicGrass->getNumChildren() > 0)
		{
			/* ���춯̬�ݵ����ݽṹDataStru */
			if(CreateLakeDataStruAndBtg(allDynamicGrass))
			{
				/* �����ļ�����*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\D%d.btg.gz", findResultAC[idx].destPath, tile); //����Ϊ D+tileID.btg.gz
				sprintf(fn, "D%d.btg.gz", tile);
			
				/* дbtg�ļ�������D+tile.btg.gzд�����stg�� */
				if(WriteDataStructToBTG(fullPathFile))
				{
					/* ��̬�ݵ�btg�ļ�����д��stg���� */
					//if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "���ļ�ʧ�� --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}

		/* ����Ƕ�̬�ݣ��Ѷ�̬�ݵ�������һ��Btg�ļ�������Ϊ D+tileID.btg */
		if(allIslands->getNumChildren() > 0)
		{
			/* ���쵺�����ݽṹDataStru */
			if(CreateLakeDataStruAndBtg(allIslands))
			{
				/* �����ļ�����*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\I%d.btg.gz", findResultAC[idx].destPath, tile); //����Ϊ D+tileID.btg.gz
				sprintf(fn, "I%d.btg.gz", tile);
			
				/* дbtg�ļ�������I+tile.btg.gzд�����stg�� */
				if(WriteDataStructToBTG(fullPathFile))
				{
					/* ��̬�ݵ�btg�ļ�����д��stg���� */
					//if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "���ļ�ʧ�� --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}

		/* �����·����·��������һ��Btg�ļ����������и����㣬����Ϊ R+tileID.btg */
		if(allRoadSkirt->getNumChildren() > 0)
		{
			/* �����·�����ݽṹDataStru */
			if(CreateRoadDataStruAndBtg(allRoadSkirt))
			{
				/* �����ļ�����*/
				char msg[STRING_STANDARD_LEN], fn[STRING_STANDARD_LEN], fullPathFile[_MAX_PATH];
				sprintf(fullPathFile, "%s\\R%d.btg.gz", findResultAC[idx].destPath, tile); //����Ϊ tileID+S.btg.gz
				sprintf(fn, "R%d.btg.gz", tile);
			
				/* дbtg�ļ�������R+tile.btg.gzд�����stg�� */
				if(WriteDataStructToBTG(fullPathFile))
				{
					/* ·��btg�ļ���д��stg���� */
					if(!WriteAllStgs(fn, 1)) return 0;
				}else{
					sprintf(msg, "���ļ�ʧ�� --- %s", fullPathFile);
					MessageToDialog(msg);
					return 0;
				}
			}
		}

RSAMB_Continue_1:

		if (ifHasIslandInCurrentTerrain) findResultAC[idx].isOK		= 2;
		else	findResultAC[idx].isOK		= 1;

		/* ������һ�������������ļ� */
		fio = _findnext(ffio, &ff);
	}
	_findclose(ffio);
	
	/* ����ԭ�й���Ŀ¼ */
	io = chdir(path);
	if(io) 
	{
		MessageToDialog("���ܷ���ԭ����Ŀ¼��");
		return 0;
	}
	
	return 1; 
}

/* �ָ�һ���ؾ��ļ�֮ǰ����Ҫ���������ڵ�Ŀ¼ */
int GoIntoDirectoryBeforeRestoreATerrainFile(int lon, int lat)
{
	int io;
	char fullPath[_MAX_PATH];

	/* �ҵ���ǰĿ¼ */
	Bucket Bt;
	char basePath[STRING_STANDARD_LEN];
	set_bucket(&Bt, lon, lat);
	gen_base_path(&Bt, basePath);
	sprintf(fullPath, "%s\\%s", iveroot, basePath);
	
	/* ��ȡһ��ive�ļ����������ĳ��ȡ� */
	io = chdir(fullPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", fullPath);
		MessageToDialog(msg);
		return 0;
	}
	return 1;
}

/* �����Ļ���Ϣ֮ǰ��Ҫ�Ȼָ���allHTI�����г������еؾ��ļ� */
int RestoreTerrainBeforeDoCulture()
{	
	int i, io;
	char ivefname[STRING_STANDARD_LEN], msg[STRING_STANDARD_LEN * 2];
	char path[_MAX_PATH];

	MessageToDialog("�����Ļ���Ϣ֮ǰ����ʼ�ָ�ԭʼ�ؾ��ļ���");
		
	/* ��ȡ��ǰĿ¼ */
	_getcwd(path, _MAX_PATH);

	/* ��ʼ��*/
	allNeedCut.clear();
	
	/* �����Ѿ����õĵؾ����֣��ҵ��ļ����������ļ��������ݻָ���*/
	for(i=0; i<numHTI; i++)
	{
		if(allHTI[i].isOK == 1)
		{
			/* ��tileID�����ҳ����Ⱥ�γ��*/
			int tmp = allHTI[i].tileID / 10000;
			int lon = tmp / 100;
			int lat = tmp - lon * 100;
			InsertItemIntoCutList(lon, lat);

			/* ����ɾ�����Ŀ¼������ǰ�����������Ļ���Ϣ */
			DeleteAllPrevCultureFiles(lon, lat);

			/* ��Ҫ��������Ŀ¼*/
			GoIntoDirectoryBeforeRestoreATerrainFile(lon, lat);

			///* ����ļ�������û�б��ı䣬û�б��ı�Ͳ���Ҫ�ָ���*/
			if(CheckFileLengthBeforeRestore(allHTI[i].tileID))
			{
				/* �Ȼָ��û��߷ֱ��ʵĵؾ� */
				sprintf(ivefname, "%d.ive", allHTI[i].tileID);
				int ret = RestoreFileFromTemp(ivefname);
				if( ret != 0 ) 
				{
					sprintf(msg, "���ڻָ��ļ��� --- %s", ivefname);
					MessageToDialog(msg);
					if( ret == 2)	BackUpFileToTemp(ivefname, ""); 	//����ָ�Ŀ¼û������ļ������ȱ���һ�¡������ڵؾ������ļ����Ȳ�Ϊ701K�������
				}
			}
		}
	}
	
	/* �ٻָ��ͷֱ��ʵĵؾ����� */
	for(std::vector<CutList>::iterator p = allNeedCut.begin(); p != allNeedCut.end(); p++)
	{
		/* ��Ҫ��������Ŀ¼*/
		GoIntoDirectoryBeforeRestoreATerrainFile(p->lon, p->lat);
		
		for(i=0; i<4; i++)
		{
			sprintf(ivefname, "%d%d_%d.ive", p->lon, p->lat, i);
			int ret = RestoreFileFromTemp(ivefname);
			if( ret != 0 ) 
			{
				sprintf(msg, "���ڻָ��ļ��� --- %s", ivefname);
				MessageToDialog(msg);
				if( ret == 2)	BackUpFileToTemp(ivefname, "");	//����ָ�Ŀ¼û������ļ������ȱ���һ�¡������ڵؾ������ļ����Ȳ�Ϊ701K�������
			}
		}
	}

	/* ������ʱ·����*/
	io = chdir(tempDirectory);
	if(!io)
	{
		/* ɾ��������ʱ�ļ� */
		char fullPathFile[_MAX_PATH];
		sprintf(fullPathFile, "%s\\AptLight.ive", tempDirectory);  
		/* ��ɾ��������ļ��� */
		DeleteFile(fullPathFile);
	}


	/* ����ԭ�й���Ŀ¼ */
	io = chdir(path);
	if(io) 
	{
		MessageToDialog("���ܷ���ԭ����Ŀ¼��");
		return 0;
	}
	return 1;
}


/**************************************************
 *
 *   ��ʼ�����Ļ���Ϣ
 *
 **************************************************/

int BeginCreateCulture()
{
	MessageTo = 1;	ifHasIslandInCurrentTerrain = 0;	fRtoCutIsland.clear();

	/* ��·������һ�� */
	UpdateDialog01();
	
	/* ����ȫ�ֱ��� iveroot�� hgtroot */
	SetGlobalVars();
	
	/* ����ǰ�û�������Ϣд��Log�ļ� */
	WriteCultureCfgIntoLogFile();

	/* �����Ļ���Ϣ֮ǰ��Ҫ�Ȼָ���allHTI�����г������еؾ��ļ� */
	RestoreTerrainBeforeDoCulture();

	/* ȫ�������ʼ�� */
	allNeedCut.clear();	allShps.clear();

	/* ���ƽ���������б� */
	csmList.clear();
	
	/* ���ϸ�������ռ��б� */
	allDetailedTxt.clear();

	/* ��һ���֣���ʼ�����Ļ���Ϣ��Ŀǰ֧��ac�ļ���shp�ļ����ָ�ʽ�� */
	if(1)
	{
		MessageToDialog("��ʼ�����Ļ���Ϣ��");

		/* ����·���ҵ�����.ac�ļ� */
		TravelDirectoryToReadAc(baseTextDir);
	
		/* ����·���ҵ�����.shp�ļ� */
		TravelDirectoryToReadShp(baseTextDir);
	
		/* ���ݾ�γ���ҵ��Ļ���Ϣ��������, �����γ�ȶ���0���Ǿ��������еġ�*/
		MessageToDialog("�����Ļ���ϢAC�ļ���");
		if(!DoAllAc())
		{
			MessageToDialog("ACת��BTGʧ�ܡ�");
			goto BCC_Continue_1;
		}
	
		MessageToDialog("�����Ļ���ϢSHP�ļ���");
		if(!DoAllShp())
		{
			MessageToDialog("SHPת��BTGʧ�ܡ�");
			goto BCC_Continue_1;
		}
	}
BCC_Continue_1:		

	/* �ڶ����֣����ݲ�������ͨ�û�����Ȼ������ͨ�û�����������������Ϊ��ѡ� */
	if(CurrentCreateGeneralAirport)
	{
		MessageToDialog("��ʼ����ͨ�û�����");
		if(!DoAllGA())
		{
			MessageToDialog("����ͨ�û���ʧ�ܡ�");
			goto BCC_Continue_2;
		}
	}
BCC_Continue_2:	

	/* �������֣������û���λ���flt����������������� */
	if(1)
	{
		char userPath[_MAX_PATH];
		int retVal;
		
		MessageToDialog("��ʼ��λ�û�������");
	
		/* ֱ�Ӳ����û�������ά����Ŀ¼���Բ����û���ά�������� */
		sprintf(userPath, "%s\\data\\Models\\User", fsData);
	
		/* ������ά����Ŀ¼�������û�������ά���� */
		retVal = TravelDirectoryToReadFlt(userPath);
		if(!retVal)
		{
			MessageToDialog("������ά����Ŀ¼��ȡ��ά������Ϣʧ��\n");
			goto BCC_Continue_3;
		}
		
		/* ��λ�û�������ά���塣*/
		if(!DoAllFlt())
		{
			MessageToDialog("��λ�û�������ά����ʧ�ܡ�");
			return 0;
			//goto BCC_Continue_3;
		}
	}
BCC_Continue_3:
	/* �����õ���Ϣ���½����ݿ� */
	CheckConfig(1);

	/* ���Ĳ��֣�����ǰ�����ֵļ�������ͳһ����ؾ����ڵؾ����и���Ļ���Ϣ��ͨ�û������û��������ڵ��档*/
	if(!DoTerrainCut())
	{
		MessageToDialog("�ؾ��ڶ�ʧ��");
		return 0;
	}

	/* �����ǰ�Ļ���Ϣ�е���ȥ������*/
	if(ifHasIslandInCurrentTerrain)
	{
		if(!DeleteIslandFromCultureBTG()) return 0;
	}

	MessageToDialog("�и�ؾ���ɡ�");
	return 1;	
}


/*============================================================================================
 *  ������ͨ�û������������̡�
 *  ����������ͨ��������ѡ��ģ�壬Ȼ���ģ�����ʵ���������ƽ�ƣ��õ�����btg�ļ���
 *============================================================================================ */

/**************************************************
 *
 *   ��ʼ����ͨ�û���
 *
 **************************************************/
 
/* ��ȡģ��BTG��ת�������дĿ��BTG */
int XchgModelBtgToDestination()
{
	int itmp;
	char templatefile[STRING_STANDARD_LEN];
	char fullpath[STRING_STANDARD_LEN], srcname[STRING_STANDARD_LEN], dstname[STRING_STANDARD_LEN];
	Geodesy cur;
	Bucket Bt;
	char basepath[32];

	sprintf(fullpath, "%s\\generalapt\\Gadata", workDirectory);
	
	/* �Ȼ�ȡģ���ļ��� */
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
	
	/* ��ȡģ��btg�ļ� */
	ReadBTGToDataStruct(srcname);	
	
	/* ת������ */
	cur._lon = currGA.dLongtitude; cur._lat = currGA.dLatitude; cur._elevation = currGA.fAltitude;
	ChangeVertexInDataStru(&cur);
	
	/* дĿ��Btg�ļ�*/
	/* ���ҳ�Ŀ��·����*/
	set_bucket(&Bt, currGA.dLongtitude, currGA.dLatitude);
	gen_base_path(&Bt, basepath);
	sprintf(dstname, "%s\\Scenery\\Terrain\\%s\\%s.btg.gz", currGA.cLocation, basepath, currGA.cICAO);
	WriteDataStructToBTG(dstname);	  
	
	return 1;
}

/*  ����ȹ�ߵ㲢�����һ����Χ�ߵ�osg::Geode�����ڵ������ڳ����� */
int ReadAptOSGAndOutputLineLoopIve()
{
	char msg[STRING_STANDARD_LEN];
	int ilon, ilat;
	Geodesy cur;
	Bucket Bt;
	char basepath[32];
	char dstname[STRING_STANDARD_LEN];

	/* �����ҳ�ȹ�ߵ㼯�� */
	ilon = (int)currGA.dLongtitude; ilat = (int)currGA.dLatitude;
	cur._lon = currGA.dLongtitude; cur._lat = currGA.dLatitude; cur._elevation = currGA.fAltitude;
	FindSkirtVertex(cur._elevation);
	
	/* ��DataStruȥ��ȹ�ߵ��� */
	RemoveSkirtsFromDataStru();

	/* ����ȹ���������Χ�� */
	osg::ref_ptr<osg::Geode> gd = GetAllLineLoopFromSkirtEdge();
		
	/* дΪive�ļ�*/
	/* ���ҳ�Ŀ��·����*/
	set_bucket(&Bt, currGA.dLongtitude, currGA.dLatitude);
	gen_base_path(&Bt, basepath);
	sprintf(dstname, "%s\\Scenery\\Terrain\\%s\\%s.ive", currGA.cLocation, basepath, currGA.cICAO);
	osgDB::writeNodeFile(*gd, dstname);

	/* ���뵽�ڶ��б��У�������ͳһ�ڶ���*/
	if(!InsertItemIntoCutList(ilon, ilat))
	{   
		sprintf(msg, "������ %s �ڶ�ʧ�ܡ�", currGA.cICAO);
		MessageToDialog(msg);
		return 0;
	}   
	return 1;
}

/* ���� MAKEMODELBTG �Σ���ģ������ɹ��̣�ͨ����ȡ.dat�ļ���������Ӧ��.btg�ļ���
   ����ģ���Ѿ�׼����ϣ���ι��̾Ͳ����ˡ� */

//#define MAKEMODELBTG
#ifdef MAKEMODELBTG
int CreateModelHgtdata(char *workDirectory);
int DoGenaptModel(char *fn)
{
	char arg[30][STRING_STANDARD_LEN];
	char *argptr[30];
	int  io;
	char path[_MAX_PATH];

	//��ȡ��ǰĿ¼
	_getcwd(path, _MAX_PATH);
	
	//�������·����Ϣ
	memset(arg, 0, sizeof(arg)); memset(argptr, 0, sizeof(argptr)); 

	//str.Format("@echo off\ngenapts --input=data/airports/ZBAA_UseForVpb.dat"
	//����������
	sprintf(arg[ 0], "geneapt.exe"); 												argptr[ 0] = (char *)&arg[ 0];
	sprintf(arg[ 1], "--input=generalapt\\Ta\\%s"										, fn	); 	argptr[ 1] = (char *)&arg[ 1];
	sprintf(arg[ 2], "--work=generalapt\\work"	);	argptr[ 2] = (char *)&arg[ 2];
	int argIdx = 3;

	/* �����õ�ǰĿ¼Ϊ����Ŀ¼ */
	io = chdir(workDirectory);
	if(io)  return 0;
	
	/* �����ɸ߳����� */
	CreateModelHgtdata(workDirectory);
	
	/* ���ÿ�����Ĺ�������ͨ�û����ļ� */
	if(!CreateAirport(argIdx, argptr))
	{
		MessageToDialog("����ͨ�û���ʧ�ܡ�");
		return 0;
	}

	/* �ع�ԭ��Ŀ¼ */
	io = chdir(path);
	if(io)  return 0;
	
	return 1;
}

/* ģ��׼����������ȡ.dat��ת����.btg */
int ReadDataAndWriteToBTG()
{
	char fullpath[STRING_STANDARD_LEN], fullname[STRING_STANDARD_LEN], msg[STRING_STANDARD_LEN];
	int io, dio, dfio;
	struct _finddata_t sdff;

	sprintf(fullpath, "%s\\generalapt\\Ta", workDirectory);
	
	/* �������Ŀ¼��Ȼ������ͬtile�ļ����������ļ��� */
	io = chdir(fullpath);
	if(io) return 0;
	
	/* ���ҵ�ǰ·���µ������ļ��� */
	dio = _findfirst("*.dat", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		/* ���ļ� */
		sprintf(fullname, "%s", sdff.name);

		sprintf(msg, "����ͨ�û��� --- %s", sdff.name);
		MessageToDialog(msg);
		
		/* */
		DoGenaptModel(fullname);
		
		/* Ѱ����һ�� */
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);

	return 1;
} 

/* ��ȡģ��btg�ļ�������У�顣 */
int ReadModelBTGAndVerify()
{
	char fullpath[STRING_STANDARD_LEN], destpath[STRING_STANDARD_LEN], fullname[STRING_STANDARD_LEN], msg[STRING_STANDARD_LEN];
	int io, dio, dfio;
	struct _finddata_t sdff;

	sprintf(fullpath, "%s\\generalapt\\Gold", workDirectory);
	sprintf(destpath, "%s\\generalapt\\Gadata", workDirectory);
	
	/* �������Ŀ¼��Ȼ������ͬtile�ļ����������ļ��� */
	io = chdir(fullpath);
	if(io) return 0;
	
	/* ���ҵ�ǰ·���µ������ļ��� */
	dio = _findfirst("*.btg.gz", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		/* ���ļ� */
		sprintf(msg, "У��ͨ�û��� --- %s", sdff.name);
		MessageToDialog(msg);
		
		/* */
		sprintf(fullname, "%s\\%s", destpath, sdff.name);
		VerifyTheModels(sdff.name, fullname);

		/* Ѱ����һ�� */
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);

	return 1;
} 
#endif //MAKEMODELBTG


/* ����ͨ�û��������� */
int DoAllGA()
{
	PerZone mPZ;

	/* ������ȫ�ֱ�����һЩ·����*/
	strcpy(currGA.cLocation, currTi.iveDir);
	mPZ.lon = (int)currGA.dLongtitude;  mPZ.lat = (int)currGA.dLatitude;
	currTi.minLongi = (float)(mPZ.lon);
	currTi.minLati  = (float)(mPZ.lat);
	currTi.maxLongi = (float)(mPZ.lon + 1);
	currTi.maxLati  = (float)(mPZ.lat + 1);
	currPZ = &mPZ;

	if(rebuildTerrainForGeneapt)
	{
		/* �������ɻ������� */
		if(!CreatePerDegreeTerrain(0))
		{
			return 0;
		}
	}

	/* ���жϻ������θ߶ȣ�����߶�Ϊ0.0����ôȡ�þ�γ�ȵ����ʵ�߶ȡ�*/	
	if(fabs(currGA.fAltitude - 0.0) < 0.05)
		currGA.fAltitude = GetElevation(currGA.dLongtitude, currGA.dLatitude, 0.0);

	/* ���������У�����ͨ�û��� */
	if(!XchgModelBtgToDestination())
	{
		return 0;
	}

	/* ���ݼ���õ�osg���һ����Χ�ߵ�ive��*/
	if(!ReadAptOSGAndOutputLineLoopIve())
	{
		return 0;
	}
	
	/* дSTG�ļ� */
	if(!WriteAllStgs(currGA.cICAO, 1))
	{
		MessageToDialog("����STG�ļ�ʧ��\n");
		return 0;
	}		

	return 1;
}


/* ͨ�û��������̣����ڰ�ť������ڡ� */
int CreateGeneralAirport()
{

#ifdef MAKEMODELBTG
	//ReadDataAndWriteToBTG();
	ReadModelBTGAndVerify();
	return 1;
#endif //MAKEMODELBTG

	if(!BeginCreateCulture())
	{
		MessageToDialog("�и�ؾ�ʧ�ܣ�����");
	}
	return 1;
}

/*============================================================================================
 *  �������û��Զ��������λ���������̡�Ŀǰ�û�����ֻ֧��flt��ʽ��
 *  ���������Ƕ�ȡflt�ļ����������flt�ļ���������Ȼ����ݾ�γ�ȶ�λ��Ϣ��������ƽ�Ƶ�Ŀ��λ�ã�
 *  ����ڵؾ����ڳ����������
 *============================================================================================ */

/**********************************************************
 *
 *  ����ģ��Ŀ¼����Ѱ.flt�ļ�����ȡ��Ϣ�������ݿ�
 *       
 **********************************************************/
/* ��totalUserBld������ҵ�ǰ�ļ����ڲ��ڣ�*/ 
int CheckFltInUserBld(char *fltname, char *curPath, int *idx)
{
	int i, retval = 0;

	for(i=0; i<bldNum; i++)
	{
		if((!strcmp(fltname, totalUserBld[i].fltFname)) && (!strcmp(curPath, totalUserBld[i].fltPname)))
		{ retval = 1; break; }
	}

	/* �����ǰ��ά�����ļ����ڣ����Ƕ�λ�ļ������ڣ�����2 */
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
 
/* ����Ŀ¼��������flt�ļ� */
/* ������ǰĿ¼�µ�������Ŀ¼���Բ��� flt �ļ���ʹ�õݹ鷽ʽ�� */
char Tpath[_MAX_PATH];
char locname[32];
int TravelSubDirectoryForFlt(char *dir)
{
	int dio, dfio, dirIo, i;
	struct _finddata_t sdff;

	//���ҵ�ǰ·���µ���Ŀ¼��
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		//���������Ŀ¼Ҳ�����ļ�
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR_FLT;

		//�������Ŀ¼�����뵱ǰ��Ŀ¼
		if(sdff.attrib == _A_SUBDIR)
		{

			//������Ŀ¼
			dirIo = chdir(sdff.name);
			_getcwd(Tpath, _MAX_PATH);
			
			//������Ŀ¼
			TravelSubDirectoryForFlt(Tpath);
			
			//�˻ص���һ��Ŀ¼
			dirIo = chdir("..");
		}
		
		//������ļ�
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//�ҵ�Tif�ļ�
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='F')||(sdff.name[i-3]=='f')) && ((sdff.name[i-2]=='L')||(sdff.name[i-2]=='l')) && ((sdff.name[i-1]=='T')||(sdff.name[i-1]=='t')))  
			{
				int idx;
				int  len = strlen(sdff.name);
				FILE *floc;
				Bucket bt;
				char basepath[16];
				LocateBuilding lb;

				/* ��ȡ��ǰĿ¼ */
				_getcwd(Tpath, _MAX_PATH);
				
				/* �Ȳ������ļ��Ƿ��Ѿ��ڿ������ˡ��Ӻ���ǰ��*/
				if(!CheckFltInUserBld(sdff.name, Tpath, &idx))
				{					
					/* ���totalUserBld�������� */
					strcpy(locname, sdff.name); locname[len-3] = 'l'; locname[len-2] = 'o'; locname[len-1] = 'c';
					
					/* �򿪶�λ�ļ�����ȡ��ά����������Ϣ */
					floc = fopen(locname, "rt");
					if(floc != NULL)
					{
						/*	���û��ͬ��loc�ļ���˵����flt��һ���ǻ����ļ������Բ���ӡ�
						/* ��ʼ֧��һ������������ص�� 2012-01-13 */
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
				
				/* ��������ά�����Ѿ����б��У����ǻ�û�б���λ���򲹳䶨λ��Ϣ */
				if(CheckFltInUserBld(sdff.name, Tpath, &idx) == 2)
				{
					/* �򿪶�λ�ļ�����ȡ��ά����������Ϣ */
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
						
						/* ���һ��������λ�ദ�ط� */
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
		//Ѱ����һ��
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* ������ǰĿ¼�����Ҹ�Ŀ¼���Լ�������Ŀ¼����չ��Ϊ.flt���ļ�����ȡ���������Ϣ�� */
int TravelDirectoryToReadFlt(char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForFlt(dPath))
	{
		MessageToDialog("�Ҳ�����������Ϣ�������ļ�");
		return 0;
	}
	CheckConfig(needUpdate);
	
	return 1;
} 

/**************************************************
 *
 *   ��ʼ��λ�û�����
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

/* �ҳ�һ��geode���ж���ľ�γ������ļ�ֵ */
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



/* ���stg�ļ������������Ƿ�����Ҫ�޸ĵĻ�����Ŀ������У����޸ġ�*/
int CheckSTGAndReplaceFltItem(int idx, char *fn)
{
	FILE *fi, *fo;
	std::vector<LocateBuilding> currStg;
	
	/* ���ı�����ʽ�� ��ǰ�� stg �ļ�*/
	fi = fopen(fn, "rt");
	if(fi == NULL)
	{
		MessageToDialog("���ܴ�����stg�ļ���\n");
		return 0;
	}

	/* ����һ���ļ����ݣ��Ƿ�������ȫ�ֻ��������ݡ�ͨ���ļ����ıȽ���ȷ����*/
	currStg.clear();
	while(!feof(fi))
	{
		LocateBuilding lb;
		
		/* ��ȡstg�ļ���ÿһ�� */
		ReadALineToArray(fi);
		strcpy(lb.locFname, cf[0]);
		strcpy(lb.fltPname, cf[1]);
		
		/* �� ·���ļ��� ���浥����ȡ �ļ��� */
		std::string fltPath = lb.fltPname;
		splitPathName(fltPath);
		strcpy(lb.fltFname, materialParts[materialParts.size() - 1].c_str());
		
		/* �Ƚ��ļ����뵱ǰidx��Ӧ�ļ����Ƿ�һ�� */
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
	
	/* ���ı�д�ķ�ʽ�򿪵�ǰ�ļ�����д ���stg */
	fo = fopen(fn, "wt");
	if(fo == NULL)
	{
		MessageToDialog("���ܴ�����stg�ļ���\n");
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


/* ����Ŀ¼��������stg�ļ� */
/* ������ǰĿ¼�µ�������Ŀ¼���Բ��� stg �ļ���ʹ�õݹ鷽ʽ�� */
int TravelSubDirectoryForStg(int idx, char *dir)
{
	int dio, dfio, dirIo, i;
	struct _finddata_t sdff;

	//���ҵ�ǰ·���µ���Ŀ¼��
	dio = _findfirst("*.*", &sdff);
	if((!dio)||(dio==-1)) return 0;
	dfio = 0;

	while(!dfio)
	{
		//���������Ŀ¼Ҳ�����ļ�
		if((sdff.name[0] == '.')) goto GET_NEXT_SUBDIR_STG;

		//�������Ŀ¼�����뵱ǰ��Ŀ¼
		if(sdff.attrib == _A_SUBDIR)
		{

			//������Ŀ¼
			dirIo = chdir(sdff.name);
			_getcwd(Tpath, _MAX_PATH);
			
			//������Ŀ¼
			TravelSubDirectoryForStg(idx, Tpath);
			
			//�˻ص���һ��Ŀ¼
			dirIo = chdir("..");
		}
		
		//������ļ�
		if((sdff.attrib != _A_SUBDIR)&&(sdff.name[0] != '.'))
		{
			//�ҵ�Stg�ļ�
			i = strlen(sdff.name);
			if(((sdff.name[i-3]=='S')||(sdff.name[i-3]=='s')) && ((sdff.name[i-2]=='T')||(sdff.name[i-2]=='t')) && ((sdff.name[i-1]=='G')||(sdff.name[i-1]=='g')))
			{
				/* �Ȳ������ļ������Ƿ����������*/
				if(!CheckSTGAndReplaceFltItem(idx, sdff.name))
				{					
						return 0;
				}
			}
		}

GET_NEXT_SUBDIR_STG:
		//Ѱ����һ��
		dfio = _findnext(dio, &sdff);
	}
	_findclose(dio);		

	return 1;
}

/* ������ǰĿ¼�����Ҹ�Ŀ¼���Լ�������Ŀ¼����չ��Ϊ.Stg���ļ�����ȡ���������Ϣ�� */
int TravelDirectoryToReadStg(int idx, char *dPath)
{
	int io;

	io = chdir(dPath);
	if(io) return 0;

	if(!TravelSubDirectoryForStg(idx, dPath))
	{
		MessageToDialog("�Ҳ���STG�ļ�");
		return 0;
	}

	return 1;
} 


/* BigEndian �� LittleEndian ����������8Bytes��*/
void endian_swap8byte(char *buff)
{
	std::swap(buff[0], buff[7]);	std::swap(buff[1], buff[6]);	std::swap(buff[2], buff[5]);	std::swap(buff[3], buff[4]);
}

/* ����ת�������ݲ������Ӿֲ�����ϵת������������ϵ */
int ConvertCoordination(int idx, osg::ref_ptr<osg::Geode> gnode, std::string fullName)
{
	Geodesy geod;
	Carton  cart;

#if 0		//�ݲ����Ƕ�ȡȫ����������⡣��Ϊ��ʹ�ҵ��˾�γ�ȣ�Ҳ�Ҳ����߶Ⱥͳ���
	/* �����ö����Ʒ�����ȡflt�ļ����鿴�Ƿ���OriginLon��OriginLat������У�����ȫ������ϵ��flt����֮���û�У����Ǿֲ�����ϵ��flt��*/
	/* ��ȡflt�ļ������ԭ�㾭γ�����꣬��ʹ�õĶ������굥λ��*/
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
		/* ����׼λ�ô�Geodesy����תΪCarton����*/
		geod._lon = originLon; geod._lat = originLat; geod._elevation = 0.0;
		GeodToCart(&geod, &cart);
		osg::Vec3d pos; pos.set(cart.x, cart.y, cart.z);

		/* ����Ŀ��ת������ */
		osg::Quat hlOr = fromLonLatRad(originLon, originLat);
		osg::Matrix obj_pos; obj_pos.set(hlOr);
		osg::Quat flip(0.0, 1.0, 0.0, 0.0);
		obj_pos.preMult(osg::Matrix(flip));
		obj_pos.setTrans(pos);

		double hdg = 0.0;
		obj_pos.preMult(osg::Matrix::rotate(hdg * SGD_DEGREES_TO_RADIANS, 0.0, 0.0, 1.0));

		/* ����ȫ��ת���������*/
		coordMatrix = obj_pos;
		rotateMatrix = osg::Matrix::rotate(hdg * SGD_DEGREES_TO_RADIANS, 0.0, 0.0, 1.0);	//��֪��ȫ������ϵ�£�heading��ֵ��ô�õ�������ȱʡ����������0.0

		/* ���ж���ȫ��ת�� */
		ChangeCoordSystem(gnode, obj_pos);
		
		/* Ҫ�޸� Objects Ŀ¼�� ��� stg ��������ֵ��*/
		findResultFlt[idx].fLongi = originLon;	findResultFlt[idx].fLati = originLat;	findResultFlt[idx].fElevation = 0.0;	findResultFlt[idx].fHeading = hdg;
		if(strlen(objroot) == 0) sprintf(objroot, "%s\\Scenery\\Objects", currTi.iveDir);
		TravelDirectoryToReadStg(idx, objroot);
	}
	else
#endif
	{
		/* ����׼λ�ô�Geodesy����תΪCarton����*/
		geod._lon = findResultFlt[idx].fLongi; geod._lat = findResultFlt[idx].fLati; geod._elevation = findResultFlt[idx].fElevation;
		GeodToCart(&geod, &cart);
		osg::Vec3d pos; pos.set(cart.x, cart.y, cart.z);
	
		/* ����Ŀ��ת������ */
		osg::Quat hlOr = fromLonLatRad(findResultFlt[idx].fLongi, findResultFlt[idx].fLati);
		osg::Matrix obj_pos; obj_pos.set(hlOr);
		osg::Quat flip(0.0, 1.0, 0.0, 0.0);
		obj_pos.preMult(osg::Matrix(flip));
		obj_pos.setTrans(pos);
	
		double hdg = findResultFlt[idx].fHeading;
		obj_pos.preMult(osg::Matrix::rotate(hdg * SGD_DEGREES_TO_RADIANS, 0.0, 0.0, 1.0));
			
		/* ����ȫ��ת���������*/
		coordMatrix = obj_pos;
		rotateMatrix = osg::Matrix::rotate(hdg * SGD_DEGREES_TO_RADIANS, 0.0, 0.0, 1.0);
	
		/* ���ж���ȫ��ת�� */
		ChangeCoordSystem(gnode, obj_pos);
		
		/* �ҳ����ж���ľ�γ�ȼ�ֵ�����ڿ������������ʵ�����塣Ȩ�ҷ�������� */
		//InitMMLLEx(&currFlt, findResultFlt[idx].fltFname);
		//FindMMLLInGeode(gnode);
	}

	return 1;
}

/* ����CorrectionSpaceMesh */
void AddCorrectionSpaceMesh(osg::ref_ptr<osg::Group> gp)
{
	osg::ref_ptr<CorrectionSpaceMesh> csm = new CorrectionSpaceMesh( gp, 1.0 );
	csmList.push_back( csm );
}
 
 
/* ����idxָ����findResultFlt��������ȡ�û���ά����flt�ļ����������������Χ��Geometry, ���뵽gnode���档*/
int ReadFltAndComputePerimeter(int idx, int lon, int lat)
{
	int io, ffio;
	char path[_MAX_PATH], dPath[_MAX_PATH], fullName[_MAX_PATH], iveName[24], currTerPath[STRING_STANDARD_LEN * 2], basepath[32], fullTrueName[_MAX_PATH];
	struct _finddata_t fltff;
	
	/* ȷ����ǰҪ�����ڵؾ����е�·�� terPath */
	Bucket Bt;
	double fLon = lon + 0.0000001;
	double fLat = lat + 0.0000001;
	set_bucket(&Bt, fLon, fLat);
	gen_base_path(&Bt, basepath);
	if(strlen(iveroot) == 0) sprintf(iveroot, "%s\\Scenery\\Terrain", currTi.iveDir);
	sprintf(currTerPath, "%s\\%s", iveroot, basepath);

	/* ��ȡ�����浱ǰĿ¼ */
	_getcwd(path, _MAX_PATH);

	/* ����flt�ļ����ڵ�Ŀ¼ */
	sprintf(dPath, "%s", findResultFlt[idx].fltPname);
	io = chdir(dPath);
	if(io) 
	{
		char msg[STRING_STANDARD_LEN];
		sprintf(msg, "���ܽ���Ŀ��Ŀ¼ -- %s", dPath);
		MessageToDialog(msg);
		return 0;
	}

	/* �������Ŀ¼�µ�flt�ļ��� */
	ffio = _findfirst(findResultFlt[idx].fltFname, &fltff);	
	if((!ffio)||(ffio == -1))
	{
		MessageToDialog("�Ҳ���flt�ļ���");
		return 0;
	}

	/* �Ȳ��ҵ�ǰ·�����Ƿ���Around.flt�ļ�������У���ʹ��Around.flt��Ϊ�����û�������Ե�������ļ���*/
	ffio = _findfirst("Around.flt", &fltff);	
	if((!ffio)||(ffio == -1))
	{
		/* ƴװȫ·���ļ��� */
		sprintf(fullName, "%s\\%s", findResultFlt[idx].fltPname, findResultFlt[idx].fltFname);
	}else{
		/* ƴװȫ·���ļ��� Around.flt */
		sprintf(fullName, "%s\\Around.flt", findResultFlt[idx].fltPname);
	}
	_findclose(ffio);

	/* ����һ��gnode, ����Ӹ�����ά�����perimeter�� */
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	gd->setName("River");
	
	/* �����û�����flt�ı�Ե���������ڵؾ����ڶ������������Ե����ΪGeometry�ӵ�gnode���� */
	if(!ReadFltAndVistIt(fullName, gd))
	{
		if(!foundAround)
			MessageToDialog("����flt����û��Around�飬���ܼ����Ե��");
		else
			MessageToDialog("������Եû�з�գ�������ȷ�и");
		MessageToDialog("�����û���ά������Χ��ʧ�ܡ�");
		return 0;
	}

	/* �����û�����flt�ĵƹ���Ϣ����д��ƹ���ʱ�ļ� AptLight.ive*/
	/* ����ƴװȫ·���ļ��� */
	sprintf(fullTrueName, "%s\\%s", findResultFlt[idx].fltPname, findResultFlt[idx].fltFname);
	if(!ReadFltAndCreateAptLight(fullTrueName))
	{
		char msg[_MAX_PATH];
		sprintf(msg, "��ǰ����û�еƹ���Ϣ --- %s", fullTrueName);
		MessageToDialog(msg);
	}

	
	/* ת������ */
	if(!ConvertCoordination(idx, gd, fullName))
	{
		MessageToDialog("ת������ʧ�ܡ�");
		return 0;
	}

	/* ת��������󣬾Ϳɼ���CorrectionSpaceMesh�� */
	osg::ref_ptr<osg::Group> gp = new osg::Group;
	gp->addChild(gd);
#if 1
	AddCorrectionSpaceMesh(gp);
#endif

	/* ȷ����ǰ�������� */
	if(findResultFlt[idx].iCode == 0)
	{
		SetCurrFltAPCode(idx);
	}

	/* �����gnode���������ͨ�û���������ive�ļ� �ļ���Ϊ GUR?.ive*/
	sprintf(iveName, "GUR%c.ive",  findResultFlt[idx].iCode);
	sprintf(fullName, "%s\\%s", currTerPath, iveName);
	osgDB::writeNodeFile(*gd, fullName);
		
	/* д�����ƹ�BTG */
	sprintf(fullName, "%s\\%s", findResultFlt[idx].fltPname, findResultFlt[idx].fltFname);
	if(ReadFltAndCreateLightsDataStru(fullName))
	{
		char fn[STRING_STANDARD_LEN],  msg[_MAX_PATH];
		char fullPathFile[_MAX_PATH];

		sprintf(iveName, "GUR%c.btg.gz",  findResultFlt[idx].iCode);
		sprintf(fullPathFile, "%s\\%s", currTerPath, iveName);
		strcpy(fn, iveName);
	
		/* �Ȱ�tile.btg.gzд�����stg�� */
		if(WriteDataStructToBTG(fullPathFile))
		{
			if(!WriteAllStgs(fn, 1)) return 0;
		}else{
			sprintf(msg, "д������ƹ��ļ�ʧ�� --- %s", fullPathFile);
			MessageToDialog(msg);
		}		
	}

	return 1;
}

/* �������е��û���ά�������� */
int DoAllFlt()
{
	int i, j, k;
	char msg[STRING_STANDARD_LEN];
	int currDoing[16];
	int numCD;
	PerZone mPZ[4];
	int numPZ, iFlag;

	/* �ҳ���û�ж�λ�û���ά�������ڵĵؾ�Ŀ¼���� */
	if(FindOutDirsNoBldDone() == 0)
	{
		MessageToDialog("û����Ҫ��λ����ά����.");
		return 1;	
	}

	/* �ȴ�dirsNoCultDone���𣬼���ЩĿ¼���涼��û���û���ά����Ķ�λ */
	for(i=0; i<numDNCD; i++)
	{
		/* ���dirsNoCultDone����ÿһ���totalUserBld��������ͬTerrainĿ¼������ϵ�findResultFlt*/
		memset(currDoing, 0, sizeof(currDoing));  numCD = 0;
		memset(mPZ, 0, sizeof(mPZ));  numPZ = 0;
		if(FindAllFltItemsWithSameDir(i) == 0) continue;
		strcpy(findJpgDirectory, findResultFlt[0].terPname);
		
		/* ���findResultFlt����ÿһ��ȼ���Ƿ������û���ά���嶨λ����û��������currDoing���¼��*/
		for(j=0; j<numFRFLT; j++)
		{
			if(findResultFlt[j].isOK == 0)
			{	
				currDoing[numCD++] = j;
			}
		}
		
		/* �Ӹո����õ����п������ҳ����в�ͬ�Ŀ� */
		for(j=0; j<numCD; j++)
		{
			for(int n=0; n<4; n++)
			{
				/* ȷ���Ի�������Ϊԭ�㣬һ����Χ�ڷ��������ĸ��ǵľ�γ��ֵ */
				double cornerLon, cornerLat;
				switch(n)
				{
					case 0:	{	cornerLon = findResultFlt[currDoing[j]].fLongi - 0.05; cornerLat = findResultFlt[currDoing[j]].fLati - 0.025; break;	}
					case 1:	{	cornerLon = findResultFlt[currDoing[j]].fLongi - 0.05; cornerLat = findResultFlt[currDoing[j]].fLati + 0.025; break;	}
					case 2:	{	cornerLon = findResultFlt[currDoing[j]].fLongi + 0.05; cornerLat = findResultFlt[currDoing[j]].fLati - 0.025; break;	}
					case 3:	{	cornerLon = findResultFlt[currDoing[j]].fLongi + 0.05; cornerLat = findResultFlt[currDoing[j]].fLati + 0.025; break;	}
				}				
				
				/* �ҳ����������ռ�ݵ�������γ������ĸ��� */
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

		/* ����ÿ���飬���������ؾ����ټӶ�λ��ά���� */
		for(j=0; j<numPZ; j++)
		{
			/* ������ʼ��ȫ�ֱ��� */
			currTi.minLongi = (float)(mPZ[j].lon);
			currTi.minLati  = (float)(mPZ[j].lat);
			currTi.maxLongi = (float)(mPZ[j].lon + 1);
			currTi.maxLati  = (float)(mPZ[j].lat + 1);
			currPZ = &mPZ[j];
			

			/* ����currDoing�������Ϣ�������û���ά����ĵ���Ǳ�Ե */
			for(k=0; k<numCD; k++)
			{
				if(findResultFlt[currDoing[k]].isOK == 0)
				{
					sprintf(msg, "���ڼ����Ե--%s...", findResultFlt[currDoing[k]].fltFname);
					MessageToDialog(msg);
					if(!ReadFltAndComputePerimeter(currDoing[k], mPZ[j].lon, mPZ[j].lat))
					{
						sprintf(msg, "�����Ե--%sʧ�ܡ�", findResultFlt[currDoing[k]].fltFname);
						MessageToDialog(msg);
						return 0;
					}
				}
			}   

			/* ����ͬ���ڵ�����btg���һ��Geode, ���뵽�ڶ��б��У�������ͳһ�ڶ��� */
			if(!InsertItemIntoCutList(mPZ[j].lon, mPZ[j].lat))
			{   
				sprintf(msg, "ת��IVEʧ�ܡ�");
				MessageToDialog(msg);
				return 0;
			}
		}

		/* ������ɲ������������isOK��״̬ */
		for(k=0; k<numCD; k++)
		{
			if(findResultFlt[currDoing[k]].isOK == 0)
			{
				findResultFlt[currDoing[k]].isOK = 1;
			}
		}

		/* ���õ��Ļ���Ϣ������������ݿ� */
		UpdateTotalUserBldWithNewFindResult();
	}
	
	return 1;
}

 
/* ��λ�û����������� */
int LocateUserAirport()
{
	return 1;
}
 
