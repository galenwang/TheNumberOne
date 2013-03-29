#include "bucket.h"

/* 初始化检查表*/
typedef struct _INITIALIZATIONCHECK_
{
	bool tempDirectory;			//临时路径
	bool backDirectory;			//原始库备份路径
	bool usrConfigFile;			//用户配置文件
	bool pathOfTexture;			//纹理素材路径
	bool pathOfHgt;					//高程素材路径
	bool pathOfTarget;			//生成结果路径
}InitializationCheck;

/* 纹理图信息结点 */
typedef struct _MAPNODE_
{
	double minLong;
	double maxLong;
	double minLati;				
	double maxLati;
	char   fileName[32];
	char   pathName[128];
	unsigned short image_width;
	unsigned short image_height;
	int	   ID;
}MapNode;


typedef struct _HIGHNODE_
{
	double minLong;
	double maxLong;
	double minLati;				
	double maxLati;
	char   fileName[32];
	char   pathName[128];
	unsigned short image_width;
	unsigned short image_height;
	int    highLevel;
	int	   ID;
}HighNode;


typedef struct _CULTNODE_
{
	double minLong;
	double maxLong;
	double minLati;				
	double maxLati;
	char   fileName[32];
	char   pathName[128];
	int	   destTileID;
	int    destLongi;
	int    destLati;
	char   destFile[32];
	char   destPath[128];	
	int	   ID;
	int	   isOK;
}CultNode;


typedef struct _TERRAININPUT_
{
	char  txtDir[64];
	char  hgtDir[64];
	char  iveDir[64];
	char  resolution[8];
	double minLongi;
	double maxLongi;
	double minLati;
	double maxLati;
	int   lod;
	float scale;
}TerrainInput;

typedef struct _GENERALAIRPORT_
{
	short 	shGenericAirportState;		//通用机场状态		[0，1]
	short 	shRunwayID;              	//通用机场跑道编号		[01，36]
	short 	shRunwayWidth;           	//通用机场跑道宽度		[0，2]
	short 	shRunwayCondition;       	//通用机场跑道条件		[0，3]
	short 	shRunwayLength;          	//通用机场跑道长度		[0，7]
	short 	shBuildingPosition;      	//通用机场建筑物位置		[0，1]
	short 	shBuildingTypes;         	//通用机场建筑物类型		[0，1]
	short 	shAppLightTypes;         	//通用机场跑道进近灯光类型		[0，2]
	short 	shVLATypes;              	//通用机场跑道助降灯光类型		[0，2]
	double	dLatitude;               	//通用机场跑道纬度	deg	[-90.0，90.0]
	double	dLongtitude;             	//通用机场跑道经度	deg	[-180.0，180.0]  
	float 	fAltitude;               	//通用机场跑道高度	ft	[-500.0，50000.0]
	float 	fHeading;                	//通用机场跑道航向角	deg	[0.0，360.0]     
	float 	fVLAAngles;					//通用机场跑道助降灯光角度	deg	[2.4，4.0] 
	float 	fRunwaySlope;            	//通用机场跑道倾斜角度	deg	[-0.5，0.5]  	
	char    cICAO[5];					//通用机场虚拟ICAO代码
	char    cLocation[64];				//FlightGear所在位置
	char    cOrgLoc[64];				//用来保存配置参数或者文件中读取的cLocation相对路径;
}GeneralAirport;

#define MAXDIRECTORIES 1000

/* 批量单ive制作地形所需变量及数据结构 */
typedef struct _PERZONE_
{
	int lon;
	int lat;
	int isOK;
}PerZone;
typedef struct _TEXTDIRINFO_
{
	char  	dirName[64];
	PerZone	zone[4];
}TextDirInfo;

//管理高精度地景块制作列表
typedef struct _HIGHTERRAININFO_
{
	char  	texFname[256];
	char	hgtFname[256];
	float	startLongi;
	float	startLati;
	float	offsetLongi;
	float	offsetLati;
	int		highLevel;			//精度级别
	int		tileID;
	int     ID;
	int		isOK;
}HighTerrainInfo;


//管理高精度文化信息块制作列表
typedef struct _HIGHCULTUREINFO_
{
	char  	texFname[256];
	float	startLongi;
	float	startLati;
	float	offsetLongi;
	float	offsetLati;
	int		tileID;
	int     ID;
	int		isOK;
}HighCultureInfo;


/* 管理用户三维数据定位的数据结构 */
#define MAXBUILDING  200
typedef struct _LOCATEBUILDING_
{
	char    fltPname[256];		//三维物体的路径名。
	char    fltFname[32];		//三维物体的文件名， 指定flt格式
	char    locFname[32];		//定位文件的文件名， 指定扩展名为loc， 指定同三维物体文件在同一个路径下面。
	char    terPname[256];		//所在地景区域的路径名。
	double  fLongi;
	double  fLati;
	double  fElevation;
	double  fHeading;
	char	iCode;				//机场代号，为 A, B, C, D......
	int     ID;
	int		isOK;
}LocateBuilding;

#define MAX_HGT_SIZE 1201	//3601
typedef struct _TGhgt_
{
	int curLong, curLati;
  double originx, originy;		// coordinates (in arc seconds) of south west corner
  int cols, rows;							// number of columns and rows
  double col_step, row_step;	// Distance between column and row data points (in arc seconds)
  FILE *fd;										// file pointer for input
  int hgt_resolution;
  short int data[MAX_HGT_SIZE][MAX_HGT_SIZE];	// pointers to the actual grid data allocated here
}TGhgt;


int CheckDirectories();

int CheckNodes();

int TravelSubDirectoryForMap(char *dir);

int TravelDirectoryToReadMap(char *dPath);

int CheckConfig(int rewrite);

int SetGlobalVars();

int CalcTilesPics();

int MergePictures();

int convertTextureGeotiff();

bool GetHgtFileList(double leftlongi, double downlati, double rightlongi, double upperlati);

bool MergeHgt();

int convertHgtGeotiff();

int MakeSTGFile(Bucket *Bt, int main);

int DoVpb(Bucket *b, int type);

int GetCurrNodeByName(MapNode *pnode, char *fn);

int FindJpgsFromBucket(Bucket &b);

int SortfindResult();

int BatchCreateTerrain();

int BeginCreateTerrain();

int initCreateTerrain();

int BeginCreateCulture();

int BeginCreateNeededCulture();

int CreateGeneralAirport();

int LocateUserAirport();
