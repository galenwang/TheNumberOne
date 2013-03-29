#include "bucket.h"

/* ��ʼ������*/
typedef struct _INITIALIZATIONCHECK_
{
	bool tempDirectory;			//��ʱ·��
	bool backDirectory;			//ԭʼ�ⱸ��·��
	bool usrConfigFile;			//�û������ļ�
	bool pathOfTexture;			//�����ز�·��
	bool pathOfHgt;					//�߳��ز�·��
	bool pathOfTarget;			//���ɽ��·��
}InitializationCheck;

/* ����ͼ��Ϣ��� */
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
	short 	shGenericAirportState;		//ͨ�û���״̬		[0��1]
	short 	shRunwayID;              	//ͨ�û����ܵ����		[01��36]
	short 	shRunwayWidth;           	//ͨ�û����ܵ����		[0��2]
	short 	shRunwayCondition;       	//ͨ�û����ܵ�����		[0��3]
	short 	shRunwayLength;          	//ͨ�û����ܵ�����		[0��7]
	short 	shBuildingPosition;      	//ͨ�û���������λ��		[0��1]
	short 	shBuildingTypes;         	//ͨ�û�������������		[0��1]
	short 	shAppLightTypes;         	//ͨ�û����ܵ������ƹ�����		[0��2]
	short 	shVLATypes;              	//ͨ�û����ܵ������ƹ�����		[0��2]
	double	dLatitude;               	//ͨ�û����ܵ�γ��	deg	[-90.0��90.0]
	double	dLongtitude;             	//ͨ�û����ܵ�����	deg	[-180.0��180.0]  
	float 	fAltitude;               	//ͨ�û����ܵ��߶�	ft	[-500.0��50000.0]
	float 	fHeading;                	//ͨ�û����ܵ������	deg	[0.0��360.0]     
	float 	fVLAAngles;					//ͨ�û����ܵ������ƹ�Ƕ�	deg	[2.4��4.0] 
	float 	fRunwaySlope;            	//ͨ�û����ܵ���б�Ƕ�	deg	[-0.5��0.5]  	
	char    cICAO[5];					//ͨ�û�������ICAO����
	char    cLocation[64];				//FlightGear����λ��
	char    cOrgLoc[64];				//�����������ò��������ļ��ж�ȡ��cLocation���·��;
}GeneralAirport;

#define MAXDIRECTORIES 1000

/* ������ive��������������������ݽṹ */
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

//����߾��ȵؾ��������б�
typedef struct _HIGHTERRAININFO_
{
	char  	texFname[256];
	char	hgtFname[256];
	float	startLongi;
	float	startLati;
	float	offsetLongi;
	float	offsetLati;
	int		highLevel;			//���ȼ���
	int		tileID;
	int     ID;
	int		isOK;
}HighTerrainInfo;


//����߾����Ļ���Ϣ�������б�
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


/* �����û���ά���ݶ�λ�����ݽṹ */
#define MAXBUILDING  200
typedef struct _LOCATEBUILDING_
{
	char    fltPname[256];		//��ά�����·������
	char    fltFname[32];		//��ά������ļ����� ָ��flt��ʽ
	char    locFname[32];		//��λ�ļ����ļ����� ָ����չ��Ϊloc�� ָ��ͬ��ά�����ļ���ͬһ��·�����档
	char    terPname[256];		//���ڵؾ������·������
	double  fLongi;
	double  fLati;
	double  fElevation;
	double  fHeading;
	char	iCode;				//�������ţ�Ϊ A, B, C, D......
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
