
typedef struct tagXY
{
	//ProEnv
	int XSum;
	int YSum;
	int X;
	int Y;
}XY;



typedef struct tagProEnv
{
	//ProEnv
	char DEMPathName[128];
	char TexPathName[128];
	char ivePathName[128];
	int Lod;
	char Type[64];
	int DateSum;
	float Scalse;  
	float lon;
	float lat;
	float LonRage;
	float LatRage;
}PROENV;

typedef struct tagRunwayProp
{
	//RunwayProp
	int Lenght;
	int Width;
	float Heading;
	float Longitude;
	float Latitude;
	int m_RunwayAlt;
	int DispThr1;
	int DispThr2;
	int Stopway1;
	int Stopway2;
	char RunwayLight1[16];
	char RunwayLight2[16];
	char ApproachLight1[16];
	char ApproachLight2[16];
	char GlideslopeLight1[16];
	char GlideslopeLight2[16];
	char SurfaceType[16];
	char Marking[16];
	char Should[16];
	float Roughness;
	BOOL Signage; 
}RUNWAYPROP;

typedef struct tagLatLon
{
	//LatLon
	char ACPathName[128];
	float LongitudeMAX;
	float LongitudeMIN;
	float LatitudeMAX;
	float LatitudeMIN;  
}LATLON;